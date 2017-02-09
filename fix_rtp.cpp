/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <string.h>
#include <stdlib.h>
#include "fix_rtp.h"
#include "atom.h"
#include "atom_masks.h"
#include "accelerator_kokkos.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

FixRtp::FixRtp(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
//  fprintf(screen,"Narg %d",narg);
  if (narg != 9) error->all(FLERR,"Illegal fix rtp command");
  //narg is the number of arguement, start from 0

  dynamic_group_allow = 1;
  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;

// check whereas the force is a constant or a variable to be computed;
// need to use a constant here;
// there are 3 values fx fy fz
// we need just the magnitude, f, and an integer, the length of the molecule

//  xstr = ystr = zstr = NULL;
//  fvalue = force->numeric(FLERR,arg[3]); // magnitude of the force
//  mol_length = force->numeric(FLERR,arg[4]); // apply the force to particles with index n*mol_length+1

  fvalue = force->numeric(FLERR,arg[3]);  // transform string to float
  ftvalue = force->numeric(FLERR,arg[4]); // force applied to impose the torque
  tau= force->numeric(FLERR,arg[5]);// transform string to inverse rate
  t_tumble = force->numeric(FLERR,arg[6]);//transform string to tumble duration time
  mol_length = force->inumeric(FLERR,arg[7]); // transform string to integer
  seed = force->inumeric(FLERR,arg[8]); // seed of random number generator

  random = new RanMars(lmp,seed + comm->me);

  // optional args
  nevery = 1;
  iregion = -1;
  idregion = NULL;
  estr = NULL;

// check for error related to the number of arguments
// I change iarg to 7, if there is no arguments behind, might no problem.
  int iarg = 9;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix Rtp command");
      nevery = atoi(arg[iarg+1]);
      if (nevery <= 0) error->all(FLERR,"Illegal fix Rtp command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix Rtp command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR,"Region ID for fix Rtp does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"energy") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix Rtp command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
        int n = strlen(&arg[iarg+1][2]) + 1;
        estr = new char[n];
        strcpy(estr,&arg[iarg+1][2]);
      } else error->all(FLERR,"Illegal fix Rtp command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix Rtp command");
  }

  force_flag = 0;
  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;

  // I don't understand what is sforce
  maxatom = 1;
  memory->create(sforce,maxatom,4,"rtp:sforce");

  
}

/* ---------------------------------------------------------------------- */

FixRtp::~FixRtp()
{
//  delete [] xstr;
//  delete [] ystr;
//  delete [] zstr;
  delete [] estr;
  delete [] idregion;
  memory->destroy(sforce);
}

/* ---------------------------------------------------------------------- */

int FixRtp::setmask()
{
  datamask_read = datamask_modify = 0;

  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixRtp::init()
{
  // check variables

  // check whereas the string for the energy is ok
  if (estr) {
    evar = input->variable->find(estr);
    if (evar < 0)
      error->all(FLERR,"Variable name for fix Rtp does not exist");
    if (input->variable->atomstyle(evar)) estyle = ATOM;
    else error->all(FLERR,"Variable for fix Rtp is invalid style");
  } else estyle = NONE;

  // set index and check validity of region

  if (iregion >= 0) {
    iregion = domain->find_region(idregion);
    if (iregion == -1)
      error->all(FLERR,"Region ID for fix Rtp does not exist");
  }

  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixRtp::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixRtp::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixRtp::post_force(int vflag)
{
  //printf("Rtppost force\n"); getchar();
  double **x = atom->x;
  double **f = atom->f;
  double **v = atom->v;
  double **Ft = atom->Ft; 
  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;
  int *state= atom->state;// not sure about here
  int *t_elapsed= atom->t_elapsed;// initial should be zero
  int *t_transition= atom->t_transition;// initial should be zero

  if (update->ntimestep % nevery) return;

  if (lmp->kokkos)
    atom->sync_modify(Host, (unsigned int) (F_MASK | MASK_MASK),
                      (unsigned int) F_MASK);

  // update region if necessary

  Region *region = NULL;
  if (iregion >= 0) {
    region = domain->regions[iregion];
    region->prematch();
  }

  // foriginal[0] = "potential energy" for added force
  // foriginal[123] = force on atoms before extra force added

  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;
  force_flag = 0;

  // constant force
  // potential energy = - x dot f in unwrapped coords

  // this is what we have to modify; i loops over the local;
  // we need to:
  // a) check if the particle is part of a molecule
  // b) if it is the first of the molecule, get the coordinate of the second,
  // the direction, and add the force
  // check ilist and list in pair_repuslive_q.cpp

  double unwrap[3];
  double unwrap2[3];
  double Ftorque[3];
  double dx,dy,dz,nx,ny,nz, dist;
  int j;
  int gid;
  bool new_torque;
  for(int i = 0; i < nlocal; i++){ // i is the local index, not the id
    if(mask[i] & groupbit){ // if the particle is in the right group
      if(region && !region->match(x[i][0],x[i][1],x[i][2])) continue;
      gid=atom->tag[i]; // global id of atom
      if( gid-mol_length*int(gid/mol_length) == 1){    // this is the first particle of a molecule
        domain->unmap(x[i],image[i],unwrap); // x of the head
        j = atom->map(gid+1); 
        domain->unmap(x[j],image[j],unwrap2); // x of the following particle
        // compute the unnormalized versor
        dx = unwrap[0]-unwrap2[0];
        dy = unwrap[1]-unwrap2[1];
        dz = unwrap[2]-unwrap2[2];
        new_torque = false;
        if(t_transition[i] < 0){
          if(random->uniform() < tau/(tau+t_tumble) ) state[i] = 0; else state[i] = 1;
        }
        if(t_elapsed[i]>t_transition[i]){// if change state, count the timesteps from 0
          if( ( (state[i]==0) && (t_tumble > 0)) || // from running to tumbling state
              (  (state[i]==1) && (tau == 0)  ) ){ // from tumbling to tumbling
            state[i]=1;
            t_transition[i] = t_tumble;
            // random axis for the force used to impose the torque
            nx = random->uniform()-0.5;
            ny = random->uniform()-0.5;
            nz = random->uniform()-0.5;
             // random direction ortogonal to the molecule axis
            Ft[i][0] = dy*nz - dz*ny;
            Ft[i][1] = dz*nx - dx*nz;
            Ft[i][2] = dx*ny - dy*nx;
            if(domain->dimension == 2) Ft[i][2] = 0;
            // get the component of the force to be applied to fix the torque
            dist = sqrt(Ft[i][0]*Ft[i][0] + Ft[i][1]*Ft[i][1] + Ft[i][2]*Ft[i][2]);
            Ft[i][0] = ftvalue*Ft[i][0]/dist;
            Ft[i][1] = ftvalue*Ft[i][1]/dist;
            Ft[i][2] = ftvalue*Ft[i][2]/dist;
            // describe the force as a vector starting from the first particle of the molecule
            Ft[i][0] = Ft[i][0]+unwrap[0];
            Ft[i][1] = Ft[i][1]+unwrap[1];
            Ft[i][2] = Ft[i][2]+unwrap[2];
            new_torque = true;
          }else if( ( (state[i]==1) && (tau > 0)) ||  // from tumbling to running
                    ( (state[i]==0) && (t_tumble == 0)) ) { // from running to running
            state[i]=0;
            t_transition[i]= int(-1.0/tau*log(random->uniform()));
          }
          t_elapsed[i] = 0; // reset the time
        }
        if(state[i]== 0 ){ // if the molecule is in state 0, which means running
          // normalize the versor
          dist = sqrt(dx*dx + dy*dy + dz*dz);
          dx /= dist;
          dy /= dist;
          dz /= dist;
          // now we need to find the position of particles i and j, and bla bla
          foriginal[0] -= xvalue*unwrap[0] + yvalue*unwrap[1] + zvalue*unwrap[2]; // energy change
          foriginal[1] += fvalue*dx;
          foriginal[2] += fvalue*dy;
          foriginal[3] += fvalue*dz;
          f[i][0] += fvalue*dx;
          f[i][1] += fvalue*dy;
          f[i][2] += fvalue*dz;
//          if(comm->me == 0){ fprintf(screen,"Add linear force %f %f %f - %f %f %f - %f \n",dx,dy,dz,f[i][0],f[i][1],f[i][2], fvalue); getchar();}
        }else {
          if(!new_torque){// the torque change every new tumble state , and move a little bit in the tumble state
            Ft[i][0] = Ft[i][0] + v[i][0]*update->dt; // move the endpoint of the force applied to impose the torque as the first molecule of the rigid rod
            Ft[i][1] = Ft[i][1] + v[i][1]*update->dt;
            Ft[i][2] = Ft[i][2] + v[i][2]*update->dt;
          }
          Ftorque[0] = Ft[i][0]-unwrap[0];
          Ftorque[1] = Ft[i][1]-unwrap[1];
          Ftorque[2] = Ft[i][2]-unwrap[2];
          if(domain->dimension == 2){ Ftorque[2] = 0;}

          // now we need to find the position of particles i and j, and bla bla
          foriginal[0] -= xvalue*unwrap[0] + yvalue*unwrap[1] + zvalue*unwrap[2]; // energy change
          foriginal[1] += Ftorque[0];
          foriginal[2] += Ftorque[1];
          foriginal[3] += Ftorque[2];
          f[i][0] += Ftorque[0];
          f[i][1] += Ftorque[1];
          f[i][2] += Ftorque[2];
          j = atom->map(gid+mol_length-1); // local index of the atom whose global index marks the end of the molecule
          f[j][0] -= Ftorque[0];
          f[j][1] -= Ftorque[1];
          f[j][2] -= Ftorque[2];
 //         if(comm->me == 0) fprintf(screen,"Add torque\n");
        }
        t_elapsed[i]+=1; // when this fix rtp run for one time, times+1
      }
    }
  }



}

/* ---------------------------------------------------------------------- */

void FixRtp::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixRtp::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   potential energy of added force
------------------------------------------------------------------------- */

double FixRtp::compute_scalar()
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return foriginal_all[0];
}

/* ----------------------------------------------------------------------
   return components of total force on fix group before force was changed
------------------------------------------------------------------------- */

double FixRtp::compute_vector(int n)
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return foriginal_all[n+1];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixRtp::memory_usage()
{
  double bytes = 0.0;
  if (varflag == ATOM) bytes = maxatom*4 * sizeof(double);
  return bytes;
}
