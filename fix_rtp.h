/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(rtp,FixRtp)

#else

#ifndef LMP_FIX_RTP_H
#define LMP_FIX_RTP_H

#include "fix.h"
#include "random_mars.h"

namespace LAMMPS_NS {

class FixRtp: public Fix {
 public:
  FixRtp(class LAMMPS *, int, char **);
  ~FixRtp();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double compute_scalar();
  double compute_vector(int);
  double memory_usage();

 private:
  double xvalue,yvalue,zvalue;
  double fvalue; // value of the pushing force; we then compute the components x,y,z
  double ftvalue; // value of the force used to impose the torque; we add +tvalue and -tvalu to the extreme of the molecule
  double tau, t_tumble; // characteristic running and tumbling times
  int mol_length;
  int varflag,iregion;
  char *estr;
  char *idregion;
  int xvar,yvar,zvar,evar,xstyle,ystyle,zstyle,estyle;
  double foriginal[4],foriginal_all[4];
  int force_flag;
  int nlevels_respa;

  int maxatom;
  double **sforce;

  protected:
  class RanMars *random;
  int seed;

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Region ID for fix Rtp does not exist

Self-explanatory.

E: Variable name for fix Rtp does not exist

Self-explanatory.

E: Variable for fix Rtp is invalid style

Self-explanatory.

E: Cannot use variable energy with constant force in fix torque

This is because for constant force, LAMMPS can compute the change
in energy directly.

E: Must use variable energy with fix torque

Must define an energy vartiable when applyting a dynamic
force during minimization.

*/
