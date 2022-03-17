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

FixStyle(ElGRAV,FixElGRAV)

#else

#ifndef LMP_FIX_ElGRAV_H
#define LMP_FIX_ElGRAV_H

#include "fix.h"

namespace LAMMPS_NS {

class FixElGRAV : public Fix {
  friend class FixPour;

 public:
  FixElGRAV(class LAMMPS *, int, char **);
  ~FixElGRAV();
  int setmask();
  double int_ellipsoid_V (double, double, double);
  void init();
  void setup(int);
  virtual void post_force(int);
  virtual void post_force_respa(int, int, int);
  double compute_scalar();
  void solve_cubic (double, double, double, double*, double*, double*);

 protected:
  int style;
  double magnitude;
  double vert,phi,theta;
  double xdir,ydir,zdir,wcb,acb,bcb;
  double xgrav,ygrav,zgrav,xacc,yacc,zacc;
  double degree2rad;
  double Dt,now;
  int ilevel_respa;
  int time_origin;
  int eflag;
  double egrav,egrav_all;

  int varflag;
  int mstyle,vstyle,pstyle,tstyle,xstyle,ystyle,zstyle,wstyle,astyle,bstyle;
  int mvar,vvar,pvar,tvar,xvar,yvar,zvar;
  char *mstr,*vstr,*pstr,*tstr,*xstr,*ystr,*zstr,*wstr,*astr,*bstr;

const double grav_const = 6.67408E-11;//Gravitational constant
double sma_ax;//semi-major axis in m (ax>=ay>=az)
double sma_ay;//#######################
double sma_az;//#######################
double a_nd;
double rho_cb;//cb densiy in kg/m3

  void set_acceleration(double, double, double, double, double, double);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Variable name for fix ElGRAV does not exist

Self-explanatory.

E: Variable for fix ElGRAV is invalid style

Only equal-style variables can be used.

*/
