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

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "fix_ElGRAV.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "respa.h"
#include "modify.h"
#include "input.h"
#include "variable.h"
#include "math_const.h"
#include "error.h"
#include "force.h"
#include "/.../GSL/E_gsl/include/gsl/gsl_poly.h" // ---> /.../ = path to gsl library
#include "/.../GSL/E_gsl/include/gsl/gsl_sf.h"
#include "/.../GSL/E_gsl/include/gsl/gsl_integration.h"
//#include "config.h"
#include "math.h"
#include "/.../GSL/E_gsl/include/gsl/gsl_math.h"


gsl_mode_t mode = GSL_PREC_DOUBLE;

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum{CHUTE,SPHERICAL,VECTOR};
enum{CONSTANT,EQUAL};

//double int_ellipsoid_V (double x, double y, double z);

/* ---------------------------------------------------------------------- */

FixElGRAV::FixElGRAV(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  mstr(NULL), vstr(NULL), pstr(NULL), tstr(NULL), xstr(NULL), ystr(NULL), zstr(NULL)
{
  if (narg < 5) error->all(FLERR,"Illegal fix gravity command");

  dynamic_group_allow = 1;
  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;
  respa_level_support = 1;
  ilevel_respa = 0;

  mstr = vstr = pstr = tstr = xstr = ystr = zstr = NULL;
  mstyle = vstyle = pstyle = tstyle = xstyle = ystyle = zstyle = CONSTANT;

  if (strstr(arg[3],"v_") == arg[3]) {
    int n = strlen(&arg[3][2]) + 1;
    mstr = new char[n];
    strcpy(mstr,&arg[3][2]);
    mstyle = EQUAL;
  } else {
    magnitude = force->numeric(FLERR,arg[3]);
    mstyle = CONSTANT;
  }

  if (strcmp(arg[4],"chute") == 0) {
    if (narg != 6) error->all(FLERR,"Illegal fix gravity command");
    style = CHUTE;
    if (strstr(arg[5],"v_") == arg[5]) {
      int n = strlen(&arg[5][2]) + 1;
      vstr = new char[n];
      strcpy(vstr,&arg[5][2]);
      vstyle = EQUAL;
    } else {
      vert = force->numeric(FLERR,arg[5]);
      vstyle = CONSTANT;
    }

  } else if (strcmp(arg[4],"spherical") == 0) {
    if (narg != 7) error->all(FLERR,"Illegal fix gravity command");
    style = SPHERICAL;
    if (strstr(arg[5],"v_") == arg[5]) {
      int n = strlen(&arg[5][2]) + 1;
      pstr = new char[n];
      strcpy(pstr,&arg[5][2]);
      pstyle = EQUAL;
    } else {
      phi = force->numeric(FLERR,arg[5]);
      pstyle = CONSTANT;
    }
    if (strstr(arg[6],"v_") == arg[6]) {
      int n = strlen(&arg[6][2]) + 1;
      tstr = new char[n];
      strcpy(tstr,&arg[6][2]);
      tstyle = EQUAL;
    } else {
      theta = force->numeric(FLERR,arg[6]);
      tstyle = CONSTANT;
    }

  } else if (strcmp(arg[4],"vector") == 0) {
    if (narg != 11) error->all(FLERR,"Illegal fix gravity command");
    style = VECTOR;
    if (strstr(arg[5],"v_") == arg[5]) {
      int n = strlen(&arg[5][2]) + 1;
      xstr = new char[n];
      strcpy(xstr,&arg[5][2]);
      xstyle = EQUAL;
    } else {
      xdir = force->numeric(FLERR,arg[5]);
      xstyle = CONSTANT;
    }
    if (strstr(arg[6],"v_") == arg[6]) {
      int n = strlen(&arg[6][2]) + 1;
      ystr = new char[n];
      strcpy(ystr,&arg[6][2]);
      ystyle = EQUAL;
    } else {
      ydir = force->numeric(FLERR,arg[6]);
      ystyle = CONSTANT;
    }
    if (strstr(arg[7],"v_") == arg[7]) {
      int n = strlen(&arg[7][2]) + 1;
      zstr = new char[n];
      strcpy(zstr,&arg[7][2]);
      zstyle = EQUAL;
    } else {
      zdir = force->numeric(FLERR,arg[7]);
      zstyle = CONSTANT;
    }
    if (strstr(arg[8],"v_") == arg[8]) {
	      int n = strlen(&arg[8][2]) + 1;
	      wstr = new char[n];
	      strcpy(wstr,&arg[8][2]);
	      wstyle = EQUAL;
	    } else {
	      wcb = force->numeric(FLERR,arg[8]);
	      wstyle = CONSTANT;
	    }
if (strstr(arg[9],"v_") == arg[9]) {
	      int n = strlen(&arg[9][2]) + 1;
	      astr = new char[n];
	      strcpy(astr,&arg[9][2]);
	      astyle = EQUAL;
	    } else {
	      acb = force->numeric(FLERR,arg[9]);
	      astyle = CONSTANT;
	    }
if (strstr(arg[10],"v_") == arg[10]) {
	      int n = strlen(&arg[10][2]) + 1;
	      bstr = new char[n];
	      strcpy(bstr,&arg[10][2]);
	      bstyle = EQUAL;
	    } else {
	      bcb = force->numeric(FLERR,arg[10]);
	      bstyle = CONSTANT;
	    }

  } else error->all(FLERR,"Illegal fix gravity command");

  degree2rad = MY_PI/180.0;
  time_origin = update->ntimestep;
  Dt=update->dt;
  now = time_origin*Dt;

  eflag = 0;
  egrav = 0.0;
}

/* ---------------------------------------------------------------------- */

FixElGRAV::~FixElGRAV()
{
  delete [] mstr;
  delete [] vstr;
  delete [] pstr;
  delete [] tstr;
  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
}

/* ---------------------------------------------------------------------- */

int FixElGRAV::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixElGRAV::init()
{
  if (strstr(update->integrate_style,"respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }

  // check variables

  if (mstr) {
    mvar = input->variable->find(mstr);
    if (mvar < 0)
      error->all(FLERR,"Variable name for fix gravity does not exist");
    if (!input->variable->equalstyle(mvar))
      error->all(FLERR,"Variable for fix gravity is invalid style");
  }
  if (vstr) {
    vvar = input->variable->find(vstr);
    if (vvar < 0)
      error->all(FLERR,"Variable name for fix gravity does not exist");
    if (!input->variable->equalstyle(vvar))
      error->all(FLERR,"Variable for fix gravity is invalid style");
  }
  if (pstr) {
    pvar = input->variable->find(pstr);
    if (pvar < 0)
      error->all(FLERR,"Variable name for fix gravity does not exist");
    if (!input->variable->equalstyle(pvar))
      error->all(FLERR,"Variable for fix gravity is invalid style");
  }
  if (tstr) {
    tvar = input->variable->find(tstr);
    if (tvar < 0)
      error->all(FLERR,"Variable name for fix gravity does not exist");
    if (!input->variable->equalstyle(tvar))
      error->all(FLERR,"Variable for fix gravity is invalid style");
  }
  if (xstr) {
    xvar = input->variable->find(xstr);
    if (xvar < 0)
      error->all(FLERR,"Variable name for fix gravity does not exist");
    if (!input->variable->equalstyle(xvar))
      error->all(FLERR,"Variable for fix gravity is invalid style");
  }
  if (ystr) {
    yvar = input->variable->find(ystr);
    if (yvar < 0)
      error->all(FLERR,"Variable name for fix gravity does not exist");
    if (!input->variable->equalstyle(yvar))
      error->all(FLERR,"Variable for fix gravity is invalid style");
  }
  if (zstr) {
    zvar = input->variable->find(zstr);
    if (zvar < 0)
      error->all(FLERR,"Variable name for fix gravity does not exist");
    if (!input->variable->equalstyle(zvar))
      error->all(FLERR,"Variable for fix gravity is invalid style");
  }

  varflag = CONSTANT;
  if (mstyle != CONSTANT || vstyle != CONSTANT || pstyle != CONSTANT ||
      tstyle != CONSTANT || xstyle != CONSTANT || ystyle != CONSTANT ||
      zstyle != CONSTANT) varflag = EQUAL;

  // set gravity components once and for all

  //if (varflag == CONSTANT) set_acceleration();
}

/* ---------------------------------------------------------------------- */

void FixElGRAV::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixElGRAV::post_force(int vflag)
{
  // update gravity due to variables

  if (varflag != CONSTANT) {
    modify->clearstep_compute();
    if (mstyle == EQUAL) magnitude = input->variable->compute_equal(mvar);
    if (vstyle == EQUAL) vert = input->variable->compute_equal(vvar);
    if (pstyle == EQUAL) phi = input->variable->compute_equal(pvar);
    if (tstyle == EQUAL) theta = input->variable->compute_equal(tvar);
    if (xstyle == EQUAL) xdir = input->variable->compute_equal(xvar);
    if (ystyle == EQUAL) ydir = input->variable->compute_equal(yvar);
    if (zstyle == EQUAL) zdir = input->variable->compute_equal(zvar);
    modify->addstep_compute(update->ntimestep + 1);

    //set_acceleration();
  }

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double massone;
  double dx;

  dx=0.0001;

  eflag = 0;
  egrav = 0.0;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	set_acceleration(x[i][0],x[i][1],x[i][2],v[i][0],v[i][1],v[i][2]);
        massone = rmass[i];
	//central difference for calculating acceleration
        f[i][0] -= massone*((int_ellipsoid_V (x[i][0]+dx,x[i][1],x[i][2])-int_ellipsoid_V (x[i][0]-dx,x[i][1],x[i][2]))/(2*dx))-xacc*massone;
        f[i][1] -= massone*((int_ellipsoid_V (x[i][0],x[i][1]+dx,x[i][2])-int_ellipsoid_V (x[i][0],x[i][1]-dx,x[i][2]))/(2*dx))-yacc*massone;
        f[i][2] -= massone*((int_ellipsoid_V (x[i][0],x[i][1],x[i][2]+dx)-int_ellipsoid_V (x[i][0],x[i][1],x[i][2]-dx))/(2*dx))-zacc*massone;
        egrav -= massone * (int_ellipsoid_V (x[i][0],x[i][1],x[i][2]));
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	set_acceleration(x[i][0],x[i][1],x[i][2],v[i][0],v[i][1],v[i][2]);
        massone = mass[type[i]];
        f[i][0] -= massone*((int_ellipsoid_V (x[i][0]+dx,x[i][1],x[i][2])-int_ellipsoid_V (x[i][0]-dx,x[i][1],x[i][2]))/(2*dx))-xacc*massone;
        f[i][1] -= massone*((int_ellipsoid_V (x[i][0],x[i][1]+dx,x[i][2])-int_ellipsoid_V (x[i][0],x[i][1]-dx,x[i][2]))/(2*dx))-yacc*massone;
        f[i][2] -= massone*((int_ellipsoid_V (x[i][0],x[i][1],x[i][2]+dx)-int_ellipsoid_V (x[i][0],x[i][1],x[i][2]-dx))/(2*dx))-yacc*massone;
        egrav -= massone * (int_ellipsoid_V (x[i][0],x[i][1],x[i][2]));
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixElGRAV::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixElGRAV::set_acceleration(double xp, double yp, double zp, double vxp, double vyp, double vzp)
{

  xacc =  -(wcb*acb*yp - acb*acb*xp - bcb*bcb*xp + wcb*bcb*zp) - 2*(acb*vzp - bcb*vyp);
  yacc =  +(wcb*wcb*yp - wcb*acb*xp - bcb*acb*zp + bcb*bcb*yp) + 2*(wcb*vzp -bcb*vxp);
  zacc =  -(bcb*wcb*xp - wcb*wcb*zp - acb*acb*zp + acb*bcb*yp) - 2*(wcb*vyp - acb*vxp);
  
  
}

/* ----------------------------------------------------------------------
   potential energy in gravity field
------------------------------------------------------------------------- */

double FixElGRAV::compute_scalar()
{
  // only sum across procs one time

  if (eflag == 0) {
    MPI_Allreduce(&egrav,&egrav_all,1,MPI_DOUBLE,MPI_SUM,world);
    eflag = 1;
  }
  return egrav_all;
}

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//To calculate potential due to central body
double FixElGRAV::int_ellipsoid_V (double x, double y, double z) 
{

	//non-dimesionalized parameters: to avoid numerical issues with 
	//operations on large numbers
	double x0, y0, z0;
	double ax0, ay0, az0;

	
	double p, q, r;  //cubic eqn's coeff.
	double	sol_cubic[3];
	double	lambda=0.;
	//~ double result;
	//~ double error;
	//~ size_t limit0=10000;
	int i;
	
	sma_ax=xdir;
	sma_ay=ydir;
	sma_az=zdir;
	rho_cb=magnitude;

	a_nd=sma_ax;

	x0=x/a_nd;
	y0=y/a_nd;
	z0=z/a_nd;
	ax0=sma_ax/a_nd;
	ay0=sma_ay/a_nd;
	az0=sma_az/a_nd;


	//solving cubic
	p = (ax0*ax0 + ay0*ay0 + az0*az0) - (x0*x0 + y0*y0 + z0*z0);
	q = ( (ax0*ax0)*(ay0*ay0) + (ay0*ay0)*(az0*az0) + (az0*az0)*(ax0*ax0) - (x0*x0)*(ay0*ay0 + az0*az0) - (y0*y0)*(ax0*ax0 + az0*az0) - (z0*z0)*(ax0*ax0 + ay0*ay0) );
	r = (ax0*ax0)*(ay0*ay0)*(az0*az0) - (x0*x0)*(ay0*ay0)*(az0*az0) - (ax0*ax0)*(y0*y0)*(az0*az0) - (ax0*ax0)*(ay0*ay0)*(z0*z0);

	for(i=1;i<3;i++) sol_cubic[i]=0.;

	solve_cubic(p,q,r, &sol_cubic[0],&sol_cubic[1],&sol_cubic[2]);

	for(i=0;i<=2;i++) {
		if (lambda < sol_cubic[i]) lambda = sol_cubic[i];
	}

	//solving semi-infinite integral


	double ax02_eff = ax0*ax0 + lambda;
	double ay02_eff = ay0*ay0 + lambda;
	double az02_eff = az0*az0 + lambda;
	double temp = (2./3.)*ax0*ay0*az0;

	double V_ell1 = x0*x0*temp*gsl_sf_ellint_RD(az02_eff, ay02_eff, ax02_eff, mode);
	double V_ell2 = y0*y0*temp*gsl_sf_ellint_RD(az02_eff, ax02_eff, ay02_eff, mode);
	double V_ell3 = z0*z0*temp*gsl_sf_ellint_RD(ax02_eff, ay02_eff, az02_eff, mode);
	double V_ell4 = -3.*temp*gsl_sf_ellint_RF(ax02_eff, ay02_eff, az02_eff, mode);

	double V_ell = M_PI*grav_const*rho_cb*a_nd*a_nd*(V_ell1 + V_ell2 + V_ell3 + V_ell4);

return V_ell;
}



//#define SWAP(a,b) do { double tmp = b ; b = a ; a = tmp ; } while(0)

void FixElGRAV::solve_cubic (double a, double b, double c, double *x0, double *x1, double *x2)
{
  double q = (a * a - 3 * b);
  double r = (2 * a * a * a - 9 * a * b + 27 * c);

  double Q = q / 9;
  double R = r / 54;

  double Q3 = Q * Q * Q;
  double R2 = R * R;

  double CR2 = 729 * r * r;
  double CQ3 = 2916 * q * q * q;

  if (R == 0 && Q == 0)
    {
      *x0 = - a / 3 ;
      *x1 = - a / 3 ;
      *x2 = - a / 3 ;

    }
  else if (CR2 == CQ3) 
    {
      /* this test is actually R2 == Q3, written in a form suitable
         for exact computation with integers */

      /* Due to finite precision some double roots may be missed, and
         considered to be a pair of complex roots z = x +/- epsilon i
         close to the real axis. */

      double sqrtQ = sqrt (Q);

      if (R > 0)
        {
          *x0 = -2 * sqrtQ  - a / 3;
          *x1 = sqrtQ - a / 3;
          *x2 = sqrtQ - a / 3;
        }
      else
        {
          *x0 = - sqrtQ  - a / 3;
          *x1 = - sqrtQ - a / 3;
          *x2 = 2 * sqrtQ - a / 3;
        }

    }
  else if (R2 < Q3)
    {
      double sgnR = (R >= 0 ? 1 : -1);
      double ratio = sgnR * sqrt (R2 / Q3);
      double theta = acos (ratio);
      double norm = -2 * sqrt (Q);
      *x0 = norm * cos (theta / 3) - a / 3;
      *x1 = norm * cos ((theta + 2.0 * M_PI) / 3) - a / 3;
      *x2 = norm * cos ((theta - 2.0 * M_PI) / 3) - a / 3;
      
      /* Sort *x0, *x1, *x2 into increasing order 

      if (*x0 > *x1)
        SWAP(*x0, *x1) ;
      
      if (*x1 > *x2)
        {
          SWAP(*x1, *x2) ;
          
          if (*x0 > *x1)
            SWAP(*x0, *x1) ;
        }*/
      

    }
  else
    {
      double sgnR = (R >= 0 ? 1 : -1);
      double A = -sgnR * pow (fabs (R) + sqrt (R2 - Q3), 1.0/3.0);
      double B = Q / A ;
      *x0 = A + B - a / 3;

    }
}

