#ifdef FIX_CLASS

FixStyle(rotate,FixRotate)

#else

#ifndef LMP_FIX_ROTATE_H_
#define LMP_FIX_ROTATE_H_

#include <stdio.h>
#include "fix.h"
#include "select_bond.h"
#include "mkdssp.h"
#define MAX_ATTEMPTS 1000

namespace LAMMPS_NS {

class FixRotate : public Fix {
 public:
  FixRotate(class LAMMPS *, int, char **);
  ~FixRotate();
  int setmask();
  void init();
  void initial_integrate(int);
  void final_integrate();
  void initial_integrate_respa(int, int, int);
  void final_integrate_respa(int, int);
  int setParameters();

  double memory_usage();
  void write_restart(FILE *);
  void restart(char *);
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  void set_arrays(int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  int pack_restart(int, double *);
  void unpack_restart(int, int);
  int maxsize_restart();
  int size_restart(int);
  void reset_dt();
  void check_energy();
  double energy_full();

 private:
  double **xprevious; //coordinates of last (accepted) step
  int attempts;
  int atom1,atom2;
  double elast,eprev;
  double **x;
  char *xvarstr,*yvarstr,*zvarstr,*vxvarstr,*vyvarstr,*vzvarstr;
  int mstyle;
  int vxflag,vyflag,vzflag,axflag,ayflag,azflag;
  double vx,vy,vz,ax,ay,az;
  double period,omega_rotate;
  double point[3],axis[3],runit[3]; //origin, axis, and normalized axis of rotation
  double dt,dtv,dtf;
  int xvar,yvar,zvar,vxvar,vyvar,vzvar;
  int xvarstyle,yvarstyle,zvarstyle,vxvarstyle,vyvarstyle,vzvarstyle;
  int extra_flag,omega_flag,angmom_flag;
  int radius_flag,ellipsoid_flag,line_flag,tri_flag,body_flag;
  int theta_flag,quat_flag;
  int nlevels_respa,nrestart;
  int time_origin;
  int count;
  int firststruct;
 // double *xatom1,*xatom2;
  double **xoriginal;			// unmapped coords of atoms (NO PBC)
  double xatom1[3],xatom2[3];	// unmapped coords of atoms 1 and 2 (used for calculations)
 // double **xprevious;			// coordinates of previous potential step
  double *toriginal;				// original theta of atoms
  double **qoriginal;			// original quat of atoms
  int displaceflag,velocityflag;
  int maxatom;
  double **displace,**velocity;
  class AtomVecEllipsoid *avec_ellipsoid;
  class AtomVecLine *avec_line;
  class AtomVecTri *avec_tri;
  class AtomVecBody *avec_body;
  class RanMars *random;
  class SelectBond *s ;
  class Compute *compute;
  class Compute *c_pe;
  class mkdssp *dssp;
};

}

#endif
#endif
