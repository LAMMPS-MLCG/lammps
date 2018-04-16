#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "fix_rotate.h"
#include "atom.h"
#include "group.h"
#include "update.h"
#include "modify.h"
#include "force.h"
#include "domain.h"
#include "lattice.h"
#include "comm.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "atom_vec_ellipsoid.h"
#include "atom_vec_line.h"
#include "atom_vec_tri.h"
#include "atom_vec_body.h"
#include "math_extra.h"
#include "memory.h"
#include "error.h"
#include "select_bond.h"
#include "compute.h"
#include "random_mars.h"
#include "compute_pe.h"
#include "neighbor.h"
#include "atom_vec.h"
#include "atom_vec_hybrid.h"
#include "molecule.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "math_const.h"
#include "math_extra.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum{LINEAR,WIGGLE,ROTATE,VARIABLE};
enum{EQUAL,ATOM};

#define INERTIA 0.2          // moment of inertia prefactor for ellipsoid

/* ---------------------------------------------------------------------- */

FixRotate::FixRotate(LAMMPS *lmp, int narg, char **arg) :
Fix(lmp, narg, arg),
xvarstr(NULL), yvarstr(NULL), zvarstr(NULL), vxvarstr(NULL),
vyvarstr(NULL), vzvarstr(NULL),
xoriginal(NULL), toriginal(NULL), qoriginal(NULL),
displace(NULL), velocity(NULL),random(NULL),s(NULL),xprevious(NULL)
{
	if (narg < 4) error->all(FLERR,"Illegal fix rotate command");
	restart_global = 1;
	restart_peratom = 1;
	peratom_flag = 1;
	size_peratom_cols = 3;
	peratom_freq = 1;
	time_integrate = 1;
	create_attribute = 1;
	displaceflag = 0;
	velocityflag = 0;
	maxatom = 0;
	int iarg;

	random = new RanMars(lmp,2347924 + comm->me);
	if (strcmp(arg[3],"rotate") == 0){
		if (narg < 7) error->all(FLERR,"Illegal fix rotate command");
		iarg = 7;
		mstyle = ROTATE;
		s = new SelectBond(lmp,arg[4]);
	}
	int scaleflag = 1;

	while (iarg < narg) {
		if (strcmp(arg[iarg],"units") == 0) {
			if (iarg+2 > narg) error->all(FLERR,"Illegal fix rotate command");
			if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
			else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
			else error->all(FLERR,"Illegal fix rotate command");
			iarg += 2;
		} else error->all(FLERR,"Illegal fix rotate command");
	}

	// error checks and warnings

	if (domain->dimension == 2) {

		if (mstyle == ROTATE && (axis[0] != 0.0 || axis[1] != 0.0))
			error->all(FLERR,
					"Fix move cannot rotate around non z-axis for 2d problem");

	}


	if ((mstyle == ROTATE) &&
			scaleflag) {
		double xscale = domain->lattice->xlattice;
		double yscale = domain->lattice->ylattice;
		double zscale = domain->lattice->zlattice;

		if (mstyle == ROTATE) {
			point[0] *= xscale;
			point[1] *= yscale;
			point[2] *= zscale;
		}
	}

	// set omega_rotate from period

	if (mstyle == ROTATE) omega_rotate = MY_2PI / period;

	// runit = unit vector along rotation axis

	//	if (mstyle == ROTATE) {
	//		double len = sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
	//		if (len == 0.0)
	//			error->all(FLERR,"Zero length rotation vector with fix move");
	//		runit[0] = axis[0]/len;
	//		runit[1] = axis[1]/len;
	//		runit[2] = axis[2]/len;
	//	}

	// set flags for extra attributes particles may store
	// relevant extra attributes = omega, angmom, theta, quat

	omega_flag = atom->omega_flag;
	angmom_flag = atom->angmom_flag;

	radius_flag = atom->radius_flag;
	ellipsoid_flag = atom->ellipsoid_flag;
	line_flag = atom->line_flag;
	tri_flag = atom->tri_flag;
	body_flag = atom->body_flag;

	theta_flag = quat_flag = 0;
	if (line_flag) theta_flag = 1;
	if (ellipsoid_flag || tri_flag || body_flag) quat_flag = 1;

	extra_flag = 0;
	if (omega_flag || angmom_flag || theta_flag || quat_flag) extra_flag = 1;

	// perform initial allocation of atom-based array
	// register with Atom class

	grow_arrays(atom->nmax);
	//cout<<"atom->nmax "<<atom->nmax;
	atom->add_callback(0);

	displace = velocity = NULL;

	// AtomVec pointers to retrieve per-atom storage of extra quantities

	avec_ellipsoid = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
	avec_line = (AtomVecLine *) atom->style_match("line");
	avec_tri = (AtomVecTri *) atom->style_match("tri");
	avec_body = (AtomVecBody *) atom->style_match("body");

	// xoriginal = initial unwrapped positions of atoms
	// toriginal = initial theta of lines
	// qoriginal = initial quat of extended particles

	double **x = atom->x;
	imageint *image = atom->image;
	int *ellipsoid = atom->ellipsoid;
	int *line = atom->line;
	int *tri = atom->tri;
	int *body = atom->body;
	int *mask = atom->mask;
	int nlocal = atom->nlocal;

	for(int i=0;i<nlocal;i++){
		xprevious[i][0] = x[i][0];
		xprevious[i][1] = x[i][1];
		xprevious[i][2] = x[i][2];
	}

	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit){
			domain->unmap(xprevious[i],image[i],xoriginal[i]);
		}
		else xoriginal[i][0] = xoriginal[i][1] = xoriginal[i][2] = 0.0;
	}

	if (theta_flag) {
		for (int i = 0; i < nlocal; i++) {
			if ((mask[i] & groupbit) && line[i] >= 0)
				toriginal[i] = avec_line->bonus[line[i]].theta;
			else toriginal[i] = 0.0;
		}
	}

	if (quat_flag) {
		double *quat;
		for (int i = 0; i < nlocal; i++) {
			quat = NULL;
			if (mask[i] & groupbit) {
				if (ellipsoid_flag && ellipsoid[i] >= 0)
					quat = avec_ellipsoid->bonus[ellipsoid[i]].quat;
				else if (tri_flag && tri[i] >= 0)
					quat = avec_tri->bonus[tri[i]].quat;
				else if (body_flag && body[i] >= 0)
					quat = avec_body->bonus[body[i]].quat;
			}
			if (quat) {
				qoriginal[i][0] = quat[0];
				qoriginal[i][1] = quat[1];
				qoriginal[i][2] = quat[2];
				qoriginal[i][3] = quat[3];
			} else qoriginal[i][0] = qoriginal[i][1] =
					qoriginal[i][2] = qoriginal[i][3] = 0.0;
		}
	}
	// nrestart = size of per-atom restart data
	// nrestart = 1 + xorig + torig + qorig

	nrestart = 4;
	if (theta_flag) nrestart++;
	if (quat_flag) nrestart += 4;

	// time origin for movement = current timestep

	time_origin = update->ntimestep;

	int icompute = lmp->modify->find_compute("thermo_pe");
	compute = lmp->modify->compute[icompute];
	elast = compute->scalar;
	cout<<"elast "<<elast<<endl;
}

/* ---------------------------------------------------------------------- */

FixRotate::~FixRotate()
{
	// unregister callbacks to this fix from Atom class

	atom->delete_callback(id,0);
	atom->delete_callback(id,1);

	// delete locally stored arrays

	memory->destroy(xoriginal);
	memory->destroy(toriginal);
	memory->destroy(qoriginal);
	memory->destroy(displace);
	memory->destroy(velocity);
	//memory->destroy(xprevious);

	delete [] xvarstr;
	delete [] yvarstr;
	delete [] zvarstr;
	delete [] vxvarstr;
	delete [] vyvarstr;
	delete [] vzvarstr;
}

/* ---------------------------------------------------------------------- */

int FixRotate::setmask()
{
	int mask = 0;
	mask |= INITIAL_INTEGRATE;
	mask |= INITIAL_INTEGRATE_RESPA;
	mask |= FINAL_INTEGRATE;
	mask |= FINAL_INTEGRATE_RESPA;
	return mask;
}

/* ---------------------------------------------------------------------- */

void FixRotate::init()
{
	dt = update->dt;
	dtv = update->dt;
	dtf = 0.5 * update->dt * force->ftm2v;

	// set indices and style of all variables

	displaceflag = velocityflag = 0;

	maxatom = atom->nmax;
	memory->destroy(displace);
	memory->destroy(velocity);
	if (displaceflag) memory->create(displace,maxatom,3,"move:displace");
	else displace = NULL;
	if (velocityflag) memory->create(velocity,maxatom,3,"move:velocity");
	else velocity = NULL;

	//	if (strstr(update->integrate_style,"respa"))
	//		nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ----------------------------------------------------------------------
   set x,v of particles
------------------------------------------------------------------------- */
int FixRotate::setParameters() {
	int *mask = atom->mask;
	int nlocal = atom->nlocal;
	imageint *image = atom->image;
	s->selectBondPairAtoms();
	s->selectDirection();
	s->makeGroup();
	atom1 = s->bondpair[0];
	atom2 = s->bondpair[1];
	return s->sizeGroup;

}

void FixRotate::initial_integrate(int vflag){
	double enext;
	//double delta = (update->ntimestep - time_origin) * dt;
	int flag;
	double ddotr,dx,dy,dz;
	double dtfm,theta_new;
	double xold[3],a[3],b[3],c[3],d[3],disp[3],w[3],ex[3],ey[3],ez[3];	//xold is used only to remap x into the PBC box
	double inertia_ellipsoid[3],qrotate[4];
	double *quat,*inertia,*shape;

	double **x = atom->x;
	//double **v = atom->v;
	//double **f = atom->f;
	imageint *image = atom->image;
	double **omega = atom->omega;
	double **angmom = atom->angmom;
	double *radius = atom->radius;
	double *rmass = atom->rmass;
	double *mass = atom->mass;
	int *type = atom->type;
	int *ellipsoid = atom->ellipsoid;
	int *line = atom->line;
	int *tri = atom->tri;
	int *body = atom->body;
	int *mask = atom->mask;
	int nlocal = atom->nlocal;
	int nmax = atom->nmax;
	memory->grow(xoriginal,nmax,3,"rotate:xoriginal");

	attempts=MAX_ATTEMPTS;
	while(attempts--){

		setParameters();
		cout<<atom1<<'\t'<<atom2<<endl;
		if(atom1 == atom2)
			error->all(FLERR, "Atom1 and Atom2 can not be the same");
		x = atom->x;
		for (int i = 0; i < nlocal; i++) {
			if (mask[i] & groupbit) domain->unmap(xprevious[i],image[i],xoriginal[i]);
			else xoriginal[i][0] = xoriginal[i][1] = xoriginal[i][2] = 0.0;
		}
		domain->unmap(xprevious[atom1],image[atom1],xatom1);
		domain->unmap(xprevious[atom2],image[atom2],xatom2);

		point[0] = (xatom1[0] + xatom2[0])/2;
		point[1] = (xatom1[1] + xatom2[1])/2;
		point[2] = (xatom1[2] + xatom2[2])/2;
		axis[0] = (xatom1[0] - xatom2[0]);
		axis[1] = (xatom1[1] - xatom2[1]);
		axis[2] = (xatom1[2] - xatom2[2]);
		double len = sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
		if (len == 0.0)
			error->all(FLERR,"Zero length rotation vector with fix move");
		runit[0] = axis[0]/len;
		runit[1] = axis[1]/len;
		runit[2] = axis[2]/len;

		double degree = fmod(((random->uniform() * 180.0) - 90 + 360),360); //a value between [-90, 90] (normalized [0,360])
		//double degree = random->uniform() * 20;
		period = 360.0 / degree;
		cout<<"degree   "<<degree<<" period   "<<period<<endl;
		omega_rotate = MY_2PI / period;
		//cout<<xprevious[atom->tag[0]][0]<<'\t'<<xprevious[atom->tag[0]][1]<<'\t'<<xprevious[atom->tag[0]][2]<<endl;
		//cout<<xprevious[atom->tag[2587]][0]<<'\t'<<xprevious[atom->tag[2587]][1]<<'\t'<<xprevious[atom->tag[2587]][2]<<endl;
		//cout<<"atom1 "<<atom1<<' '<<xatom1[0]<<' '<<xatom1[1]<<' '<<xatom1[2]<<endl;

		//cout<<"atom2 "<<atom2<<' '<<xatom2[0]<<' '<<xatom2[1]<<' '<<xatom2[2]<<endl;
		// for rotate by right-hand rule around omega:
		// P = point = vector = point of rotation
		// R = vector = axis of rotation
		// w = omega of rotation (from period)
		// X0 = xoriginal = initial coord of atom
		// R0 = runit = unit vector for R
		// D = X0 - P = vector from P to X0
		// C = (D dot R0) R0 = projection of atom coord onto R line
		// A = D - C = vector from R line to X0
		// B = R0 cross A = vector perp to A in plane of rotation
		// A,B define plane of circular rotation around R line
		// X = P + C + A cos(w*dt) + B sin(w*dt)
		// V = w R0 cross (A cos(w*dt) + B sin(w*dt))

		if (mstyle == ROTATE) {
			double arg = omega_rotate * 1; // simplified legacy code
			double cosine = cos(arg);
			double sine = sin(arg);
			double qcosine = cos(0.5*arg);
			double qsine = sin(0.5*arg);
			qrotate[0] = qcosine;
			qrotate[1] = runit[0]*qsine;
			qrotate[2] = runit[1]*qsine;
			qrotate[3] = runit[2]*qsine;

			for (int i = 0; i < nlocal; i++) {
				if (mask[i] & groupbit) {
					xold[0] = xoriginal[i][0];
					xold[1] = xoriginal[i][1];
					xold[2] = xoriginal[i][2];

					d[0] = xoriginal[i][0] - point[0];
					d[1] = xoriginal[i][1] - point[1];
					d[2] = xoriginal[i][2] - point[2];
					ddotr = d[0]*runit[0] + d[1]*runit[1] + d[2]*runit[2];
					c[0] = ddotr*runit[0];
					c[1] = ddotr*runit[1];
					c[2] = ddotr*runit[2];
					a[0] = d[0] - c[0];
					a[1] = d[1] - c[1];
					a[2] = d[2] - c[2];
					b[0] = runit[1]*a[2] - runit[2]*a[1];
					b[1] = runit[2]*a[0] - runit[0]*a[2];
					b[2] = runit[0]*a[1] - runit[1]*a[0];
					disp[0] = a[0]*cosine  + b[0]*sine;
					disp[1] = a[1]*cosine  + b[1]*sine;
					disp[2] = a[2]*cosine  + b[2]*sine;

					x[i][0] = point[0] + c[0] + disp[0];
					x[i][1] = point[1] + c[1] + disp[1];
					x[i][2] = point[2] + c[2] + disp[2];

					domain->remap_near(x[i],xold);
				}
			}
			// Place holder for extra calculations (in case we need them again)
		}
		enext = energy_full();

		cout<<"attempt # "<<(MAX_ATTEMPTS - attempts)<<endl;
		fprintf(logfile,"attempt # %d\n",(MAX_ATTEMPTS - attempts));
		double delta = enext - elast;
		cout<<"enext: "<<enext<<"\telast: "<<elast<<"\tdelta: "<<delta<<endl;
		if (logfile)
			fprintf(logfile,"enext: %f elast: %f delta: %f\n",enext,elast,delta);

		if(attempts == 0) error->all(FLERR,"SIMULATION QUIT BECAUSE ATTEMPTS REACHED MAX LIMIT");

		if(delta <= 0){
			cout<<"======>\t\t\t structure accepted due to LOW energy"<<endl;
			if (logfile)
				fprintf(logfile,"structure accepted due to low energy\n");
			break;
		} else{
			double factor = exp(-delta/force->boltz/300);
			double randno = random->uniform();
			cout<<"randno \t"<<randno<<" factor \t"<<factor<<endl;
			if (logfile)
				fprintf(logfile,"random no: %f exponential factor: %f\n",randno,factor);
			if (randno <= factor){
				cout<<"======>\t\t\t structure accepted due to RANDOM"<<endl;
				if (logfile)
					fprintf(logfile,"structure accepted due to random\n");
				break;
			}
		}
		//at the end of each rejected attempt, restore atom coordinates.
		for(int i=0;i< nlocal;i++){
			x[i][0] = xprevious[i][0];
			x[i][1] = xprevious[i][1];
			x[i][2] = xprevious[i][2];
		}

		//FIXME test purpose only
		//break;

		//cout<<atom->map(1)<<' '<<xprevious[atom->map(1)][0]<<'\t'<<xprevious[atom->map(1)][1]<<'\t'<<xprevious[atom->map(1)][2]<<endl;
		//cout<<atom->map(2587)<<' '<<xprevious[atom->map(2587)][0]<<'\t'<<xprevious[atom->map(2587)][1]<<'\t'<<xprevious[atom->map(2587)][2]<<endl;

		//TODO use someway to declare that the loop didn't find anything

	} //attempts loop ends here
	elast = enext;
	//keep a copy of the current atom coordinates, as previous
	for(int i=0;i<atom->nlocal;i++) {
		xprevious[i][0] = x[i][0];
		xprevious[i][1] = x[i][1];
		xprevious[i][2] = x[i][2];
	}
}



double FixRotate::energy_full()
{
	int imolecule;

	// if (triclinic) domain->x2lamda(atom->nlocal);
	domain->pbc();
	comm->exchange();
	atom->nghost = 0;
	comm->borders();
	//if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
	if (modify->n_pre_neighbor) modify->pre_neighbor();
	neighbor->build(1);
	int eflag = 1;
	int vflag = 0;

	// if overlap check requested, if overlap,
	// return signal value for energy

	// clear forces so they don't accumulate over multiple
	char *id_pe = (char *) "thermo_pe";
	int ipe = modify->find_compute(id_pe);
	c_pe = modify->compute[ipe];
	size_t nbytes = sizeof(double) * (atom->nlocal + atom->nghost);
	if (nbytes) memset(&atom->f[0][0],0,3*nbytes);

	if (modify->n_pre_force) modify->pre_force(vflag);

	if (force->pair) force->pair->compute(eflag,vflag);

	if (atom->molecular) {
		if (force->bond) force->bond->compute(eflag,vflag);
		if (force->angle) force->angle->compute(eflag,vflag);
		if (force->dihedral) force->dihedral->compute(eflag,vflag);
		if (force->improper) force->improper->compute(eflag,vflag);
	}

	if (force->kspace) force->kspace->compute(eflag,vflag);

	// unlike Verlet, not performing a reverse_comm() or forces here
	// b/c MC does not care about forces
	// don't think it will mess up energy due to any post_force() fixes

	if (modify->n_post_force) modify->post_force(vflag);
	if (modify->n_end_of_step) modify->end_of_step();

	// NOTE: all fixes with THERMO_ENERGY mask set and which
	//   operate at pre_force() or post_force() or end_of_step()
	//   and which user has enable via fix_modify thermo yes,
	//   will contribute to total MC energy via pe->compute_scalar()

	update->eflag_global = update->ntimestep;
	double total_energy = c_pe->compute_scalar();

	return total_energy;
}

/* ----------------------------------------------------------------------
   final NVE of particles with NULL components
------------------------------------------------------------------------- */

void FixRotate::final_integrate()
{
	double dtfm;
	int xflag = 1;
	/*if (mstyle == LINEAR && vxflag) xflag = 0;
	else if (mstyle == WIGGLE && axflag) xflag = 0;
	else if (mstyle == ROTATE) xflag = 0;
	else if (mstyle == VARIABLE && (xvarstr || vxvarstr)) xflag = 0;*/

	int yflag = 1;
	/*if (mstyle == LINEAR && vyflag) yflag = 0;
	else if (mstyle == WIGGLE && ayflag) yflag = 0;
	else if (mstyle == ROTATE) yflag = 0;
	else if (mstyle == VARIABLE && (yvarstr || vyvarstr)) yflag = 0;*/

	int zflag = 1;
	/*if (mstyle == LINEAR && vzflag) zflag = 0;
	else if (mstyle == WIGGLE && azflag) zflag = 0;
	else if (mstyle == ROTATE) zflag = 0;
	else if (mstyle == VARIABLE && (zvarstr || vzvarstr)) zflag = 0; */

	if (mstyle == ROTATE) xflag = 0;
	if (mstyle == ROTATE) yflag = 0;
	if (mstyle == ROTATE) zflag = 0;

	double **v = atom->v;
	double **f = atom->f;
	double *rmass = atom->rmass;
	double *mass = atom->mass;
	int *type = atom->type;
	int *mask = atom->mask;
	int nlocal = atom->nlocal;
	/*	int gid = group->find("rotategroup");
	if (gid > 0) {
		char *cmd[2];
		cmd[0] = "rotategroup";
		cmd[1] = (char *)"clear";
		group->assign(2,cmd);
	}*/
	//	for (int i = 0; i < nlocal; i++) {
	//		if (mask[i] & groupbit) {
	//			if (xflag) {
	//				if (rmass) {
	//					dtfm = dtf / rmass[i];
	//					v[i][0] += dtfm * f[i][0];
	//				} else {
	//					dtfm = dtf / mass[type[i]];
	//					v[i][0] += dtfm * f[i][0];
	//				}
	//			}
	//
	//			if (yflag) {
	//				if (rmass) {
	//					dtfm = dtf / rmass[i];
	//					v[i][1] += dtfm * f[i][1];
	//				} else {
	//					dtfm = dtf / mass[type[i]];
	//					v[i][1] += dtfm * f[i][1];
	//				}
	//			}
	//
	//			if (zflag) {
	//				if (rmass) {
	//					dtfm = dtf / rmass[i];
	//					v[i][2] += dtfm * f[i][2];
	//				} else {
	//					dtfm = dtf / mass[type[i]];
	//					v[i][2] += dtfm * f[i][2];
	//				}
	//			}
	//		}
	//	}
}

/* ---------------------------------------------------------------------- */

void FixRotate::initial_integrate_respa(int vflag, int ilevel, int iloop)
{
	// outermost level - update v and x
	// all other levels - nothing

	if (ilevel == nlevels_respa-1) initial_integrate(vflag);
}

/* ---------------------------------------------------------------------- */

void FixRotate::final_integrate_respa(int ilevel, int iloop)
{
	if (ilevel == nlevels_respa-1) final_integrate();
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixRotate::memory_usage()
{
	//cout<<"flags "<<theta_flag<<"  "<<quat_flag<<"  "<<displaceflag<<"  "<<velocityflag<<endl;
	double bytes = atom->nmax*3 * sizeof(double);
	if (theta_flag) bytes += atom->nmax * sizeof(double);
	if (quat_flag) bytes += atom->nmax*4 * sizeof(double);
	if (displaceflag) bytes += atom->nmax*3 * sizeof(double);
	if (velocityflag) bytes += atom->nmax*3 * sizeof(double);
	return bytes;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixRotate::write_restart(FILE *fp)
{
	int n = 0;
	double list[1];
	list[n++] = time_origin;

	if (comm->me == 0) {
		int size = n * sizeof(double);
		fwrite(&size,sizeof(int),1,fp);
		fwrite(list,sizeof(double),n,fp);
	}
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixRotate::restart(char *buf)
{
	int n = 0;
	double *list = (double *) buf;

	time_origin = static_cast<int> (list[n++]);
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixRotate::grow_arrays(int nmax)
{
	memory->grow(xprevious,atom->nmax,3,"rotate:xprevious");
	memory->grow(xoriginal,nmax,3,"rotate:xoriginal");
	if (theta_flag) memory->grow(toriginal,nmax,"rotate:toriginal");
	if (quat_flag) memory->grow(qoriginal,nmax,4,"rotate:qoriginal");
	array_atom = xprevious;
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixRotate::copy_arrays(int i, int j, int delflag)
{
	xprevious[j][0] = xprevious[i][0];
	xprevious[j][1] = xprevious[i][1];
	xprevious[j][2] = xprevious[i][2];
	if (theta_flag) toriginal[j] = toriginal[i];
	if (quat_flag) {
		qoriginal[j][0] = qoriginal[i][0];
		qoriginal[j][1] = qoriginal[i][1];
		qoriginal[j][2] = qoriginal[i][2];
		qoriginal[j][3] = qoriginal[i][3];
	}
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixRotate::set_arrays(int i)
{
	double theta;
	double *quat;

	double **x = atom->x;
	imageint *image = atom->image;
	int *ellipsoid = atom->ellipsoid;
	int *line = atom->line;
	int *tri = atom->tri;
	int *body = atom->body;
	int *mask = atom->mask;

	// particle not in group

	if (!(mask[i] & groupbit)) {
		xoriginal[i][0] = xoriginal[i][1] = xoriginal[i][2] = 0.0;
		return;
	}

	// current time still equal fix creation time

	if (update->ntimestep == time_origin) {
		domain->unmap(x[i],image[i],xoriginal[i]);
		return;
	}

	// backup particle to time_origin
	/*
	if (mstyle == VARIABLE)
		error->all(FLERR,"Cannot add atoms to fix move variable");
	 */
	domain->unmap(x[i],image[i],xoriginal[i]);
	double delta = (update->ntimestep - time_origin) * update->dt;
	/*
	if (mstyle == LINEAR) {
		if (vxflag) xoriginal[i][0] -= vx * delta;
		if (vyflag) xoriginal[i][1] -= vy * delta;
		if (vzflag) xoriginal[i][2] -= vz * delta;
	} else if (mstyle == WIGGLE) {
		double arg = omega_rotate * delta;
		double sine = sin(arg);
		if (axflag) xoriginal[i][0] -= ax*sine;
		if (ayflag) xoriginal[i][1] -= ay*sine;
		if (azflag) xoriginal[i][2] -= az*sine;
	} else if (mstyle == ROTATE) {*/
	//	if (mstyle == ROTATE) {
	//		double a[3],b[3],c[3],d[3],disp[3],ddotr;
	//		double arg = - omega_rotate * delta;
	//		double sine = sin(arg);
	//		double cosine = cos(arg);
	//		d[0] = x[i][0] - point[0];
	//		d[1] = x[i][1] - point[1];
	//		d[2] = x[i][2] - point[2];
	//		ddotr = d[0]*runit[0] + d[1]*runit[1] + d[2]*runit[2];
	//		c[0] = ddotr*runit[0];
	//		c[1] = ddotr*runit[1];
	//		c[2] = ddotr*runit[2];
	//
	//		a[0] = d[0] - c[0];
	//		a[1] = d[1] - c[1];
	//		a[2] = d[2] - c[2];
	//		b[0] = runit[1]*a[2] - runit[2]*a[1];
	//		b[1] = runit[2]*a[0] - runit[0]*a[2];
	//		b[2] = runit[0]*a[1] - runit[1]*a[0];
	//		disp[0] = a[0]*cosine  + b[0]*sine;
	//		disp[1] = a[1]*cosine  + b[1]*sine;
	//		disp[2] = a[2]*cosine  + b[2]*sine;
	//
	//		xoriginal[i][0] = point[0] + c[0] + disp[0];
	//		xoriginal[i][1] = point[1] + c[1] + disp[1];
	//		xoriginal[i][2] = point[2] + c[2] + disp[2];
	//
	//		// set theta and quat extra attributes affected by rotation
	//
	//		if (extra_flag) {
	//
	//			// theta for lines
	//
	//			if (theta_flag && line[i] >= 0.0) {
	//				theta = avec_line->bonus[atom->line[i]].theta;
	//				toriginal[i] = theta - 0.0;  // NOTE: edit this line
	//			}
	//
	//			// quats for ellipsoids, tris, and bodies
	//
	//			if (quat_flag) {
	//				quat = NULL;
	//				if (ellipsoid_flag && ellipsoid[i] >= 0)
	//					quat = avec_ellipsoid->bonus[ellipsoid[i]].quat;
	//				else if (tri_flag && tri[i] >= 0)
	//					quat = avec_tri->bonus[tri[i]].quat;
	//				else if (body_flag && body[i] >= 0)
	//					quat = avec_body->bonus[body[i]].quat;
	//				if (quat) {
	//					// qoriginal = f(quat,-delta);   // NOTE: edit this line
	//				}
	//			}
	//		}
	//	}
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixRotate::pack_exchange(int i, double *buf)
{
	int n = 0;
	buf[n++] = xprevious[i][0];
	buf[n++] = xprevious[i][1];
	buf[n++] = xprevious[i][2];
	if (theta_flag) buf[n++] = toriginal[i];
	if (quat_flag) {
		buf[n++] = qoriginal[i][0];
		buf[n++] = qoriginal[i][1];
		buf[n++] = qoriginal[i][2];
		buf[n++] = qoriginal[i][3];
	}
	return n;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixRotate::unpack_exchange(int nlocal, double *buf)
{
	int n = 0;
	xprevious[nlocal][0] = buf[n++];
	xprevious[nlocal][1] = buf[n++];
	xprevious[nlocal][2] = buf[n++];
	if (theta_flag) toriginal[nlocal] = buf[n++];
	if (quat_flag) {
		qoriginal[nlocal][0] = buf[n++];
		qoriginal[nlocal][1] = buf[n++];
		qoriginal[nlocal][2] = buf[n++];
		qoriginal[nlocal][3] = buf[n++];
	}
	return n;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixRotate::pack_restart(int i, double *buf)
{
	int n = 1;
	buf[n++] = xprevious[i][0];
	buf[n++] = xprevious[i][1];
	buf[n++] = xprevious[i][2];
	if (theta_flag) buf[n++] = toriginal[i];
	if (quat_flag) {
		buf[n++] = qoriginal[i][0];
		buf[n++] = qoriginal[i][1];
		buf[n++] = qoriginal[i][2];
		buf[n++] = qoriginal[i][3];
	}
	buf[0] = n;
	return n;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixRotate::unpack_restart(int nlocal, int nth)
{
	double **extra = atom->extra;

	// skip to Nth set of extra values

	int m = 0;
	for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
	m++;

	xprevious[nlocal][0] = extra[nlocal][m++];
	xprevious[nlocal][1] = extra[nlocal][m++];
	xprevious[nlocal][2] = extra[nlocal][m++];
	if (theta_flag) toriginal[nlocal] = extra[nlocal][m++];
	if (quat_flag) {
		qoriginal[nlocal][0] = extra[nlocal][m++];
		qoriginal[nlocal][1] = extra[nlocal][m++];
		qoriginal[nlocal][2] = extra[nlocal][m++];
		qoriginal[nlocal][3] = extra[nlocal][m++];
	}
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixRotate::maxsize_restart()
{
	return nrestart;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixRotate::size_restart(int nlocal)
{
	return nrestart;
}

/* ---------------------------------------------------------------------- */

void FixRotate::reset_dt()
{
	error->all(FLERR,"Resetting timestep size is not allowed with fix move");
}
