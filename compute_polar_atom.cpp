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

/* ----------------------------------------------------------------------
   Contributing author: Sai Jayaraman (University of Notre Dame)
------------------------------------------------------------------------- */
#include <math.h>
#include <string.h>
#include "atom.h"
#include "update.h"
#include "group.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "fix_store.h"
#include "memory.h"
#include "error.h"
#include "input.h"
#include "variable.h"
#include "compute_polar_atom.h"


using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputePolarAtom::ComputePolarAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute polar/atom command");

  peratom_flag = 1;
  size_peratom_cols = 4;//
  create_attribute = 1;
  atom->q_flag = 1;

  // create a new fix STORE style
  // id = compute-ID + COMPUTE_STORE, fix group = compute group

  int n = strlen(id) + strlen("_COMPUTE_STORE") + 1;
  id_fix = new char[n];
  strcpy(id_fix,id);
  strcat(id_fix,"_COMPUTE_STORE");

  char **newarg = new char*[6];
  newarg[0] = id_fix;
  newarg[1] = group->names[igroup];
  newarg[2] = (char *) "STORE";
  newarg[3] = (char *) "peratom";
  newarg[4] = (char *) "1";
  newarg[5] = (char *) "3";
  modify->add_fix(6,newarg);
  fix = (FixStore *) modify->fix[modify->nfix-1];
  delete [] newarg;

  // calculate xu,yu,zu for fix store array
  // skip if reset from restart file

  if (fix->restart_reset) fix->restart_reset = 0;
  else {
    double **xoriginal = fix->astore;

    double **x = atom->x;
    int *mask = atom->mask;
    double *q = atom->q;
    imageint *image = atom->image;
    int nlocal = atom->nlocal;
    int natoms = atom->natoms;//
    double volume = domain->xprd * domain->yprd * domain->zprd;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) domain->unmap(x[i],image[i],xoriginal[i]);
      else xoriginal[i][0] = xoriginal[i][1] = xoriginal[i][2] = 0.0;
  }

  // per-atom polar array

  nmax = 0;
  polaratom = NULL;
}

/* ---------------------------------------------------------------------- */

ComputePolarAtom::~ComputePolarAtom()
{
  // check nfix in case all fixes have already been deleted

  if (modify->nfix) modify->delete_fix(id_fix);

  delete [] id_fix;
  memory->destroy(polaratom);
}

/* ---------------------------------------------------------------------- */

void ComputePolarAtom::init()
{
  // set fix which stores original atom coords

  int ifix = modify->find_fix(id_fix);
  if (ifix < 0) error->all(FLERR,"Could not find compute polar/atom fix ID");
  fix = (FixStore *) modify->fix[ifix];

//bigint ncd = group->count(igroup);
}

/* ---------------------------------------------------------------------- */


void ComputePolarAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow local charge displacement array if necessary

  if (atom->nlocal > nmax) {
    memory->destroy(polaratom);
    nmax = atom->nmax;
    memory->create(polaratom,nmax,4,"polar/atom:polar_atom");
    array_atom = polaratom;
  }
/*memory->create(mx,2*nlocal,1000,"setup:mx");

This will allocate a 2d array and store the ptr for it
in mx, wherever you have defined mx.*/

//v
  // dx,dy,dz = displacement of atom from original position
  // original unwrapped position is stored by fix
  // for triclinic, need to unwrap current atom coord via h matrix

  double **xoriginal = fix->astore;

  double **x = atom->x;
  double *q = atom->q;
  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;
  int natoms = atom->natoms;//

  double *h = domain->h;
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double volume = xprd * yprd * zprd;
  int xbox,ybox,zbox;
  double dx,dy,dz;

  if (domain->triclinic == 0) {
    for (int i = 0; i < nlocal; i++){
      if (mask[i] & groupbit) {
        xbox = (image[i] & IMGMASK) - IMGMAX;
        ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
        zbox = (image[i] >> IMG2BITS) - IMGMAX;
        dx = x[i][0] + xbox*xprd - xoriginal[i][0];
        dy = x[i][1] + ybox*yprd - xoriginal[i][1];
        dz = x[i][2] + zbox*zprd - xoriginal[i][2];
        polaratom[i][0] = 1.60217662 * 10 * dx * q[i] / volume;//(volume * natoms);
        polaratom[i][1] = 1.60217662 * 10 * dy * q[i] / volume;//(volume * natoms);
        polaratom[i][2] = 1.60217662 * 10 * dz * q[i] / volume;//(volume * natoms);
       // polaratom[i][3] = sqrt(dx*dx + dy*dy + dz*dz) * fabs(q[i])*\
                                sgn(polaratom[i][0])*sgn(polaratom[i][1])*sgn(polaratom[i][2]);

      } else polaratom[i][0] = polaratom[i][1] =
	     polaratom[i][2] = polaratom[i][3] = 0.0;
}

  } else {
    for (int i = 0; i < nlocal; i++){
      if (mask[i] & groupbit) {
        xbox = (image[i] & IMGMASK) - IMGMAX;
        ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
        zbox = (image[i] >> IMG2BITS) - IMGMAX;
        dx = x[i][0] + h[0]*xbox + h[5]*ybox + h[4]*zbox - xoriginal[i][0];
        dy = x[i][1] + h[1]*ybox + h[3]*zbox - xoriginal[i][1];
        dz = x[i][2] + h[2]*zbox - xoriginal[i][2];
        polaratom[i][0] = 1.60217662 * 10 * dx * q[i] / volume;//(volume * natoms);
        polaratom[i][1] = 1.60217662 * 10 * dy * q[i] / volume;//(volume * natoms);
        polaratom[i][2] = 1.60217662 * 10 * dz * q[i] / volume;//(volume * natoms);
        //polaratom[i][3] = sqrt(dx*dx + dy*dy + dz*dz) * fabs(q[i])*\
                                sgn(polaratom[i][0])*sgn(polaratom[i][1])*sgn(polaratom[i][2]);

      } else polaratom[i][0] = polaratom[i][1] =
	     polaratom[i][2] = polaratom[i][3] = 0.0;
}
}
        
        //charge_displace[i][3] = sqrt(dx*dx + dy*dy + dz*dz)*fabs(q[i]);

}

/* ----------------------------------------------------------------------
   initialize one atom's storage values, called when atom is created
------------------------------------------------------------------------- */

void ComputePolarAtom::set_arrays(int i)
{
  double **xoriginal = fix->astore;
  double **x = atom->x;
  xoriginal[i][0] = x[i][0];
  xoriginal[i][1] = x[i][1];
  xoriginal[i][2] = x[i][2];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */
//V
double ComputePolarAtom::memory_usage()
{
  double bytes = 4 * sizeof(double);
  return bytes;
} 
