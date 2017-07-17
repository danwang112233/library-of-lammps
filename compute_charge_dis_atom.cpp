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

#include <math.h>
#include <string.h>
#include "atom.h"
#include "update.h"
#include "group.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "fix_store.h"
#include "compute_displace_atom.h"
#include "memory.h"
#include "error.h"
#include "compute_charge_dis_atom.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeChargeDisAtom::ComputeChargeDisAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), id_position(NULL)
{
  if (narg != 3) error->all(FLERR,"Illegal compute chargedis/atom command");


  peratom_flag = 1;
  size_peratom_cols = 4;
  create_attribute = 1;
  atom->q_flag = 1;

  // create a new fix STORE style
  // id = compute-ID + COMPUTE_STORE, fix group = compute group

  int n = strlen(id) + strlen("_COMPUTE_STORE") + 1;
  id_position = new char[n];
  strcpy(id_position,id);
  strcat(id_position,"_COMPUTE_STORE");

  /*char **newarg = new char*[6];
  newarg[0] = id_fix;
  newarg[1] = group->names[igroup];
  newarg[2] = (char *) "STORE";
  newarg[3] = (char *) "peratom";
  newarg[4] = (char *) "1";
  newarg[5] = (char *) "3";
  modify->add_fix(6,newarg);
  fix = (FixStore *) modify->fix[modify->nfix-1];
  delete [] newarg;*/

  // calculate xu,yu,zu for fix store array
  // skip if reset from restart file

  /*if (fix->restart_reset) fix->restart_reset = 0;
  else {*/
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
 /* }*/

  // per-atom charge displacement array

  nmax = 0;
  charge_dis = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeChargeDisAtom::~ComputeChargeDisAtom()
{
  // check nfix in case all fixes have already been deleted

  if (modify->nfix) modify->delete_fix(id_position);

  delete [] id_position;
  memory->destroy(charge_dis);
}

/* ---------------------------------------------------------------------- */

void ComputeChargeDisAtom::init()
{
  // set fix which stores original atom coords

  //int ifix = modify->find_fix(id_fix);
  //if (ifix < 0) error->all(FLERR,"Could not find compute chargedis/atom fix ID");
  //fix = (FixStore *) modify->fix[ifix];

  double icompute = modify->find_compute(id_position);
  if (icompute < 0)
    error->all(FLERR,"compute displace/atom compute does not exist for "
               "compute chargedis/atom");

  ddisplace = (ComputeDisplaceAtom **) modify->compute[icompute];
  //if (strcmp(displace->style,"displace/atom") != 0)//
    //error->all(FLERR,"Compute chargedis/atom does not use displace/atom compute");

//bigint ncd = group->count(igroup);
}

/* ---------------------------------------------------------------------- */

void ComputeChargeDisAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow local charge displacement array if necessary

  if (atom->nlocal > nmax) {
    memory->destroy(charge_dis);
    nmax = atom->nmax;
    memory->create(charge_dis,nmax,4,"chargedis/atom:charge_dis");
    array_atom = charge_dis;
  }
/*memory->create(mx,2*nlocal,1000,"setup:mx");

This will allocate a 2d array and store the ptr for it
in mx, wherever you have defined mx.*/

//v
  // dx,dy,dz = displacement of atom from original position
  // original unwrapped position is stored by fix
  // for triclinic, need to unwrap current atom coord via h matrix

  double **xoriginal = fix->astore;
  double **ddisplace = id+

  double **x = atom->x;
  double *q = atom->q;
  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;
  int natoms = atom->natoms;//

  //double *h = domain->h;
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double volume = xprd * yprd * zprd;
  //int xbox,ybox,zbox;
  //double dx,dy,dz;

  if (domain->triclinic == 0) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
       double **ddisplace =  
        /*xbox = (image[i] & IMGMASK) - IMGMAX;
        ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
        zbox = (image[i] >> IMG2BITS) - IMGMAX;
        dx = x[i][0] + xbox*xprd - xoriginal[i][0];
        dy = x[i][1] + ybox*yprd - xoriginal[i][1];
        dz = x[i][2] + zbox*zprd - xoriginal[i][2];
        charge_dis[i][0] = dx * q[i];
        charge_dis[i][1] = dy * q[i];
        charge_dis[i][2] = dz * q[i];
        charge_dis[i][3] = sqrt(dx*dx + dy*dy + dz*dz)*fabs(q[i]);*/

	charge_dis[i][0] = displace[i][0] * q[i] / volume;//
        charge_dis[i][1] = displace[i][1] * q[i] / volume;//
        charge_dis[i][2] = displace[i][2] * q[i] / volume;//
        charge_dis[i][3] = displace[i][3] * q[i] / volume;//
      } else charge_dis[i][0] = charge_dis[i][1] =
	     charge_dis[i][2] = charge_dis[i][3] = 0.0;

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        /*xbox = (image[i] & IMGMASK) - IMGMAX;
        ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
        zbox = (image[i] >> IMG2BITS) - IMGMAX;
        dx = x[i][0] + xbox*xprd - xoriginal[i][0];
        dy = x[i][1] + ybox*yprd - xoriginal[i][1];
        dz = x[i][2] + zbox*zprd - xoriginal[i][2];
        charge_dis[i][0] = dx * q[i];
        charge_dis[i][1] = dy * q[i];
        charge_dis[i][2] = dz * q[i];
        charge_dis[i][3] = sqrt(dx*dx + dy*dy + dz*dz)*fabs(q[i]);*/

	charge_dis[i][0] = displace[i][0] * q[i] / volume;//
        charge_dis[i][1] = displace[i][1] * q[i] / volume;//
        charge_dis[i][2] = displace[i][2] * q[i] / volume;//
        charge_dis[i][3] = displace[i][3] * q[i] / volume;//
      } else charge_dis[i][0] = charge_dis[i][1] =
	     charge_dis[i][2] = charge_dis[i][3] = 0.0;

  }
        
        //charge_displace[i][3] = sqrt(dx*dx + dy*dy + dz*dz)*fabs(q[i]);

}

/* ----------------------------------------------------------------------
   initialize one atom's storage values, called when atom is created
------------------------------------------------------------------------- */

void ComputeChargeDisAtom::set_arrays(int i)
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
double ComputeChargeDisAtom::memory_usage()
{
  double bytes = nmax*4 * sizeof(double);
  return bytes;
}
