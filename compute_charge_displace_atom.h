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

#ifdef COMPUTE_CLASS

ComputeStyle(chargedisplace/atom,ComputeChargeDisplaceAtom)

#else

#ifndef LMP_COMPUTE_CHARGE_DISPLACE_ATOM_H
#define LMP_COMPUTE_CHARGE_DISPLACE_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeChargeDisplaceAtom : public Compute {
 public:
  ComputeChargeDisplaceAtom(class LAMMPS *, int, char **);
  ~ComputeChargeDisplaceAtom();
  void init();
  void compute_peratom();
  void set_arrays(int);
  double memory_usage();

 private:
  int nmax;
  //double *q;
  double **charge_displace;
  char *id_fix;
  class FixStore *fix;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Could not find compute chargedisplace/atom fix ID

Self-explanatory.

*/
