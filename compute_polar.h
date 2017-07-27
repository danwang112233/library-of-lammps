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

ComputeStyle(polar,ComputePolar)

#else

#ifndef LMP_COMPUTE_POLAR_H
#define LMP_COMPUTE_POLAR_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputePolar : public Compute {
 public:
  ComputePolar(class LAMMPS *, int, char **);
  ~ComputePolar();
  void init();
  void compute_peratom();
  void set_arrays(int);
  double memory_usage();
  /* ---------------------------------------------------------------------- */
  int sgn(double d){
      if(d<0) return -1;
      else if (d==0) return 0;
      else return 1;}
/* ---------------------------------------------------------------------- */

 private:
  int nmax;
  //double *q;
  double **polaratom;
  double *polar;
  //double **charge_displacement;
  char *id_fix;
  class ComputeDisplaceAtom **displace;
  class FixStore *fix;
};

}

#endif
#endif
