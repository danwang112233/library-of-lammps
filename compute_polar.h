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
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <limits>


using namespace std;
namespace LAMMPS_NS {

class ComputePolar : public Compute {
 public:
  ComputePolar(class LAMMPS *, int, char **);
  ~ComputePolar();
  void init();
  void compute_vector();
  void set_arrays(int);
 // double memory_usage();

  /* ---------------------------------------------------------------------- */
  int sgn(double d){
      if(d<0) return -1;
      else if (d==0) return 0;
      else return 1;}
/* ---------------------------------------------------------------------- */

/*----------------------------------------------------------------------*/
int CountLines(char *filename)  
{  
    ifstream ReadDump;  
    int n=0;  
    string tmp;  
    ReadDump.open(filename, ios::in);
    if(ReadDump.fail())
    {  
        return 0;  
    }  
    else  
    {  
        while(getline(ReadDump,tmp,'\n'))  
        {  
            n++;  
        }  
       
        ReadDump.close();  
        return n;  
    }  
}  
/*----------------------------------------------------------------------*/

 private:
  //bigint natoms = atom->natoms;//change to bigint in the lab
  double *q;
  double **data;
  double polar[4];
  char *id_fix;
  class FixStore *fix;
};

}

#endif
#endif
