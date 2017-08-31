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
#include "compute_polar.h"


using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputePolar::ComputePolar(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute polar command");

  vector_flag = 1;
  size_vector = 5;
  extvector = 1;
  create_attribute = 1;
  dynamic_group_allow = 0;

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
  else 
{

    double **xoriginal = fix->astore;

    double **x = atom->x;
    double *q = atom->q;
    bigint natoms = atom->natoms;//
    double polaratom[natoms][4];/**/
    double volume = domain->xprd * domain->yprd * domain->zprd;
//read number of lines of Displdump

    char filename[]="displdump";  
    int n = CountLines(filename);

//Displdump data
   data = new double* [natoms];
    for(int i = 0; i < natoms; i++)
    {
      data[i] = new double[5];
    }

  
    std::ifstream ReadDump;  
    ReadDump.open(filename,ios::in);  
    
    
    if(ReadDump.fail())  
    {  
        cout << "No Displdump file, please check again." << endl;
    }  
    
    for(int i = 0; i < n-natoms; ++i)
    {
       ReadDump.ignore (numeric_limits<streamsize>::max(),'\n' );
    }

    while (!ReadDump.eof())
   {  
         for(int j = 0; j < natoms; j++)
         {    
             for(int k = 0; k < 5; k++)
             {
              ReadDump >>  data[j][k];
              
             }
         }
         
    
   }
    ReadDump.close();  
    for(int i = 0; i < natoms; i++)
    {
      delete []data[i];
      
    }
    delete []data;


}
}

/* ---------------------------------------------------------------------- */

ComputePolar::~ComputePolar(){
  // check nfix in case all fixes have already been deleted

  if (modify->nfix) modify->delete_fix(id_fix);

  delete [] id_fix;


}

/* ---------------------------------------------------------------------- */

void ComputePolar::init()
{
  // set fix which stores original atom coords

  int ifix = modify->find_fix(id_fix);
  if (ifix < 0) error->all(FLERR,"Could not find compute Polar fix ID");
  fix = (FixStore *) modify->fix[ifix];


}

/* ---------------------------------------------------------------------- */


void ComputePolar::compute_vector(){
  invoked_vector = update->ntimestep;

  double *q = atom->q;
  bigint natoms = atom->natoms;//
  double polaratom[natoms][4];/**/
  double volume = domain->xprd * domain->yprd * domain->zprd;

//read number of lines of Displdump

    char filename[]="displdump";  
    int n = CountLines(filename);

//Displdump data
   data = new double* [natoms];
    for(int i = 0; i < natoms; i++)
    {
      data[i] = new double[5];
    }

  
    std::ifstream ReadDump;  
    ReadDump.open(filename,ios::in);  
    
    
    if(ReadDump.fail())  
    {  
        cout << "No Displdump file, please check again." << endl;
        
    }  
    
    for(int i = 0; i < n-natoms; ++i)
    {
       ReadDump.ignore (numeric_limits<streamsize>::max(),'\n' );
    }

    while (ReadDump.eof())
   {  
         for(int j = 0; j < natoms; j++)
         {    
             for(int k = 0; k < 5; k++)
             {
              ReadDump >>  data[j][k];
              
             }
         }
         
    
   }
    ReadDump.close();  
        
  if (domain->triclinic == 0) 
{
    for (int i = 0; i < natoms; i++)
   {
      for(int j = 0; j < 5; j++) 
      {

        polaratom[i][0] = 1.60217662 * 10 * data[i][3] * q[i] / (volume * natoms);
        polaratom[i][1] = 1.60217662 * 10 * data[i][4] / (volume * natoms);
        polaratom[i][2] = 1.60217662 * 10 * data[i][5] / (volume * natoms);
        polaratom[i][3] = sqrt(data[i][3]*data[i][3] + data[i][4]*data[i][4] + data[i][5]*data[i][5]) * fabs(q[i])*\
                                sgn(polaratom[i][0])*sgn(polaratom[i][1])*sgn(polaratom[i][2]);

      } 
       polar[0] +=polaratom[i][0];
       polar[1] +=polaratom[i][1];
       polar[2] +=polaratom[i][2];
       polar[3] +=polaratom[i][3];

  } 
}
else 
{
    for (int i = 0; i < natoms; i++)
   {
      for(int j = 0; j < 5; j++) 
      {

        polaratom[i][0] = 1.60217662 * 10 * data[i][3] * q[i] / (volume * natoms);
        polaratom[i][1] = 1.60217662 * 10 * data[i][4] / (volume * natoms);
        polaratom[i][2] = 1.60217662 * 10 * data[i][5] / (volume * natoms);
        polaratom[i][3] = sqrt(data[i][3]*data[i][3] + data[i][4]*data[i][4] + data[i][5]*data[i][5]) * fabs(q[i])*\
                                sgn(polaratom[i][0])*sgn(polaratom[i][1])*sgn(polaratom[i][2]);

      } 
       polar[0] +=polaratom[i][0];
       polar[1] +=polaratom[i][1];
       polar[2] +=polaratom[i][2];
       polar[3] +=polaratom[i][3];

   } 
        

}

}

/* ----------------------------------------------------------------------
   initialize one atom's storage values, called when atom is created
------------------------------------------------------------------------- */

void ComputePolar::set_arrays(int i)
{
  double **xoriginal = fix->astore;
  double **x = atom->x;
  xoriginal[i][0] = x[i][0];
  xoriginal[i][1] = x[i][1];
  xoriginal[i][2] = x[i][2];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array

//V
double ComputePolar::memory_usage()
{
  double bytes = 4 * sizeof(double);
  return bytes;
} 
------------------------------------------------------------------------- */
