#include <fstream>
#include <iostream>
#include <string.h>
#include <sstream>
#include <limits>
#include <math.h>
#include <time.h>
#include <stdlib.h>

using namespace std;

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

int sgn(double d)
{
      if(d<0) return -1;
      else if (d==0) return 0;
      else return 1;
}

//
int ReadLine(char *filename, int line, double nnn, int natoms, const char *output)//, const char *file)  
{  
   
   double data[natoms][5];  
   double polar[5];

    ifstream ReadDump;  
    ReadDump.open(filename,ios::in);  
    
    
    if(!ReadDump)  
    {  
        cout << "No displdump file, please check again." << endl;
        return 1;  
    }  
    
    for(int i = 0; i < line-natoms; ++i)
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
//cout << data[j][k] << " " << endl; 
             }
         }
         
   }
    ReadDump.close();
polar[0] = nnn;
//cout<< polar[0] <<endl;
      for(int k = 2; k < 5; k++)
         {   double sum =0.0;
             for(int j = 0; j < natoms; j++)
             {
              sum += data[j][k];
              polar[k - 1] = sum;
		//cout << sum << " " << endl; 

             }
             //cout << polar[k-1] << " " << endl;
         }  
           polar[4] = sqrt(polar[2] * polar[2] + polar[2] * polar[2] + polar[3] * polar[3]) * sgn(polar[1]) * sgn(polar[2]) * sgn(polar[3]);

    ofstream fout(output, ios::app|ios::out);
    
            fout << "["<< polar[0] << ", " <<  polar[1] << ", " << polar[2] << ", " << polar[3] << ", " << polar[4] << "] " << endl;
          fout.close(); 
   
    return 0;  
}  

int main(int argc,char **argv)  
{  
    
    clock_t t1, t2;
    t1 = clock();
	if(argc != 4) 
		{
		cout << "Usage: ./readpol nz ne_min ne_max delta.\n" << "nz: the number of unit cells in z-direction\n" << "ne_min and ne_max are the min and max number your loop used in infile of lammps\n" << "delta is the number multipled the min and the max number to get the final loop of electric strength." << endl;
		}
	else
			{	
			int nz = atoi(argv[1]);
        		//int ne_min = atoi(argv[2]);
			char *b = argv[2];
			char *d = argv[3];
			char *c;
			c = new char[strlen("polarization_com_")+sizeof(nz)];
   	 		//char *d;
			//d = new char[strlen("polarization_fin_")+sizeof(nz)];
			sprintf(c, "%s%d", "polarization_com_", nz);
			//sprintf(d, "%s%d", "polarization_fin_", nz);
				//for(int i = ne_min; i < ne_max +1; ++i)
				//{
				//double e = double(i) * delta;
    				char *a, *f;
        			a = new char[strlen("displdump")+sizeof(b)];
				f = new char[strlen("displdump")+sizeof(d)];
        			//char *b;
        			//b = new char[strlen("poldump")+sizeof(i)];
    				sprintf(a, "%s%s", "displdump", b);
    				sprintf(f, "%s%s", "displdump", d);
    				char *filename = a; 
    				char *finalname = f;

    				int n = CountLines(filename);
    				int nf = CountLines(finalname);
    				int natoms = nz * nz * nz * 10 ;//10 * 10 * 10 * nz;//10 * nx * ny * nz;    
    				cout<<"number of lines："<< n <<"\n"<<endl;  
   	 			cout<<"number of lines："<< nf <<"\n"<<endl;
   	 			
    

  	 			ReadLine(filename, n, 0.0, natoms, c);
 				ReadLine(finalname, nf, 0.0, natoms, c);
  				t2 = clock();
  				cout<<"Time: "<< double (t2-t1)/CLOCKS_PER_SEC << " s."<<endl;
				delete []a;
				delete []f;
        			//delete []b;
				//}
			delete []c;
			//delete []d;
    			}
	
}  

