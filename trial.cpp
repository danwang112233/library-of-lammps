#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <limits>
//#include <vector>

#include <time.h>

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

//
/*int ReadLine(char *filename,int line, int natoms)  
{  
   
  // double data[natoms][5];  
    double **data;
    data = new double* [natoms];
    for(int i = 0; i < natoms; i++)
    {
      data[i] = new double[5];
    }

  
    ifstream ReadDump;  
    ReadDump.open(filename,ios::in);  
    
    
    if(ReadDump.fail())  
    {  
        cout << "No Displdump file, please check again." << endl;
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
              
             }
         }
         //for(int j = 0; j < natoms; j++)
         //{    
             for(int k = 0; k < 5; k++)
             {
              cout << data[9][k] << " " << endl;
             }
         //}
    
   }

    ReadDump.close();  

    for(int i; i < natoms; i++)
    {
      delete []data[i];
      
    }
    delete []data;

    return 0;  
}  */

int main()  
{  
    
    clock_t t1, t2;
    t1 = clock();
    char filename[]="displdump";  
    int n = CountLines(filename);
    int natoms = 10;
    cout<<"number of lines："<< n <<"\n"<<endl;  
   
    double **data;
    data = new double* [natoms];
    for(int i = 0; i < natoms; i++)
    {
      data[i] = new double[5];
    }

  
    ifstream ReadDump;  
    ReadDump.open(filename,ios::in);  
    
    
    if(ReadDump.fail())  
    {  
        cout << "No Displdump file, please check again." << endl;
        return 1;  
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
                for(int j = 0; j < natoms; j++)
         {    
             for(int k = 0; k < 5; k++)
             {
              cout << data[j][k] << " " << endl;
             }
         }
         
    
   }
    ReadDump.close();  

    for(int i = 0; i < natoms; i++)
    {
      delete []data[i];
      
    }
    delete []data;

     
    t2 = clock();
    cout<<"Time: "<< double (t2-t1)/CLOCKS_PER_SEC << " s."<<endl;

    return 0; 
}  

