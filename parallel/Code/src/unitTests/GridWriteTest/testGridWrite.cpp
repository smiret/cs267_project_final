#include <cstdio>   
#include <cmath>    
#include <cstring>  
#include <cassert>  
#include <vector>
#include <fstream>
#include "GridWrite.H"

int main(int argc, char** argv)
{
    
  int nNodes = 1331;
  int nCells = 1000;
  vector<vector<double>> nodes(nNodes, vector<double>(3));
  vector<vector<int>> conn(nCells, vector<int>(8));
  vector<vector<double>> solution(nNodes, vector<double>(3));
  
  int nSol = 11;
  int nGroup = nNodes/nSol;
  for(int i=0; i<nSol; i++)
    {
      for(int j=0; j<nGroup; j++)
        {
          solution[j+(i*nGroup)][0] = 0.;
          solution[j+(i*nGroup)][1] = i*9.93130240774;
          solution[j+(i*nGroup)][2] = 0.;
        }
    }
  int nSide = 11;
  for(int i=0; i<nSide; i++)
    {
      for(int j=0; j<nSide; j++)
        {
          for(int k=0; k<nSide; k++)
            {
              nodes[k+(j*nSide)+(i*nSide*nSide)][0] = 0.1*k;
              nodes[k+(j*nSide)+(i*nSide*nSide)][1] = 0.1*j;
              nodes[k+(j*nSide)+(i*nSide*nSide)][2] = 0.1*i;
            }
        }
    }
  int nSet = 10;
  int count = 0;
  for(int i=0; i<nSet; i++)
    {
      for(int k=0; k<nSet; k++)
        {
          for(int j=0; j<nSet; j++)
            {
              conn[j+(k*nSet)+(i*nSet*nSet)][0] = count;
              conn[j+(k*nSet)+(i*nSet*nSet)][1] = count+1;
              conn[j+(k*nSet)+(i*nSet*nSet)][2] = count+nSet+2;
              conn[j+(k*nSet)+(i*nSet*nSet)][3] = count+nSet+1;
              conn[j+(k*nSet)+(i*nSet*nSet)][4] = count+((nSet+1)*(nSet+1));
              conn[j+(k*nSet)+(i*nSet*nSet)][5] = count+((nSet+1)*(nSet+1))+1;
              conn[j+(k*nSet)+(i*nSet*nSet)][6] = count+((nSet+1)*(nSet+1))+nSet+2;
              conn[j+(k*nSet)+(i*nSet*nSet)][7] = count+((nSet+1)*(nSet+1))+nSet+1;
              count++;
            }
          count++;
        }
      count+=nSet+1;
    }
  
  GridWrite(nodes,conn,solution,"test.vtu");

  return 0;
}
