#include <cstdio>   
#include <cmath>    
#include <cstring>  
#include <cassert>  
#include <vector>
#include <fstream>
#include "GridWrite.H"
#include "mesher.H"

int main(int argc, char** argv)
{
    
  int N = 10;
  Mesher mesh(N);
  mesh.createmesh();
  vector<vector<double>> node_table; //Node table initialized
  vector<vector<int>> conn_table; //Connectivity table initialized
  node_table = mesh.getNode();
  conn_table = mesh.getConn();
  int nNodes = node_table.size();
  vector<vector<double>> solution(nNodes, vector<double>(3));
  int nGroup = nNodes/N;
  for(int i=0; i<N; i++)
    {
      for(int j=0; j<nGroup; j++)
        {
          solution[j+(i*nGroup)][0] = 0.;
          solution[j+(i*nGroup)][1] = i*9.93130240774;
          solution[j+(i*nGroup)][2] = 0.;
        }
    }
  
  GridWrite(node_table,conn_table,solution,"testmesh.vtu");
  return 0;
}
