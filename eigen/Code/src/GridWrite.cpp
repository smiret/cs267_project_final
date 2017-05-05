#include <cstdio>   
#include <cmath>    
#include <cstring>  
#include <cassert>  
#include <vector>
#include <fstream>
#include "GridWrite.H"

void GridWrite(const vector<vector<double>> &a_nodes,const vector<vector<int>> &a_conn,const vector<vector<double>> &a_solution, const char* a_filename)
{
  int nNodes = a_nodes.size();
  int dim = a_nodes[0].size();
  int nCells = a_conn.size();
  int conDim = a_conn[0].size();
  
  FILE *fp = NULL;
  fp = fopen(a_filename, "w+");
  fprintf(fp,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
  fprintf(fp,"  <UnstructuredGrid>\n");
  fprintf(fp,"    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",nNodes,nCells);
  fprintf(fp,"     <PointData Scalars=\"scalars\">\n");
  fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Displacement\" NumberOfComponents=\"%d\" format=\"ascii\">\n",dim);
  for(int i=0; i<nNodes; i++)
    {
      fprintf(fp,"     ");
      for(int j=0; j<dim;j++)
        {
          fprintf(fp,"%14.10f ",a_solution[i][j]);
        }
      fprintf(fp,"\n");
    }
  fprintf(fp,"      </DataArray>\n");
  fprintf(fp,"    </PointData>\n");
  fprintf(fp,"    <Points>\n");
  fprintf(fp,"      <DataArray type=\"Float32\" NumberOfComponents=\"%d\" format=\"ascii\">\n",dim);
  for(int i=0; i<nNodes; i++)
    {
      fprintf(fp,"     ");
      for(int j=0; j<dim;j++)
        {
          fprintf(fp,"%14.10f ",a_nodes[i][j]);
        }
      fprintf(fp,"\n");
    }
  fprintf(fp,"      </DataArray>\n");
  fprintf(fp,"     </Points>\n");
  fprintf(fp,"    <Cells>\n");
  fprintf(fp,"	  <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
  for(int i=0; i<nCells; i++)
    {
      for(int j=0; j<conDim;j++)
        {
          fprintf(fp,"%d ",a_conn[i][j]);
        }
      fprintf(fp,"\n");
    }
  fprintf(fp,"      </DataArray>\n");
  fprintf(fp,"	      <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
  for(int i=0; i<nCells; i++)
    {
      fprintf(fp,"    ");
      fprintf(fp,"%d",(i+1)*8);
      fprintf(fp,"    ");
    }
  fprintf(fp,"\n      </DataArray>\n");
  fprintf(fp,"	      <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for(int i=0; i<nCells; i++)
    {
      fprintf(fp,"    ");
      fprintf(fp,"12");
      fprintf(fp,"    ");
    }
  fprintf(fp,"\n      </DataArray>\n");
  fprintf(fp,"    </Cells>\n");
  fprintf(fp,"    </Piece>\n");
  fprintf(fp,"  </UnstructuredGrid>\n");
  fprintf(fp,"</VTKFile>");
  
  fclose(fp);
  fp = NULL;
}
