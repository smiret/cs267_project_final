#include <iostream>
#include <cassert>
#include <cmath>
#include <math.h>
#include <vector>
#include "mesher.H"

using namespace std;

Mesher::Mesher(){};

Mesher::Mesher(const int& a_N)
{
  m_N = a_N;
  m_nNodes = m_N*m_N*m_N; //Number of nodes
  m_nElem = (m_N-1)*(m_N-1)*(m_N-1); //number of elements
  vector<vector<double>> node_table(m_nNodes,vector<double>(DIM)); //Node table initialized
  vector<vector<int>> conn_table(m_nElem,vector<int>(8)); //Connectivity table initialized
  m_nodes = node_table;
  m_conn = conn_table;
};
  

void Mesher::createmesh()
{
    
  int nSet = m_N-1;
  int count = 0;
  for(int i=0; i<nSet; i++)
    {
      for(int k=0; k<nSet; k++)
        {
          for(int j=0; j<nSet; j++)
            {
              m_conn[j+(k*nSet)+(i*nSet*nSet)][0] = count;
              m_conn[j+(k*nSet)+(i*nSet*nSet)][1] = count+1;
              m_conn[j+(k*nSet)+(i*nSet*nSet)][2] = count+nSet+2;
              m_conn[j+(k*nSet)+(i*nSet*nSet)][3] = count+nSet+1;
              m_conn[j+(k*nSet)+(i*nSet*nSet)][4] = count+((nSet+1)*(nSet+1));
              m_conn[j+(k*nSet)+(i*nSet*nSet)][5] = count+((nSet+1)*(nSet+1))+1;
              m_conn[j+(k*nSet)+(i*nSet*nSet)][6] = count+((nSet+1)*(nSet+1))+nSet+2;
              m_conn[j+(k*nSet)+(i*nSet*nSet)][7] = count+((nSet+1)*(nSet+1))+nSet+1;
              count++;
            }
          count++;
        }
      count+=nSet+1;
    }
  //The first row of the connectivity table
  /*
  m_conn[0][0] = 1;
  m_conn[0][1] = 2;
  m_conn[0][2] = m_N + 2;
  m_conn[0][3] = m_N + 1;
  m_conn[0][4] = m_N*m_N + 1;
  m_conn[0][5] = m_N*m_N + 2;
  m_conn[0][6] = m_N*m_N + m_N + 2;
  m_conn[0][7] = m_N*m_N + m_N + 1;


  int counter = 0; //counter for loops
  int jj = 0; //loop variable
  int kk = 0; //loop variable
  
  for (int ii = 0; ii < (m_N - 1); ii++)
    {
      if (ii < (m_N - 1) )
        {
          for (int qq = 0; qq < (m_N - 1); ii++)
            {
              m_conn[counter][qq] = m_conn[counter - (m_N - 1)*(m_N - 1)][qq] + m_N*m_N;
            }
          m_conn[0][0] = 1;
          m_conn[0][1] = 2;
          m_conn[0][2] = m_N + 2;
          m_conn[0][3] = m_N + 1;
          m_conn[0][4] = m_N*m_N + 1;
          m_conn[0][5] = m_N*m_N + 2;
          m_conn[0][6] = m_N*m_N + m_N + 2;
          m_conn[0][7] = m_N*m_N + m_N + 1;
          jj = 0;
        }
      else
        {
          jj = 0;
        }

      for (int jj = 0; jj < (m_N - 1); jj++)
        {
          if (jj < (m_N - 2))
            {
              for (int qq = 0; qq < (m_N - 1); ii++)
                {
                  m_conn[counter][qq] = m_conn[counter - 1][qq] + 1;
                }
              kk = 0;
            }
	      else
	        {
              kk = 0;
            
              for (int kk = 0; kk < (m_N - 2); kk++)
                {
                  for (int qq = 0; qq < (m_N - 1); ii++)
                    {
                      m_conn[counter][qq] = m_conn[counter - 1][qq] + 2;
                    }
                  counter  += 1;

                  for (int jj = 0; jj < (m_N - 2); kk++)
                    {
                      for (int qq = 0; qq < (m_N - 1); ii++)
                        {
                          m_conn[counter][qq] = m_conn[counter - 1][qq] + 1;
                        }
                      counter  += 1;
                    }
                }
	        }
	    }
    }

  //Adjust to zero base
  for (int ww = 0; ww < m_nElem ; ww++)
    {
      for (int aa = 0; aa < 4; aa++)
        {
          m_conn[ww][aa] -= 1;
        }
    }
*/
  //Create the node table

  //Equally spaced vector between 0 and 1

  vector <float> space(m_N);
  float step = 1.0/(m_N - 1);

  for (int qq = 0; qq < m_N; qq++)
    {
      space[qq] = qq*step;
    }

  //Actual Node table

  int cc1 = 0; //index variable

  for (int zz = 0; zz < m_N; zz++)
    {
      m_nodes[cc1][2] = space[zz];
      m_nodes[cc1][1] = 0;
      m_nodes[cc1][0] = 0;
      m_nodes[0][0] = 0;
      m_nodes[0][1] = 0;
      m_nodes[0][2] = 0;

      for (int yy = 0; yy < m_N; yy++)
        {
          m_nodes[cc1][2] = space[zz];
          m_nodes[cc1][1] = space[yy];

          for (int xx = 0; xx < m_N; xx++)
            {
              m_nodes[cc1][2] = space[zz];
              m_nodes[cc1][1] = space[yy];
              m_nodes[cc1][0] = space[xx];
              cc1 += 1;
            }
        }
    }
};

vector<vector<double>> Mesher::getNode()
{
  return(m_nodes);
};

vector<vector<int>> Mesher::getConn()
{
  return(m_conn);
};