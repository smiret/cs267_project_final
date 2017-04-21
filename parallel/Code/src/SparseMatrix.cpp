#include "SparseMatrix.H"
#include <iostream>
using namespace std;
SparseMatrix::SparseMatrix(){};
SparseMatrix::SparseMatrix(int a_M, int a_N)
{
  m_m = a_M;
  m_n = a_N;
  m_zero = 0.;
  vector<vector<double> > data(m_m,vector<double>(0));
  vector<vector<int> > colIndex(m_m,vector<int>(0));
  m_data = data;
  m_colIndex = colIndex;
};
vector<double> SparseMatrix::operator*(const vector<double>& a_v) const
{
  vector<double> prod(m_m);
  const SparseMatrix A = *this;
  for(int i = 0; i < m_m; i++)
    {
      double value = 0.;
      for(int j = 0; j < m_n; j++)
        {
          const double matrixElement = A[{{i,j}}];
          const double vectorElement = a_v[j];
          const double temp = matrixElement*vectorElement;
          value += temp;
        }
      prod[i] = value;
    }
  return prod;
};
double& SparseMatrix::operator[](array<int, 2> a_index)
{
  int l = m_colIndex[a_index[0]].size();
  bool test = 0;
  int index;
  if((a_index[0] > m_m) || (a_index[1] > m_n))
    {
      cout << "Indices out of bounds" <<endl;
      exit(-1);
    }
  for(int i = 0; i < l; i++)
    {
      if(m_colIndex[a_index[0]][i] == a_index[1])
        {
          test = 1;
          index = i;
        }
    }
  if(test)
    {
      return m_data[a_index[0]][index];
    }
  else
    {
      m_data[a_index[0]].push_back(m_zero);
      m_colIndex[a_index[0]].push_back(a_index[1]);
      return m_data[a_index[0]][l];
    }
};
const double& SparseMatrix::operator[](array<int, 2> a_index) const
{
  int l = m_colIndex[a_index[0]].size();
  bool test = 0;
  int index;
  if((a_index[0] > m_m) || (a_index[1] > m_n))
    {
      cout << "Indices out of bounds" <<endl;
      exit(-1);
    }
  for(int i = 0; i < l; i++)
    {
      if(m_colIndex[a_index[0]][i] == a_index[1])
        {
          test = 1;
          index = i;
        }
    }
  if(test)
    {
      return m_data[a_index[0]][index];
    }
  else
    {
      return m_zero;
    }
};
void SparseMatrix::zero()
{
  vector<vector<double> > data(m_m,vector<double>(0));
  vector<vector<int> > colIndex(m_m,vector<int>(0));
  
  m_data = data;
  m_colIndex = colIndex;
};
SparseMatrix SparseMatrix::transpose() const
{
  int m = m_n;
  int n = m_m;
  
  SparseMatrix trans(m,n);
  for(int i = 0; i < m_m; i++)
    {
      for(int j = 0; j < m_colIndex[i].size(); j++)
        {
          //trans[{{m_colIndex[i][j],i}}] = m_data[i][m_colIndex[i][j]];
          trans[{{m_colIndex[i][j],i}}] = m_data[i][j];
        }
    }
  return trans;
};
unsigned int SparseMatrix::M() const
{
  return m_m;
};
unsigned int SparseMatrix::N() const
{
  return m_n;
};
bool SparseMatrix::symmetric() const
{
  bool sym = 1;
  const SparseMatrix trans = transpose();
  sym = sym && (m_m == trans.m_m) && (m_n == trans.m_n);
  const SparseMatrix A = *this;
  if(sym)
    {
      for(int i = 0; i < m_m; i++)
        {
          for(int j = 0; j < m_n; j++)
            {
              const double a = A[{{i,j}}];
              const double b = trans[{{i,j}}];
              sym = sym && (a == b);
            }
        }
    }
  return sym;
};
void SparseMatrix::print() const
{
  cout << "Data:" << endl;
  for (int i = 0; i < m_m; i++) 
    {
      int l = m_data[i].size();
      for (int j = 0; j < l; j++)
        {
          const double data = m_data[i][j];
          cout << data;
          if (j < (m_data[i].size()-1)) 
            {
              cout << ",";
            }
        }
      cout << endl;
    }
  cout << "Indices:" << endl;
  for (int i = 0; i < m_m; i++) 
    {
      int l = m_colIndex[i].size();
      for (int j = 0; j < l; j++)
        {
          cout << m_colIndex[i][j];
          if (j < (m_colIndex[i].size()-1)) 
            {
              cout << ",";
            }
        }
      cout << endl;
    }
  cout << "M: " ;
  cout << m_m;
  cout << endl;
  cout << "N: " ;
  cout << m_n;
  cout << endl;
};