#ifndef _SPARSEMATRIX_H_
#define _SPARSEMATRIX_H_

#include <vector>
#include <cassert>
#include <cmath>
#include <array>
#include <iostream>
using namespace std;
class SparseMatrix
{
public:
  /// Set up an M rows and N columns sparse matrix with all values of zero (no non-zero elements)
  SparseMatrix();
  SparseMatrix(int a_M, int a_N);

  /// Matrix Vector multiply.  a_v.size()==N, returns vector of size M
  vector<double> operator*(const vector<double>& a_v) const;

  /// Accessor functions for get and set operations of matrix elements. 
  double& operator[](array<int, 2> a_index);

  /// Accessor function just to get a value
  const double& operator[](array<int, 2> a_index) const;

  /// Zero out all the elements, but leave the sparse structure in place.
  void zero();

  /// Build and return a new SparseMatrix that is the transpose of the input matrix.
  SparseMatrix transpose() const;

  /// In place sparse matrix matrix multiply
  //void operator*(SparseMatrix M) const;

  /// Get first dimension size.
  unsigned int M() const;
  /// Get second dimension size.
  unsigned int N() const;

  /// Return true if symmetric.
  bool symmetric() const;

  /// Print function for SparseMatrix type.
  void print() const;

private:
  unsigned int m_m, m_n;
  double m_zero;
  /* 
  m_data[p][q]) contains the 
  (p,colIndex[p][q]) element of the array.
  */ 
  vector<vector<double> > m_data;
  vector<vector<int> >   m_colIndex;
};

#endif
