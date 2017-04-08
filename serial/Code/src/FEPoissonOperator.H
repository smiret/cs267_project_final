#ifndef _FEPOISSONOPERATOR_H_
#define _FEPOISSONOPERATOR_H_
#include <cstdio> 
#include <cmath> 
#include <cassert>
#include <cmath>
#include <vector>
#include "Node.H"
#include "Element.H"
#include "FEGrid.H"
#include "SparseMatrix.H"
using namespace std;
class FEPoissonOperator
{
 public:
  /// default constructor
  FEPoissonOperator();
  /// constructs a FEPoissonOperator object with given FEGrid
  FEPoissonOperator(const FEGrid& a_grid);
  /// builds the right hand side for Poisson operator
  void makeRHS(vector<double> & a_rhsAtNodes, const vector<double> & a_FCentroids) const;
  /// returns FEGrid in m_grid
  const FEGrid& getFEGrid() const;
  ///returns SparseMatrix in m_matrix
  const SparseMatrix& matrix() const;

 private:
  SparseMatrix m_matrix;
  FEGrid m_grid;
};
#endif
