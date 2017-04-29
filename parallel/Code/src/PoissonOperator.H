#ifndef _FEPOISSONOPERATOR_H_
#ifndef _FEPOISSONOPERATOR_H_
#define _FEPOISSONOPERATOR_H_
#include <cstdio> 
#include <cmath> 
#include <cassert>
#include <cmath>
#include <vector>
using namespace std;

#include "Node.H"
#include "Element.H"
#include "FEGrid.H"
#include "SparseMatrix.H"

class FEPoissonOperator
{
 public:
  FEPoissonOperator();
  FEPoissonOperator(const FEGrid& a_grid);
  void applyOperator(Vector<float> & a_LOfPhi, const Vector<float> & a_phi);
  const FEGrid& getFEGrid() const;
  ~FEPoissonOperator();

 private:
  SparseMatrix m_matrix;
  FEGrid m_grid;
}
#endif  
 
