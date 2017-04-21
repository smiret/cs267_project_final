#ifndef _CGSOLVER_H_
#define _CGSOLVER_H_

#include "SparseMatrix.H"
#include <vector>

class CGSolver
{
public:
  /** Conjugate gradient solver. Solves until max norm of residual is less than tolerance, 
  or reaches passed in a_iter.  Returns final residual. **/
  double solve(const SparseMatrix& a_A, const vector<double>& a_rhs, double a_tolerance, int a_iter, vector<double>& a_phi);
};

#endif


