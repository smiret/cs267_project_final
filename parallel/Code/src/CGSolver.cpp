#include "CGSolver.H"

double CGSolver::solve(const SparseMatrix& a_A, const vector<double>& a_rhs, double a_tolerance, int a_iter, vector<double>& a_phi)
{
  // const double a_rhs[], int a_rhs_size,
  cout << "Solver" << endl;
  int length = a_rhs.size();
  vector<double> temp(length);
  a_phi = temp;
  double normrhs = 0;
  double maxres = 0;
  vector<double> res(length);
  vector<double> product = a_A * a_phi;

  #pragma omp parallel for reduction(max: normrhs, maxres)
  for(int j = 0; j < length; j++)
  {
    if(a_rhs[j] > normrhs)
      normrhs = a_rhs[j];

    res[j] = a_rhs[j] - product[j];

    if(res[j] > maxres)
      maxres = res[j];
  }

  vector<double> p(length);
  p = res;
  int iters = 0;
  for(int i = 0; (i < a_iter) && ((maxres/normrhs) > a_tolerance); i++)
    {
      double num = 0;
      double den = 0;
      vector<double> ap = a_A * p;

      #pragma parallel for reduction(+: num, den)
      for(int j = 0; j < length; j++)
      {
          num += p[j] * res[j];
          den += p[j] * ap[j];
      }

      double alpha = num/den;

      #pragma parallel for
      for(int j = 0; j < length; j++)
      {
          a_phi[j] += alpha * p[j];
      }
      product = a_A * a_phi;

      double num2 = 0;

      #pragma omp parallel for reduction(-:num2)
      for(int j = 0; j < length; j++)
      {
          res[j] = a_rhs[j] - product[j];
          num2 -= res[j] * ap[j];
      }
      double beta = num2/den;

      maxres = 0;

      #pragma omp parallel for reduction(max: maxres)
      for(int j = 0; j < length; j++)
      {
          p[j] = res[j] + (beta * p[j]);
          if(res[j] > maxres)
            maxres = res[j];
      }


      iters = i+1;
      cout << "Iteration " << iters << endl;
      cout << "Residual " << maxres/normrhs << endl;
    }
  
  cout << "norm(residual)/norm(rhs) = " << maxres/normrhs << endl;
  cout << "Number of Iterations = " << iters << endl;
  return maxres;
};