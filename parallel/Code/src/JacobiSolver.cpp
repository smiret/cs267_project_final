#include "JacobiSolver.H"
double JacobiSolver::solve(const SparseMatrix& a_A, const vector<double>& a_rhs, double a_tolerance, int a_iter, vector<double>& a_phi)
{
  cout << "Solver" << endl;
  int length = a_rhs.size();
  vector<double> temp(length);
  a_phi = temp;
  double maxL = 0;
  int size = a_A.M();
  for(int k = 0; k < size; k++)
    {
      const double Lk = a_A[{{k,k}}];
      if(Lk > maxL)
        {
          maxL = Lk;
        }
    }
  double normrhs = 0;
  for(int j = 0; j < length; j++)
    {
      const double element = a_rhs[j];
      if(element > normrhs)
        {
          normrhs = element;
        }
    }
  cout << "norm(rhs) = " << normrhs << endl;
  double alpha = 0.85;
  double lambda = alpha/maxL;
  vector<double> res(length);
  vector<double> product = a_A * a_phi;
  for(int j = 0; j < length; j++)
    {
      res[j] = a_rhs[j] - product[j];
    }
  double maxres = 0;
  for(int j = 0; j < length; j++)
    {
      const double element = res[j];
      if(element > maxres)
        {
          maxres = element;
        }   
    }
  int iters = 0;
  for(int i = 0; (i < a_iter) && ((maxres/normrhs) > a_tolerance); i++)
    {
      product = a_A * a_phi;
      for(int j = 0; j < length; j++)
        {
          a_phi[j] = a_phi[j] + lambda * (a_rhs[j] - product[j]);
        }
      //a_phi = a_phi + lambda * (a_rhs - (a_A * a_phi));
      product = a_A * a_phi;
      for(int j = 0; j < length; j++)
        {
          res[j] = a_rhs[j] - product[j];
        }
      maxres = 0;
      for(int j = 0; j < length; j++)
        {
          const double element = res[j];
          if(element > maxres)
            {
              maxres = element;
            }   
        }
      iters = i+1;
      cout << "Iteration " << iters << endl;
      cout << "Residual " << maxres/normrhs << endl;
    }
  cout << "norm(residual)/norm(rhs) = " << maxres/normrhs << endl;
  cout << "Number of Iterations = " << iters << endl;
  return maxres;
};