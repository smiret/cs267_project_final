#include "CGSolver.H"

double CGSolver::solve(const SparseMatrix& a_A, const vector<double>& a_rhs, double a_tolerance, int a_iter, vector<double>& a_phi)
{
  // const double a_rhs[], int a_rhs_size,
  cout << "Solver" << endl;
  int length = a_rhs.size();
  vector<double> temp(length);
  a_phi = temp;
  double normrhs = 0;
  for(int j = 0; j < length; j++)
    {
      const double element = a_rhs[j];
      if(element > normrhs)
        {
          normrhs = element;
        }
    }
  vector<double> res(length);
  vector<double> product = a_A * a_phi;

  #pragma omp parallel for
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
  vector<double> p(length);
  p = res;
  int iters = 0;
  for(int i = 0; (i < a_iter) && ((maxres/normrhs) > a_tolerance); i++)
    {
      /*double num = 0;
      for(int j = 0; j < length; j++)
      {
          num += res[j] * res[j];
      }
      vector<double> ap(length);
      ap = a_A * p;
      double den = 0;
      for(int j = 0; j < length; j++)
      {
          den += p[j] * ap[j];
      }
      double alpha = num/den;
      for(int j = 0; j < length; j++)
      {
          a_phi[j] += alpha * p[j];
          res[j] += -alpha * ap[j];
      }
      double num2 = 0;
      for(int j = 0; j < length; j++)
      {
          num2 += res[j] * res[j];
      }
      double beta = num2/num;
      for(int j = 0; j < length; j++)
      {
          p[j] = res[j] + (beta * p[j]);
      }
      maxres = 0;
      for(int j = 0; j < length; j++)
        {
          const double element = res[j];
          if(element > maxres)
            {
              maxres = element;
            }   
        }*/
      double num = 0;
      for(int j = 0; j < length; j++)
      {
          num += p[j] * res[j];
      }
      vector<double> ap(length);
      ap = a_A * p;
      double den = 0;
      for(int j = 0; j < length; j++)
      {
          den += p[j] * ap[j];
      }
      double alpha = num/den;
      #pragma omp parallel for
      for(int j = 0; j < length; j++)
      {
          a_phi[j] += alpha * p[j];
      }
      product = a_A * a_phi;
      #pragma omp parallel for
      for(int j = 0; j < length; j++)
      {
          res[j] = a_rhs[j] - product[j];
      }
      double num2 = 0;
      for(int j = 0; j < length; j++)
      {
          num2 -= res[j] * ap[j];
      }
      double beta = num2/den;
    #pragma omp parallel for
    for(int j = 0; j < length; j++)
      {
          p[j] = res[j] + (beta * p[j]);
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