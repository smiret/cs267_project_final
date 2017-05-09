#include "femfunctions.H"
#include <vector>
using namespace std;


//vector<vector<float> > m_matrix;
//vector<vector<float> > m_matrix_inv;
femfunctions::femfunctions(){};
//Dense Matrix Multiplication

void femfunctions::mat_mult(const vector<vector<double> > A_matrix,
    const vector<vector<double> > B_matrix,
    vector<vector<double> >& C_matrix)
{
int row_a = A_matrix.size();
int col_a = A_matrix[0].size();
int row_b = B_matrix.size();
int col_b = B_matrix[0].size();
for (int jj = 0; jj < col_b; jj++)
      {
        for(int ii = 0; ii < row_a; ii++)
          {
            C_matrix[ii][jj] = 0;
          }
      }
for (int ii = 0; ii < row_a; ii++)
  {
    for (int jj = 0; jj < col_b; jj++)
      {
        for (int X = 0; X < col_a; X++)
          {
            C_matrix[ii][jj] += A_matrix[ii][X]*B_matrix[X][jj];
          }
      }
  }
};


//Dense Matrix Multiplication to Yield Vector
void femfunctions::mat_mult_vec(const vector<vector<double> > A_matrix,
			    const vector<double> B_vector,
			    vector<double>& C_vector)

{
int row_a = A_matrix.size();
int col_a = A_matrix[0].size();
int row_b = B_vector.size();
for (int ii = 0; ii < row_a; ii++)
  {
    C_vector[ii] = 0;
  }
for (int ii = 0; ii < row_a; ii++)
  {
    for (int jj = 0; jj < row_b; jj++)
      {
        C_vector[ii] += A_matrix[ii][jj]*B_vector[jj];
      }
  }


};

void femfunctions::vec_mult_mat(vector<double>& A_vector,
			    vector<vector<double>>& B_matrix,
			    vector<double>& C_vector)

{
int row_a = A_vector.size();
int row_b = B_matrix.size();
int col_b = B_matrix[0].size();
for (int ii = 0; ii < col_b; ii++)
  {
    C_vector[ii] = 0;
  }
for (int ii = 0; ii < col_b; ii++)
  {
    for (int jj = 0; jj < row_a; jj++)
      {
        C_vector[ii] += A_vector[jj]*B_matrix[jj][ii];
      }
  }
};

//Matrix Tranpose

void femfunctions::mat_transpose(vector<vector<double> >& A_matrix,
         vector<vector<double> >& A_matrix_trans)
{
int row_a = A_matrix.size();
int col_a = A_matrix[0].size();

for (int ii = 0; ii < row_a; ii++)
  {
    for (int jj = 0; jj < col_a; jj++)
      {
        A_matrix_trans[jj][ii] = A_matrix[ii][jj];
      }
  }

};

    


// computes the inverse of a matrix m

void femfunctions::inverse_mat(vector<vector<double> >& a_matrix,
       vector<vector<double> >& a_matrix_inv,
       vector<vector<double> >& a_matrix_inv_trans)
{
double det = a_matrix[0][0] * (a_matrix[1][1] * a_matrix[2][2] - a_matrix[2][1] * a_matrix[1][2]) -
         a_matrix[0][1] * (a_matrix[1][0] * a_matrix[2][2] - a_matrix[1][2] * a_matrix[2][0]) +
         a_matrix[0][2] * (a_matrix[1][0] * a_matrix[2][1] - a_matrix[1][1] * a_matrix[2][0]);

double invdet = 1.0 / det;

//Matrix33d minv; // inverse of matrix m
a_matrix_inv[0][0] = (a_matrix[1][1] * a_matrix[2][2] - a_matrix[2][1] * a_matrix[1][2]) * invdet;
a_matrix_inv[0][1] = (a_matrix[0][2] * a_matrix[2][1] - a_matrix[0][1] * a_matrix[2][2]) * invdet;
a_matrix_inv[0][2] = (a_matrix[0][1] * a_matrix[1][2] - a_matrix[0][2] * a_matrix[1][1]) * invdet;
a_matrix_inv[1][0] = (a_matrix[1][2] * a_matrix[2][0] - a_matrix[1][0] * a_matrix[2][2]) * invdet;
a_matrix_inv[1][1] = (a_matrix[0][0] * a_matrix[2][2] - a_matrix[0][2] * a_matrix[2][0]) * invdet;
a_matrix_inv[1][2] = (a_matrix[1][0] * a_matrix[0][2] - a_matrix[0][0] * a_matrix[1][2]) * invdet;
a_matrix_inv[2][0] = (a_matrix[1][0] * a_matrix[2][1] - a_matrix[2][0] * a_matrix[1][1]) * invdet;
a_matrix_inv[2][1] = (a_matrix[2][0] * a_matrix[0][1] - a_matrix[0][0] * a_matrix[2][1]) * invdet;
a_matrix_inv[2][2] = (a_matrix[0][0] * a_matrix[1][1] - a_matrix[1][0] * a_matrix[0][1]) * invdet;



a_matrix_inv_trans[0][0] = a_matrix_inv[0][0];
a_matrix_inv_trans[0][1] = a_matrix_inv[1][0];
a_matrix_inv_trans[0][2] = a_matrix_inv[2][0];
a_matrix_inv_trans[1][0] = a_matrix_inv[0][1];
a_matrix_inv_trans[1][1] = a_matrix_inv[1][1];
a_matrix_inv_trans[1][2] = a_matrix_inv[2][1];
a_matrix_inv_trans[2][0] = a_matrix_inv[0][2];
a_matrix_inv_trans[2][1] = a_matrix_inv[1][2];
a_matrix_inv_trans[2][2] = a_matrix_inv[2][2];


};

//Jacobian

double femfunctions::Jacobian(vector<vector<double> >& a_matrix)
{
return ( a_matrix[0][0] * (a_matrix[1][1] * a_matrix[2][2] - a_matrix[2][1] * a_matrix[1][2]) -
         a_matrix[0][1] * (a_matrix[1][0] * a_matrix[2][2] - a_matrix[1][2] * a_matrix[2][0]) +
         a_matrix[0][2] * (a_matrix[1][0] * a_matrix[2][1] - a_matrix[1][1] * a_matrix[2][0]));

};


//Mapping Functions
double femfunctions::phi0(double z1, double z2, double z3)
{
return ((1.0/8.0)*(1.0 - z1)*(1.0 - z2)*(1.0 - z3));
};

double femfunctions::phi1(double z1, double z2, double z3)
{
return ((1.0/8.0)*(1.0 + z1)*(1.0 - z2)*(1.0 - z3));
};

double femfunctions::phi2(double z1, double z2, double z3)
{
return ((1.0/8.0)*(1.0 + z1)*(1.0 + z2)*(1.0 - z3));
};

double femfunctions::phi3(double z1, double z2, double z3)
{
return ((1.0/8.0)*(1.0 - z1)*(1.0 + z2)*(1.0 - z3));
};

double femfunctions::phi4(double z1, double z2, double z3)
{
return ((1.0/8.0)*(1.0 - z1)*(1.0 - z2)*(1.0 + z3));
};

double femfunctions::phi5(double z1, double z2, double z3)
{
return ((1.0/8.0)*(1.0 + z1)*(1.0 - z2)*(1.0 + z3));
};

double femfunctions::phi6(double z1, double z2, double z3)
{
return ((1.0/8.0)*(1.0 + z1)*(1.0 + z2)*(1.0 + z3));
};

double femfunctions::phi7(double z1, double z2, double z3)
{
return ((1.0/8.0)*(1.0 - z1)*(1.0 + z2)*(1.0 + z3));
};

//Derivatives of Mapping Functions with respect to z1

double femfunctions::dphi0dz1(double z1, double z2, double z3)
{
return ((1.0/8.0)*(-1.0)*(1.0 - z2)*(1.0 - z3));
};

double femfunctions::dphi1dz1(double z1, double z2, double z3)
{
return ((1.0/8.0)*(1.0 - z2)*(1.0 - z3));
};

double femfunctions::dphi2dz1(double z1, double z2, double z3)
{
return ((1.0/8.0)*(1.0 + z2)*(1.0 - z3));
};

double femfunctions::dphi3dz1(double z1, double z2, double z3)
{
return ((1.0/8.0)*(-1.0)*(1.0 + z2)*(1.0 - z3));
};

double femfunctions::dphi4dz1(double z1, double z2, double z3)
{
return ((1.0/8.0)*(-1.0)*(1.0 - z2)*(1.0 + z3));
};

double femfunctions::dphi5dz1(double z1, double z2, double z3)
{
return ((1.0/8.0)*(1.0 - z2)*(1.0 + z3));
};

double femfunctions::dphi6dz1(double z1, double z2, double z3)
{
return ((1.0/8.0)*(1.0 + z2)*(1.0 + z3));
};

double femfunctions::dphi7dz1(double z1, double z2, double z3)
{
return ((1.0/8.0)*(-1.0)*(1.0 + z2)*(1.0 + z3));
};


//Derivatives of Mapping Functions with respect to z2

double femfunctions::dphi0dz2(double z1, double z2, double z3)
{
return ((1.0/8.0)*(1.0 - z1)*(-1.0)*(1.0 - z3));
};

double femfunctions::dphi1dz2(double z1, double z2, double z3)
{
return ((1.0/8.0)*(1.0 + z1)*(-1.0)*(1.0 - z3));
};

double femfunctions::dphi2dz2(double z1, double z2, double z3)
{
return ((1.0/8.0)*(1.0 + z1)*(1.0 - z3));
};

double femfunctions::dphi3dz2(double z1, double z2, double z3)
{
return ((1.0/8.0)*(1.0 - z1)*(1.0 - z3));
};

double femfunctions::dphi4dz2(double z1, double z2, double z3)
{
return ((1.0/8.0)*(1.0 - z1)*(-1.0)*(1.0 + z3));
};

double femfunctions::dphi5dz2(double z1, double z2, double z3)
{
return ((1.0/8.0)*(1.0 + z1)*(-1.0)*(1.0 + z3));
};

double femfunctions::dphi6dz2(double z1, double z2, double z3)
{
return ((1.0/8.0)*(1.0 + z1)*(1.0 + z3));
};

double femfunctions::dphi7dz2(double z1, double z2, double z3)
{
return ((1.0/8.0)*(1.0 - z1)*(1.0 + z3));
};



//Derivatives of Mapping Functions with respect to z3

double femfunctions::dphi0dz3(double z1, double z2, double z3)
{
return ((1.0/8.0)*(1.0 - z1)*(1.0 - z2)*(-1.0));
};

double femfunctions::dphi1dz3(double z1, double z2, double z3)
{
return ((1.0/8.0)*(1.0 + z1)*(1.0 - z2)*(-1.0));
};

double femfunctions::dphi2dz3(double z1, double z2, double z3)
{
return ((1.0/8.0)*(1.0 + z1)*(1.0 + z2)*(-1.0));
};

double femfunctions::dphi3dz3(double z1, double z2, double z3)
{
return ((1.0/8.0)*(1.0 - z1)*(1.0 + z2)*(-1.0));
};

double femfunctions::dphi4dz3(double z1, double z2, double z3)
{
return ((1.0/8.0)*(1.0 - z1)*(1.0 - z2));
};

double femfunctions::dphi5dz3(double z1, double z2, double z3)
{
return ((1.0/8.0)*(1.0 + z1)*(1.0 - z2));
};

double femfunctions::dphi6dz3(double z1, double z2, double z3)
{
return ((1.0/8.0)*(1.0 + z1)*(1.0 + z2));
};

double femfunctions::dphi7dz3(double z1, double z2, double z3)
{
return ((1.0/8.0)*(1.0 - z1)*(1.0 + z2));
};