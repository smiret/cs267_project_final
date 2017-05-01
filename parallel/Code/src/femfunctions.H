#ifndef _FEMFUNCTIONS_H_
#define _FEMFUNCTIONS_H_
#include <iostream>
#include <cassert>
#include <cmath>
#include <math.h>
#include <vector>



using namespace std;

class femfunctions
{
public:

  //vector<vector<float> > m_matrix;
  //vector<vector<float> > m_matrix_inv;
  ///Constructor
  femfunctions();
  ///Dense Matrix Multiplication
  void mat_mult(const vector<vector<double> > A_matrix,
                const vector<vector<double> > B_matrix,
                vector<vector<double> >& C_matrix);

  ///Dense Matrix Multiplication with Vector
  void mat_mult_vec(const vector<vector<double> > A_matrix,
                    const vector<double> B_vector,
                    vector<double>& C_vector);
                
 ///Dense Matrix Multiplication with Vector
  void vec_mult_mat(vector<double>& A_vector,
			    vector<vector<double>>& B_matrix,
			    vector<double>& C_vector);
 
  ///Matrix Tranpose
  void mat_transpose(vector<vector<double> >& A_matrix,
		     vector<vector<double> >& A_matrix_trans);

  ///Computes the inverse of a matrix m
  void inverse_mat(vector<vector<double> >& a_matrix,
		   vector<vector<double> >& a_matrix_inv,
		   vector<vector<double> >& a_matrix_inv_trans);

  ///Jacobian of matrix
  double Jacobian(vector<vector<double> >& a_matrix);
  

  ///Mapping Functions
  double phi0(double z1, double z2, double z3);
  ///Mapping Functions
  double phi1(double z1, double z2, double z3);
  ///Mapping Functions
  double phi2(double z1, double z2, double z3);
  ///Mapping Functions
  double phi3(double z1, double z2, double z3);
  ///Mapping Functions
  double phi4(double z1, double z2, double z3);
  ///Mapping Functions
  double phi5(double z1, double z2, double z3);
  ///Mapping Functions
  double phi6(double z1, double z2, double z3);
  ///Mapping Functions
  double phi7(double z1, double z2, double z3);

  ///Derivatives of Mapping Functions with respect to z1
  double dphi0dz1(double z1, double z2, double z3);
  ///Derivatives of Mapping Functions with respect to z1
  double dphi1dz1(double z1, double z2, double z3);
  ///Derivatives of Mapping Functions with respect to z1
  double dphi2dz1(double z1, double z2, double z3);
  ///Derivatives of Mapping Functions with respect to z1
  double dphi3dz1(double z1, double z2, double z3);
  ///Derivatives of Mapping Functions with respect to z1
  double dphi4dz1(double z1, double z2, double z3);
  ///Derivatives of Mapping Functions with respect to z1
  double dphi5dz1(double z1, double z2, double z3);
  ///Derivatives of Mapping Functions with respect to z1
  double dphi6dz1(double z1, double z2, double z3);
  ///Derivatives of Mapping Functions with respect to z1
  double dphi7dz1(double z1, double z2, double z3);

  ///Derivatives of Mapping Functions with respect to z2
  double dphi0dz2(double z1, double z2, double z3);
  ///Derivatives of Mapping Functions with respect to z2
  double dphi1dz2(double z1, double z2, double z3);
  ///Derivatives of Mapping Functions with respect to z2
  double dphi2dz2(double z1, double z2, double z3);
   ///Derivatives of Mapping Functions with respect to z2 
  double dphi3dz2(double z1, double z2, double z3);
  ///Derivatives of Mapping Functions with respect to z2
  double dphi4dz2(double z1, double z2, double z3);
  ///Derivatives of Mapping Functions with respect to z2
  double dphi5dz2(double z1, double z2, double z3);
  ///Derivatives of Mapping Functions with respect to z2
  double dphi6dz2(double z1, double z2, double z3);
  ///Derivatives of Mapping Functions with respect to z2
  double dphi7dz2(double z1, double z2, double z3);

  ///Derivatives of Mapping Functions with respect to z3
  double dphi0dz3(double z1, double z2, double z3);
  ///Derivatives of Mapping Functions with respect to z3
  double dphi1dz3(double z1, double z2, double z3);
  ///Derivatives of Mapping Functions with respect to z3
  double dphi2dz3(double z1, double z2, double z3);
  ///Derivatives of Mapping Functions with respect to z3
  double dphi3dz3(double z1, double z2, double z3);
  ///Derivatives of Mapping Functions with respect to z3
  double dphi4dz3(double z1, double z2, double z3);
  ///Derivatives of Mapping Functions with respect to z3
  double dphi5dz3(double z1, double z2, double z3);
  ///Derivatives of Mapping Functions with respect to z3
  double dphi6dz3(double z1, double z2, double z3);
  ///Derivatives of Mapping Functions with respect to z3
  double dphi7dz3(double z1, double z2, double z3);
 
};

#endif
