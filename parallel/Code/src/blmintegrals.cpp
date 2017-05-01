#include <iostream>
#include <cmath>
#include <math.h>
#include <vector>
#include "blmintegrals.H"

using namespace std;

//vector<vector<double> >&
void blmintegrals :: integral1(double (&a_matrix)[24][24],
                               double (&a_x1)[8],
                               double (&a_x2)[8],
                               double (&a_x3)[8],
		 const vector<vector<double> > Etensor)

{
  //Gauss Points for numerical integration

  vector<double> gauss5(5);
  vector<double> w5(5);
  double F11, F12, F13, F21, F22, F23, F31, F32, F33;
  double Jacobian;
  vector <double> dphidz1vec(8), dphidz2vec(8), dphidz3vec(8);
  vector<vector<double> > Fmat(3,vector <double>(3));
  vector < vector <double> > Fmat_inv(3,vector <double> (3));
  vector < vector <double> > Fmat_inv_trans(3,vector <double> (3));
  //vector<vector<double> > Tphi1(6,8), Tphi2(6,8), Tphi3(6,8);
  vector<vector<double> > Tphi(6, vector <double> (24));
  vector<vector<double> > Tphi_trans(24, vector <double> (6));
  vector<vector<double> > inter1(24,vector <double> (6));
  vector<vector<double> > inter2(24,vector <double> (24));
  vector<vector<double> > fS1(24,vector <double>(24));
  femfunctions femfunc;

    //cout << "Hello from " << omp_get_thread_num() << endl;

  gauss5 = {0.000000000000000, 0.538469310105683, -0.538469310105683,0.906179845938664, -0.906179845938664};

  w5 = {0.568888888888889, 0.478628670499366, 0.478628670499366,0.236926885056189, 0.236926885056189};

  for (int qq = 0; qq < 24; qq++)
    {
      for (int ww = 0; ww < 24; ww++)
        {
          fS1[qq][ww] = 0;
        }
    }

//#pragma parallel for
  for (int ii = 0; ii < 5; ii++)
    {
      for (int jj = 0; jj < 5; jj++)
        {
          for (int kk = 0; kk < 5; kk++)
            {
              F11 = (a_x1[0]*femfunc.dphi0dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[1]*femfunc.dphi1dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[2]*femfunc.dphi2dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[3]*femfunc.dphi3dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[4]*femfunc.dphi4dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[5]*femfunc.dphi5dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[6]*femfunc.dphi6dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[7]*femfunc.dphi7dz1(gauss5[ii],gauss5[jj],gauss5[kk]) );

              F12 = (a_x1[0]*femfunc.dphi0dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[1]*femfunc.dphi1dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[2]*femfunc.dphi2dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[3]*femfunc.dphi3dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[4]*femfunc.dphi4dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[5]*femfunc.dphi5dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[6]*femfunc.dphi6dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[7]*femfunc.dphi7dz2(gauss5[ii],gauss5[jj],gauss5[kk]) );

              F13 = (a_x1[0]*femfunc.dphi0dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[1]*femfunc.dphi1dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[2]*femfunc.dphi2dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[3]*femfunc.dphi3dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[4]*femfunc.dphi4dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[5]*femfunc.dphi5dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[6]*femfunc.dphi6dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[7]*femfunc.dphi7dz3(gauss5[ii],gauss5[jj],gauss5[kk]) );

              F21 = (a_x2[0]*femfunc.dphi0dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[1]*femfunc.dphi1dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[2]*femfunc.dphi2dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[3]*femfunc.dphi3dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[4]*femfunc.dphi4dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[5]*femfunc.dphi5dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[6]*femfunc.dphi6dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[7]*femfunc.dphi7dz1(gauss5[ii],gauss5[jj],gauss5[kk]) );

              F22 = (a_x2[0]*femfunc.dphi0dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[1]*femfunc.dphi1dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[2]*femfunc.dphi2dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[3]*femfunc.dphi3dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[4]*femfunc.dphi4dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[5]*femfunc.dphi5dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[6]*femfunc.dphi6dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[7]*femfunc.dphi7dz2(gauss5[ii],gauss5[jj],gauss5[kk]) );

              F23 = (a_x2[0]*femfunc.dphi0dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[1]*femfunc.dphi1dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[2]*femfunc.dphi2dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[3]*femfunc.dphi3dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[4]*femfunc.dphi4dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[5]*femfunc.dphi5dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[6]*femfunc.dphi6dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[7]*femfunc.dphi7dz3(gauss5[ii],gauss5[jj],gauss5[kk]) );

              F31 = (a_x3[0]*femfunc.dphi0dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[1]*femfunc.dphi1dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[2]*femfunc.dphi2dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[3]*femfunc.dphi3dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[4]*femfunc.dphi4dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[5]*femfunc.dphi5dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[6]*femfunc.dphi6dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[7]*femfunc.dphi7dz1(gauss5[ii],gauss5[jj],gauss5[kk]) );

              F32 = (a_x3[0]*femfunc.dphi0dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[1]*femfunc.dphi1dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[2]*femfunc.dphi2dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[3]*femfunc.dphi3dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[4]*femfunc.dphi4dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[5]*femfunc.dphi5dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[6]*femfunc.dphi6dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[7]*femfunc.dphi7dz2(gauss5[ii],gauss5[jj],gauss5[kk]) );

              F33 = (a_x3[0]*femfunc.dphi0dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[1]*femfunc.dphi1dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[2]*femfunc.dphi2dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[3]*femfunc.dphi3dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[4]*femfunc.dphi4dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[5]*femfunc.dphi5dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[6]*femfunc.dphi6dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[7]*femfunc.dphi7dz3(gauss5[ii],gauss5[jj],gauss5[kk]) );

              Fmat[0][0] = F11; Fmat[0][1] = F12; Fmat[0][2] = F13;
              Fmat[1][0] = F21; Fmat[1][1] = F22; Fmat[1][2] = F23;
              Fmat[2][0] = F31; Fmat[2][1] = F32; Fmat[2][2] = F33;

              femfunc.inverse_mat(Fmat, Fmat_inv, Fmat_inv_trans);

              dphidz1vec = {femfunc.dphi0dz1(gauss5[ii],gauss5[jj],gauss5[kk]) , \
              femfunc.dphi1dz1(gauss5[ii],gauss5[jj],gauss5[kk]) , \
              femfunc.dphi2dz1(gauss5[ii],gauss5[jj],gauss5[kk]) , \
              femfunc.dphi3dz1(gauss5[ii],gauss5[jj],gauss5[kk]) , \
              femfunc.dphi4dz1(gauss5[ii],gauss5[jj],gauss5[kk]) , \
              femfunc.dphi5dz1(gauss5[ii],gauss5[jj],gauss5[kk]) , \
              femfunc.dphi6dz1(gauss5[ii],gauss5[jj],gauss5[kk]) , \
              femfunc.dphi7dz1(gauss5[ii],gauss5[jj],gauss5[kk]) };

              dphidz2vec = {femfunc.dphi0dz2(gauss5[ii],gauss5[jj],gauss5[kk]) , \
              femfunc.dphi1dz2(gauss5[ii],gauss5[jj],gauss5[kk]) , \
              femfunc.dphi2dz2(gauss5[ii],gauss5[jj],gauss5[kk]) , \
              femfunc.dphi3dz2(gauss5[ii],gauss5[jj],gauss5[kk]) , \
              femfunc.dphi4dz2(gauss5[ii],gauss5[jj],gauss5[kk]) , \
              femfunc.dphi5dz2(gauss5[ii],gauss5[jj],gauss5[kk]) , \
              femfunc.dphi6dz2(gauss5[ii],gauss5[jj],gauss5[kk]) , \
              femfunc.dphi7dz2(gauss5[ii],gauss5[jj],gauss5[kk]) };

              dphidz3vec = {femfunc.dphi0dz3(gauss5[ii],gauss5[jj],gauss5[kk]) , \
              femfunc.dphi1dz3(gauss5[ii],gauss5[jj],gauss5[kk]) , \
              femfunc.dphi2dz3(gauss5[ii],gauss5[jj],gauss5[kk]) , \
              femfunc.dphi3dz3(gauss5[ii],gauss5[jj],gauss5[kk]) , \
              femfunc.dphi4dz3(gauss5[ii],gauss5[jj],gauss5[kk]) , \
              femfunc.dphi5dz3(gauss5[ii],gauss5[jj],gauss5[kk]) , \
              femfunc.dphi6dz3(gauss5[ii],gauss5[jj],gauss5[kk]) , \
              femfunc.dphi7dz3(gauss5[ii],gauss5[jj],gauss5[kk]) };


              //Tphi1 - normal index
              
              for (int X = 0; X < 8; X++)
                {
                  Tphi[0][X] = dphidz1vec[X]*Fmat_inv[0][0] + dphidz2vec[X]*Fmat_inv[1][0] + dphidz3vec[X]*Fmat_inv[2][0];
                }

              for (int X = 0; X < 8; X++)
                {
                  Tphi[3][X] = dphidz1vec[X]*Fmat_inv[0][1] + dphidz2vec[X]*Fmat_inv[1][1] + dphidz3vec[X]*Fmat_inv[2][1];
                }

              for (int X = 0; X < 8; X++)
                {
                  Tphi[5][X] = dphidz1vec[X]*Fmat_inv[0][2] + dphidz2vec[X]*Fmat_inv[1][2] + dphidz3vec[X]*Fmat_inv[2][2];
                }


              //Tphi2 - normal index + 8
              
              for (int X = 0; X < 8; X++)
                {
                  Tphi[1][X+8] = dphidz1vec[X]*Fmat_inv[0][1] + dphidz2vec[X]*Fmat_inv[1][1] + dphidz3vec[X]*Fmat_inv[2][1];
                }

              for (int X = 0; X < 8; X++)
                {
                  Tphi[3][X+8] = dphidz1vec[X]*Fmat_inv[0][0] + dphidz2vec[X]*Fmat_inv[1][0] + dphidz3vec[X]*Fmat_inv[2][0];
                }

              for (int X = 0; X < 8; X++)
                {
                  Tphi[4][X+8] = dphidz1vec[X]*Fmat_inv[0][2] + dphidz2vec[X]*Fmat_inv[1][2] + dphidz3vec[X]*Fmat_inv[2][2];
                }

              //Tphi3 - normal index + 16
              
              for (int X = 0; X < 8; X++)
                {
                  Tphi[2][X+16] = dphidz1vec[X]*Fmat_inv[0][2] + dphidz2vec[X]*Fmat_inv[1][2] + dphidz3vec[X]*Fmat_inv[2][2];
                }

              for (int X = 0; X < 8; X++)
                {
                  Tphi[4][X+16] = dphidz1vec[X]*Fmat_inv[0][1] + dphidz2vec[X]*Fmat_inv[1][1] + dphidz3vec[X]*Fmat_inv[2][1];
                }

              for (int X = 0; X < 8; X++)
                {
                  Tphi[5][X+16] = dphidz1vec[X]*Fmat_inv[0][0] + dphidz2vec[X]*Fmat_inv[1][0] + dphidz3vec[X]*Fmat_inv[2][0];
                }


              femfunc.mat_transpose(Tphi, Tphi_trans);

              femfunc.mat_mult(Tphi_trans, Etensor, inter1);

              femfunc.mat_mult(inter1, Tphi, inter2);

              Jacobian = femfunc.Jacobian(Fmat);

              for (int qq = 0; qq < 24; qq++)
                {
                  for (int ww = 0; ww < 24; ww++)
                    {
                      fS1[qq][ww] += w5[ii]*w5[jj]*w5[kk]*inter2[qq][ww]*Jacobian;
                    }
                }
            }
        }
    }
//#pragma omp parallel for
  for (int qq = 0; qq < 24; qq++)
    {
      for (int ww = 0; ww < 24; ww++)
        {
          //#pragma omp atomic
          a_matrix[qq][ww] = fS1[qq][ww];
        }
    }
};


void blmintegrals :: integral2(double (&a_vector)[24],
           double (&a_x1)[8],
           double (&a_x2)[8],
           double (&a_x3)[8],
           vector<double>& a_force)

{
//Gauss Points for numerical integration

  vector<double> gauss5(5);
  vector<double> w5(5);
  double F11, F12, F13, F21, F22, F23, F31, F32, F33;
  double Jacobian;
  vector <double> phivec(8);
  vector<vector<double> > phimat(3,vector <double>(24));
  vector<vector<double> > phimat_trans(24,vector <double> (3));
  vector<vector<double> > Fmat(3,vector <double>(3));
  vector < vector <double> > Fmat_inv(3,vector <double> (3));
  vector < vector <double> > Fmat_inv_trans(3,vector <double>(3));
  vector<vector<double> > Tphi(6,vector <double>(24));
  vector<vector<double> > Tphi_trans(24,vector <double>(6));
  vector<double> inter1(24);
  //vector<vector<double> > inter2(24,1);
  vector<double> fR1(24);
  femfunctions femfunc;
  

  gauss5 = {0.000000000000000, 0.538469310105683, -0.538469310105683,0.906179845938664, -0.906179845938664};

  w5 = {0.568888888888889, 0.478628670499366, 0.478628670499366,0.236926885056189, 0.236926885056189};

  for (int qq = 0; qq < 24; qq++)
    {
      fR1[qq] = 0;
    }
    
  for (int ii = 0; ii < 5; ii++)
    {
      for (int jj = 0; jj < 5; jj++)
        {
          for (int kk = 0; kk < 5; kk++)
            {
              F11 = (a_x1[0]*femfunc.dphi0dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[1]*femfunc.dphi1dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[2]*femfunc.dphi2dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[3]*femfunc.dphi3dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[4]*femfunc.dphi4dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[5]*femfunc.dphi5dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[6]*femfunc.dphi6dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[7]*femfunc.dphi7dz1(gauss5[ii],gauss5[jj],gauss5[kk]) );

              F12 = (a_x1[0]*femfunc.dphi0dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[1]*femfunc.dphi1dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[2]*femfunc.dphi2dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[3]*femfunc.dphi3dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[4]*femfunc.dphi4dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[5]*femfunc.dphi5dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[6]*femfunc.dphi6dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[7]*femfunc.dphi7dz2(gauss5[ii],gauss5[jj],gauss5[kk]) );

              F13 = (a_x1[0]*femfunc.dphi0dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[1]*femfunc.dphi1dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[2]*femfunc.dphi2dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[3]*femfunc.dphi3dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[4]*femfunc.dphi4dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[5]*femfunc.dphi5dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[6]*femfunc.dphi6dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x1[7]*femfunc.dphi7dz3(gauss5[ii],gauss5[jj],gauss5[kk]) );

              F21 = (a_x2[0]*femfunc.dphi0dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[1]*femfunc.dphi1dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[2]*femfunc.dphi2dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[3]*femfunc.dphi3dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[4]*femfunc.dphi4dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[5]*femfunc.dphi5dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[6]*femfunc.dphi6dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[7]*femfunc.dphi7dz1(gauss5[ii],gauss5[jj],gauss5[kk]) );

              F22 = (a_x2[0]*femfunc.dphi0dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[1]*femfunc.dphi1dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[2]*femfunc.dphi2dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[3]*femfunc.dphi3dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[4]*femfunc.dphi4dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[5]*femfunc.dphi5dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[6]*femfunc.dphi6dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[7]*femfunc.dphi7dz2(gauss5[ii],gauss5[jj],gauss5[kk]) );

              F23 = (a_x2[0]*femfunc.dphi0dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[1]*femfunc.dphi1dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[2]*femfunc.dphi2dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[3]*femfunc.dphi3dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[4]*femfunc.dphi4dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[5]*femfunc.dphi5dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[6]*femfunc.dphi6dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x2[7]*femfunc.dphi7dz3(gauss5[ii],gauss5[jj],gauss5[kk]) );

              F31 = (a_x3[0]*femfunc.dphi0dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[1]*femfunc.dphi1dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[2]*femfunc.dphi2dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[3]*femfunc.dphi3dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[4]*femfunc.dphi4dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[5]*femfunc.dphi5dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[6]*femfunc.dphi6dz1(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[7]*femfunc.dphi7dz1(gauss5[ii],gauss5[jj],gauss5[kk]) );

              F32 = (a_x3[0]*femfunc.dphi0dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[1]*femfunc.dphi1dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[2]*femfunc.dphi2dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[3]*femfunc.dphi3dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[4]*femfunc.dphi4dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[5]*femfunc.dphi5dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[6]*femfunc.dphi6dz2(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[7]*femfunc.dphi7dz2(gauss5[ii],gauss5[jj],gauss5[kk]) );

              F33 = (a_x3[0]*femfunc.dphi0dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[1]*femfunc.dphi1dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[2]*femfunc.dphi2dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[3]*femfunc.dphi3dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[4]*femfunc.dphi4dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[5]*femfunc.dphi5dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[6]*femfunc.dphi6dz3(gauss5[ii],gauss5[jj],gauss5[kk]) +
                 a_x3[7]*femfunc.dphi7dz3(gauss5[ii],gauss5[jj],gauss5[kk]) );

              Fmat[0][0] = F11; Fmat[0][1] = F12; Fmat[0][2] = F13;
              Fmat[1][0] = F21; Fmat[1][1] = F22; Fmat[1][2] = F23;
              Fmat[2][0] = F31; Fmat[2][1] = F32; Fmat[2][2] = F33;

              femfunc.inverse_mat(Fmat, Fmat_inv, Fmat_inv_trans);

              phivec = {femfunc.phi0(gauss5[ii],gauss5[jj],gauss5[kk]) , \
              femfunc.phi1(gauss5[ii],gauss5[jj],gauss5[kk]) , \
              femfunc.phi2(gauss5[ii],gauss5[jj],gauss5[kk]) , \
              femfunc.phi3(gauss5[ii],gauss5[jj],gauss5[kk]) , \
              femfunc.phi4(gauss5[ii],gauss5[jj],gauss5[kk]) , \
              femfunc.phi5(gauss5[ii],gauss5[jj],gauss5[kk]) , \
              femfunc.phi6(gauss5[ii],gauss5[jj],gauss5[kk]) , \
              femfunc.phi7(gauss5[ii],gauss5[jj],gauss5[kk])};

              for (int X = 0; X < 8; X++)
            {
              phimat[0][X] = phivec[X];
              phimat[1][X+8] = phivec[X];
              phimat[2][X+16] = phivec[X];

            }

              femfunc.mat_transpose(phimat, phimat_trans);

              femfunc.mat_mult_vec(phimat_trans, a_force, inter1);

              Jacobian = femfunc.Jacobian(Fmat);

              for (int qq = 0; qq < 24; qq++)
            {
              
                  fR1[qq] += w5[ii]*w5[jj]*w5[kk]*inter1[qq]*Jacobian;
            }

            }
        }
    }
    for (int qq = 0; qq < 24; qq++)
      {
        a_vector[qq] = fR1[qq];
      }
};
		
	      
void blmintegrals :: integral3xx(double (&a_vector)[24],
                                 double (&a_x1)[8],
                                 double (&a_x2)[8],
                                 double (&a_x3)[8],
		 vector<double>& a_zvector,
		 vector<double>& a_traction)

{

  //Gauss Points for numerical integration

  vector<double> gauss5(5);
  vector<double> w5(5);
  double F11, F12, F13, F21, F22, F23, F31, F32, F33;
  double Jacobian, znormal;
  vector <double> phivec(8);
  vector<vector<double> > phimat(3,vector<double> (24));
  vector<vector<double> > phimat_trans(24,vector <double> (3));
  vector<vector<double> > Fmat(3,vector <double> (3));
  vector<vector<double> > Fmat_inv(3,vector <double> (3));
  vector<vector<double> > Fmat_inv_trans(3,vector <double> (3));
  vector<vector<double> > Tphi(6,vector <double> (24));
  vector<vector<double> > Tphi_trans(24,vector <double> (6));
  vector<double> inter21(3);
  vector<double> inter22(3);
  vector<double> inter3(24);
  double inter2;
  vector<double> fR1(24);
  femfunctions femfunc;

  gauss5 = {0.000000000000000, 0.538469310105683, -0.538469310105683,0.906179845938664, -0.906179845938664};

  w5 = {0.568888888888889, 0.478628670499366, 0.478628670499366,0.236926885056189, 0.236926885056189};

  znormal = a_zvector[0];

  for (int qq = 0; qq < 24; qq++)
    {
      fR1[qq] = 0;
    }
  
  for (int jj = 0; jj < 5; jj++)
    {
      for (int kk = 0; kk < 5; kk++)
        {
          F11 = (a_x1[0]*femfunc.dphi0dz1(znormal,gauss5[jj],gauss5[kk]) +
             a_x1[1]*femfunc.dphi1dz1(znormal,gauss5[jj],gauss5[kk]) +
             a_x1[2]*femfunc.dphi2dz1(znormal,gauss5[jj],gauss5[kk]) +
             a_x1[3]*femfunc.dphi3dz1(znormal,gauss5[jj],gauss5[kk]) +
             a_x1[4]*femfunc.dphi4dz1(znormal,gauss5[jj],gauss5[kk]) +
             a_x1[5]*femfunc.dphi5dz1(znormal,gauss5[jj],gauss5[kk]) +
             a_x1[6]*femfunc.dphi6dz1(znormal,gauss5[jj],gauss5[kk]) +
             a_x1[7]*femfunc.dphi7dz1(znormal,gauss5[jj],gauss5[kk]) );
          
          F12 = (a_x1[0]*femfunc.dphi0dz2(znormal,gauss5[jj],gauss5[kk]) +
             a_x1[1]*femfunc.dphi1dz2(znormal,gauss5[jj],gauss5[kk]) +
             a_x1[2]*femfunc.dphi2dz2(znormal,gauss5[jj],gauss5[kk]) +
             a_x1[3]*femfunc.dphi3dz2(znormal,gauss5[jj],gauss5[kk]) +
             a_x1[4]*femfunc.dphi4dz2(znormal,gauss5[jj],gauss5[kk]) +
             a_x1[5]*femfunc.dphi5dz2(znormal,gauss5[jj],gauss5[kk]) +
             a_x1[6]*femfunc.dphi6dz2(znormal,gauss5[jj],gauss5[kk]) +
             a_x1[7]*femfunc.dphi7dz2(znormal,gauss5[jj],gauss5[kk]) );

          F13 = (a_x1[0]*femfunc.dphi0dz3(znormal,gauss5[jj],gauss5[kk]) +
             a_x1[1]*femfunc.dphi1dz3(znormal,gauss5[jj],gauss5[kk]) +
             a_x1[2]*femfunc.dphi2dz3(znormal,gauss5[jj],gauss5[kk]) +
             a_x1[3]*femfunc.dphi3dz3(znormal,gauss5[jj],gauss5[kk]) +
             a_x1[4]*femfunc.dphi4dz3(znormal,gauss5[jj],gauss5[kk]) +
             a_x1[5]*femfunc.dphi5dz3(znormal,gauss5[jj],gauss5[kk]) +
             a_x1[6]*femfunc.dphi6dz3(znormal,gauss5[jj],gauss5[kk]) +
             a_x1[7]*femfunc.dphi7dz3(znormal,gauss5[jj],gauss5[kk]) );

          F21 = (a_x2[0]*femfunc.dphi0dz1(znormal,gauss5[jj],gauss5[kk]) +
             a_x2[1]*femfunc.dphi1dz1(znormal,gauss5[jj],gauss5[kk]) +
             a_x2[2]*femfunc.dphi2dz1(znormal,gauss5[jj],gauss5[kk]) +
             a_x2[3]*femfunc.dphi3dz1(znormal,gauss5[jj],gauss5[kk]) +
             a_x2[4]*femfunc.dphi4dz1(znormal,gauss5[jj],gauss5[kk]) +
             a_x2[5]*femfunc.dphi5dz1(znormal,gauss5[jj],gauss5[kk]) +
             a_x2[6]*femfunc.dphi6dz1(znormal,gauss5[jj],gauss5[kk]) +
             a_x2[7]*femfunc.dphi7dz1(znormal,gauss5[jj],gauss5[kk]) );

          F22 = (a_x2[0]*femfunc.dphi0dz2(znormal,gauss5[jj],gauss5[kk]) +
             a_x2[1]*femfunc.dphi1dz2(znormal,gauss5[jj],gauss5[kk]) +
             a_x2[2]*femfunc.dphi2dz2(znormal,gauss5[jj],gauss5[kk]) +
             a_x2[3]*femfunc.dphi3dz2(znormal,gauss5[jj],gauss5[kk]) +
             a_x2[4]*femfunc.dphi4dz2(znormal,gauss5[jj],gauss5[kk]) +
             a_x2[5]*femfunc.dphi5dz2(znormal,gauss5[jj],gauss5[kk]) +
             a_x2[6]*femfunc.dphi6dz2(znormal,gauss5[jj],gauss5[kk]) +
             a_x2[7]*femfunc.dphi7dz2(znormal,gauss5[jj],gauss5[kk]) );

          F23 = (a_x2[0]*femfunc.dphi0dz3(znormal,gauss5[jj],gauss5[kk]) +
             a_x2[1]*femfunc.dphi1dz3(znormal,gauss5[jj],gauss5[kk]) +
             a_x2[2]*femfunc.dphi2dz3(znormal,gauss5[jj],gauss5[kk]) +
             a_x2[3]*femfunc.dphi3dz3(znormal,gauss5[jj],gauss5[kk]) +
             a_x2[4]*femfunc.dphi4dz3(znormal,gauss5[jj],gauss5[kk]) +
             a_x2[5]*femfunc.dphi5dz3(znormal,gauss5[jj],gauss5[kk]) +
             a_x2[6]*femfunc.dphi6dz3(znormal,gauss5[jj],gauss5[kk]) +
             a_x2[7]*femfunc.dphi7dz3(znormal,gauss5[jj],gauss5[kk]) );

          F31 = (a_x3[0]*femfunc.dphi0dz1(znormal,gauss5[jj],gauss5[kk]) +
             a_x3[1]*femfunc.dphi1dz1(znormal,gauss5[jj],gauss5[kk]) +
             a_x3[2]*femfunc.dphi2dz1(znormal,gauss5[jj],gauss5[kk]) +
             a_x3[3]*femfunc.dphi3dz1(znormal,gauss5[jj],gauss5[kk]) +
             a_x3[4]*femfunc.dphi4dz1(znormal,gauss5[jj],gauss5[kk]) +
             a_x3[5]*femfunc.dphi5dz1(znormal,gauss5[jj],gauss5[kk]) +
             a_x3[6]*femfunc.dphi6dz1(znormal,gauss5[jj],gauss5[kk]) +
             a_x3[7]*femfunc.dphi7dz1(znormal,gauss5[jj],gauss5[kk]) );

          F32 = (a_x3[0]*femfunc.dphi0dz2(znormal,gauss5[jj],gauss5[kk]) +
             a_x3[1]*femfunc.dphi1dz2(znormal,gauss5[jj],gauss5[kk]) +
             a_x3[2]*femfunc.dphi2dz2(znormal,gauss5[jj],gauss5[kk]) +
             a_x3[3]*femfunc.dphi3dz2(znormal,gauss5[jj],gauss5[kk]) +
             a_x3[4]*femfunc.dphi4dz2(znormal,gauss5[jj],gauss5[kk]) +
             a_x3[5]*femfunc.dphi5dz2(znormal,gauss5[jj],gauss5[kk]) +
             a_x3[6]*femfunc.dphi6dz2(znormal,gauss5[jj],gauss5[kk]) +
             a_x3[7]*femfunc.dphi7dz2(znormal,gauss5[jj],gauss5[kk]) );

          F33 = (a_x3[0]*femfunc.dphi0dz3(znormal,gauss5[jj],gauss5[kk]) +
             a_x3[1]*femfunc.dphi1dz3(znormal,gauss5[jj],gauss5[kk]) +
             a_x3[2]*femfunc.dphi2dz3(znormal,gauss5[jj],gauss5[kk]) +
             a_x3[3]*femfunc.dphi3dz3(znormal,gauss5[jj],gauss5[kk]) +
             a_x3[4]*femfunc.dphi4dz3(znormal,gauss5[jj],gauss5[kk]) +
             a_x3[5]*femfunc.dphi5dz3(znormal,gauss5[jj],gauss5[kk]) +
             a_x3[6]*femfunc.dphi6dz3(znormal,gauss5[jj],gauss5[kk]) +
             a_x3[7]*femfunc.dphi7dz3(znormal,gauss5[jj],gauss5[kk]) );

          Fmat[0][0] = F11; Fmat[0][1] = F12; Fmat[0][2] = F13;
          Fmat[1][0] = F21; Fmat[1][1] = F22; Fmat[1][2] = F23;
          Fmat[2][0] = F31; Fmat[2][1] = F32; Fmat[2][2] = F33;

          femfunc.inverse_mat(Fmat, Fmat_inv, Fmat_inv_trans);


          phivec = {femfunc.phi0(znormal,gauss5[jj],gauss5[kk]) , \
              femfunc.phi1(znormal,gauss5[jj],gauss5[kk]) , \
              femfunc.phi2(znormal,gauss5[jj],gauss5[kk]) , \
              femfunc.phi3(znormal,gauss5[jj],gauss5[kk]) , \
              femfunc.phi4(znormal,gauss5[jj],gauss5[kk]) , \
              femfunc.phi5(znormal,gauss5[jj],gauss5[kk]) , \
              femfunc.phi6(znormal,gauss5[jj],gauss5[kk]) , \
              femfunc.phi7(znormal,gauss5[jj],gauss5[kk])};

              for (int X = 0; X < 8; X++)
            {
              phimat[0][X] = phivec[X];
              phimat[1][X+8] = phivec[X];
              phimat[2][X+16] = phivec[X];

            }

              femfunc.mat_transpose(phimat, phimat_trans);

              femfunc.vec_mult_mat(a_zvector, Fmat_inv, inter21);

              femfunc.vec_mult_mat(inter21, Fmat_inv_trans, inter22);

              inter2 = sqrt(inter22[0]*a_zvector[0] + inter22[1]*a_zvector[1] + inter22[2]*a_zvector[2]);

              femfunc.mat_mult_vec(phimat_trans, a_traction, inter3);

              Jacobian = femfunc.Jacobian(Fmat);

              for (int qq = 0; qq < 24; qq++)
            {
              
                  fR1[qq] += w5[jj]*w5[kk]*inter3[qq]*Jacobian*inter2;
            }
        }
    }
    for (int qq = 0; qq < 24; qq++)
      {
        a_vector[qq] = fR1[qq];
      }
};

void blmintegrals :: integral3yy(double (&a_vector)[24],
                                 double (&a_x1)[8],
                                 double (&a_x2)[8],
                                 double (&a_x3)[8],
		 vector<double>& a_zvector,
		 vector<double>& a_traction)

{

  //Gauss Points for numerical integration

  vector<double> gauss5(5);
  vector<double> w5(5);
  double F11, F12, F13, F21, F22, F23, F31, F32, F33;
  double Jacobian, znormal;
  vector <double> phivec(8);
  vector<vector<double> > phimat(3,vector<double> (24));
  vector<vector<double> > phimat_trans(24,vector <double> (3));
  vector<vector<double> > Fmat(3,vector <double> (3));
  vector<vector<double> > Fmat_inv(3,vector <double> (3));
  vector<vector<double> > Fmat_inv_trans(3,vector <double> (3));
  vector<vector<double> > Tphi(6,vector <double> (24));
  vector<vector<double> > Tphi_trans(24,vector <double> (6));
  vector<double> inter21(3);
  vector<double> inter22(3);
  vector<double> inter3(24);
  double inter2;
  vector<double> fR1(24);
  femfunctions femfunc;

  gauss5 = {0.000000000000000, 0.538469310105683, -0.538469310105683,0.906179845938664, -0.906179845938664};

  w5 = {0.568888888888889, 0.478628670499366, 0.478628670499366,0.236926885056189, 0.236926885056189};

  znormal = a_zvector[1];
  for (int qq = 0; qq < 24; qq++)
    {
      fR1[qq] = 0;
    }
  
  for (int ii = 0; ii < 5; ii++)
    {
      for (int kk = 0; kk < 5; kk++)
        {
          F11 = (a_x1[0]*femfunc.dphi0dz1(gauss5[ii],znormal,gauss5[kk]) +
             a_x1[1]*femfunc.dphi1dz1(gauss5[ii],znormal,gauss5[kk]) +
             a_x1[2]*femfunc.dphi2dz1(gauss5[ii],znormal,gauss5[kk]) +
             a_x1[3]*femfunc.dphi3dz1(gauss5[ii],znormal,gauss5[kk]) +
             a_x1[4]*femfunc.dphi4dz1(gauss5[ii],znormal,gauss5[kk]) +
             a_x1[5]*femfunc.dphi5dz1(gauss5[ii],znormal,gauss5[kk]) +
             a_x1[6]*femfunc.dphi6dz1(gauss5[ii],znormal,gauss5[kk]) +
             a_x1[7]*femfunc.dphi7dz1(gauss5[ii],znormal,gauss5[kk]) );
          
          F12 = (a_x1[0]*femfunc.dphi0dz2(gauss5[ii],znormal,gauss5[kk]) +
             a_x1[1]*femfunc.dphi1dz2(gauss5[ii],znormal,gauss5[kk]) +
             a_x1[2]*femfunc.dphi2dz2(gauss5[ii],znormal,gauss5[kk]) +
             a_x1[3]*femfunc.dphi3dz2(gauss5[ii],znormal,gauss5[kk]) +
             a_x1[4]*femfunc.dphi4dz2(gauss5[ii],znormal,gauss5[kk]) +
             a_x1[5]*femfunc.dphi5dz2(gauss5[ii],znormal,gauss5[kk]) +
             a_x1[6]*femfunc.dphi6dz2(gauss5[ii],znormal,gauss5[kk]) +
             a_x1[7]*femfunc.dphi7dz2(gauss5[ii],znormal,gauss5[kk]) );

          F13 = (a_x1[0]*femfunc.dphi0dz3(gauss5[ii],znormal,gauss5[kk]) +
             a_x1[1]*femfunc.dphi1dz3(gauss5[ii],znormal,gauss5[kk]) +
             a_x1[2]*femfunc.dphi2dz3(gauss5[ii],znormal,gauss5[kk]) +
             a_x1[3]*femfunc.dphi3dz3(gauss5[ii],znormal,gauss5[kk]) +
             a_x1[4]*femfunc.dphi4dz3(gauss5[ii],znormal,gauss5[kk]) +
             a_x1[5]*femfunc.dphi5dz3(gauss5[ii],znormal,gauss5[kk]) +
             a_x1[6]*femfunc.dphi6dz3(gauss5[ii],znormal,gauss5[kk]) +
             a_x1[7]*femfunc.dphi7dz3(gauss5[ii],znormal,gauss5[kk]) );

          F21 = (a_x2[0]*femfunc.dphi0dz1(gauss5[ii],znormal,gauss5[kk]) +
             a_x2[1]*femfunc.dphi1dz1(gauss5[ii],znormal,gauss5[kk]) +
             a_x2[2]*femfunc.dphi2dz1(gauss5[ii],znormal,gauss5[kk]) +
             a_x2[3]*femfunc.dphi3dz1(gauss5[ii],znormal,gauss5[kk]) +
             a_x2[4]*femfunc.dphi4dz1(gauss5[ii],znormal,gauss5[kk]) +
             a_x2[5]*femfunc.dphi5dz1(gauss5[ii],znormal,gauss5[kk]) +
             a_x2[6]*femfunc.dphi6dz1(gauss5[ii],znormal,gauss5[kk]) +
             a_x2[7]*femfunc.dphi7dz1(gauss5[ii],znormal,gauss5[kk]) );

          F22 = (a_x2[0]*femfunc.dphi0dz2(gauss5[ii],znormal,gauss5[kk]) +
             a_x2[1]*femfunc.dphi1dz2(gauss5[ii],znormal,gauss5[kk]) +
             a_x2[2]*femfunc.dphi2dz2(gauss5[ii],znormal,gauss5[kk]) +
             a_x2[3]*femfunc.dphi3dz2(gauss5[ii],znormal,gauss5[kk]) +
             a_x2[4]*femfunc.dphi4dz2(gauss5[ii],znormal,gauss5[kk]) +
             a_x2[5]*femfunc.dphi5dz2(gauss5[ii],znormal,gauss5[kk]) +
             a_x2[6]*femfunc.dphi6dz2(gauss5[ii],znormal,gauss5[kk]) +
             a_x2[7]*femfunc.dphi7dz2(gauss5[ii],znormal,gauss5[kk]) );

          F23 = (a_x2[0]*femfunc.dphi0dz3(gauss5[ii],znormal,gauss5[kk]) +
             a_x2[1]*femfunc.dphi1dz3(gauss5[ii],znormal,gauss5[kk]) +
             a_x2[2]*femfunc.dphi2dz3(gauss5[ii],znormal,gauss5[kk]) +
             a_x2[3]*femfunc.dphi3dz3(gauss5[ii],znormal,gauss5[kk]) +
             a_x2[4]*femfunc.dphi4dz3(gauss5[ii],znormal,gauss5[kk]) +
             a_x2[5]*femfunc.dphi5dz3(gauss5[ii],znormal,gauss5[kk]) +
             a_x2[6]*femfunc.dphi6dz3(gauss5[ii],znormal,gauss5[kk]) +
             a_x2[7]*femfunc.dphi7dz3(gauss5[ii],znormal,gauss5[kk]) );

          F31 = (a_x3[0]*femfunc.dphi0dz1(gauss5[ii],znormal,gauss5[kk]) +
             a_x3[1]*femfunc.dphi1dz1(gauss5[ii],znormal,gauss5[kk]) +
             a_x3[2]*femfunc.dphi2dz1(gauss5[ii],znormal,gauss5[kk]) +
             a_x3[3]*femfunc.dphi3dz1(gauss5[ii],znormal,gauss5[kk]) +
             a_x3[4]*femfunc.dphi4dz1(gauss5[ii],znormal,gauss5[kk]) +
             a_x3[5]*femfunc.dphi5dz1(gauss5[ii],znormal,gauss5[kk]) +
             a_x3[6]*femfunc.dphi6dz1(gauss5[ii],znormal,gauss5[kk]) +
             a_x3[7]*femfunc.dphi7dz1(gauss5[ii],znormal,gauss5[kk]) );

          F32 = (a_x3[0]*femfunc.dphi0dz2(gauss5[ii],znormal,gauss5[kk]) +
             a_x3[1]*femfunc.dphi1dz2(gauss5[ii],znormal,gauss5[kk]) +
             a_x3[2]*femfunc.dphi2dz2(gauss5[ii],znormal,gauss5[kk]) +
             a_x3[3]*femfunc.dphi3dz2(gauss5[ii],znormal,gauss5[kk]) +
             a_x3[4]*femfunc.dphi4dz2(gauss5[ii],znormal,gauss5[kk]) +
             a_x3[5]*femfunc.dphi5dz2(gauss5[ii],znormal,gauss5[kk]) +
             a_x3[6]*femfunc.dphi6dz2(gauss5[ii],znormal,gauss5[kk]) +
             a_x3[7]*femfunc.dphi7dz2(gauss5[ii],znormal,gauss5[kk]) );

          F33 = (a_x3[0]*femfunc.dphi0dz3(gauss5[ii],znormal,gauss5[kk]) +
             a_x3[1]*femfunc.dphi1dz3(gauss5[ii],znormal,gauss5[kk]) +
             a_x3[2]*femfunc.dphi2dz3(gauss5[ii],znormal,gauss5[kk]) +
             a_x3[3]*femfunc.dphi3dz3(gauss5[ii],znormal,gauss5[kk]) +
             a_x3[4]*femfunc.dphi4dz3(gauss5[ii],znormal,gauss5[kk]) +
             a_x3[5]*femfunc.dphi5dz3(gauss5[ii],znormal,gauss5[kk]) +
             a_x3[6]*femfunc.dphi6dz3(gauss5[ii],znormal,gauss5[kk]) +
             a_x3[7]*femfunc.dphi7dz3(gauss5[ii],znormal,gauss5[kk]) );

          Fmat[0][0] = F11; Fmat[0][1] = F12; Fmat[0][2] = F13;
          Fmat[1][0] = F21; Fmat[1][1] = F22; Fmat[1][2] = F23;
          Fmat[2][0] = F31; Fmat[2][1] = F32; Fmat[2][2] = F33;

          femfunc.inverse_mat(Fmat, Fmat_inv, Fmat_inv_trans);


          phivec = {femfunc.phi0(gauss5[ii],znormal,gauss5[kk]) , \
              femfunc.phi1(gauss5[ii],znormal,gauss5[kk]) , \
              femfunc.phi2(gauss5[ii],znormal,gauss5[kk]) , \
              femfunc.phi3(gauss5[ii],znormal,gauss5[kk]) , \
              femfunc.phi4(gauss5[ii],znormal,gauss5[kk]) , \
              femfunc.phi5(gauss5[ii],znormal,gauss5[kk]) , \
              femfunc.phi6(gauss5[ii],znormal,gauss5[kk]) , \
              femfunc.phi7(gauss5[ii],znormal,gauss5[kk])};

              for (int X = 0; X < 8; X++)
            {
              phimat[0][X] = phivec[X];
              phimat[1][X+8] = phivec[X];
              phimat[2][X+16] = phivec[X];

            }

              femfunc.mat_transpose(phimat, phimat_trans);

              femfunc.vec_mult_mat(a_zvector, Fmat_inv, inter21);

              femfunc.vec_mult_mat(inter21, Fmat_inv_trans, inter22);

              inter2 = sqrt(inter22[0]*a_zvector[0] + inter22[1]*a_zvector[1] + inter22[2]*a_zvector[2]);

              femfunc.mat_mult_vec(phimat_trans, a_traction, inter3);

              Jacobian = femfunc.Jacobian(Fmat);

              for (int qq = 0; qq < 24; qq++)
            {
              
                  fR1[qq] += w5[ii]*w5[kk]*inter3[qq]*Jacobian*inter2;
            }
        }
    }
    for (int qq = 0; qq < 24; qq++)
      {
        a_vector[qq] = fR1[qq];
      }
};

void blmintegrals :: integral3zz(double (&a_vector)[24],
                                 double (&a_x1)[8],
                                 double (&a_x2)[8],
                                 double (&a_x3)[8],
		 vector<double>& a_zvector,
		 vector<double>& a_traction)

{

  //Gauss Points for numerical integration

  vector<double> gauss5(5);
  vector<double> w5(5);
  double F11, F12, F13, F21, F22, F23, F31, F32, F33;
  double Jacobian, znormal;
  vector <double> phivec(8);
  vector<vector<double> > phimat(3,vector<double> (24));
  vector<vector<double> > phimat_trans(24,vector <double> (3));
  vector<vector<double> > Fmat(3,vector <double> (3));
  vector<vector<double> > Fmat_inv(3,vector <double> (3));
  vector<vector<double> > Fmat_inv_trans(3,vector <double> (3));
  vector<vector<double> > Tphi(6,vector <double> (24));
  vector<vector<double> > Tphi_trans(24,vector <double> (6));
  vector<double> inter21(3);
  vector<double> inter22(3);
  vector<double> inter3(24);
  double inter2;
  vector<double> fR1(24);
  femfunctions femfunc;

  gauss5 = {0.000000000000000, 0.538469310105683, -0.538469310105683,0.906179845938664, -0.906179845938664};

  w5 = {0.568888888888889, 0.478628670499366, 0.478628670499366,0.236926885056189, 0.236926885056189};

  znormal = a_zvector[2];

  for (int qq = 0; qq < 24; qq++)
    {
      fR1[qq] = 0;
    }
  
  for (int ii = 0; ii < 5; ii++)
    {
      for (int jj = 0; jj < 5; jj++)
        {
          F11 = (a_x1[0]*femfunc.dphi0dz1(gauss5[ii],gauss5[jj],znormal) +
             a_x1[1]*femfunc.dphi1dz1(gauss5[ii],gauss5[jj],znormal) +
             a_x1[2]*femfunc.dphi2dz1(gauss5[ii],gauss5[jj],znormal) +
             a_x1[3]*femfunc.dphi3dz1(gauss5[ii],gauss5[jj],znormal) +
             a_x1[4]*femfunc.dphi4dz1(gauss5[ii],gauss5[jj],znormal) +
             a_x1[5]*femfunc.dphi5dz1(gauss5[ii],gauss5[jj],znormal) +
             a_x1[6]*femfunc.dphi6dz1(gauss5[ii],gauss5[jj],znormal) +
             a_x1[7]*femfunc.dphi7dz1(gauss5[ii],gauss5[jj],znormal) );
          
          F12 = (a_x1[0]*femfunc.dphi0dz2(gauss5[ii],gauss5[jj],znormal) +
             a_x1[1]*femfunc.dphi1dz2(gauss5[ii],gauss5[jj],znormal) +
             a_x1[2]*femfunc.dphi2dz2(gauss5[ii],gauss5[jj],znormal) +
             a_x1[3]*femfunc.dphi3dz2(gauss5[ii],gauss5[jj],znormal) +
             a_x1[4]*femfunc.dphi4dz2(gauss5[ii],gauss5[jj],znormal) +
             a_x1[5]*femfunc.dphi5dz2(gauss5[ii],gauss5[jj],znormal) +
             a_x1[6]*femfunc.dphi6dz2(gauss5[ii],gauss5[jj],znormal) +
             a_x1[7]*femfunc.dphi7dz2(gauss5[ii],gauss5[jj],znormal) );

          F13 = (a_x1[0]*femfunc.dphi0dz3(gauss5[ii],gauss5[jj],znormal) +
             a_x1[1]*femfunc.dphi1dz3(gauss5[ii],gauss5[jj],znormal) +
             a_x1[2]*femfunc.dphi2dz3(gauss5[ii],gauss5[jj],znormal) +
             a_x1[3]*femfunc.dphi3dz3(gauss5[ii],gauss5[jj],znormal) +
             a_x1[4]*femfunc.dphi4dz3(gauss5[ii],gauss5[jj],znormal) +
             a_x1[5]*femfunc.dphi5dz3(gauss5[ii],gauss5[jj],znormal) +
             a_x1[6]*femfunc.dphi6dz3(gauss5[ii],gauss5[jj],znormal) +
             a_x1[7]*femfunc.dphi7dz3(gauss5[ii],gauss5[jj],znormal) );

          F21 = (a_x2[0]*femfunc.dphi0dz1(gauss5[ii],gauss5[jj],znormal) +
             a_x2[1]*femfunc.dphi1dz1(gauss5[ii],gauss5[jj],znormal) +
             a_x2[2]*femfunc.dphi2dz1(gauss5[ii],gauss5[jj],znormal) +
             a_x2[3]*femfunc.dphi3dz1(gauss5[ii],gauss5[jj],znormal) +
             a_x2[4]*femfunc.dphi4dz1(gauss5[ii],gauss5[jj],znormal) +
             a_x2[5]*femfunc.dphi5dz1(gauss5[ii],gauss5[jj],znormal) +
             a_x2[6]*femfunc.dphi6dz1(gauss5[ii],gauss5[jj],znormal) +
             a_x2[7]*femfunc.dphi7dz1(gauss5[ii],gauss5[jj],znormal) );

          F22 = (a_x2[0]*femfunc.dphi0dz2(gauss5[ii],gauss5[jj],znormal) +
             a_x2[1]*femfunc.dphi1dz2(gauss5[ii],gauss5[jj],znormal) +
             a_x2[2]*femfunc.dphi2dz2(gauss5[ii],gauss5[jj],znormal) +
             a_x2[3]*femfunc.dphi3dz2(gauss5[ii],gauss5[jj],znormal) +
             a_x2[4]*femfunc.dphi4dz2(gauss5[ii],gauss5[jj],znormal) +
             a_x2[5]*femfunc.dphi5dz2(gauss5[ii],gauss5[jj],znormal) +
             a_x2[6]*femfunc.dphi6dz2(gauss5[ii],gauss5[jj],znormal) +
             a_x2[7]*femfunc.dphi7dz2(gauss5[ii],gauss5[jj],znormal) );

          F23 = (a_x2[0]*femfunc.dphi0dz3(gauss5[ii],gauss5[jj],znormal) +
             a_x2[1]*femfunc.dphi1dz3(gauss5[ii],gauss5[jj],znormal) +
             a_x2[2]*femfunc.dphi2dz3(gauss5[ii],gauss5[jj],znormal) +
             a_x2[3]*femfunc.dphi3dz3(gauss5[ii],gauss5[jj],znormal) +
             a_x2[4]*femfunc.dphi4dz3(gauss5[ii],gauss5[jj],znormal) +
             a_x2[5]*femfunc.dphi5dz3(gauss5[ii],gauss5[jj],znormal) +
             a_x2[6]*femfunc.dphi6dz3(gauss5[ii],gauss5[jj],znormal) +
             a_x2[7]*femfunc.dphi7dz3(gauss5[ii],gauss5[jj],znormal) );

          F31 = (a_x3[0]*femfunc.dphi0dz1(gauss5[ii],gauss5[jj],znormal) +
             a_x3[1]*femfunc.dphi1dz1(gauss5[ii],gauss5[jj],znormal) +
             a_x3[2]*femfunc.dphi2dz1(gauss5[ii],gauss5[jj],znormal) +
             a_x3[3]*femfunc.dphi3dz1(gauss5[ii],gauss5[jj],znormal) +
             a_x3[4]*femfunc.dphi4dz1(gauss5[ii],gauss5[jj],znormal) +
             a_x3[5]*femfunc.dphi5dz1(gauss5[ii],gauss5[jj],znormal) +
             a_x3[6]*femfunc.dphi6dz1(gauss5[ii],gauss5[jj],znormal) +
             a_x3[7]*femfunc.dphi7dz1(gauss5[ii],gauss5[jj],znormal) );

          F32 = (a_x3[0]*femfunc.dphi0dz2(gauss5[ii],gauss5[jj],znormal) +
             a_x3[1]*femfunc.dphi1dz2(gauss5[ii],gauss5[jj],znormal) +
             a_x3[2]*femfunc.dphi2dz2(gauss5[ii],gauss5[jj],znormal) +
             a_x3[3]*femfunc.dphi3dz2(gauss5[ii],gauss5[jj],znormal) +
             a_x3[4]*femfunc.dphi4dz2(gauss5[ii],gauss5[jj],znormal) +
             a_x3[5]*femfunc.dphi5dz2(gauss5[ii],gauss5[jj],znormal) +
             a_x3[6]*femfunc.dphi6dz2(gauss5[ii],gauss5[jj],znormal) +
             a_x3[7]*femfunc.dphi7dz2(gauss5[ii],gauss5[jj],znormal) );

          F33 = (a_x3[0]*femfunc.dphi0dz3(gauss5[ii],gauss5[jj],znormal) +
             a_x3[1]*femfunc.dphi1dz3(gauss5[ii],gauss5[jj],znormal) +
             a_x3[2]*femfunc.dphi2dz3(gauss5[ii],gauss5[jj],znormal) +
             a_x3[3]*femfunc.dphi3dz3(gauss5[ii],gauss5[jj],znormal) +
             a_x3[4]*femfunc.dphi4dz3(gauss5[ii],gauss5[jj],znormal) +
             a_x3[5]*femfunc.dphi5dz3(gauss5[ii],gauss5[jj],znormal) +
             a_x3[6]*femfunc.dphi6dz3(gauss5[ii],gauss5[jj],znormal) +
             a_x3[7]*femfunc.dphi7dz3(gauss5[ii],gauss5[jj],znormal) );

          Fmat[0][0] = F11; Fmat[0][1] = F12; Fmat[0][2] = F13;
          Fmat[1][0] = F21; Fmat[1][1] = F22; Fmat[1][2] = F23;
          Fmat[2][0] = F31; Fmat[2][1] = F32; Fmat[2][2] = F33;

          femfunc.inverse_mat(Fmat, Fmat_inv, Fmat_inv_trans);


          phivec = {femfunc.phi0(gauss5[ii],gauss5[jj],znormal) , \
              femfunc.phi1(gauss5[ii],gauss5[jj],znormal) , \
              femfunc.phi2(gauss5[ii],gauss5[jj],znormal) , \
              femfunc.phi3(gauss5[ii],gauss5[jj],znormal) , \
              femfunc.phi4(gauss5[ii],gauss5[jj],znormal) , \
              femfunc.phi5(gauss5[ii],gauss5[jj],znormal) , \
              femfunc.phi6(gauss5[ii],gauss5[jj],znormal) , \
              femfunc.phi7(gauss5[ii],gauss5[jj],znormal)};

              for (int X = 0; X < 8; X++)
            {
              phimat[0][X] = phivec[X];
              phimat[1][X+8] = phivec[X];
              phimat[2][X+16] = phivec[X];

            }

              femfunc.mat_transpose(phimat, phimat_trans);

              femfunc.vec_mult_mat(a_zvector, Fmat_inv, inter21);

              femfunc.vec_mult_mat(inter21, Fmat_inv_trans, inter22);

              inter2 = sqrt(inter22[0]*a_zvector[0] + inter22[1]*a_zvector[1] + inter22[2]*a_zvector[2]);

              femfunc.mat_mult_vec(phimat_trans, a_traction, inter3);

              Jacobian = femfunc.Jacobian(Fmat);

              for (int qq = 0; qq < 24; qq++)
                {
                  fR1[qq] += w5[ii]*w5[jj]*inter3[qq]*Jacobian*inter2;
                }
        }
    }
    for (int qq = 0; qq < 24; qq++)
      {
        a_vector[qq] = fR1[qq];
      }
};