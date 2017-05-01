#include <iostream>
#include <vector>
#include <fstream>
#include "blmintegrals.H"
#include "femfunctions.H"
#include "mesher.H"
#include "hsbounds.H"
#include "GridWrite.H"
#include "CGSolver.H"
#include "SparseMatrix.H"
#include "JacobiSolver.H"
#include <omp.h>

using namespace std;
int main(int argc, char* argv[])
{
  //Step 1: Generate the mesh
    /*
  if(argc != 2)
    {
      cout << "this program takes one argument that is the input file\n";
      return 1;
    }
     */
    char* fileName = NULL;
    int numThreads;

  for (int i=1; i < argc; i++)
  {
      if (string(argv[i]) == "-i") {
          fileName = argv[i + 1];
      }
      else if (string(argv[i]) == "-n") {
          numThreads = atoi(argv[i+1]);
          omp_set_num_threads(numThreads);
      }
    }

  if (!fileName)
    {
        cout << "Provide an input file!" << endl;
        return 1;
    }

  string inFile(fileName);
  string lines[8];
  string words[20];
  ifstream myfile(inFile.c_str());
  myfile >> words[0] >> words[1] >> words[2] >> words[3] >> words[4]  >> words[5];
  int N;
  myfile >> N;
  myfile >> words[6] >> words[7];
  double forceIn[3];
  myfile >> forceIn[0] >> forceIn[1] >> forceIn[2];
  myfile >> words[8] >> words[9];
  double tractionIn[3];
  myfile >> tractionIn[0] >> tractionIn[1] >> tractionIn[2];
  myfile >> words[10] >> words[11] >> words[12] >> words[13] >> words[14];
  double k1, k2, u1, u2, v2;
  myfile >> k1 >> u1 >> k2 >> u2 >> v2;

  // Create clock
  double start;
  double startTotal = omp_get_wtime();

  //int N = 10; //Number of nodes in each direction
  int n_nodes = N*N*N; //Number of nodes
  int n_elem = (N-1)*(N-1)*(N-1); //number of elements
  //FILE *fp = NULL;
  //FILE *fp2 = NULL;
  vector<vector<double>> node_table; //Node table initialized
  vector<vector<int>> conn_table; //Connectivity table initialized
  Mesher mesh(N);
  mesh.createmesh();
  node_table = mesh.getNode();
  conn_table = mesh.getConn();
  //Extend the connectivity to include the 3 dimension of the elasticity solution
  vector<vector<int>> conn_table2(n_elem,vector<int>(8));
  vector<vector<int>> conn_table3(n_elem,vector<int>(8));


  //#pragma omp parallel for //Not worth it?
  for (int qq = 0; qq < n_elem; qq++)
    {
      for (int ww = 0; ww < 8; ww++)
        {
          conn_table2[qq][ww] = conn_table[qq][ww] + n_nodes;
          conn_table3[qq][ww] = conn_table2[qq][ww] + n_nodes;
        }
    }

  //Step 2 - Define the inputs and initialize the required arrays
  vector<double> traction_elem(24); //traction vector per element

  //Define the conditions
  vector<double> force(3), traction(3); //body force and traction vector
  //force = {0.0,0.0,0.0};
  for(int i=0;i<3;i++)
    {
      force[i]=forceIn[i];
    }
  double kstar, ustar; //material properties
  //k1 = 80; k2 = 100; u1 = 60; u2 = 50; v2 = 0.3;

  //Call the hsbounds to get the effective properties
  hsbounds hsbound;
  hsbound.HS_bounds(k1,k2,u1,u2,v2,kstar,ustar);

  vector<vector<double>> Etensor(6,vector<double>(6)); //elasticity tensor
  //Construct the elasticity tensor
  //#pragma omp parallel for
  for (int qq = 0; qq < 6; qq++)
    {
      for (int ww = 0; ww < 6; ww++)
        {
          Etensor[qq][ww] = 0;
        }
    }

  Etensor[0][0] = kstar + (4.0/3.0)*ustar;
  Etensor[0][1] = kstar - (2.0/3.0)*ustar;
  Etensor[0][2] = kstar - (2.0/3.0)*ustar;
  Etensor[1][0] = kstar - (2.0/3.0)*ustar;
  Etensor[1][1] = kstar + (4.0/3.0)*ustar;
  Etensor[1][2] = kstar - (2.0/3.0)*ustar;
  Etensor[2][0] = kstar - (2.0/3.0)*ustar;
  Etensor[2][1] = kstar - (2.0/3.0)*ustar;
  Etensor[2][2] = kstar + (4.0/3.0)*ustar;
  Etensor[3][3] = ustar;
  Etensor[4][4] = ustar;
  Etensor[5][5] = ustar;
  
  //Step 3 - Start the compution by looping over the elements
    start = omp_get_wtime();
    ////cout << "Starting matrix construction " << endl;
  blmintegrals blm_integrate;
  SparseMatrix stiffness_matrix(3*n_nodes,3*n_nodes); //full stiffness matrix - Sparse matrix
  vector<double> load_vector(3*n_nodes); //full load vector

  double x1[8], x2[8], x3[8]; //space coordinates
  double stiff_elem[24][24];

  #pragma omp parallel for private(x1, x2, x3, stiff_elem)
  for(int e = 0; e < n_elem; e++)
    {
      ////cout << "Hello from thread " << omp_get_thread_num() << endl;

      for (int X = 0; X < 8; X++) {
          x1[X] = node_table[conn_table[e][X]][0];
          x2[X] = node_table[conn_table[e][X]][1];
          x3[X] = node_table[conn_table[e][X]][2];
      }
    blm_integrate.integral1(stiff_elem, x1, x2, x3, Etensor);
      for (int ii = 0; ii < 8; ii++) {
          for (int jj = 0; jj < 8; jj++) {
              double stf = stiff_elem[ii][jj];;
              if (stf > 0.00001 || stf < -0.00001) {
                  stiffness_matrix[{{conn_table[e][ii], conn_table[e][jj]}}] += stiff_elem[ii][jj]; //first component
              }
              stf = stiff_elem[ii + 8][jj + 8];;
              if (stf > 0.00001 || stf < -0.00001) {
                  stiffness_matrix[{{conn_table2[e][ii], conn_table2[e][jj]}}] += stiff_elem[ii + 8][jj + 8]; //second component
              }
              stf = stiff_elem[ii + 16][jj + 16];;
              if (stf > 0.00001 || stf < -0.00001) {
                  stiffness_matrix[{{conn_table3[e][ii], conn_table3[e][jj]}}] += stiff_elem[ii + 16][jj +
                                                                                                      16]; //second component
              }
          }
          }
      double load_elem[24]; //load vector per element
          blm_integrate.integral2(load_elem, x1, x2, x3, force);

      for (int ii = 0; ii < 8; ii++)
        {
          load_vector[conn_table[e][ii]] += load_elem[ii]; //first component
          load_vector[conn_table2[e][ii]] += load_elem[ii+8]; //second component
          load_vector[conn_table3[e][ii]] += load_elem[ii+16]; //third component
        }
        ////cout << "Thread " << omp_get_thread_num() << " stiff_elem[0][0] = " << stiff_elem[0][1] << endl;
    }
    //return 0;

    //cout << "Matrix construction, " << (omp_get_wtime() - start) << endl;

  //Integral 3 - for boundary conditions

    ////cout << "Setting first boundary conditions construction " << endl;
    start = omp_get_wtime();
  int p = stiffness_matrix.N();
  const SparseMatrix constStiff = stiffness_matrix;

  #pragma omp parallel for
  for(int e = 0; e < ((N-1)*(N-1)); e++)
    {
          //First Component
          for (int aa = 0; aa < p; aa++) {
              const double stf = constStiff[{{conn_table[e][0], aa}}];
              if (stf != 0.) {
                  #pragma omp critical (zero)
                  stiffness_matrix[{{conn_table[e][0], aa}}] = 0.0;
              }
          }
          for (int aa = 0; aa < p; aa++) {
              const double stf = constStiff[{{conn_table[e][1], aa}}];
              if (stf != 0.) {
                  #pragma omp critical (one)
                  stiffness_matrix[{{conn_table[e][1], aa}}] = 0.0;
              }
          }
          for (int aa = 0; aa < p; aa++) {
              const double stf = constStiff[{{conn_table[e][2], aa}}];
              if (stf != 0.) {
                  #pragma omp critical (two)
                  stiffness_matrix[{{conn_table[e][2], aa}}] = 0.0;
              }
          }
          for (int aa = 0; aa < p; aa++) {
              const double stf = constStiff[{{conn_table[e][3], aa}}];
              if (stf != 0.) {
                  #pragma omp critical (three)
                  stiffness_matrix[{{conn_table[e][3], aa}}] = 0.0;
              }
          }

      #pragma omp critical
      {
          stiffness_matrix[{{conn_table[e][0], conn_table[e][0]}}] = 1.0;
          stiffness_matrix[{{conn_table[e][1], conn_table[e][1]}}] = 1.0;
          stiffness_matrix[{{conn_table[e][2], conn_table[e][2]}}] = 1.0;
          stiffness_matrix[{{conn_table[e][3], conn_table[e][3]}}] = 1.0;
      }
          //Second Component
          for (int aa = 0; aa < p; aa++) {
              const double stf = constStiff[{{conn_table2[e][0], aa}}];
              if (stf != 0.) {
                  #pragma omp critical (zero2)
                  stiffness_matrix[{{conn_table2[e][0], aa}}] = 0.0;
              }
          }
          for (int aa = 0; aa < p; aa++) {
              const double stf = constStiff[{{conn_table2[e][1], aa}}];
              if (stf != 0.) {
                  #pragma omp critical (one2)
                  stiffness_matrix[{{conn_table2[e][1], aa}}] = 0.0;
              }
          }
          for (int aa = 0; aa < p; aa++) {
              const double stf = constStiff[{{conn_table2[e][2], aa}}];
              if (stf != 0.) {
                  #pragma omp critical (two2)
                  stiffness_matrix[{{conn_table2[e][2], aa}}] = 0.0;
              }
          }
          for (int aa = 0; aa < p; aa++) {
              const double stf = constStiff[{{conn_table2[e][3], aa}}];
              if (stf != 0.) {
                  #pragma omp critical (three2)
                  stiffness_matrix[{{conn_table2[e][3], aa}}] = 0.0;
              }
          }
      #pragma omp critical
      {
          stiffness_matrix[{{conn_table2[e][0], conn_table2[e][0]}}] = 1.0;
          stiffness_matrix[{{conn_table2[e][1], conn_table2[e][1]}}] = 1.0;
          stiffness_matrix[{{conn_table2[e][2], conn_table2[e][2]}}] = 1.0;
          stiffness_matrix[{{conn_table2[e][3], conn_table2[e][3]}}] = 1.0;
      }
          //Third Component
          for (int aa = 0; aa < p; aa++) {
              const double stf = constStiff[{{conn_table3[e][0], aa}}];
              if (stf != 0.) {
                  #pragma omp critical (zero3)
                  stiffness_matrix[{{conn_table3[e][0], aa}}] = 0.0;
              }
          }
          for (int aa = 0; aa < p; aa++) {
              const double stf = constStiff[{{conn_table3[e][1], aa}}];
              if (stf != 0.) {
                  #pragma omp critical (one3)
                  stiffness_matrix[{{conn_table3[e][1], aa}}] = 0.0;
              }
          }
          for (int aa = 0; aa < p; aa++) {
              const double stf = constStiff[{{conn_table3[e][2], aa}}];
              if (stf != 0.) {
                  #pragma omp critical (two3)
                  stiffness_matrix[{{conn_table3[e][2], aa}}] = 0.0;
              }
          }
          for (int aa = 0; aa < p; aa++) {
              const double stf = constStiff[{{conn_table3[e][3], aa}}];
              if (stf != 0.) {
                  #pragma omp critical (three3)
                  stiffness_matrix[{{conn_table3[e][3], aa}}] = 0.0;
              }
          }
      #pragma omp critical
      {
          stiffness_matrix[{{conn_table3[e][0], conn_table3[e][0]}}] = 1.0;
          stiffness_matrix[{{conn_table3[e][1], conn_table3[e][1]}}] = 1.0;
          stiffness_matrix[{{conn_table3[e][2], conn_table3[e][2]}}] = 1.0;
          stiffness_matrix[{{conn_table3[e][3], conn_table3[e][3]}}] = 1.0;
      }
    }

    //cout << "Dirichlet BC, " << (omp_get_wtime() - start) << endl;
    start = omp_get_wtime();

  //Traction boundary condition
  //traction = {0.0,0.0,10.0};
  for(int i=0;i<3;i++)
    {
      traction[i]=tractionIn[i];
    }
  vector<double> zvector(3); //direction for the surface boundary condition
  zvector = {0.0,0.0,1.0};

  #pragma omp for
  for (int e = (n_elem - (N-1)*(N-1)); e < n_elem; e++)
    {
        double x1[8], x2[8], x3[8]; //space coordinates
      for (int X = 0; X < 8; X++)
        {
          x1[X] = node_table[conn_table[e][X]][0];
          x2[X] = node_table[conn_table[e][X]][1];
          x3[X] = node_table[conn_table[e][X]][2];
        }
        double surface_elem[24]; //load vector per element for surface boundary condition

      blm_integrate.integral3zz(surface_elem, x1, x2, x3, zvector, traction);

      for (int ii = 0; ii < 8; ii++)
        {
          load_vector[conn_table[e][ii]] += surface_elem[ii]; //first component
          load_vector[conn_table2[e][ii]] += surface_elem[ii+8]; //second component
          load_vector[conn_table3[e][ii]] += surface_elem[ii+16]; //third component  
        }
    }
    //cout << "Neumann BC, " << (omp_get_wtime() - start) << endl;
    //cout << "Total Construction, " << (omp_get_wtime() - startTotal) << endl;

  start = omp_get_wtime();
  //Step 4 - Solve the system of equations
  vector<double> solution_vector(3*n_nodes);
  vector<vector<double>> ufull(n_nodes,vector<double>(3));

  CGSolver solver;
  //JacobiSolver solver;
  int iterations = 10000;
  float residual = solver.solve(stiffness_matrix,load_vector,1E-6,iterations,solution_vector);

  //cout << "CG Solver, " << (omp_get_wtime() - start) << endl;

  vector<vector<double>> u_full(n_nodes,vector<double>(3));
  //Step 4.5 - Separate the components of the solution
  for(int ii = 0; ii < n_nodes; ii++)
    {
      u_full[ii][0] = solution_vector[ii];
      u_full[ii][1] = solution_vector[ii+n_nodes];
      u_full[ii][2] = solution_vector[ii+2*n_nodes];
    }

  //Step 5 - create the graphical output files
  GridWrite(node_table,conn_table,u_full,"projectSolution.vtu");
    //cout << "Total Runtime, " << (omp_get_wtime() - startTotal) << endl;
};
