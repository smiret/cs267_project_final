#include <iostream>
#include <vector>
#include <list>
#include <fstream>
#include "blmintegrals.H"
#include "mesher.H"
#include "hsbounds.H"
#include "GridWrite.H"
#include <omp.h>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

using namespace std;

// Eigen typedefs
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> triplet;

int main(int argc, char** argv)
{
    // Initialize Eigen threads
    Eigen::initParallel();

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

    cout << "Threads, " << Eigen::nbThreads() << endl;

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

    // Timing variables
    double start;
    double startTotal = omp_get_wtime();

    int n_nodes = N*N*N; //Number of nodes
    int n_elem = (N-1)*(N-1)*(N-1); //number of elements
    vector<vector<double>> node_table; //Node table initialized
    vector<vector<int>> conn_table; //Connectivity table initialized
    Mesher mesh(N);
    mesh.createmesh();
    node_table = mesh.getNode();
    conn_table = mesh.getConn();

    //Extend the connectivity to include the 3 dimension of the elasticity solution
    vector<vector<int>> conn_table2(n_elem,vector<int>(8));
    vector<vector<int>> conn_table3(n_elem,vector<int>(8));

    for (int qq = 0; qq < n_elem; qq++)
    {
        for (int ww = 0; ww < 8; ww++)
        {
            conn_table2[qq][ww] = conn_table[qq][ww] + n_nodes;
            conn_table3[qq][ww] = conn_table2[qq][ww] + n_nodes;
        }
    }

    //Step 2 - Define the inputs and initialize the required arrays
    //Define the conditions
    vector<double> force(3), traction(3); //body force and traction vector
    //force = {0.0,0.0,0.0};
    for(int i=0;i<3;i++)
    {
        force[i]=forceIn[i];
    }
    double kstar, ustar; //material properties

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

    blmintegrals blm_integrate;
    list<triplet> stiffness_tripletList;
    SpMat stiffness_matrix(3*n_nodes,3*n_nodes); //full stiffness matrix - Sparse matrix
    Eigen::VectorXd load_vector(3*n_nodes); //full load vector

    double x1[8], x2[8], x3[8]; //space coordinates
    double stiff_elem[24][24];
    double load_elem[24];

    #pragma omp parallel for private(x1, x2, x3, stiff_elem, load_elem)
    for(int e = 0; e < ((N-1)*(N-1)); e++)
    {
        for(int X = 0; X < 8; X++)
        {
            x1[X] = node_table[conn_table[e][X]][0];
            x2[X] = node_table[conn_table[e][X]][1];
            x3[X] = node_table[conn_table[e][X]][2];
        }

        //Integral 1
        blm_integrate.integral1(stiff_elem, x1, x2, x3, Etensor);

        for (int ii = 4; ii < 8; ii++)
        {
            for (int jj = 4; jj < 8; jj++)
            {
                double stf = stiff_elem[ii][jj];;
                if(stf > 0.00001 || stf < -0.00001)
                {
                    //first component
#pragma omp critical
                    stiffness_tripletList.push_back(triplet(conn_table[e][ii],conn_table[e][jj], stiff_elem[ii][jj]));
                }
                stf = stiff_elem[ii+8][jj+8];;
                if(stf > 0.00001 || stf < -0.00001)
                {
                    //second component
#pragma omp critical
                    stiffness_tripletList.push_back(triplet(conn_table2[e][ii],conn_table2[e][jj], stiff_elem[ii+8][jj+8]));
                }
                stf = stiff_elem[ii+16][jj+16];;
                if(stf > 0.00001 || stf < -0.00001)
                {
                    //third component
#pragma omp critical
                    stiffness_tripletList.push_back(triplet(conn_table3[e][ii],conn_table3[e][jj], stiff_elem[ii+16][jj+16]));
                }
            }
        }

        //Integral 2
        blm_integrate.integral2(load_elem, x1, x2, x3, force);
        for (int ii = 4; ii < 8; ii++)
        {
            load_vector[conn_table[e][ii]] += load_elem[ii]; //first component
            load_vector[conn_table2[e][ii]] += load_elem[ii+8]; //second component
            load_vector[conn_table3[e][ii]] += load_elem[ii+16]; //third component
        }
    }

    #pragma omp parallel for private(x1, x2, x3, stiff_elem, load_elem)
    for(int e = ((N-1)*(N-1)); e < n_elem; e++)
    {
        for(int X = 0; X < 8; X++)
        {
            x1[X] = node_table[conn_table[e][X]][0];
            x2[X] = node_table[conn_table[e][X]][1];
            x3[X] = node_table[conn_table[e][X]][2];
        }

        //Integral 1
        blm_integrate.integral1(stiff_elem, x1, x2, x3, Etensor);
        for (int ii = 0; ii < 8; ii++)
        {
            for (int jj = 4; jj < 8; jj++)
            {
                double stf = stiff_elem[ii][jj];;
                if(stf > 0.00001 || stf < -0.00001)
                {
                    //first component
#pragma omp critical
                    stiffness_tripletList.push_back(triplet(conn_table[e][ii],conn_table[e][jj], stiff_elem[ii][jj]));
                }
                stf = stiff_elem[ii+8][jj+8];;
                if(stf > 0.00001 || stf < -0.00001)
                {
                    //second component
#pragma omp critical
                    stiffness_tripletList.push_back(triplet(conn_table2[e][ii],conn_table2[e][jj], stiff_elem[ii+8][jj+8]));
                }
                stf = stiff_elem[ii+16][jj+16];;
                if(stf > 0.00001 || stf < -0.00001)
                {
                    //third component
#pragma omp critical
                    stiffness_tripletList.push_back(triplet(conn_table3[e][ii],conn_table3[e][jj], stiff_elem[ii+16][jj+16]));
                }
            }
        }

        //Integral 2
        blm_integrate.integral2(load_elem, x1, x2, x3, force);
        for (int ii = 0; ii < 8; ii++)
        {
            load_vector[conn_table[e][ii]] += load_elem[ii]; //first component
            load_vector[conn_table2[e][ii]] += load_elem[ii+8]; //second component
            load_vector[conn_table3[e][ii]] += load_elem[ii+16]; //third component
        }
    }

    cout << "Matrix construction, " << (omp_get_wtime() - start) << endl;
    start = omp_get_wtime();

    //Integral 3 - for boundary conditions
//#pragma omp parallel for
    for(int e = 0; e < ((N-1)*(N-1)); e++)
    {
        //First Component
        stiffness_tripletList.push_back(triplet(conn_table[e][0],conn_table[e][0], 1.0));
        stiffness_tripletList.push_back(triplet(conn_table[e][1],conn_table[e][1], 1.0));
        stiffness_tripletList.push_back(triplet(conn_table[e][2],conn_table[e][2], 1.0));
        stiffness_tripletList.push_back(triplet(conn_table[e][3],conn_table[e][3], 1.0));

        load_vector[conn_table[e][0]] = 0.0;
        load_vector[conn_table[e][1]] = 0.0;
        load_vector[conn_table[e][2]] = 0.0;
        load_vector[conn_table[e][3]] = 0.0;

        //Second Component
        stiffness_tripletList.push_back(triplet(conn_table2[e][0],conn_table2[e][0], 1.0));
        stiffness_tripletList.push_back(triplet(conn_table2[e][1],conn_table2[e][1], 1.0));
        stiffness_tripletList.push_back(triplet(conn_table2[e][2],conn_table2[e][2], 1.0));
        stiffness_tripletList.push_back(triplet(conn_table2[e][3],conn_table2[e][3], 1.0));

        load_vector[conn_table2[e][0]] = 0.0;
        load_vector[conn_table2[e][1]] = 0.0;
        load_vector[conn_table2[e][2]] = 0.0;
        load_vector[conn_table2[e][3]] = 0.0;

        //Third Component
        stiffness_tripletList.push_back(triplet(conn_table3[e][0],conn_table3[e][0], 1.0));
        stiffness_tripletList.push_back(triplet(conn_table3[e][1],conn_table3[e][1], 1.0));
        stiffness_tripletList.push_back(triplet(conn_table3[e][2],conn_table3[e][2], 1.0));
        stiffness_tripletList.push_back(triplet(conn_table3[e][3],conn_table3[e][3], 1.0));

        load_vector[conn_table3[e][0]] = 0.0;
        load_vector[conn_table3[e][1]] = 0.0;
        load_vector[conn_table3[e][2]] = 0.0;
        load_vector[conn_table3[e][3]] = 0.0;
    }

    stiffness_matrix.setFromTriplets(stiffness_tripletList.begin(), stiffness_tripletList.end());

    cout << "Dirichlet BC, " << (omp_get_wtime() - start) << endl;
    start = omp_get_wtime();

    for(int i=0;i<3;i++)
    {
        traction[i]=tractionIn[i];
    }
    vector<double> zvector(3); //direction for the surface boundary condition
    zvector = {0.0,0.0,1.0};

//#pragma omp for
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

    cout << "Neumann BC, " << (omp_get_wtime() - start) << endl;
    cout << "Total Construction, " << (omp_get_wtime() - startTotal) << endl;


    start = omp_get_wtime();
    //Step 4 - Solve the system of equations

    //vector<double> solution_vector(3*n_nodes);
    Eigen::VectorXd solution_vector(3*n_nodes);
    vector<vector<double>> ufull(n_nodes,vector<double>(3));

//    for (auto it = stiffness_tripletList.begin(); it!=stiffness_tripletList.end(); ++it)
//        cout << (*it).row() << " " << (*it).col() << " " << (*it).value() << endl;

//    for (int i = 0; i < load_vector.size(); ++i)
//        cout << load_vector[i] << endl;

    //return 0;

    //CGSolver solver;
    Eigen::ConjugateGradient<SpMat, Eigen::Lower|Eigen::Upper> solver;
    solver.compute(stiffness_matrix);
    solver.setMaxIterations(59);
    solution_vector = solver.solve(load_vector);
//    std::cout << "#iterations:     " << solver.iterations() << std::endl;
//    std::cout << "estimated error: " << solver.error()      << std::endl;

    cout << "CG Solver, " << (omp_get_wtime() - start) << endl;

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
    cout << "Total Runtime, " << (omp_get_wtime() - startTotal) << endl;
};
