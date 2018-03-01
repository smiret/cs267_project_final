# cs267_project_final
Parallelized C++ Finite Element Solver

This code includes a C++ based finite element solver for the partial differential equation governing the balance of linear momemtum (BLM).
The BLM equation gives information about the mechanical properties and behavior of a given material system.

In this framework, the BLM equation is solved using a three-dimensional finite element code written in C++. 
The framework includes the setup of the finite element method with all mathematical manipulations to the equations required.
The resulting matrix system of equations is solved using a Conjugate Gradient Solver. 
The solution being written to an xml file that can be opened by a graphics software, such as Paraview.

The repository includes a serial implementation of the code, as well as a parallelized version.
In the parallelized version, some operations are performed in parallel to enhance the wall-clock time of the computation.

