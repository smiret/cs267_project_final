#ifndef _BLMINTEGRALS_H_
#define _BLMINTEGRALS_H_
#include <iostream>
#include <cassert>
#include <cmath>
#include <math.h>
#include <vector>
#include "femfunctions.H"




using namespace std;

class blmintegrals

{
public:
  /** Perform Integral 1 for balance of linear momentum. First argument is matrix to return into.
  Next 3 are input vectors from node table. Last is input matrix of elasticity tensor.
  **/
  void integral1(double (&a_matrix)[24][24],
		 double (&a_x1)[8],
         double (&a_x2)[8],
         double (&a_x3)[8],
		 const vector<vector<double> > Etensor);
  /** Perform Integral 2 for balance of linear momentum. First argument is vector to return into.
  Next 3 are input vectors from node table. Last is input vector of force components.
  **/
  void integral2(double (&a_vector)[24],
         double (&a_x1)[8],
         double (&a_x2)[8],
         double (&a_x3)[8],
		 vector<double>& a_force);
  /** Perform surface Integral 3 for balance of linear momentum on x surface. First argument is vector to return into.
  Next 3 are input vectors from node table. Last is input vector of traction components.
  **/
  void integral3xx(double (&a_vector)[24],
       double (&a_x1)[8],
       double (&a_x2)[8],
       double (&a_x3)[8],
		 vector<double>& a_zvector,
		 vector<double>& a_traction);
  /** Perform surface Integral 3 for balance of linear momentum on y surface. First argument is vector to return into.
  Next 3 are input vectors from node table. Last is input vector of traction components.
  **/
  void integral3yy(double (&a_vector)[24],
       double (&a_x1)[8],
       double (&a_x2)[8],
       double (&a_x3)[8],
		 vector<double>& a_zvector,
		 vector<double>& a_traction);
  /** Perform surface Integral 3 for balance of linear momentum on z surface. First argument is vector to return into.
  Next 3 are input vectors from node table. Last is input vector of traction components.
  **/
  void integral3zz(double (&a_vector)[24],
       double (&a_x1)[8],
       double (&a_x2)[8],
       double (&a_x3)[8],
		 vector<double>& a_zvector,
		 vector<double>& a_traction);


  







};





#endif

