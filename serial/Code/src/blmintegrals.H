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
  void integral1(vector<vector<double> >& a_matrix,
		 vector<double>& a_x1,
		 vector<double>& a_x2,
		 vector<double>& a_x3,
		 vector<vector<double> >& Etensor);
  /** Perform Integral 2 for balance of linear momentum. First argument is vector to return into.
  Next 3 are input vectors from node table. Last is input vector of force components.
  **/
  void integral2(vector<double>& a_vector,
		 vector<double>& a_x1,
		 vector<double>& a_x2,
		 vector<double>& a_x3,
		 vector<double>& a_force);
  /** Perform surface Integral 3 for balance of linear momentum on x surface. First argument is vector to return into.
  Next 3 are input vectors from node table. Last is input vector of traction components.
  **/
  void integral3xx(vector<double>& a_vector,
		 vector<double>& a_x1,
		 vector<double>& a_x2,
		 vector<double>& a_x3,
		 vector<double>& a_zvector,
		 vector<double>& a_traction);
  /** Perform surface Integral 3 for balance of linear momentum on y surface. First argument is vector to return into.
  Next 3 are input vectors from node table. Last is input vector of traction components.
  **/
  void integral3yy(vector<double>& a_vector,
		 vector<double>& a_x1,
		 vector<double>& a_x2,
		 vector<double>& a_x3,
		 vector<double>& a_zvector,
		 vector<double>& a_traction);
  /** Perform surface Integral 3 for balance of linear momentum on z surface. First argument is vector to return into.
  Next 3 are input vectors from node table. Last is input vector of traction components.
  **/
  void integral3zz(vector<double>& a_vector,
		 vector<double>& a_x1,
		 vector<double>& a_x2,
		 vector<double>& a_x3,
		 vector<double>& a_zvector,
		 vector<double>& a_traction);


  







};





#endif

