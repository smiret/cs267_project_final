#ifndef _MESHER_H_
#define _MESHER_H_
#include <iostream>
#include <cassert>
#include <cmath>
#include <math.h>
#include <vector>


using namespace std;

class Mesher

{
public:
  ///Default Constructor
  Mesher();
  ///Constructor: just require int number of nodes in each direction
  Mesher(const int& a_N);
  ///Mesh Creator
  void createmesh();
  ///Get node table (vector of vector of double)
  vector<vector<double>> getNode();
  ///Get connectivity table (vector of vector of double)
  vector<vector<int>> getConn();
  
  //Member Data
private:
  vector<vector<double>> m_nodes;
  vector<vector<int>> m_conn;
  int m_N;
  int m_nNodes;
  int m_nElem;
};





#endif
