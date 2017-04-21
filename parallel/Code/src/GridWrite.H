#ifndef _GRIDWRITE_H_
#define _GRIDWRITE_H_

#include <cstdio>   
#include <cmath>    
#include <cassert>
#include <cstring> 
#include <vector> 
#include <string>

using namespace std;
/// Writes VTU file. Takes node table (3 x # of nodes), connectivity table (8 x # of elements), solution table (3 x # of nodes), filename
void GridWrite(const vector<vector<double>> &a_nodes,const vector<vector<int>> &a_conn,const vector<vector<double>> &a_solution,const char* a_filename);
#endif
