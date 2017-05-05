#include <cassert>
#include <cmath>
#include <math.h>
#include <vector>
#include "femfunctions.H"




using namespace std;

class hsbounds

{
public:
  /** Determines effective bulk (kstar) and shear (ustar) moduli from input bulk (k1,k2) 
  and shear (u1,u2) moduli and proportion (v2)**/
  void HS_bounds(double& k1,
		 double& k2,
		 double& u1,
		 double& u2,
		 double& v2,
		 double& kstar,
		 double& ustar);
};
