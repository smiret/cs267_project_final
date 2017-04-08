#include <iostream>
#include <cassert>
#include <cmath>
#include <math.h>
#include <vector>
#include "hsbounds.H"
using namespace std;

void hsbounds::HS_bounds(double& k1,
		 double& k2,
		 double& u1,
		 double& u2,
		 double& v2,
		 double& kstar,
		 double& ustar)
{

  double k_hs_u, k_hs_l, u_hs_u, u_hs_l;

  k_hs_l = k1 + v2/((1/(k2-k1)) + (3*(1-v2)/(3*k1+4*u1)));
  k_hs_u = k2 + (1-v2)/((1/(k1-k2)) + (3*v2/(3*k2+4*u2)));

  u_hs_l = u1 + v2/((1/(u2-u1)) + (6*(1-v2)*(k1+2*u1))/(5*u1*(3*k1+4*u1)));
  u_hs_u = u2 + (1-v2)/((1/(u1-u2)) + (6*v2*(k2+2*u2))/(5*u2*(3*k2+4*u2)));

  kstar = 0.5*(k_hs_l + k_hs_u);
  ustar = 0.5*(u_hs_l + u_hs_u);


};
