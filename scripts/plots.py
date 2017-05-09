#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  1 19:26:48 2017

@author: lbluque
"""

import matplotlib.pyplot as plt
import numpy as np

length = np.array([10,	14,	18,	22,	26,	30,	40,	50,	65,	85])
serial = np.array([[0.588945,	1.53378,	3.25455,	5.60245,	10.0611,	15.6578,	37.6734,	73.2446,	163.665,	370.796],
[0.043577,	0.20035,	0.786935,	2.02905,	5.14299,	10.6629,	45.8899,	141.107,	531.008,	2065.17],
[0.00103,	0.00174,	0.00319,	0.004693,	0.006984,	0.009264,	0.016735,	0.027568,	0.047738,	0.077398],
[0.634364,	1.73789,	4.04892,	7.64478,	15.2205,	26.3566,	83.6304,	214.417,	694.884,	2436.33]])

parallel = np.array([[0.241697,	0.511107,	1.09222,	1.8114,	3.05505,	4.7293,	11.3565,	22.1571,	50.4189,	112.529],
[0.0127487,	0.067322,	0.194912,	0.399103,	0.815835,	1.57779,	5.81777,	18.1444,	66.9255,	255.76],
[0.00442908,	0.0128053,	0.0185936,	0.0183287,	0.0358848,	0.0375735,	0.0617036,	0.203903,	0.191541,	0.396372],
[0.259623,	0.592083,	1.3076,	2.23252,	3.91178,	6.35612,	17.2598,	40.5581,	117.686,	368.904]])



plt.figure()
plt.plot(length, serial[3], '-v', length, parallel[3], '-o')
plt.xlabel(r'System length (arb)')
plt.ylabel(r'Total stiffness matrix construction time (s)')
plt.legend(['Serial', 'Parallel'])
plt.tight_layout()
plt.savefig('total_time')

plt.figure()
plt.loglog(length, serial[3], '-v', length, parallel[3], '-o')
plt.xlabel(r'System length (arb)')
plt.ylabel(r'Total stiffness matrix construction time (s)')
plt.legend(['Serial', 'Parallel'])
plt.tight_layout()
plt.savefig('total_time')     

plt.figure()
plt.plot(length, serial[2], '-v', length, parallel[2], '-o')
plt.xlabel(r'System length (arb)')
plt.ylabel(r'Neumann boundary condition setup time (s)')
plt.legend(['Serial', 'Parallel'])
plt.tight_layout()
plt.savefig('neumann_time')

plt.figure()
plt.plot(length, serial[1], '-v', length, parallel[1], '-o')
plt.xlabel(r'System length (arb)')
plt.ylabel(r'Dirichlet boundary condition setup  (s)')
plt.legend(['Serial', 'Parallel'])
plt.tight_layout()
plt.savefig('dirichlet_time')

plt.figure()
plt.plot(length, serial[0], '-v', length, parallel[0], '-o')
plt.xlabel(r'System length (arb)')
plt.ylabel(r'Stiffness matrix construction time (s)')
plt.legend(['Serial', 'Parallel'])
plt.tight_layout()
plt.savefig('mat_time')


size = np.array([10,	14, 18,	22,	26])
cg_serial = np.array([11.057,	121.488,	593.162,	2671.36,	8449.53])
cg_parallel = np.array([8.67854,	82.0936,	464.75,	2056.06,	6075.91])
eigen = np.array([])

plt.figure()
plt.plot(size, cg_serial, '-v', size, cg_parallel, '-o')
plt.xlabel(r'System length (arb)')
plt.ylabel(r'Total CG Solver time (s)')
plt.legend(['Serial', 'Parallel'])
plt.tight_layout()
plt.savefig('total_cg_time')

plt.figure()
plt.loglog(size, cg_serial, '-v', size, cg_parallel, '-o')
plt.xlabel(r'System length (arb)')
plt.ylabel(r'Total CG Solver time (s)')
plt.legend(['Serial', 'Parallel'])
plt.tight_layout()
plt.savefig('total_cg_time_log')