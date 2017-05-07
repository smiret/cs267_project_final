#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 29 20:46:13 2017

@author: lbluque
"""

import numpy as np
import scipy.sparse as sp
from itertools import takewhile
from matplotlib.pyplot import spy

filename = 'datos.txt'

with open(filename, 'r') as f:
    f_string = f.readlines()
    data = [list(map(float, line.split(','))) for line in takewhile(lambda x: 'Indices:' not in x, f_string[1:])]
    indices = [list(map(int, line.split(','))) for line in takewhile(lambda x: 'M:' not in x,  f_string[len(data) + 2:])]


indptr = [len(col) for col in indices]
indptr = [0] + list(np.cumsum(indptr))

data, indices = sum(data, []), sum(indices, [])

stiffmat = sp.csr_matrix((data, indices, indptr))

spy(stiffmat)    