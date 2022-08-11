# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 19:59:14 2016

@author: rharnish
"""

#%%

import numpy as np

#%%

def circular_kernel(radius=1):
    print radius
    r = radius
    kw = 2*r + 1
    c  = (kw+1) / 2
    kern = np.ones([kw,kw])
    for i in range(kw):
        for j in range(kw):
                if np.sqrt( (i-c+1)**2 + (j-c+1)**2 ) > r:
                    kern[i][j] = 0
    return kern