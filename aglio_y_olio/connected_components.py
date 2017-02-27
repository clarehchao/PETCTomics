# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 19:22:54 2016

@author: rharnish
"""

#%%

import numpy as np
from skimage.measure import label, regionprops

#%%

def get_connected_components_larger_than_min_vol(MASK,min_vol=0):
    dim      = len(MASK.shape)
    LABEL    = label(MASK,connectivity=dim,background=0)
    regions  = regionprops(LABEL)
    LARGECC  = np.zeros_like(LABEL)     
    regions_to_include = []
    for r in regions:
        if r.area > min_vol:
            regions_to_include.append(r)  
    for r in regions_to_include:
        COORDS = r.coords
        for C in COORDS:
            if dim == 3:
                LARGECC[C[0],C[1],C[2]] = 1
            if dim == 2:
                LARGECC[C[0],C[1]] = 1
    return LARGECC
    
def get_largest_connected_component(MASK):
    dim = len(MASK.shape)
    LABEL    = label(MASK,connectivity=3,background=0)
    regions  = regionprops(LABEL)
    LARGECC  = np.zeros_like(LABEL)     
    max_area = -1
    for r in regions:
        if r.area > max_area:
            max_area = r.area
            max_r = r 
    COORDS = max_r.coords
    for C in COORDS:
        if dim == 3:
            LARGECC[C[0],C[1],C[2]] = 1
        if dim == 2:
            LARGECC[C[0],C[1]] = 1
    return LARGECC