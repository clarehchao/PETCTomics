# -*- coding: utf-8 -*-
"""
Created on Sun Nov  1 22:34:00 2015

@author: rharnish
"""

import os
import socket

#%%

class MachineSpecificPaths:
    
    def __init__(self,*args):

        self._db_dir       = None
        self._test_dir     = None
        self._data_dir     = None
        self._incoming_dir = None
        
        self.determine_paths_for_machine()
    
    # ---------------------------------------------------------------------    
    # Properties
    # --------------------------------------------------------------------- 
    @property
    def db_dir(self):
        return self._db_dir  
        
    @property
    def test_dir(self):
        return self._test_dir  
    
    @property
    def data_dir(self):
        return self._data_dir 
    
    @property
    def incoming_dir(self):
        return self._incoming_dir 
        
    def determine_paths_for_machine(self):
        
        hostName = socket.gethostname()
        print 'host: {}'.format(hostName)

        # Roy's Mac
        RoysMac = 'cbl-mbp-3369'
        if hostName == RoysMac:
    
            print 'RoysMac'
            self._test_dir     = os.path.join('/Users/rharnish/Projects/Franc/PETCT','TEST')
            self._db_dir       = '/data/db'
            self._data_dir     = '/Volumes/Untitled 1/data'
            self._incoming_dir = '/Volumes/Untitled 1/data'
            
        # Lassen
#        Lassen = 'lassen.radiology.ucsf.edu'
#        if hostName == Lassen:
        if hostName.endswith('radiology.ucsf.edu'):
            
            print '{} -- using RRCS linux paths'.format(hostName)
            self._test_dir     = os.path.join('/data/francgrp1/PETCT','TEST')
            self._db_dir       = '/data/francgrp1/PETCT/DB'
#            self._data_dir     = '/data/francgrp1/incoming_breast_v2_processing'
            self._data_dir     = '/data/francgrp1/breast_radiomics/her2/PETCT'
            self._incoming_dir = '/data/francgrp1/incoming'

#%%

if __name__ == "__main__":
    paths = MachineSpecificPaths()
    