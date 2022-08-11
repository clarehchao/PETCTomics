# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 22:30:04 2015

@author: rharnish
"""

import subprocess

def run_command(command, verbose=True):
    process = subprocess.Popen(command, stdout=subprocess.PIPE)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            if verbose == True:
                print output.strip()
                
    rc = process.poll()
    return rc