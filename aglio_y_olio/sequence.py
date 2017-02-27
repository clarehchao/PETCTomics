# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 19:39:27 2016

@author: rharnish
"""

#%%

import collections
#http://stackoverflow.com/questions/8664708/in-python-how-does-one-efficiently-find-the-largest-consecutive-set-of-numbers
def longest_consecutive_sequence(sequence):
    # map starting values to largest ending value so far
    map = collections.OrderedDict()

    for i in sequence:
        found = False
        for k, v in map.iteritems():
            if i == v:
                map[k] += 1
                found = True

        if not found and i not in map:
            map[i] = i + 1

    return xrange(*max(map.iteritems(), key=lambda i: i[1] - i[0]))