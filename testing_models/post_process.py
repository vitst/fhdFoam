#!/usr/bin/env python3

import numpy as np

import sys

res_file = sys.argv[1]

with open(res_file, 'r') as f:
    f.readline()
    times = []
    for line in f:
        times.append(line.split()[0])
        dens = line.split()[1:]
        dens = np.asarray(dens, dtype=float)

        dens2 = np.power(dens, 2)
        
        fl = np.average(dens2) - np.average(dens)**2 

        print("{}   {}".format(line.split()[0], fl))
