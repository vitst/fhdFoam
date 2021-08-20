#!/usr/bin/env python3

import numpy as np
import os.path

total_time = 100000
write_each = 10
print_each = 1000

Npart = 100000
Ncels = 100

# resulting file name
rfn_0 = "density1D"
count = 0
rfn = rfn_0
while os.path.isfile(rfn):
    count+=1
    rfn = "{}_{}".format(rfn_0, count)

with open(rfn, 'a') as f:
    f.write("Total time: {}   Npart: {}   Ncells: {}\n".format(total_time, Npart, Ncels))

# each particle coord in 1D
x = np.zeros(Npart, dtype=int)

# initial condition
for i,coord in enumerate(x):
    x[i]=i%Ncels

#number_per_cell, bin_edges = np.histogram(x)
#print("Initial distribution")
#print(number_per_cell)

# start time
for t in range(total_time):

    if t%print_each==0:
        print('time: {}    final: {}\r'.format(t, total_time), end='')

    # generate random array
    rnd_ar = np.random.randint(2, size=Npart)
    rnd_ar = 2 * rnd_ar - 1

    # move particles
    x = ( x + rnd_ar )%Ncels
    # apply periodic bc
    #x[x==-1]=Ncels-1
    #x[x==Ncels]=0

    # write res to file
    if t%write_each==0:
        number_per_cell, bin_edges = np.histogram(x, bins=Ncels)
        with open(rfn, 'a') as f:
            f.write("{}   ".format(t))
            for npc in number_per_cell:
                f.write("{:.5f}   ".format(npc/float(Npart)))
            f.write("\n")


