#!/usr/bin/env python
###############################################################################
##                         III-V Calculator Database                         ##
###############################################################################
"""Calculates the lattice constant 'a' for the entire parameter space which
satisfies the III-V semidonductor composition equation found on the III-V Calc
information page. <http://ahrenkiel.sdsmt.edu/III_V_Calc/info/>

WARNING: For small fractional_resolution values the code takes a long time to 
         run and creates a large output file. 

Returns:
    Creates a compressed .npz file containing the composition information and
    the calculated lattice constants.
"""
import argparse
import numpy


parser = argparse.ArgumentParser(description="Calculates the lattice constant 'a' for the entire parameter space which satisfies the III-V semidonductor composition equation found on the III-V Calc information page. <http://ahrenkiel.sdsmt.edu/III_V_Calc/info/>")
parser.add_argument("resolution", type=int, help="The resolution of the step size in composition in percent where (1 = 1 percent).")
parser.add_argument("output_file", type=str, help="Name of the compressed npz file where the array of composition and corresponding lattice constant is saved.")
args = parser.parse_args()

#define all reference lattice constants
a_AlP = 5.4510
a_GaP = 5.4505
a_InP = 5.8686
a_AlAs = 5.6605
a_GaAs = 5.6533
a_InAs = 6.0584
a_AlSb = 6.1355
a_GaSb = 6.0950
a_InSb = 6.4794

#This sets the resolution of the steps in compsition in percent (1 = 1%)
fraction_resolution = args.resolution

#This code counts the number of iterations in order to pre-allocate a numpy array
i = 0
for y_P in numpy.arange(0, 100 + fraction_resolution, fraction_resolution):
    for y_As in numpy.arange(0, 100 + fraction_resolution - y_P, fraction_resolution):
        for x_Al in numpy.arange(0, 100 + fraction_resolution, fraction_resolution):
            for x_Ga in numpy.arange(0, 100 + fraction_resolution - x_Al, fraction_resolution):
                i+=1

lst = numpy.zeros((i, 7), dtype=numpy.float32) #Pre allocate a numpy array of float32, which is enough precision, size is the number of iterations the below code does

i = 0
for y_P in numpy.arange(0, 100 + fraction_resolution, fraction_resolution):
    for y_As in numpy.arange(0, 100 + fraction_resolution - y_P, fraction_resolution):
        for x_Al in numpy.arange(0, 100 + fraction_resolution, fraction_resolution):
            for x_Ga in numpy.arange(0, 100 + fraction_resolution - x_Al, fraction_resolution):
                y_Sb = 100 - y_P - y_As
                x_In = 100 - x_Al - x_Ga
                a = (y_P*(x_Al*a_AlP + x_Ga*a_GaP + x_In*a_InP) + \
                y_As*(x_Al*a_AlAs + x_Ga*a_GaAs + x_In*a_InAs) + \
                y_Sb*(x_Al*a_AlSb + x_Ga*a_GaSb + x_In*a_InSb))/10000.0
                lst[i] = [x_Al/100.0, x_Ga/100.0, x_In/100.0, y_P/100.0, y_As/100.0, y_Sb/100.0, a]
                i+=1


numpy.savez_compressed(args.output_file, lst)