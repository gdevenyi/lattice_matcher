#!/usr/bin/env python
###############################################################################
##                         III-V Calculator Database                         ##
###############################################################################
"""Calculates the lattice constant 'a' for the entire parameter space which
satisfies the III-V semidonductor composition equation found on the III-V Calc
information page. <http://ahrenkiel.sdsmt.edu/III_V_Calc/info/>

Returns:
    Creates a compressed .npz file containing the composition information and
    the calculated lattice constants.
"""

import numpy

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

lst = numpy.zeros((26532801, 7), dtype=numpy.float32) #Pre allocate a numpy array of float32, which is enough precision, size is the number of iterations the below code does
i = 0
for y_P in numpy.arange(0, 101, 1):
    for y_As in numpy.arange(0, 101 - y_P, 1):
        for x_Al in numpy.arange(0, 101, 1):
            for x_Ga in numpy.arange(0, 101 - x_Al, 1):
                y_Sb = 100 - y_P - y_As
                x_In = 100 - x_Al - x_Ga

                a = (y_P*(x_Al*a_AlP + x_Ga*a_GaP + x_In*a_InP) + \
                y_As*(x_Al*a_AlAs + x_Ga*a_GaAs + x_In*a_InAs) + \
                y_Sb*(x_Al*a_AlSb + x_Ga*a_GaSb + x_In*a_InSb))/10000.0
		lst[i] = [x_Al/100.0, x_Ga/100.0, x_In/100.0, y_P/100.0, y_As/100.0, y_Sb/100.0, a]
		i = i+1
numpy.savez_compressed("lattice_constants", lst)
