#!/usr/bin/env python
###############################################################################
##                         III-V Calculator Database                         ##
###############################################################################
"""Calculates the lattice constant 'a' for the entire parameter space which
satisfies the III-V semidonductor composition equation found on the III-V Calc
information page. <http://ahrenkiel.sdsmt.edu/III_V_Calc/info/>
"""

import numpy

#define all reference lattice constants
a_AlP = 5.451
a_GaP = 5.4505
a_InP = 5.8686
a_AlAs = 5.6605
a_GaAs = 5.6533
a_InAs = 6.0584
a_AlSb = 6.1355
a_GaSb = 6.095
a_InSb = 6.4794

results_file = open("lattice_constant_list.txt", 'w')
results_file.write("#y_P\ty_As\ty_Sb\tx_Al\tx_Ga\tx_In\ta\n")

for y_P in numpy.arange(0, 1.1, 0.1):
    y = 1. - y_P
    for y_As in numpy.arange(0, y+0.1, 0.1):
        y_Sb = 1. - y_P - y_As
        for x_Al in numpy.arange(0, 1.1, 0.1):
            x = 1. - x_Al
            for x_Ga in numpy.arange(0, x+0.1, 0.1):
                x_In = 1. - x_Al - x_Ga
                a = y_P*(x_Al*a_AlP + x_Ga*a_GaP + x_In*a_InP) + \
                    y_As*(x_Al*a_AlAs + x_Ga*a_GaAs + x_In*a_InAs) + \
                    y_Sb*(x_Al*a_AlSb + x_Ga*a_GaSb + x_In*a_InSb)
                results_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(y_P, 
                                   y_As, y_Sb, x_Al, x_Ga, x_In, a))
                
results_file.close()