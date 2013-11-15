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
a_AlP = 5.451
a_GaP = 5.4505
a_InP = 5.8686
a_AlAs = 5.6605
a_GaAs = 5.6533
a_InAs = 6.0584
a_AlSb = 6.1355
a_GaSb = 6.095
a_InSb = 6.4794

lst = []
for y_P in numpy.arange(0., 1.1, 0.1):
    for y_As in numpy.arange(0., 1.1 - numpy.around(y_P, 1), 0.1):
        for x_Al in numpy.arange(0., 1.1, 0.1):
            for x_Ga in numpy.arange(0., 1.1 - numpy.around(x_Al, 1), 0.1):
                y_Sb = numpy.around(1. - y_P - y_As, 1)
                x_In = numpy.around(1. - x_Al - x_Ga, 1)
                if (x_Al < 0.) or (x_Ga < 0.) or (x_In < 0.) or (y_P < 0.) or (y_As < 0.) or (y_Sb < 0.):
                    pass
                else:
                    a = numpy.around((y_P*(x_Al*a_AlP + x_Ga*a_GaP + x_In*a_InP) + \
                        y_As*(x_Al*a_AlAs + x_Ga*a_GaAs + x_In*a_InAs) + \
                        y_Sb*(x_Al*a_AlSb + x_Ga*a_GaSb + x_In*a_InSb)), 6)
                    entry = [x_Al, x_Ga, x_In, y_P, y_As, y_Sb, a]
                lst.append(entry)
numpy.savez_compressed("lattice_constants", lst)