#!/usr/bin/env python
###############################################################################
##                         III-V Calculator Database                         ##
###############################################################################
"""Calculates the lattice constant 'a' for the entire parameter space which
satisfies the III-V semidonductor composition equation found on the III-V Calc
information page. <http://ahrenkiel.sdsmt.edu/III_V_Calc/info/>

Returns:
    Creates a compressed pickle file containing the composition information and
    the calculated lattice constants.
"""

import numpy
import pickle
import gzip

#define all reference lattice constants
a_AlP = 5.4510
a_GaP = 5.4505
a_InP = 5.8686
a_AlAs = 5.6605
a_GaAs = 5.6533
a_InAs = 6.0584
a_AlSb = 6.1355
a_GaSb = 6.095
a_InSb = 6.4794

lst = []
for y_P in numpy.arange(0, 110, 1):
    for y_As in numpy.arange(0, 110 - y_P, 1):
        for x_Al in numpy.arange(0, 110, 1):
            for x_Ga in numpy.arange(0, 110 - x_Al, 1):
                y_Sb = 100 - y_P - y_As
                x_In = 100 - x_Al - x_Ga

                a = (y_P*(x_Al*a_AlP + x_Ga*a_GaP + x_In*a_InP) + \
                y_As*(x_Al*a_AlAs + x_Ga*a_GaAs + x_In*a_InAs) + \
                y_Sb*(x_Al*a_AlSb + x_Ga*a_GaSb + x_In*a_InSb))/10000.0

                entry = ['{0:.{1}f}'.format(x_Al/100.0, 2), '{0:.{1}f}'.format(x_Ga/100.0,2) , '{0:.{1}f}'.format(x_In/100.0,2), '{0:.{1}f}'.format(y_P/100.0,2) , '{0:.{1}f}'.format(y_As/100.0,2), '{0:.{1}f}'.format(y_Sb/100.0,2), numpy.around(a,4)]
                lst.append(entry)

outputfile = gzip.open('lattice_constants.pklz','wb')
pickle.dump(lst,outputfile)
outputfile.close()
