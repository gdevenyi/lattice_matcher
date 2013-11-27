#!/usr/bin/env python
###############################################################################
##                           Composition Calculator                          ##
###############################################################################

"""
Software to calculate a range of materials which may be epitaxially grown on a 
given substrate. 

CURRENT PROGESS: Reads in both a substrate database and a reference III-V 
semiconductor lattice constant .npz file, compares lattice constants and writes 
matches to a tab delimited txt file.
"""

import numpy #includes numpy.sqrt()
import argparse #for command line implementation

parser = argparse.ArgumentParser(description="Software for calculating a range of material composition for an epitaxially grown film on a given substrate.")
parser.add_argument("substrate", type=str, help="Tab-delimited txt file with substrate material data.")
parser.add_argument("lattice_constant_database", type=str, help="Compressed npz file containing composition and lattice constant information.")
args = parser.parse_args()

def check_substrate_file(sub_file, lattice_const_file, tolerance_percentage, output_file):
    """Calls functions for calculations based on the information obtained from
       a supplied database file.
    
    Args:
        sub_file: database file with substrate material information
        lattice_const_file: npz file containing lattice constants and composition
                            information
        tolerance_percentage: tolerance percentage of mismatch error represented
                              as a decimal value
        output_file: tab delimited .txt file where results are written
        
    Returns:
        A tab delimited .txt file with the maximum and minimum values related
        to the lattice constants of the substrate material. 
    """
    output_file.write("#Film Composition\tFilm Symmetry\tFlim a\tSubstrate\tSymmetry\n")
    for i, l in enumerate(sub_file):
        if sub_file[i][1] == "C":
            cubic_sub(sub_file[i][0], sub_file[i][1], sub_file[i][2], lattice_const_file, tolerance_percentage, output_file)
        if sub_file[i][1] == "T":
            tetragonal_sub(sub_file[i][0], sub_file[i][1], sub_file[i][2], sub_file[i][3], lattice_const_file, tolerance_percentage, output_file)
        if sub_file[i][1] == "H":
            hexagonal_sub(sub_file[i][0], sub_file[i][1], sub_file[i][2], sub_file[i][3], lattice_const_file, tolerance_percentage, output_file)


def cubic_sub(sub_comp, sub_sym, sub_a_val, lattice_consts, tol, result_file):
    """Calculates max/min lattice constant values for a cubic substrate.
    
    Args:
        sub_comp: substrate composition
        sub_sym: substrate symmetry
        sub_a_val: value of lattice constant 'a'
        lattice_consts: npz file of composition and lattice constants. Lattice 
                        constants must be the right-most entry in each line of 
                        the array
        tol: tolerance percentage of mismatch error represented as a decimal
        result_file: a tab delimited .txt file with the calculated maximum and
                      minimum lattice constant values for a specified tolerance
                      value.
        
    Returns:
        Writes a new line to the result_file which contains the results of the
        lattice constant comparision.    
    """
    result_string = ""    
    good_lattice_vals = (lattice_consts[:,-1] > (1. - tol)*sub_a_val) & (lattice_consts[:,-1] < (1. + tol)*sub_a_val)
    for i, line in enumerate(good_lattice_vals):
        if line:
            film_comp = cust_print(lattice_consts[i][0], lattice_consts[i][1], lattice_consts[i][2], 
                                   lattice_consts[i][3], lattice_consts[i][4], lattice_consts[i][5])
            result_string += "{}\t{}\t{}\t{}\t{}\n".format(film_comp, "C", round(lattice_consts[i][6], 4), sub_comp, sub_sym)
    good_lattice_vals45 = (lattice_consts[:,-1] > (1. - tol)*sub_a_val*numpy.sqrt(2.0)) & (lattice_consts[:,-1] < (1. + tol)*sub_a_val*numpy.sqrt(2.0))
    for i, line in enumerate(good_lattice_vals45):
        if line:
            film_comp = cust_print(lattice_consts[i][0], lattice_consts[i][1], lattice_consts[i][2], 
                                   lattice_consts[i][3], lattice_consts[i][4], lattice_consts[i][5])
            result_string += "{}\t{}\t{}\t{}\t{}\n".format(film_comp, "C (45deg)", round(lattice_consts[i][6], 4), sub_comp, sub_sym)
    result_file.write(result_string)                           
def tetragonal_sub(sub_comp, sub_sym, sub_a_val, sub_c_val, lattice_consts, tol, result_file):
    """Calculates max/min lattice constant values for a tetragonal substrate.
    
    Args:
        sub_comp: substrate composition
        sub_sym: substrate symmetry
        sub_a_val: value of lattice constant 'a'
        sub_c_val: value of lattice constant 'c'
        lattice_consts: npz file of composition and lattice constants. Lattice 
                        constants must be the right-most entry in each line of 
                        the array
        tol: tolerance percentage of mismatch error represented as a decimal
        result_file: a tab delimited .txt file with the calculated maximum and
                      minimum lattice constant values for a specified tolerance
                      value.
    Returns:
        Writes a new line to the result_file which contains the results of the
        lattice constant comparision.    
    """
    result_string = ""
    good_lattice_vals = (lattice_consts[:,-1] > (1. - tol)*sub_a_val) & (lattice_consts[:,-1] < (1. + tol)*sub_a_val)
    for i, line in enumerate(good_lattice_vals):
        if line:
            film_comp = cust_print(lattice_consts[i][0], lattice_consts[i][1], lattice_consts[i][2], 
                                   lattice_consts[i][3], lattice_consts[i][4], lattice_consts[i][5])
            result_string += "{}\t{}\t{}\t{}\t{}\n".format(film_comp, "C", round(lattice_consts[i][6], 4), sub_comp, sub_sym)
    good_lattice_vals45 = (lattice_consts[:,-1] > (1. - tol)*sub_a_val*numpy.sqrt(2.0)) & (lattice_consts[:,-1] < (1. + tol)*sub_a_val*numpy.sqrt(2.0))
    for i, line in enumerate(good_lattice_vals45):
        if line:
            film_comp = cust_print(lattice_consts[i][0], lattice_consts[i][1], lattice_consts[i][2], 
                                   lattice_consts[i][3], lattice_consts[i][4], lattice_consts[i][5])
            result_string += "{}\t{}\t{}\t{}\t{}\n".format(film_comp, "C (45deg)", round(lattice_consts[i][6], 4), sub_comp, sub_sym)
    result_file.write(result_string)            

def hexagonal_sub(sub_comp, sub_sym, sub_a_val, sub_c_val, lattice_consts, tol, result_file):
    """Calculates max/min lattice constant values for a hexagonal substrate.
    
    Args:
        sub_comp: substrate composition
        sub_sym: substrate symmetry
        sub_a_val: value of lattice constant 'a'
        sub_c_val: value of lattice constant 'c'
        lattice_consts: npz file of composition and lattice constants. Lattice 
                        constants must be the right-most entry in each line of 
                        the array
        tol: tolerance percentage of mismatch error represented as a decimal
        result_file: a tab delimited .txt file with the calculated maximum and
                      minimum lattice constant values for a specified tolerance
                      value.
    Returns:
        Writes a new line to the result_file which contains the results of the
        lattice constant comparision.    
    """
    result_string = ""
    good_lattice_vals = (lattice_consts[:,-1] > (1. - tol)*sub_a_val*numpy.sqrt(2.0)) & (lattice_consts[:,-1] < (1. + tol)*sub_a_val*numpy.sqrt(2.0))
    for i, line in enumerate(good_lattice_vals):
        if line:
            film_comp = cust_print(lattice_consts[i][0], lattice_consts[i][1], lattice_consts[i][2], 
                                   lattice_consts[i][3], lattice_consts[i][4], lattice_consts[i][5])
            result_string += "{}\t{}\t{}\t{}\t{}\n".format(film_comp, "C (111)", round(lattice_consts[i][6], 4), sub_comp, sub_sym)
    result_file.write(result_string)
    
def cust_print(x_Al, x_Ga, x_In, y_P, y_As, y_Sb):
    """Custom print function to simplify output of composition. Input arguments
    agree with the following rules:
            x_Al + x_Ga + x_In = 1
            y_P  + y_As + y_Sb = 1
            
    Args:
        x_Al: amount of Al in material 
        x_Ga: amount of Ga in material
        x_In: amount of In in material
        y_P:  amount of P in material
        y_As: amount of As in material
        y_Sb: amount of Sb in material
        
    Returns:
        composition: a string representing the composition formula 
    """
    composition = ""
    composition += "Al"+'{0:.{1}f}'.format(x_Al, 2)
    composition += "Ga"+'{0:.{1}f}'.format(x_Ga, 2)
    composition += "In"+'{0:.{1}f}'.format(x_In, 2)
    composition += "P"+'{0:.{1}f}'.format(y_P, 2)
    composition += "As"+'{0:.{1}f}'.format(y_As, 2)
    composition += "Sb"+'{0:.{1}f}'.format(y_Sb, 2)
    return composition
    

if __name__ == "__main__":
    # create a label for the matches file.
    results_file_label = "composition_matches_for_" + args.substrate[:-4] + ".txt"
    results_file = open(results_file_label, "w")
    substrate_file = numpy.genfromtxt(args.substrate, names=True, comments='#', dtype=None)
    # Tolerance is set narrow becuase there are a large number of outputs which causes the code to take a very long time to run
    tolerance = 0.005  
    npz_lattice_constant_database = numpy.load(args.lattice_constant_database)
    lattice_constants = npz_lattice_constant_database['arr_0']
    npz_lattice_constant_database.close()
    #call checker
    check_substrate_file(substrate_file, lattice_constants, tolerance, results_file)
    results_file.close()
