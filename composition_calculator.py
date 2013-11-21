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
parser.add_argument("substrate", type=str, help="File with substrate material data.")
#parser.add_argument("reference", type=str, help="File with reference lattice constant and composition data.")
#parser.add_argument("tolerance", type=float, help="Tolerance level for mismatch. Enter percentage as decimal.")
args = parser.parse_args()

def check_substrate_file(sub_file, lattice_const_file, output_file):
    """Calls functions for calculations based on the information obtained from
       a supplied database file.
    
    Args:
        sub_file: database file with substrate material information
        output_file: tab delimited .txt file where results are written
        
    Returns:
        A tab delimited .txt file with the maximum and minimum values related
        to the lattice constants of the substrate material. 
    """
    output_file.write("#Film Composition\tFilm Symmetry\tFlim a\tSubstrate\tSymmetry\n")
    for i, l in enumerate(sub_file):
        if sub_file[i][1] == "C":
            cubic_sub(sub_file[i][0], sub_file[i][1], sub_file[i][2], lattice_const_file, output_file)
        '''
        if sub_file[i][1] == "T":
            tetragonal_sub(sub_file[i][0], sub_file[i][1], sub_file[i][2], sub_file[i][3], lattice_const_file, output_file)
        if sub_file[i][1] == "H":
            hexagonal_sub(sub_file[i][0], sub_file[i][1], sub_file[i][2], sub_file[i][3], lattice_const_file, output_file)
        '''

def cubic_sub(sub_comp, sub_sym, sub_a_val, lattice_consts , result_file):
    """Calculates max/min lattice constant values for a cubic substrate.
    
    Args:
        sub_comp: substrate composition
        sub_sym: substrate symmetry
        sub_a_val: value of lattice constant 'a'
        lattice_const_array: numpy array of composition and lattice constants.
                             Lattice constants must be the right-most entry in 
                             each line of the array
        result_file: a tab delimited .txt file with the calculated maximum and
                      minimum lattice constant values for a specified tolerance
                      value.
        
    Returns:
        Writes a new line to the result_file which contains the results of the
        lattice constant comparision.    
    """    
    good_vals = (lattice_consts[:,-1] > (1.-tolerance)*sub_a_val) & (lattice_consts[:,-1] < (1. + tolerance)*sub_a_val)
    for i, line in enumerate(good_vals):
        if line:
            film_comp = cust_print(lattice_consts[i][0], lattice_consts[i][1], lattice_consts[i][2], 
                                   lattice_consts[i][3], lattice_consts[i][4], lattice_consts[i][5])
            result_file.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, "C", lattice_consts[i][6], sub_comp, sub_sym, sub_a_val))
    good_vals45 = (lattice_consts[:,-1] > (1.-tolerance)*sub_a_val*numpy.sqrt(2.0)) & (lattice_consts[:,-1] < (1. + tolerance)*sub_a_val*numpy.sqrt(2.0))
    for i, line in enumerate(good_vals45):
        if line:
            film_comp = cust_print(lattice_consts[i][0], lattice_consts[i][1], lattice_consts[i][2], 
                                   lattice_consts[i][3], lattice_consts[i][4], lattice_consts[i][5])
            result_file.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, "C (45deg)", lattice_consts[i][6], sub_comp, sub_sym, sub_a_val))

'''                            
def tetragonal_sub(sub_comp, sub_sym, sub_a_val, sub_c_val, lattice_const_array , result_file):
    """Calculates max/min lattice constant values for a tetragonal substrate.
    
    Args:
        sub_comp: substrate composition
        sub_sym: substrate symmetry
        sub_a_val: value of lattice constant 'a'
        sub_c_val: value of lattice constant 'c'
        lattice_const_array: numpy array of composition and lattice constants.
                             Lattice constants must be the right-most entry in 
                             each line of the array
        result_file: a tab delimited .txt file with the calculated maximum and
                      minimum lattice constant values for a specified tolerance
                      value.
    Returns:
        Writes a new line to the result_file which contains the results of the
        lattice constant comparision.    
    """
    accepted_lattice_const = create_tolerance_array(lattice_const_array, sub_a_val)
    if len(accepted_lattice_const) == 0:
        pass
    for i, line in enumerate(accepted_lattice_const):
        film_comp = cust_print(accepted_lattice_const[i][0], accepted_lattice_const[i][1], accepted_lattice_const[i][2], 
                               accepted_lattice_const[i][3], accepted_lattice_const[i][4], accepted_lattice_const[i][5])
        result_file.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, "C", accepted_lattice_const[i][6], sub_comp, sub_sym, sub_a_val))
    accepted_lattice_const45 = create_tolerance_array(lattice_const_array, (numpy.sqrt(2.0)*sub_a_val))
    if len(accepted_lattice_const45) == 0:
        pass
    for i, line in enumerate(accepted_lattice_const45):
        film_comp = cust_print(accepted_lattice_const45[i][0], accepted_lattice_const45[i][1], accepted_lattice_const45[i][2], 
                               accepted_lattice_const45[i][3], accepted_lattice_const45[i][4], accepted_lattice_const45[i][5])
        result_file.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, "C (45deg)", accepted_lattice_const45[i][6], sub_comp, sub_sym, sub_a_val))

def hexagonal_sub(sub_comp, sub_sym, sub_a_val, sub_c_val, lattice_const_array , result_file):
    """Calculates max/min lattice constant values for a hexagonal substrate.
    
    Args:
        sub_comp: substrate composition
        sub_sym: substrate symmetry
        sub_a_val: value of lattice constant 'a'
        sub_c_val: value of lattice constant 'c'
        lattice_const_array: numpy array of composition and lattice constants.
                             Lattice constants must be the right-most entry in 
                             each line of the array
        result_file: a tab delimited .txt file with the calculated maximum and
                      minimum lattice constant values for a specified tolerance
                      value.
    Returns:
        Writes a new line to the result_file which contains the results of the
        lattice constant comparision.    
    """
    accepted_lattice_const = create_tolerance_array(lattice_const_array, (sub_a_val*numpy.sqrt(2.0)))
    if len(accepted_lattice_const) == 0:
        pass
    for i, line in enumerate(accepted_lattice_const):
        film_comp = cust_print(accepted_lattice_const[i][0], accepted_lattice_const[i][1], accepted_lattice_const[i][2], 
                               accepted_lattice_const[i][3], accepted_lattice_const[i][4], accepted_lattice_const[i][5])
        result_file.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, "C", accepted_lattice_const[i][6], sub_comp, sub_sym, sub_a_val))
'''                               
def upper_value(lattice_constant):
    """Calculate the maximum lattice constant based on the specified tolerance.
    
    Args:
        lattice_constant: lattice constant of a material
        
    Returns:
        upper: the maximum accepted lattice constant based on tolerance level
    """
    upper = (1.0 + tolerance) * lattice_constant
    return upper

def lower_value(lattice_constant):
    """Calculate the minimum lattice constant based on the specified tolerance.
    
    Args:
        lattice_constant: lattice constant of a material
        
    Returns:
        lower: the minimum accepted lattice constant based on tolerance level
    """
    lower = (1.0 - tolerance) * lattice_constant
    return lower
    
def cust_print(x_Al, x_Ga, x_In, y_P, y_As, y_Sb):
    """Custom print function to simplify output of composition. Input arguments
    agree with the following rules:
            x_Al + x_Ga + x_In = 1
            y_P + y_As + y_Sb = 1
            
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
    if x_Al != 0.0:
        composition += "Al"+'{0:.{1}f}'.format(x_Al, 2)
    if x_Ga != 0.0:
        composition += "Ga"+'{0:.{1}f}'.format(x_Ga, 2)
    if x_In != 0.0:
        composition += "In"+'{0:.{1}f}'.format(x_In, 2)
    if y_P != 0.0:
        composition += "P"+'{0:.{1}f}'.format(y_P, 2)
    if y_As != 0.0:
        composition += "As"+'{0:.{1}f}'.format(y_As, 2)
    if y_Sb != 0.0:
        composition += "Sb"+'{0:.{1}f}'.format(y_Sb, 2)
    return composition
    
def create_tolerance_array(input_array, sub_lattice_const):
    """Creates an array of materials with lattice constants which match the 
    specified tolerance percentage.
    
    Args:
        input_array: a numpy array with the lattice constant stored in the 
                     right-most position of a line in the array
        sub_lattice_const: reference substrate lattice constant
    
    Returns:
        out_array: an array with materials with lattice constants within the 
                   tolerance level
    """
    if len(input_array) == 0:
        pass
    temp = input_array[input_array[:, -1] >= lower_value(sub_lattice_const)]
    if len(temp) == 0:
        return temp
    else:
        out_array = temp[temp[:, -1] <= upper_value(sub_lattice_const)]
        return out_array
    

if __name__ == "__main__":
    # create a label for the matches file.
    results_file_label = "composition_matches_for_" + args.substrate[:-4] + ".txt"
    results_file = open(results_file_label, "w")
    substrate_file = numpy.genfromtxt(args.substrate, comments='#', delimiter="\t", dtype=None)
    tolerance = 0.1 #high tolerance set to check if program works. Will reduce after testing
    #tolerance = args.tolerance # percent tolerance for lattice mismatch as decimal
    #load .npz lattice constant file
    initial_array = numpy.load("test_array_file.npz")
    lattice_constants = initial_array['arr_0']
    initial_array.close()
    #call checker
    check_substrate_file(substrate_file, lattice_constants, results_file)
    results_file.close()
