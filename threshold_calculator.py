#!/usr/bin/env python
###############################################################################
##                           Composition Calculator                          ##
###############################################################################

"""
Software to calculate a range of materials which may be epitaxially grown on a 
given substrate. 

CURRENT PROGESS: Calculates the max/min lattice constants for each substrate
given a threshold value.
"""

import numpy #includes numpy.sqrt()
import argparse #for command line implementation

# set up argument parser
parser = argparse.ArgumentParser(description="Software for a range of material composition for an epitaxially grown film on a given substrate.")
# add arguments
parser.add_argument("substrate", type=str, help="File with substrate material data.")
parser.add_argument("tolerance", type=float, help="Tolerance level for mismatch. Enter percentage as decimal.")
# parse all arguments
args = parser.parse_args()

def check_substrate_file(sub_file, output_file):
    """Calls functions for calculations based on the information obtained from
       a supplied database file.
    
    Args:
        sub_file: database file with substrate material information
        output_file: tab delimited .txt file where results are written
        
    Returns:
        A tab delimited .txt file with the maximum and minimum values related
        to the lattice constants of the substrate material. 
    """
    output_file.write("#Substrate\tSymmetry\ta-man\ta-min\tc-max\tc-min\n")
    for i, l in enumerate(sub_file):
        if sub_file[i][1] == "C":
            cubic_sub(sub_file[i][0], sub_file[i][1], sub_file[i][2], output_file)
        if sub_file[i][1] == "T":
            tetragonal_sub(sub_file[i][0], sub_file[i][1], sub_file[i][2], sub_file[i][3], output_file)
        if sub_file[i][1] == "H":
            hexagonal_sub(sub_file[i][0], sub_file[i][1], sub_file[i][2], sub_file[i][3], output_file)

def cubic_sub(sub_comp, sub_sym, a_val, results_file):
    """Calculates max/min lattice constant values for a cubic substrate.
    
    Args:
        sub_comp: substrate composition
        sub_sym: substrate symmetry
        a_val: value of lattice constant 'a'
        results_file: a tab delimited .txt file with the calculated maximum and
                      minimum lattice constant values for a specified tolerance
                      value.
        
    Returns:
        Writes a new line to the results_file which contains the results of the
        max/min calculations.    
    """
    a_max = upper_value(a_val)
    a_min = lower_value(a_val)
    c_max = upper_value((numpy.sqrt(2.0)*a_val))
    c_min = lower_value((numpy.sqrt(2.0)*a_val))
    results_file.write("{}\t{}\t{}\t{}\n".format(sub_comp, sub_sym, a_max, a_min))
    results_file.write("\t{} (110)\t{}\t{}\t{}\t{}\n".format(sub_sym, a_max, a_min, c_max, c_min))
    results_file.write("\t{} (111)\t{}\t{}\n".format(sub_sym, c_max, c_min))

def tetragonal_sub(sub_comp, sub_sym, a_val, c_val, results_file):
    """Calculates max/min lattice constant values for a tetragonal substrate.
    
    Args:
        sub_comp: substrate composition
        sub_sym: substrate symmetry
        a_val: value of lattice constant 'a'
        c_val: value of lattice constant 'c'
        results_file: a tab delimited .txt file with the calculated maximum and
                      minimum lattice constant values for a specified tolerance
                      value.
        
    Returns:
        Writes a new line to the results_file which contains the results of the
        max/min calculations.    
    """
    a_max = upper_value(a_val)
    a_min = lower_value(a_val)
    c_max = upper_value(c_val)
    c_min = lower_value(c_val)
    results_file.write("{}\t{}\t{}\t{}\n".format(sub_comp, sub_sym, a_max, a_min))
    results_file.write("\t{} (a-plane)\t{}\t{}\t{}\t{}\n".format(sub_sym, a_max, a_min, c_max, c_min))

def hexagonal_sub(sub_comp, sub_sym, a_val, c_val, results_file):
    """Calculates max/min lattice constant values for a hexagonal substrate.
    
    Args:
        sub_comp: substrate composition
        sub_sym: substrate symmetry
        a_val: value of lattice constant 'a'
        c_val: value of lattice constant 'c'
        results_file: a tab delimited .txt file with the calculated maximum and
                      minimum lattice constant values for a specified tolerance
                      value.
        
    Returns:
        Writes a new line to the results_file which contains the results of the
        max/min calculations.    
    """
    a_max = upper_value(a_val)
    a_min = lower_value(a_val)
    c_max = upper_value(c_val)
    c_min = lower_value(c_val)
    c_r_max = upper_value((numpy.sqrt(c_val**2 + (3.0*a_val**2))))
    c_r_min = lower_value((numpy.sqrt(c_val**2 + (3.0*a_val**2))))
    results_file.write("{}\t{}\t{}\t{}\n".format(sub_comp, sub_sym, a_max, a_min))
    results_file.write("\t{} (a-plane)\t{}\t{}\t{}\t{}\n".format(sub_sym, a_max, a_min, c_max, c_min))
    results_file.write("\t{} (r-plane)\t{}\t{}\t{}\t{}\n".format(sub_sym, a_max, a_min, c_r_max, c_r_min))

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

if __name__ == "__main__":
    # create a label for the matches file.
    minmax_file_label = "tolerance_values_for_" + args.substrate[:-4] + ".txt"
    min_max_file = open(minmax_file_label, "w")
    #read input .txt file with numpy.genfromtxt()
    substrate_file = numpy.genfromtxt(args.substrate, comments='#', delimiter="\t", dtype=None)
    tolerance = args.tolerance # percent tolerance for lattice mismatch as decimal
    check_substrate_file(substrate_file, min_max_file)
    min_max_file.close()
