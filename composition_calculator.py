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

parser = argparse.ArgumentParser(description="Software for calculating a range of material composition for an epitaxially grown film on a given substrate.")
parser.add_argument("substrate", type=str, help="File with substrate material data.")
parser.add_argument("reference", type=str, help="File with reference lattice constant and composition data.")
parser.add_argument("tolerance", type=float, help="Tolerance level for mismatch. Enter percentage as decimal.")
args = parser.parse_args()

def check_substrate_file(sub_file, composition_reference_file, output_file):
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
            cubic_sub(sub_file[i][0], sub_file[i][1], sub_file[i][2], composition_reference_file, output_file)
        if sub_file[i][1] == "T":
            tetragonal_sub(sub_file[i][0], sub_file[i][1], sub_file[i][2], sub_file[i][3], composition_reference_file, output_file)
        if sub_file[i][1] == "H":
            hexagonal_sub(sub_file[i][0], sub_file[i][1], sub_file[i][2], sub_file[i][3], composition_reference_file, output_file)

def cubic_sub(sub_comp, sub_sym, a_val, comp_ref_file, result_file):
    """Calculates max/min lattice constant values for a cubic substrate.
    
    Args:
        sub_comp: substrate composition
        sub_sym: substrate symmetry
        a_val: value of lattice constant 'a'
        result_file: a tab delimited .txt file with the calculated maximum and
                      minimum lattice constant values for a specified tolerance
                      value.
        
    Returns:
        Writes a new line to the result_file which contains the results of the
        max/min calculations.    
    """
    a_max = upper_value(a_val)
    a_min = lower_value(a_val)
    a_45_max = upper_value((numpy.sqrt(2.0)*a_val))
    a_45_min = lower_value((numpy.sqrt(2.0)*a_val))
    for i, line in enumerate(comp_ref_file):
        if comp_ref_file[i][6] > a_min and comp_ref_file[i][6] < a_max:
            result_file.write("Al%.1fGa%.1fIn%.1fP%.1fAs%.1fSb%.1f\tC\t%.6f\t%s\t%s\n" \
                              %(comp_ref_file[i][0], comp_ref_file[i][1], comp_ref_file[i][2], comp_ref_file[i][3], 
                                comp_ref_file[i][4], comp_ref_file[i][5], comp_ref_file[i][6], sub_comp, sub_sym))
        if comp_ref_file[i][6] > a_45_min and comp_ref_file[i][6] < a_45_max:
            result_file.write("Al{}Ga{}In{}P{}As{}Sb{}\t{}\t{}\t{}\t{}\n".format(comp_ref_file[i][0], comp_ref_file[i][1],
                               comp_ref_file[i][2], comp_ref_file[i][3], comp_ref_file[i][4], comp_ref_file[i][5], "C (45 deg)",
                               comp_ref_file[i][6], sub_comp, sub_sym))  
                               
def tetragonal_sub(sub_comp, sub_sym, a_val, c_val, comp_ref_file, result_file):
    """Calculates max/min lattice constant values for a tetragonal substrate.
    
    Args:
        sub_comp: substrate composition
        sub_sym: substrate symmetry
        a_val: value of lattice constant 'a'
        c_val: value of lattice constant 'c'
        result_file: a tab delimited .txt file with the calculated maximum and
                      minimum lattice constant values for a specified tolerance
                      value.
        
    Returns:
        Writes a new line to the result_file which contains the results of the
        max/min calculations.    
    """
    a_max = upper_value(a_val)
    a_min = lower_value(a_val)
    a_45_max = upper_value((numpy.sqrt(2.0)*a_val))
    a_45_min = lower_value((numpy.sqrt(2.0)*a_val))
    for i, line in enumerate(comp_ref_file):
        if comp_ref_file[i][6] > a_min and comp_ref_file[i][6] < a_max:
            result_file.write("Al{}Ga{}In{}P{}As{}Sb{}\t{}\t{}\t{}\t{}\n".format(comp_ref_file[i][0], comp_ref_file[i][1],
                               comp_ref_file[i][2], comp_ref_file[i][3], comp_ref_file[i][4], comp_ref_file[i][5], "C",
                               comp_ref_file[i][6], sub_comp, sub_sym))  
        if comp_ref_file[i][6] > a_45_min and comp_ref_file[i][6] < a_45_max:
            result_file.write("Al{}Ga{}In{}P{}As{}Sb{}\t{}\t{}\t{}\t{} (45 deg)\n".format(comp_ref_file[i][0], comp_ref_file[i][1],
                               comp_ref_file[i][2], comp_ref_file[i][3], comp_ref_file[i][4], comp_ref_file[i][5], "C (45 deg)",
                               comp_ref_file[i][6], sub_comp, sub_sym))  

def hexagonal_sub(sub_comp, sub_sym, a_val, c_val, comp_ref_file, result_file):
    """Calculates max/min lattice constant values for a hexagonal substrate.
    
    Args:
        sub_comp: substrate composition
        sub_sym: substrate symmetry
        a_val: value of lattice constant 'a'
        c_val: value of lattice constant 'c'
        result_file: a tab delimited .txt file with the calculated maximum and
                      minimum lattice constant values for a specified tolerance
                      value.
        
    Returns:
        Writes a new line to the result_file which contains the results of the
        max/min calculations.    
    """
    a_max = upper_value(a_val)
    a_min = lower_value(a_val)
    for i, line in enumerate(comp_ref_file):
        if comp_ref_file[i][6] > (numpy.sqrt(2.0)*a_min) and comp_ref_file[i][6] < (numpy.sqrt(2.0)*a_max):
            result_file.write("Al{}Ga{}In{}P{}As{}Sb{}\t{}\t{}\t{}\t{}\n".format(comp_ref_file[i][0], comp_ref_file[i][1],
                               comp_ref_file[i][2], comp_ref_file[i][3], comp_ref_file[i][4], comp_ref_file[i][5], "C (111)",
                               comp_ref_file[i][6], sub_comp, sub_sym))  
                               
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
    results_file_label = "composition_matches_for_" + args.substrate[:-4] + ".txt"
    results_file = open(results_file_label, "w")
    substrate_file = numpy.genfromtxt(args.substrate, comments='#', delimiter="\t", dtype=None)
    lattice_constant_file = numpy.genfromtxt(args.reference, comments='#', delimiter="\t", dtype=None)
    tolerance = args.tolerance # percent tolerance for lattice mismatch as decimal
    check_substrate_file(substrate_file, lattice_constant_file, results_file)
    results_file.close()
