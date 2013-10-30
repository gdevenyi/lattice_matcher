#!/usr/bin/env python
###############################################
##              Lattice Matcher              ##
###############################################

import numpy # includes numpy.sqrt()
import argparse # command line implementation


#### Command line code ###
# Set up argument parser
parser = argparse.ArgumentParser(description="Software for calculating epitaxial lattice matches considering symmetries of different crystal systems.")
# Add arguments needed
parser.add_argument("film", type=str, help="File with film material data")
parser.add_argument("substrate", type=str, help="File with substrate material data")
parser.add_argument("tolerance", type=float, help="Tolerance level for mismatch. Enter percent as a decimal")
# Parse all arguments
args = parser.parse_args()
# Create a label for the matches file. [:-4] strips last 4 characters of file name string
matches_file_label = args.film[:-4] + "_on_" + args.substrate[:-4] + ".txt" # Added ".txt" to specify type of file
matches_file = open(matches_file_label, "w")
# Read input .txt files using numpy.genfromtxt()
film_database = numpy.genfromtxt(args.film, comments="#", delimiter="\t", dtype=None)
substrate_database = numpy.genfromtxt(args.substrate, comments="#", delimiter="\t", dtype=None)
tolerance = args.tolerance # Percent tolerance for lattice mismatch as a decimal

def lattice_matcher(film_database, substrate_database):
    """Creates a file of all acceptable lattice symmetry matches for two input database files.

    Args:
        film_file: a tab delimited .txt file containing data for film materials
        substrate_file: a tab delimited .txt file containing data for substrate materials
    Returns:
        A tab delimited .txt file with all acceptable lattice matches. Lattice matches are calculated by passing
        values to a function which in turn calls other functions to perform the desired calculations and write
        the values to the results .txt file.
    """
    matches_file.write("#Film\tSymmetry\tSubstrate\tSymmetry\tMismatch\tRounded Ratio\tOriginal Ratio\tC Mismatch\tC Rounded Ratio\tC Original Ratio\n")
    for i, l in enumerate(substrate_database):
        film_substrate_comparison(film_database, substrate_database[i][0], substrate_database[i][1], substrate_database[i][2], substrate_database[i][3])

def film_substrate_comparison(film_file, substrate_composition, substrate_symmetry, sub_a, sub_c):
    """Passes values to various functions to perform mismatch calculations.

    Args:
        film_file: a tab delimited .txt file containing data for film materials
        substrate_composition: substrate material composition
        substrate_symmetry: substrate material symmetry
        sub_a: substrate lattice constant 'a' value
        sub_c: substrate lattice constant 'c' value
    Returns:
        A tab delimited .txt file with all acceptable lattice matches. Lattice matches are calculated by calling
        various functions which perform the desired calculations and write the results to the .txt file.
    """
    for i, l in enumerate(film_file):
        if film_file[i][0] == substrate_composition and film_file[i][1] == substrate_symmetry:
            pass
        elif film_file[i][1] == "C":
            cubic_film(film_file[i][0], film_file[i][1], substrate_composition, substrate_symmetry, sub_a, sub_c, film_file[i][2], film_file[i][3])
        elif film_file[i][1] == "T":
            tetragonal_film(film_file[i][0], film_file[i][1], substrate_composition, substrate_symmetry, sub_a, sub_c, film_file[i][2], film_file[i][3])
        elif film_file[i][1] == "H":
            hexagonal_film(film_file[i][0], film_file[i][1], substrate_composition, substrate_symmetry, sub_a, sub_c, film_file[i][2], film_file[i][3])

def cubic_film(film_comp, film_sym, sub_comp, sub_sym, sub_a, sub_c, film_a, film_c):
    """Performs mismatch and ratio checks for a cubic film on various substrates.

    Args:
        film_comp: film composition
        film_sym: film symmetry
        sub_comp: substrate composition
        sub_sym: substrate symmetry
        sub_a: substrate lattice constant 'a' value
        sub_c: substrate lattice constant 'c' value
        film_a: film lattice constant 'a' value
        film_c: film lattice constant 'c' value
    Returns:
        Write a new line containing match information in a tab delimited .txt file. The calculated mismatch
        values must be less than the chosen tolerance level for data to be written to the results file. An
        example output line (with column labels) is:
            Film    Symm    Sub     Symm    Mismatch    Rounded Ratio   Original Ratio   
            CdSe    C       In      T       0.069       0.5             0.538
    """
    if sub_sym == "C":
        # called if the substrate has cubic symmetry
        original_ratio_a = sub_a/film_a
        original_ratio_45 = (numpy.sqrt(2.0)*sub_a)/film_a #cubic film rotated 45 degrees with respect to the substrate
        ratio_a = round_ratio(original_ratio_a)
        ratio_45 = round_ratio(original_ratio_45)
        mismatch_a = ((sub_a - (ratio_a*film_a)) / sub_a)
        mismatch_45 = (((numpy.sqrt(2.0)*sub_a) - (ratio_45*film_a)) / (numpy.sqrt(2.0)*sub_a))
        if abs(mismatch_a) < tolerance:
            matches_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, film_sym, sub_comp, sub_sym, mismatch_a, ratio_a, original_ratio_a))
        if abs(mismatch_45) < tolerance:
            matches_file.write("{}\t{} (45 deg)\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, film_sym, sub_comp, sub_sym, mismatch_45, ratio_45, original_ratio_45))
    elif sub_sym == "T":
        # called if the substrate has tetragonal symmetry
        original_ratio_a = sub_a/film_a #ratio of a-values
        original_ratio_45 = (numpy.sqrt(2.0)*sub_a)/film_a #cubic film rotated 45 degrees with respect to the substrate square face
        original_ratio_c = sub_c/(numpy.sqrt(2.0)*film_a) #cubic (110) and a-plane tetragonal c-value ratio
        ratio_a = round_ratio(original_ratio_a)
        ratio_45 = round_ratio(original_ratio_45)
        ratio_c = round_ratio(original_ratio_c)
        mismatch_a = ((sub_a - (ratio_a*film_a)) / sub_a)
        mismatch_45 = (((numpy.sqrt(2.0)*sub_a) - (ratio_45*film_a)) / (numpy.sqrt(2.0)*sub_a))
        mismatch_c = ((sub_c - (ratio_c*numpy.sqrt(2.0)*film_a)) / sub_c)
        if abs(mismatch_a) < tolerance:
            matches_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, film_sym, sub_comp, sub_sym, mismatch_a, ratio_a, original_ratio_a))
        if abs(mismatch_45) < tolerance:
            matches_file.write("{}\t{} (45 deg)\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, film_sym, sub_comp, sub_sym, mismatch_45, ratio_45, original_ratio_45))
        if ratio_check(sub_c, sub_a) < tolerance and abs(mismatch_c) < tolerance and abs(mismatch_a) < tolerance:
            matches_file.write("{}\t{} (110)\t{}\t{} (a-plane)\t{}\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, film_sym, sub_comp, sub_sym, mismatch_a, ratio_a, original_ratio_a, mismatch_c, ratio_c, original_ratio_c))
    elif sub_sym == "H":
        # called if the substrate has hexagonal symmetry
        original_ratio_a_111 = sub_a/(numpy.sqrt(2.0)*film_a) #ratio of a-values for cubic (111) matches
        original_ratio_a = sub_a/film_a #ratio of a-values
        original_ratio_c = sub_c/(numpy.sqrt(2.0)*film_a) #cubic (110) on a-plane hexagonal c-value ratio
        original_ratio_c_r = numpy.sqrt((sub_c**2)+(3*(sub_a**2)))/(numpy.sqrt(2.0)*film_a) #cubic (110) on r-plane hexagonal c-value ratio
        ratio_a_111 = round_ratio(original_ratio_a_111)
        ratio_a = round_ratio(original_ratio_a)
        ratio_c = round_ratio(original_ratio_c)
        ratio_c_r = round_ratio(original_ratio_c_r)
        mismatch_a_111 = ((sub_a - (ratio_a_111*film_a*numpy.sqrt(2.0))) / sub_a)
        mismatch_a = ((sub_a - (ratio_a*film_a)) / sub_a)
        mismatch_c = ((sub_c - (ratio_c*film_a*numpy.sqrt(2.0))) / sub_c)
        mismatch_c_r = ((numpy.sqrt((sub_c**2)+(3*(sub_a**2))) - (ratio_c_r*film_a*numpy.sqrt(2.0))) / numpy.sqrt((sub_c**2)+(3*(sub_a**2))))
        r_plane_c = numpy.sqrt((sub_c**2)+(3*(sub_a**2))) #side length for camparison of r-plane hex side lengths
        if abs(mismatch_a_111) < tolerance:
            matches_file.write("{}\t{} (111)\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, film_sym, sub_comp, sub_sym, mismatch_a_111, ratio_a_111, original_ratio_a_111))
        if ratio_check(sub_c, sub_a) < tolerance and abs(mismatch_c) and abs(mismatch_a) < tolerance:
            matches_file.write("{}\t{} (110)\t{}\t{} (a-plane)\t{}\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, film_sym, sub_comp, sub_sym, mismatch_a, ratio_a, original_ratio_a, mismatch_c, ratio_c, original_ratio_c))
        if ratio_check(r_plane_c, sub_a) < tolerance and abs(mismatch_c_r) and abs(mismatch_a) < tolerance:
            matches_file.write("{}\t{} (110)\t{}\t{} (r-plane)\t{}\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, film_sym, sub_comp, sub_sym, mismatch_a, ratio_a, original_ratio_a, mismatch_c_r, ratio_c_r, original_ratio_c_r))

def tetragonal_film(film_comp, film_sym, sub_comp, sub_sym, sub_a, sub_c, film_a, film_c):
    """Performs mismatch and ratio checks for a tetragonal film on various substrates.

    Args:
        film_comp: film composition
        film_sym: film symmetry
        sub_comp: substrate composition
        sub_sym: substrate symmetry
        sub_a: substrate lattice constant 'a' value
        sub_c: substrate lattice constant 'c' value
        film_a: film lattice constant 'a' value
        film_c: film lattice constant 'c' value
    Returns:
        Write a new line containing match information in a tab delimited .txt file. The calculated mismatch
        values must be less than the chosen tolerance level for data to be written to the results file. An
        example output line (with column labels) is:
            Film    Symm    Sub     Symm    Mismatch    Rounded Ratio   Original Ratio   
            In      T       Ni      C       0.077       1               1.084
    """
    if sub_sym == "C":
        # called if the substrate has cubic symmetry
        original_ratio_a = sub_a/film_a #tetragonal against cubic a-values
        original_ratio_45 = (numpy.sqrt(2.0)*sub_a)/film_a #tetragonal film rotated 45 degrees with respect to the substrate
        original_ratio_c = (numpy.sqrt(2.0)*sub_a)/film_c # tetragonal c-value against cubic (110) 'c-value'
        ratio_a = round_ratio(original_ratio_a)
        ratio_45 = round_ratio(original_ratio_45)
        ratio_c = round_ratio(original_ratio_c)
        mismatch_a = ((sub_a - (ratio_a*film_a)) / sub_a)
        mismatch_45 = (((numpy.sqrt(2.0)*sub_a) - (ratio_45*film_a)) / (numpy.sqrt(2.0)*sub_a))
        mismatch_c = (((numpy.sqrt(2.0)*sub_a) - (ratio_c*film_a)) / (numpy.sqrt(2.0)*sub_a))
        if abs(mismatch_a) < tolerance:
            matches_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, film_sym, sub_comp, sub_sym, mismatch_a, ratio_a, original_ratio_a))
        if abs(mismatch_45) < tolerance:
            matches_file.write("{}\t{} (45 deg)\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, film_sym, sub_comp, sub_sym, mismatch_45, ratio_45, original_ratio_45))
        if ratio_check(film_c, film_a) < tolerance and abs(mismatch_c) < tolerance and abs(mismatch_a) < tolerance:
            matches_file.write("{}\t{} (a-plane)\t{}\t{} (110)\t{}\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, film_sym, sub_comp, sub_sym, mismatch_a, ratio_a, original_ratio_a, mismatch_c, ratio_c, original_ratio_c))
    elif sub_sym == "T":
        # called if the substrate has tetragonal symmetry
        original_ratio = sub_a/film_a
        original_ratio_45 = (numpy.sqrt(2.0)*sub_a)/film_a #tetragonal film rotated 45 degrees with respect to the substrate
        ratio = round_ratio(original_ratio)
        ratio_45 = round_ratio(original_ratio_45)
        mismatch = ((sub_a - (ratio*film_a)) / sub_a)
        mismatch_45 = (((numpy.sqrt(2.0)*sub_a) - (ratio_45*film_a)) / (numpy.sqrt(2.0)*sub_a))
        if abs(mismatch) < tolerance:
            matches_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, film_sym, sub_comp, sub_sym, mismatch, ratio, original_ratio))
        if abs(mismatch_45) < tolerance:
            matches_file.write("{}\t{} (45 deg)\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, film_sym, sub_comp, sub_sym, mismatch_45, ratio_45, original_ratio_45))
    elif sub_sym == "H":
        # called if the substrate has hexagonal symmetry
        original_ratio_a = sub_a/film_a #tetragonal against hexagonal a-values
        original_ratio_c = sub_c/film_c # tetragonal against hexagonal c-values
        original_ratio_c_r = numpy.sqrt((sub_c**2)+(3*(sub_a**2)))/film_c # uses 'c-value' for r-plane hexagonal
        ratio_a = round_ratio(original_ratio_a)
        ratio_c = round_ratio(original_ratio_c)
        ratio_c_r = round_ratio(original_ratio_c_r)
        mismatch_a = ((sub_a - ratio_a*film_a) / sub_a)
        mismatch_c = ((sub_c - ratio_c*film_c) / sub_c)
        mismatch_c_r = ((numpy.sqrt((sub_c**2)+(3*(sub_a**2))) - ratio_c_r*film_c) / numpy.sqrt((sub_c**2)+(3*(sub_a**2))))
        if abs(mismatch_a) < tolerance and abs(mismatch_c) < tolerance:
            matches_file.write("{}\t{} (a-plane)\t{}\t{} (a-plane)\t{}\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, film_sym, sub_comp, sub_sym, mismatch_a, ratio_a, original_ratio_a, mismatch_c, ratio_c, original_ratio_c))
        if abs(mismatch_a) < tolerance and abs(mismatch_c_r) < tolerance:
            matches_file.write("{}\t{} (a-plane)\t{}\t{} (r-plane)\t{}\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, film_sym, sub_comp, sub_sym, mismatch_a, ratio_a, original_ratio_a, mismatch_c_r, ratio_c_r, original_ratio_c_r))

def hexagonal_film(film_comp, film_sym, sub_comp, sub_sym, sub_a, sub_c, film_a, film_c):
    """Performs mismatch and ratio checks for a hexagonal film on various substrates.

    Args:
        film_comp: film composition
        film_sym: film symmetry
        sub_comp: substrate composition
        sub_sym: substrate symmetry
        sub_a: substrate lattice constant 'a' value
        sub_c: substrate lattice constant 'c' value
        film_a: film lattice constant 'a' value
        film_c: film lattice constant 'c' value
    Returns:
        Write a new line containing match information in a tab delimited .txt file. The calculated mismatch
        values must be less than the chosen tolerance level for data to be written to the results file. An
        example output line (with column labels) is:
            Film    Symm    Sub     Symm    Mismatch    Rounded Ratio   Original Ratio   
            Co      H       C(111)  C       0.076       2               2.012
    """
    if sub_sym == "C":
        # called if the substrate has cubic symmetry
        r_plane_c = numpy.sqrt((sub_c**2)+(3*(sub_a**2))) # 'c-value' of the r-plane hex
        original_ratio_a_111 = (numpy.sqrt(2.0)*sub_a)/film_a # ratio of a values for hex on cubic 111
        original_ratio_a = sub_a/film_a # ratio of a-values
        original_ratio_c = (numpy.sqrt(2.0)*sub_a)/film_a # ratio of c-values
        original_ratio_c_r = (numpy.sqrt(2.0)*sub_a)/r_plane_c
        ratio_a_111 = round_ratio(original_ratio_a_111)
        ratio_a = round_ratio(original_ratio_a)
        ratio_c = round_ratio(original_ratio_c)
        ratio_c_r = round_ratio(original_ratio_c_r)
        mismatch_a_111 = ((numpy.sqrt(2.0)*sub_a - ratio_a_111*film_a) / numpy.sqrt(2.0)*sub_a)
        mismatch_a = ((sub_a - ratio_a*film_a) / sub_a)
        mismatch_c = ((numpy.sqrt(2.0)*sub_a - ratio_c*film_c) / numpy.sqrt(2.0)*sub_a)
        mismatch_c_r = ((numpy.sqrt(2.0)*sub_a - ratio_c_r*r_plane_c) / numpy.sqrt(2.0)*sub_a)
        if abs(mismatch_a_111) < tolerance:
            matches_file.write("{}\t{}\t{}\t{} (111)\t{}\t{}\t{}\n".format(film_comp, film_sym, sub_comp, sub_sym, mismatch_a_111, ratio_a_111, original_ratio_a_111))
        if ratio_check(film_c, film_a) < tolerance and abs(mismatch_c) < tolerance and abs(mismatch_a) < tolerance:
            matches_file.write("{}\t{} (a-plane)\t{}\t{} (110)\t{}\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, film_sym, sub_comp, sub_sym, mismatch_a, ratio_a, original_ratio_a, mismatch_c, ratio_c, original_ratio_c))
        if ratio_check(r_plane_c, film_a) < tolerance and abs(mismatch_c) < tolerance and abs(mismatch_c) < tolerance:
            matches_file.write("{}\t{} (r-plane)\t{}\t{} (110)\t{}\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, film_sym, sub_comp, sub_sym, mismatch_a, ratio_a, original_ratio_a, mismatch_c_r, ratio_c_r, original_ratio_c_r))
    elif sub_sym == "T":
        # called if the substrate has tetragonal symmetry
        r_plane_c = numpy.sqrt((film_c**2)+(3*(film_a**2))) #side length for camparison of r-plane hex side lengths
        original_ratio_a = sub_a/film_a # ratio of a-values
        original_ratio_c = sub_c/film_c #ratio of c-values
        original_ratio_c_r = sub_c/r_plane_c # ratio using the 'c-value' side length of r-plane hex
        ratio_a = round_ratio(original_ratio_a)
        ratio_c = round_ratio(original_ratio_c)
        ratio_c_r = round_ratio(original_ratio_c_r)
        mismatch_a = ((sub_a - ratio_a*film_a) / sub_a)
        mismatch_c = ((sub_c - ratio_c*film_c) / sub_c)
        mismatch_c_r = ((sub_c - ratio_c_r*r_plane_c) / sub_c)
        if abs(mismatch_a) < tolerance and abs(mismatch_c) < tolerance:
            matches_file.write("{}\t{} (a-plane)\t{}\t{} (a-plane)\t{}\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, film_sym, sub_comp, sub_sym, mismatch_a, ratio_a, original_ratio_a, mismatch_c, ratio_c, original_ratio_c))
        if abs(mismatch_a) < tolerance and abs(mismatch_c_r) < tolerance:
            matches_file.write("{}\t{} (r-plane)\t{}\t{} (a-plane)\t{}\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, film_sym, sub_comp, sub_sym, mismatch_a, ratio_a, original_ratio_a, mismatch_c_r, ratio_c_r, original_ratio_c_r))
    elif sub_sym == "H":
        # called if the substrate has hexagonal symmetry
        original_ratio_a = sub_a/film_a #ratio of a-values
        ratio_a = round_ratio(original_ratio_a)
        mismatch_a = ((sub_a - ratio_a*film_a) / sub_a)
        if abs(mismatch_a) < tolerance:
            matches_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, film_sym, sub_comp, sub_sym, mismatch_a, ratio_a, original_ratio_a))

def round_ratio(original_ratio):
    """Rounds a ratio to an integer value.

    Args:
        original_ratio: the ratio value to be rounded.
    Returns:
        ratio: the rounded value of original_ratio.
    """
    if original_ratio < 1:
        ratio = 1.0 / round(1.0/original_ratio)
    else:
        ratio = round(original_ratio)
    return ratio

def ratio_check(c_value, a_value):
    """Checks to see if the ratio of lattice constants (a,c) is close to c = sqrt(2)*a.

    Args:
        a_value: value of lattice constant 'a'
        c_value: value of lattice constant 'c'
    Returns:
        abs(percent_off): the absolute value of the percent (as a decimal) of how far off the ratio is
    """
    percent_off = ((c_value/(numpy.sqrt(2.0)*a_value)) - 1)
    return abs(percent_off) #returns a percentage of how far off the ratio is

if __name__ == "__main__":
    # Call lattice_check to perform the check
    lattice_matcher(film_database, substrate_database)
    # Close any open files
    matches_file.close()
