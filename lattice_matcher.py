###############################################
##              Lattice Checker              ##
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
film_file = numpy.genfromtxt(args.film, comments="#", delimiter="\t", dtype=None)
substrate_file = numpy.genfromtxt(args.substrate, comments="#", delimiter="\t", dtype=None)
tolerance = args.tolerance # Percent tolerance for lattice mismatch as a decimal

def check_film_file(film_file, substrate_composition, substrate_symmetry, sub_a, sub_c):
    # Searches the film file for matches with the substrate 
    for i, l in enumerate(film_file):
        if film_file[i][1] == "C":
            cubic_film(film_file[i][0], film_file[i][1], substrate_composition, substrate_symmetry, sub_a, sub_c, film_file[i][2], film_file[i][3])
        elif film_file[i][1] == "T":
            tetragonal_film(film_file[i][0], film_file[i][1], substrate_composition, substrate_symmetry, sub_a, sub_c, film_file[i][2], film_file[i][3])
        elif film_file[i][1] == "H":
            hexagonal_film(film_file[i][0], film_file[i][1], substrate_composition, substrate_symmetry, sub_a, sub_c, film_file[i][2], film_file[i][3]) 

def lattice_check(film_file, substrate_file):
    # Reads the material in the substrate list and searches the entire film file for matches
    # Read the next substrate material and search the film file again
    matches_file.write("Film\tSymmetry\tSubstrate\tSymmetry\tMismatch\tRounded Ratio\tOriginal Ratio\tC Mismatch\tC Rounded Ratio\tC Original Ratio\n")
    for i, l in enumerate(substrate_file):
        check_film_file(film_file, substrate_file[i][0], substrate_file[i][1], substrate_file[i][2], substrate_file[i][3])

def cubic_film(film_comp, film_sym, sub_comp, sub_sym, sub_a, sub_c, film_a, film_c):
    # perform various mismatch and ratio checks for a cubic film
    # saves matches to a CSV files with title given above
    if sub_sym == "C":
        # called if the substrate has cubic symmetry
        original_ratio_a = sub_a/film_a
        ratio_a = round_ratio(original_ratio_a)
        mismatch_a = ((sub_a - (ratio_a*film_a)) / sub_a)
        if abs(mismatch_a) < tolerance:
            matches_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, film_sym, sub_comp, sub_sym, mismatch_a, ratio_a, original_ratio_a))
    elif sub_sym == "T":
        # called if the substrate has tetragonal symmetry
        original_ratio_a = sub_a/film_a #ratio of a-values
        original_ratio_c = sub_c/(numpy.sqrt(2.0)*film_a) #cubic (110) and a-plane tetragonal c-value ratio
        ratio_a = round_ratio(original_ratio_a)
        ratio_c = round_ratio(original_ratio_c)
        mismatch_a = ((sub_a - (ratio_a*film_a)) / sub_a)
        mismatch_c = ((sub_c - (ratio_c*numpy.sqrt(2.0)*film_a)) / sub_c) 
        if abs(mismatch_a) < tolerance:
            matches_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, film_sym, sub_comp, sub_sym, mismatch_a, ratio_a, original_ratio_a))
        if ratio_check(sub_c, sub_a) < tolerance and abs(mismatch_a) < tolerance: #check if ratio of sub_c/sub_a is close to sqrt(2.0). If not then cubic (110) will definitely not match
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
        if abs(mismatch_a) < tolerance and ratio_check(sub_c, sub_a) < tolerance:
            matches_file.write("{}\t{} (110)\t{}\t{} (a-plane)\t{}\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, film_sym, sub_comp, sub_sym, mismatch_a, ratio_a, original_ratio_a, mismatch_c, ratio_c, original_ratio_c))
        if abs(mismatch_a) < tolerance and ratio_check(r_plane_c, sub_a) < tolerance:
            matches_file.write("{}\t{} (110)\t{}\t{} (r-plane)\t{}\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, film_sym, sub_comp, sub_sym, mismatch_a, ratio_a, original_ratio_a, mismatch_c_r, ratio_c_r, original_ratio_c_r))

def tetragonal_film(film_comp, film_sym, sub_comp, sub_sym, sub_a, sub_c, film_a, film_c):
    # perform various mismatch and ratio checks for a tetragonal film
    if sub_sym == "C":
        # called if the substrate has cubic symmetry
        original_ratio_a = sub_a/film_a #tetragonal against cubic a-values
        original_ratio_c = (numpy.sqrt(2.0)*sub_a)/film_c # tetragonal c-value against cubic (110) 'c-value'
        ratio_a = round_ratio(original_ratio_a)
        ratio_c = round_ratio(original_ratio_c)
        mismatch_a = ((sub_a - (ratio_a*film_a)) / sub_a)
        mismatch_c = (((numpy.sqrt(2.0)*sub_a) - (ratio_c*film_a)) / (numpy.sqrt(2.0)*sub_a))
        if abs(mismatch_a) < tolerance:
            matches_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, film_sym, sub_comp, sub_sym, mismatch_a, ratio_a, original_ratio_a))
        if ratio_check(film_c, film_a) < tolerance and abs(mismatch_a) < tolerance:
            matches_file.write("{}\t{} (a-plane)\t{}\t{} (110)\t{}\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, film_sym, sub_comp, sub_sym, mismatch_a, ratio_a, original_ratio_a, mismatch_c, ratio_c, original_ratio_c))
    elif sub_sym == "T":
        # called if the substrate has tetragonal symmetry
        original_ratio = sub_a/film_a
        ratio = round_ratio(original_ratio)
        mismatch = ((sub_a - (ratio*film_a)) / sub_a)
        if abs(mismatch) < tolerance:
            matches_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, film_sym, sub_comp, sub_sym, mismatch, ratio, original_ratio))
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
    # perform various mismatch and ratio checks for a tetragonal film
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
        if abs(mismatch_a) < tolerance and ratio_check(film_c, film_a) < tolerance:
            matches_file.write("{}\t{} (a-plane)\t{}\t{} (110)\t{}\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, film_sym, sub_comp, sub_sym, mismatch_a, ratio_a, original_ratio_a, mismatch_c, ratio_c, original_ratio_c))
        if abs(mismatch_a) < tolerance and ratio_check(r_plane_c, film_a) < tolerance:
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
    # rounds the original ratio
    if original_ratio < 1:
        ratio = 1.0 / round(1.0/original_ratio)
    else:
        ratio = round(original_ratio)
    return ratio

def ratio_check(c_value, a_value):
    #check to see if ratio of a and c values for hexagonal and tetragonal symmetries is close to sqrt(2)
    percent_off = ((c_value/(numpy.sqrt(2.0)*a_value)) - 1)
    return abs(percent_off) #returns a percentage of how far off the ratio is

if __name__ == "__main__":
    # Call lattice_check to perform the check
    lattice_check(film_file, substrate_file)
    # Close any open files
    matches_file.close()
