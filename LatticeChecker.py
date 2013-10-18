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
film_file = numpy.genfromtxt(args.film, skip_header=1, delimiter="\t", dtype=None)
substrate_file = numpy.genfromtxt(args.substrate, skip_header=1, delimiter="\t", dtype=None)
tolerance = args.tolerance # Percent tolerance for lattice mismatch as a decimal

'''
### Fixed input code ###
# Output file will be a tab delimited .txt file
matches_file = open("matches_es.txt", "w")
# Read input .txt files using numpy.genfromtxt
film_file = numpy.genfromtxt("elements.txt", skip_header=1, delimiter="\t", dtype=None)
substrate_file = numpy.genfromtxt("semiconductors.txt", skip_header=1, delimiter="\t", dtype=None)
tolerance = 0.08 # Percent tolerance for lattice mismatch as a decimal
'''

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
    matches_file.write("Film\tSymmetry\tSubstrate\tSymmetry\tMismatch(%)\tRounded Ratio\tOriginal Ratio\tC Mismatch(%)\tC Rounded Ratio\tC Original Ratio\n")
    for i, l in enumerate(substrate_file):
        check_film_file(film_file, substrate_file[i][0], substrate_file[i][1], substrate_file[i][2], substrate_file[i][3])

def cubic_film(film_comp, film_sym, sub_comp, sub_sym, sub_a, sub_c, film_a, film_c):
    # perform various mismatch and ratio checks for a cubic film
    # saves matches to a CSV files with title given above
    if sub_sym == "C":
        # cubic substrate
        original_ratio = sub_a/film_a
        ratio = ratio_cal(original_ratio)
        mismatch = ((sub_a - (ratio*film_a)) / sub_a)
        if abs(mismatch) < tolerance:
            matches_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, film_sym, sub_comp, sub_sym, mismatch, ratio, original_ratio))
    elif sub_sym == "T":
        # tetragonal substrate
        original_ratio = sub_a/film_a
        ratio = ratio_cal(original_ratio)
        mismatch = ((sub_a - (ratio*film_a)) / sub_a)
        if abs(mismatch) < tolerance:
            matches_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, film_sym, sub_comp, sub_sym, mismatch, ratio, original_ratio))
    elif sub_sym == "H":
        # hexagonal substrate 
        original_ratio = sub_a/(numpy.sqrt(2.0)*film_a)
        ratio = ratio_cal(original_ratio)
        mismatch = ((sub_a - (ratio*film_a*numpy.sqrt(2.0))) / sub_a)
        if abs(mismatch) < tolerance:
            matches_file.write("{}\t{} (111)\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, film_sym, sub_comp, sub_sym, mismatch, ratio, original_ratio))

def tetragonal_film(film_comp, film_sym, sub_comp, sub_sym, sub_a, sub_c, film_a, film_c):
    # perform various mismatch and ratio checks for a tetragonal film
    if sub_sym == "C":
        # cubic substrate
        original_ratio = sub_a/film_a
        ratio = ratio_cal(original_ratio)
        mismatch = ((sub_a - (ratio*film_a)) / sub_a)
        if abs(mismatch) < tolerance:
            matches_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, film_sym, sub_comp, sub_sym, mismatch, ratio, original_ratio))
    elif sub_sym == "T":
        # tetragonal substrate
        original_ratio = sub_a/film_a
        ratio = ratio_cal(original_ratio)
        mismatch = ((sub_a - (ratio*film_a)) / sub_a)
        if abs(mismatch) < tolerance:
            matches_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, film_sym, sub_comp, sub_sym, mismatch, ratio, original_ratio))
    elif sub_sym == "H":
        # hexagonal substrate 
        pass

def hexagonal_film(film_comp, film_sym, sub_comp, sub_sym, sub_a, sub_c, film_a, film_c):
    # perform various mismatch and ratio checks for a tetragonal film
    if sub_sym == "C":
        # cubic substrate
        original_ratio = (sub_a*numpy.sqrt(2.0))/film_a
        ratio = ratio_cal(original_ratio)
        mismatch = (((numpy.sqrt(2.0)*sub_a) - (ratio*film_a)) / (numpy.sqrt(2.0)*sub_a))
        if abs(mismatch) < tolerance:
            matches_file.write("{}\t{}\t{}\t{} (111)\t{}\t{}\t{}\n".format(film_comp, film_sym, sub_comp, sub_sym, mismatch, ratio, original_ratio))
    elif sub_sym == "T":
        # tetragonal substrate
        pass
    elif sub_sym == "H":
        # hexagonal substrate
        original_ratio = sub_a/film_a
        ratio = ratio_cal(original_ratio)
        mismatch = ((sub_a - (ratio*film_a)) / sub_a)
        if abs(mismatch) < tolerance:
            matches_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, film_sym, sub_comp, sub_sym, mismatch, ratio, original_ratio))
            
def ratio_cal(original_ratio):
    # rounds the original ratio
    if original_ratio < 1:
        ratio = 1.0 / round(1.0/original_ratio)
    else:
        ratio = round(original_ratio)
    return ratio

if __name__ == "__main__":
    # Call lattice_check to perform the check
    lattice_check(film_file, substrate_file)
    # Close any open files
    matches_file.close()
