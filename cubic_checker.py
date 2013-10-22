#######################################################
##      Checks a cubic film against a substrate      ##
#######################################################

## This program is written as a test. The goal is to calculate lattice matches between
## a film material with cubic symmetry and a substrate. The substrate can be any file 
## of materials but the film file must contain ONLY materials with cubic symmetries.
import numpy

# film file must have ONLY cubic symmetries for this program
film_file = numpy.genfromtxt("cubic.txt", skip_header=1, delimiter="\t", dtype=None)
substrate_file = numpy.genfromtxt("hexagonal.txt", skip_header=1, delimiter="\t", dtype=None)
matches_file = open("cubic_matches_h.txt", "w")

tolerance = 0.08 #tolerance of 8% for mismatch

def search_film_file(film_file, sub_comp, sub_sym, sub_a, sub_c):
    # searches the film file for matches with the substrate 
    for i, l in enumerate(film_file):
        if sub_sym == "C":
            cubic_sub(film_file[i][0], film_file[i][1], sub_comp, sub_sym, sub_a, sub_c, film_file[i][2], film_file[i][3])
        elif sub_sym == "T":
            tetragonal_sub(film_file[i][0], film_file[i][1], sub_comp, sub_sym, sub_a, sub_c, film_file[i][2], film_file[i][3])
        elif sub_sym == "H":
            hexagonal_sub(film_file[i][0], film_file[i][1], sub_comp, sub_sym, sub_a, sub_c, film_file[i][2], film_file[i][3])

def lattice_check(film_file, substrate_file):
    # Reads the material in the substrate list and searches the entire film file for matches
    # Read the next substrate material and search the film file again
    matches_file.write("Film\tSymmetry\tSubstrate\tSymmetry\tMismatch(%)\tRounded Ratio\tOriginal Ratio\tC Mismatch(%)\tC Rounded Ratio\tC Original Ratio\n")
    for i, l in enumerate(substrate_file):
        search_film_file(film_file, substrate_file[i][0], substrate_file[i][1], substrate_file[i][2], substrate_file[i][3])

def round_ratio(original_ratio):
    #rounds the original ratio
    if original_ratio < 1:
        ratio = 1.0 / round(1.0/original_ratio)
    else:
        ratio = round(original_ratio)
    return ratio

def cubic_sub(film_comp, film_sym, sub_comp, sub_sym, sub_a, sub_c, film_a, film_c):
    # called if the substrate has cubic symmetry
    original_ratio = sub_a/film_a
    ratio = round_ratio(original_ratio)
    mismatch = ((sub_a - (ratio*film_a)) / sub_a)
    if abs(mismatch) < tolerance:
        matches_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(film_comp, film_sym, sub_comp, sub_sym, mismatch, ratio, original_ratio))

def tetragonal_sub(film_comp, film_sym, sub_comp, sub_sym, sub_a, sub_c, film_a, film_c):
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

def hexagonal_sub(film_comp, film_sym, sub_comp, sub_sym, sub_a, sub_c, film_a, film_c):
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

def ratio_check(c_value, a_value):
    #check to see if ratio of a and c values for hexagonal and tetragonal symmetries is close to sqrt(2)
    percent_off = ((c_value/(numpy.sqrt(2.0)*a_value)) - 1)
    return abs(percent_off) #returns a percentage of how far off the ratio is

if __name__ == "__main__":
    lattice_check(film_file, substrate_file)
    matches_file.close()
