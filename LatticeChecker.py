###############################################
##              Lattice Checker              ##
###############################################

import numpy # includes numpy.sqrt()

# Choose file for substrate and film
# Refer to documentation for details on input file
film_file = open("elements.txt", "r")
substrate_file = open("semiconductors.txt", "r")
matches_file = open("matches_es.csv", "w")

tolerance = 0.08 # percent tolerance for lattice mismatch

def check_film_file(film_file, substrate_composition, substrate_symmetry, sub_a, sub_c):
    # searches the film file for matches with the substrate 
    for line in film_file:
        column = line.split("\t") #splits the line into a list using tabs (\t) as the delimiter
        film_comp = str(column[0])
        film_sym = str(column[1])
        film_a = float(column[2])
        film_c = float(column[3])
        if film_sym == "C":
            cubic_film(film_comp, film_sym, substrate_composition, substrate_symmetry, sub_a, sub_c, film_a, film_c)
        elif film_sym == "T":
            tetragonal_film(film_comp, film_sym, substrate_composition, substrate_symmetry, sub_a, sub_c, film_a, film_c)
        elif film_sym == "H":
            hexagonal_film(film_comp, film_sym, substrate_composition, substrate_symmetry, sub_a, sub_c, film_a, film_c) 
    film_file.seek(0) # .seek(0) returns the pointer to the top of the file so Python may scan through the file again

def lattice_check(film_file, substrate_file):
    # Reads the material in the substrate list and searches the entire film file for matches
    # Read the next substrate material and search the film file again
    matches_file.write("Film,Symmetry,Substrate,Symmetry,Mismatch(%),Rounded Ratio,Original Ratio" + "\n")
    for line in substrate_file:
        column = line.split("\t") #splits the line into a list using tabs (\t) as the delimiter
        sub_comp = str(column[0])
        sub_sym = str(column[1])
        a_val = float(column[2])
        c_val = float(column[3])
        check_film_file(film_file, sub_comp, sub_sym, a_val, c_val)
    substrate_file.seek(0)

def cubic_film(film_comp, film_sym, sub_comp, sub_sym, sub_a, sub_c, film_a, film_c):
    # perform various mismatch and ratio checks for a cubic film
    # saves matches to a CSV files with title given above
    if sub_sym == "C":
        # cubic substrate
        original_ratio = sub_a/film_a
        ratio = ratio_cal(original_ratio)
        mismatch = ((sub_a - (ratio*film_a)) / sub_a)
        if abs(mismatch) < tolerance:
            matches_file.write(film_comp + "," + film_sym + "," + sub_comp + "," + sub_sym + "," + str(mismatch) + "," + str(ratio) + "," + str(original_ratio) + "\n")
    elif sub_sym == "T":
        # tetragonal substrate
        original_a_ratio = sub_a/film_a
        a_ratio = ratio_cal(original_a_ratio)
        a_mismatch = ((sub_a - (a_ratio*film_a)) / sub_a)
        if abs(a_mismatch) < tolerance:
            matches_file.write(film_comp + "," + film_sym + "," + sub_comp + "," + sub_sym + "," + str(a_mismatch) + "," + str(a_ratio) + "," + str(original_a_ratio) + "\n")
    elif sub_sym == "H":
        # hexagonal substrate 
        original_ratio = sub_a/(numpy.sqrt(2.0)*film_a)
        ratio = ratio_cal(original_ratio)
        mismatch = ((sub_a - (ratio*film_a*numpy.sqrt(2.0))) / sub_a)
        if abs(mismatch) < tolerance:
            matches_file.write(film_comp + "," + film_sym + " (111)," + sub_comp + "," + sub_sym + "," + str(mismatch) + "," + str(ratio) + "," + str(original_ratio) + "\n")

def tetragonal_film(film_comp, film_sym, sub_comp, sub_sym, sub_a, sub_c, film_a, film_c):
    # perform various mismatch and ratio checks for a tetragonal film
    if sub_sym == "C":
        # cubic substrate
        original_ratio = sub_a/film_a
        ratio = ratio_cal(original_ratio)
        mismatch = ((sub_a - (ratio*film_a)) / sub_a)
        if abs(mismatch) < tolerance:
            matches_file.write(film_comp + "," + film_sym + "," + sub_comp + "," + sub_sym + "," + str(mismatch) + "," + str(ratio) + "," + str(original_ratio) + "\n")
    elif sub_sym == "T":
        # tetragonal substrate
        original_ratio = sub_a/film_a
        ratio = ratio_cal(original_ratio)
        mismatch = ((sub_a - (ratio*film_a)) / sub_a)
        if abs(mismatch) < tolerance:
            matches_file.write(film_comp + "," + film_sym + "," + sub_comp + "," + sub_sym + "," + str(mismatch) + "," + str(ratio) + "," + str(original_ratio) + "\n")
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
            matches_file.write(film_comp + "," + film_sym + "," + sub_comp + "," + sub_sym + " (111)," + str(mismatch) + "," + str(ratio) + "," + str(original_ratio) + "\n")
    elif sub_sym == "T":
        # tetragonal substrate
        pass
    elif sub_sym == "H":
        # hexagonal substrate
        original_ratio = sub_a/film_a
        ratio = ratio_cal(original_ratio)
        mismatch = ((sub_a - (ratio*film_a)) / sub_a)
        if abs(mismatch) < tolerance:
            matches_file.write(film_comp + "," + film_sym + "," + sub_comp + "," + sub_sym + "," + str(mismatch) + "," + str(ratio) + "," + str(original_ratio) + "\n")

def ratio_cal(original_ratio):
    # rounds the original ratio
    if original_ratio < 1:
        ratio = 1.0 / round(1.0/original_ratio)
    else:
        ratio = round(original_ratio)
    return ratio

lattice_check(film_file, substrate_file) 

film_file.close()
substrate_file.close()
matches_file.close()
