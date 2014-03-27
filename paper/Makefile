.SECONDARY:

# A list of the SVG source
# https://www.gnu.org/software/make/manual/html_node/Wildcard-Function.html
SVGS = $(wildcard graphics/raw/*.svg)

# A list of PDFs for each SVG
# https://www.gnu.org/software/make/manual/html_node/Text-Functions.html
PDFS = $(patsubst graphics/raw/%.svg, graphics/%.pdf, $(SVGS))

# The first target is the default
# https://www.gnu.org/software/make/manual/html_node/Rules.html
all: lattice_matcher.pdf

# Remove all auto-generated content with:
#   $ make clean
clean:
	rm -f *.aux *.fdb_latexmk *.fls *.gz *.log *.out *.pdf *.spl graphics/*.pdf

# Convert LaTeX source to PDF
# The pattern rule syntax:
#   %.pdf: %.tex ...
#     recipe...
# builds lattice_matcher.pdf from lattice_matcher.tex (and any other
# prerequisites).
# https://www.gnu.org/software/make/manual/html_node/Pattern-Rules.html
#
# $* expands to the target stem (lattice_matcher)
# https://www.gnu.org/software/make/manual/html_node/Automatic-Variables.html
%.pdf: %.tex $(PDFS)
	latexmk -pdf -g -silent -f -pdflatex="pdflatex -interaction=nonstopmode -synctex=1 %O %S" "$*"

# Convert SVGs to PDFs
# $< expands to the first prerequisite (graphics/raw/%.svg)
# $@ expands to the target (graphics/%.pdf)
# https://www.gnu.org/software/make/manual/html_node/Automatic-Variables.html
graphics/%.pdf: graphics/raw/%.svg
	rsvg-convert -f pdf "$<" > "$@"
