# makefile: tex -> pdf via xelatex
# v 2.0
# .tex extension is automatically added
file_name = bvar_chart

$(file_name)_cropped.eps: $(file_name).pdf
	pdfcrop $(file_name).pdf $(file_name)_cropped.pdf
	-rm $(file_name)_cropped.eps
	pdftops -eps $(file_name)_cropped.pdf $(file_name)_cropped.eps

$(file_name).pdf: $(file_name).tex
	# protection against biber error
	# http://tex.stackexchange.com/questions/140814/
	rm -rf `biber --cache`

	# create pdf
	# will automatically run pdflatex/biber if necessary
	latexmk -pdf $(file_name).tex
	# add -xelatex option

	# clean auxillary files
	latexmk -c $(file_name).tex

# $(file_name).tex : $(file_name).Rnw
#	Rscript -e "library(knitr); knit('$(file_name).Rnw')"

clean:
	-rm $(file_name)_cropped.pdf
	-rm $(file_name)_cropped.eps
	-rm $(file_name).pdf
