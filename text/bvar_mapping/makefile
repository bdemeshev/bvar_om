# makefile: Rnw -> tex -> pdf
# v 2.0
# .Rnw extension is automatically added
file_name = bvar_mapping

$(file_name).pdf: $(file_name).tex bibliobase2.bib
	# protection against biber error
	# http://tex.stackexchange.com/questions/140814/
	rm -rf `biber --cache`

	# create pdf
	# will automatically run pdflatex/biber if necessary
	latexmk -pdf $(file_name).tex

	# clean auxillary files
	latexmk -c $(file_name).tex

	# manual clean
	rm $(file_name).bbl
	rm $(file_name).tdo
	rm $(file_name).run.xml


# $(file_name).tex : $(file_name).Rnw
#	Rscript -e "library(knitr); knit('$(file_name).Rnw')"

clean:
	-rm $(file_name).pdf $(file_name).bbl $(file_name).tdo $(file_name).run.xml $(file_name).idx $(file_name).aux $(file_name).out $(file_name).tex.blg $(file_name).log $(file_name).bcf $(file_name).fls $(file_name).fdb_latexmk
