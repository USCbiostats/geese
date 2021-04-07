build:
	Rscript -e 'Rcpp::compileAttributes();roxygen2::roxygenize()' && \
		cd .. && R CMD build geese/
install:
	$(MAKE) build && R CMD INSTALL ../geese_*
check:
	$(MAKE) build && R CMD check --as-cran ../geese_*

# Once running, we can set a debug point using 'break [filename].hpp:[linenumber]
# and then type 'run'
debug:
	R -d gdb switch
