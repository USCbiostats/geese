build:
	Rscript -e 'Rcpp::compileAttributes();roxygen2::roxygenize()' && \
		cd .. && R CMD build geese/
install:
	$(MAKE) build && R CMD INSTALL ../geese_*
