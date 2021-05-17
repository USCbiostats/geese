.PHONY: build
../geese.tar.gz: R/* src/*.cpp src/*.h
	$(MAKE) clean ; Rscript -e 'Rcpp::compileAttributes();roxygen2::roxygenize()' && \
		cd .. && R CMD build geese/ && mv geese_*.tar.gz geese.tar.gz
build: ../geese.tar.gz

install: build
	R CMD INSTALL ../geese.tar.gz
check: build
	cd .. && R CMD check --as-cran geese.tar.gz

# Once running, we can set a debug point using 'break [filename].hpp:[linenumber]
# and then type 'run'
debug:
	R -d gdb switch
profile: install
	R --debugger=valgrind --debugger-args='--tool=cachegrind --cachegrind-out-file=test.cache.out'

update:
	rsync -av ../barry/include/barry inst/include

.PHONY: clean
clean:
	rm -rf src/*.o; rm -rf src/*.a; rm -f ../geese.tar.gz
