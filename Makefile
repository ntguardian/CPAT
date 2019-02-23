SHELL=/bin/sh
RSCRIPT=Rscript

SMALLREPLICATIONS=20

LEVEL=0.05

.PHONY : all
all : inst/Makefile inst/package

inst/Makefile : Makefile
	cp $< $@

inst/package : package
	cp $< $@

package : R/*.R src/*.cpp
	$(RSCRIPT) exec/RemakePackage.R
	touch package
