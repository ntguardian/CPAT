SHELL=/bin/sh
RSCRIPT=Rscript

SMALLREPLICATIONS=20

POWERREPLICATIONS=5000
CONTEXTPREFIX=Context
SIMDATAPREFIX=SimData
SIMSTATPREFIX=SimTest
CONTEXTGENERATORS=$(wildcard exec/$(CONTEXTPREFIX)*.R)
SIMDATAGENERATORS=$(wildcard exec/$(SIMDATAPREFIX)*.R)
SIMSTATGENERATORS=$(wildcard exec/$(SIMSTATPREFIX)*.R)
CONTEXTDATA=$(wildcard data/$(CONTEXTPREFIX)*.Rda)
SIMDATADATA=$(wildcard data/$(SIMDATAPREFIX)*.Rda)
SIMSTATDATA=$(wildcard data/$(SIMSTATPREFIX)*.Rda)
SIMSEED=201902

SIMSNORMALRENYIRESID=data/SimsNormalRenyiResid.Rda
ALLSIMS=$(SIMSNORMALRENYIRESID)

SIMSNORMALRENYIRESIDDF=$(subst .Rda,DataFrame.Rda,$(SIMSNORMALRENYIRESID))
ALLSIMSDATAFRAME=$(subst .Rda,DataFrame.Rda,$(ALLSIMS))

LEVEL=0.05

.PHONY : all
all : inst/Makefile inst/package $(ALLSIMSDATAFRAME)

data/$(CONTEXTPREFIX)%.Rda : exec/$(CONTEXTPREFIX)%.R
	$(RSCRIPT) $< -o $@

data/$(SIMDATAPREFIX)%.Rda : exec/$(SIMDATAPREFIX)%.R
	$(RSCRIPT) $< -o $@

data/$(SIMSTATPREFIX)%.Rda : exec/$(SIMSTATPREFIX)%.R
	$(RSCRIPT) $< -o $@

$(SIMSNORMALRENYIRESID) : data/$(CONTEXTPREFIX)Main.Rda \
                          data/$(SIMDATAPREFIX)NormalXY.Rda \
                          data/$(SIMSTATPREFIX)RenyiTypeResid.Rda
	$(RSCRIPT) exec/PowerSimRegression.R -C $(word 1, $^) -S $(word 2, $^) \
		 -T $(word 3, $^) -o $@ -N $(POWERREPLICATIONS) -s $(SIMSEED)23 -v

$(SIMSNORMALRENYIRESIDDF) : $(SIMSNORMALRENYIRESID)
	$(RSCRIPT) exec/Aggregator.R -i $< -o $@ -a $(LEVEL) \
		 -T data/$(SIMSTATPREFIX)RenyiTypeResid.Rda \
		 -C data/$(CONTEXTPREFIX)Main.Rda

.PHONY : simconfig
simconfig : $(CONTEXTGENERATORS) $(SIMDATAGENERATORS) $(SIMSTATGENERATORS)
	make $(CONTEXTGENERATORS:exec/%.R=data/%.Rda)
	make $(SIMDATAGENERATORS:exec/%.R=data/%.Rda)
	make $(SIMSTATGENERATORS:exec/%.R=data/%.Rda)

inst/Makefile : Makefile
	cp $< $@

inst/package : package
	cp $< $@

package : R/*.R src/*.cpp
	$(RSCRIPT) exec/RemakePackage.R
	touch package

.PHONY : clean
clean :
	-rm $(CONTEXTDATA)
	-rm $(SIMDATADATA)
	-rm $(SIMSTATDATA)
	-rm $(ALLSIMS)

.PHONY : init
init :
	make simconfig

.PHONY : dependencies
dependencies :
	$(RSCRIPT) exec/GetPackages.R

.PHONY : small
small :
	make POWERREPLICATIONS=$(SMALLREPLICATIONS)
