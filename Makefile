SHELL=/bin/sh
RSCRIPT=Rscript

SMALLREPLICATIONS=20

POWERREPLICATIONS=20000
CONTEXTPREFIX=Context
SIMDATAPREFIX=SimData
SIMSTATPREFIX=SimTest
CONTEXTGENERATORS=$(wildcard exec/$(CONTEXTPREFIX)*.R)
SIMDATAGENERATORS=$(wildcard exec/$(SIMDATAPREFIX)*.R)
SIMSTATGENERATORS=$(wildcard exec/$(SIMSTATPREFIX)*.R)
CONTEXTDATA=$(wildcard data/$(CONTEXTPREFIX)*.Rda)
SIMDATADATA=$(wildcard data/$(SIMDATAPREFIX)*.Rda)
SIMSTATDATA=$(wildcard data/$(SIMSTATPREFIX)*.Rda)
SIMSEED=2019

LEVEL=0.05

SIMSNORMALRENYIRESID=data/SimsNormalRenyiResid.Rda
SIMSNORMALCUSUM=data/SimsNormalCUSUM.Rda
SIMSNORMALHS=data/SimsNormalHS.Rda
SIMSARMARENYIRESID=data/SimsARMARenyiResid.Rda
SIMSARMACUSUM=data/SimsARMACUSUM.Rda
SIMSARMAHS=data/SimsARMAHS.Rda
SIMSGARCHRENYIRESID=data/SimsGARCHRenyiResid.Rda
SIMSGARCHCUSUM=data/SimsGARCHCUSUM.Rda
SIMSGARCHHS=data/SimsGARCHHS.Rda
ALLSIMS=$(SIMSNORMALRENYIRESID) $(SIMSNORMALCUSUM) $(SIMSNORMALHS) \
        $(SIMSARMARENYIRESID) $(SIMSARMACUSUM) $(SIMSARMAHS) \
        $(SIMSGARCHRENYIRESID) $(SIMSGARCHCUSUM) $(SIMSGARCHHS)

SIMSNORMALRENYIRESIDDF=$(subst .Rda,DataFrame.Rda,$(SIMSNORMALRENYIRESID))
SIMSNORMALCUSUMDF=$(subst .Rda,DataFrame.Rda,$(SIMSNORMALCUSUM))
SIMSNORMALHSDF=$(subst .Rda,DataFrame.Rda,$(SIMSNORMALHS))
SIMSARMARENYIRESIDDF=$(subst .Rda,DataFrame.Rda,$(SIMSARMARENYIRESID))
SIMSARMACUSUMDF=$(subst .Rda,DataFrame.Rda,$(SIMSARMACUSUM))
SIMSARMAHSDF=$(subst .Rda,DataFrame.Rda,$(SIMSARMAHS))
SIMSGARCHRENYIRESIDDF=$(subst .Rda,DataFrame.Rda,$(SIMSGARCHRENYIRESID))
SIMSGARCHCUSUMDF=$(subst .Rda,DataFrame.Rda,$(SIMSGARCHCUSUM))
SIMSGARCHHSDF=$(subst .Rda,DataFrame.Rda,$(SIMSGARCHHS))

SIMSNORMALDFPREREQ=$(SIMSNORMALRENYIRESIDDF) $(SIMSNORMALCUSUMDF) \
                   $(SIMSNORMALHSDF)
SIMSNORMALDF=data/SimsNormal.Rda
SIMSARMADFPREREQ=$(SIMSARMARENYIRESIDDF) $(SIMSARMACUSUMDF) $(SIMSARMAHSDF)
SIMSARMADF=data/SimsARMA.Rda
SIMSGARCHDFPREREQ=$(SIMSGARCHRENYIRESIDDF) $(SIMSGARCHCUSUMDF) $(SIMSGARCHHSDF)
SIMSGARCHDF=data/SimsGARCH.Rda

ALLDISTDF=$(SIMSNORMALDF) $(SIMSARMADF) $(SIMSGARCHDF)

ALLSIMSDATAFRAME=$(subst .Rda,DataFrame.Rda,$(ALLSIMS))
POWERPLOTPREFIX=vignettes/power_plot_
POWERPLOTNORMALPREFIX=$(POWERPLOTPREFIX)norm_
POWERPLOTARMAPREFIX=$(POWERPLOTPREFIX)ARMA_
POWERPLOTGARCHPREFIX=$(POWERPLOTPREFIX)GARCH_
POWERPLOTWIDTH=3
POWERPLOTHEIGHT=2
POWERPLOTLEVELLINE=dotted
POWERPLOTS=$(wildcard $(POWERPLOTNORMALPREFIX)*.pdf) \
           $(wildcard $(POWERPLOTARMAPREFIX)*.pdf) \
           $(wildcard $(POWERPLOTGARCHPREFIX)*.pdf)

.PHONY : all
all : inst/Makefile inst/package $(POWERPLOTS) $(ALLSIMSDATAFRAME)

.PRECIOUS : $(ALLSIMS)

data/$(CONTEXTPREFIX)%.Rda : exec/$(CONTEXTPREFIX)%.R R/Utils.R
	make package
	$(RSCRIPT) $< -o $@

data/$(SIMDATAPREFIX)%.Rda : exec/$(SIMDATAPREFIX)%.R R/Utils.R
	make package
	$(RSCRIPT) $< -o $@

data/$(SIMSTATPREFIX)%.Rda : exec/$(SIMSTATPREFIX)%.R R/Utils.R \
                             R/ChangePointTests.R src/ChangePointTests.cpp \
                             src/BoostMath.cpp R/ProbabilityFunctions.R \
                             R/MathFunctions.R
	make package
	$(RSCRIPT) $< -o $@

$(SIMSNORMALRENYIRESID) : data/$(CONTEXTPREFIX)Main.Rda \
                          data/$(SIMDATAPREFIX)NormalXY.Rda \
                          data/$(SIMSTATPREFIX)RenyiTypeResid.Rda
$(SIMSNORMALCUSUM) : data/$(CONTEXTPREFIX)Main.Rda \
                     data/$(SIMDATAPREFIX)NormalXY.Rda \
                     data/$(SIMSTATPREFIX)CUSUM.Rda
$(SIMSNORMALHS) : data/$(CONTEXTPREFIX)Main.Rda \
                  data/$(SIMDATAPREFIX)NormalXY.Rda \
                  data/$(SIMSTATPREFIX)HS.Rda
$(SIMSARMARENYIRESID) : data/$(CONTEXTPREFIX)Main.Rda \
                        data/$(SIMDATAPREFIX)ARMAXY.Rda \
                        data/$(SIMSTATPREFIX)RenyiTypeResid.Rda
$(SIMSARMACUSUM) : data/$(CONTEXTPREFIX)Main.Rda \
                   data/$(SIMDATAPREFIX)ARMAXY.Rda \
                   data/$(SIMSTATPREFIX)CUSUM.Rda
$(SIMSARMAHS) : data/$(CONTEXTPREFIX)Main.Rda \
                data/$(SIMDATAPREFIX)ARMAXY.Rda \
                data/$(SIMSTATPREFIX)HS.Rda
$(SIMSGARCHRENYIRESID) : data/$(CONTEXTPREFIX)Main.Rda \
                         data/$(SIMDATAPREFIX)GARCHXY.Rda \
                         data/$(SIMSTATPREFIX)RenyiTypeResid.Rda
$(SIMSGARCHCUSUM) : data/$(CONTEXTPREFIX)Main.Rda \
                    data/$(SIMDATAPREFIX)GARCHXY.Rda \
                    data/$(SIMSTATPREFIX)CUSUM.Rda
$(SIMSGARCHHS) : data/$(CONTEXTPREFIX)Main.Rda \
                 data/$(SIMDATAPREFIX)GARCHXY.Rda \
                 data/$(SIMSTATPREFIX)HS.Rda

$(ALLSIMS) : exec/PowerSimRegression.R R/Utils.R

$(ALLSIMS) :
	make package
	$(RSCRIPT) exec/PowerSimRegression.R -C $(word 1, $^) -S $(word 2, $^) \
		 -T $(word 3, $^) -o $@ -N $(POWERREPLICATIONS) -v \
		 -s $(SIMSEED)$(shell echo $@ $^ | md5sum | grep -Eo "[[:digit:]]{3,9}" | head -n1)

$(SIMSNORMALRENYIRESIDDF) : $(SIMSNORMALRENYIRESID) \
                            data/$(SIMSTATPREFIX)RenyiTypeResid.Rda \
                            data/$(CONTEXTPREFIX)Main.Rda
$(SIMSNORMALCUSUMDF) : $(SIMSNORMALCUSUM) data/$(SIMSTATPREFIX)CUSUM.Rda \
                       data/$(CONTEXTPREFIX)Main.Rda
$(SIMSNORMALHSDF) : $(SIMSNORMALHS) data/$(SIMSTATPREFIX)HS.Rda \
                    data/$(CONTEXTPREFIX)Main.Rda
$(SIMSARMARENYIRESIDDF) : $(SIMSARMARENYIRESID) \
                          data/$(SIMSTATPREFIX)RenyiTypeResid.Rda \
                          data/$(CONTEXTPREFIX)Main.Rda
$(SIMSARMACUSUMDF) : $(SIMSARMACUSUM) data/$(SIMSTATPREFIX)CUSUM.Rda \
                     data/$(CONTEXTPREFIX)Main.Rda
$(SIMSARMAHSDF) : $(SIMSARMAHS) data/$(SIMSTATPREFIX)HS.Rda \
                  data/$(CONTEXTPREFIX)Main.Rda
$(SIMSGARCHRENYIRESIDDF) : $(SIMSGARCHRENYIRESID) \
                           data/$(SIMSTATPREFIX)RenyiTypeResid.Rda \
                           data/$(CONTEXTPREFIX)Main.Rda
$(SIMSGARCHCUSUMDF) : $(SIMSGARCHCUSUM) data/$(SIMSTATPREFIX)CUSUM.Rda \
                      data/$(CONTEXTPREFIX)Main.Rda
$(SIMSGARCHHSDF) : $(SIMSGARCHHS) data/$(SIMSTATPREFIX)HS.Rda \
                   data/$(CONTEXTPREFIX)Main.Rda

$(ALLSIMSDATAFRAME) : exec/Aggregator.R R/Utils.R

$(ALLSIMSDATAFRAME) :
	make package
	$(RSCRIPT) exec/Aggregator.R -i $(word 1, $^) -o $@ -a $(LEVEL) \
		 -T $(word 2, $^) -C $(word 3, $^)

$(SIMSNORMALDF) : $(SIMSNORMALDFPREREQ)
$(SIMSARMADF) : $(SIMSARMADFPREREQ)
$(SIMSGARCHDF) : $(SIMSGARCHDFPREREQ)
$(ALLDISTDF) : exec/Appender.R

$(ALLDISTDF) :
	make package
	$(RSCRIPT) exec/Appender.R -i $(filter-out exec/Appender.R, $^) -o $@

$(POWERPLOTNORMALPREFIX)%.pdf : $(SIMSNORMALDF) exec/UnivariatePlotter.R \
                                R/Utils.R
	make package
	$(RSCRIPT) exec/UnivariatePlotter.R $< -p $(POWERPLOTNORMALPREFIX) \
		 --width $(POWERPLOTWIDTH) --height $(POWERPLOTHEIGHT) \
		 -l $(LEVEL) --levellinetype $(POWERPLOTLEVELLINE)

$(POWERPLOTARMAPREFIX)%.pdf : $(SIMSARMADF) exec/UnivariatePlotter.R \
                              R/Utils.R
	make package
	$(RSCRIPT) exec/UnivariatePlotter.R $< -p $(POWERPLOTARMAPREFIX) \
		 --width $(POWERPLOTWIDTH) --height $(POWERPLOTHEIGHT) \
		 -l $(LEVEL) --levellinetype $(POWERPLOTLEVELLINE)

$(POWERPLOTGARCHPREFIX)%.pdf : $(SIMSGARCHDF) exec/UnivariatePlotter.R \
                               R/Utils.R
	make package
	$(RSCRIPT) exec/UnivariatePlotter.R $< -p $(POWERPLOTGARCHPREFIX) \
		 --width $(POWERPLOTWIDTH) --height $(POWERPLOTHEIGHT) \
		 -l $(LEVEL) --levellinetype $(POWERPLOTLEVELLINE)

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

.PHONY : mostlyclean
mostlyclean :
	-rm $(CONTEXTDATA)
	-rm $(SIMDATADATA)
	-rm $(SIMSTATDATA)
	-rm $(ALLSIMSDATAFRAME)
	-rm $(ALLDISTDF)
	-rm $(subst .pdf,.tex,$(POWERPLOTS))
	-rm $(POWERPLOTS)

.PHONY : clean
clean :
	make mostlyclean
	-rm $(ALLSIMS)

.PHONY : init
init :
	make simconfig
	echo "I'm empty for now" > $(POWERPLOTNORMALPREFIX)n50.pdf
	echo "I'm empty for now" > $(POWERPLOTARMAPREFIX)n50.pdf
	echo "I'm empty for now" > $(POWERPLOTGARCHPREFIX)n50.pdf

.PHONY : dependencies
dependencies :
	$(RSCRIPT) exec/GetPackages.R

.PHONY : small
small :
	make POWERREPLICATIONS=$(SMALLREPLICATIONS)
