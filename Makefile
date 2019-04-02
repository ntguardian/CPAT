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

SIMSNORMALRENYI=data/SimsNormalRenyiReg.Rda
SIMSNORMALRENYIRESID=data/SimsNormalRenyiResid.Rda
SIMSNORMALRENYINOKERN=data/SimsNormalRenyiRegNoKern.Rda
SIMSNORMALCUSUM=data/SimsNormalCUSUM.Rda
SIMSNORMALHS=data/SimsNormalHS.Rda
SIMSNORMALHETERORENYI=data/SimsNormalHeteroRenyiReg.Rda
SIMSNORMALHETERORENYINOKERN=data/SimsNormalHeteroRenyiRegNoKern.Rda
SIMSNORMALHETERORENYIRESID=data/SimsNormalHeteroRenyiResid.Rda
SIMSNORMALHETEROCUSUM=data/SimsNormalHeteroCUSUM.Rda
SIMSNORMALHETEROHS=data/SimsNormalHeteroHS.Rda
SIMSARMARENYI=data/SimsARMARenyiReg.Rda
SIMSARMARENYIRESID=data/SimsARMARenyiResid.Rda
SIMSARMACUSUM=data/SimsARMACUSUM.Rda
SIMSARMAHS=data/SimsARMAHS.Rda
SIMSARMAHETERORENYI=data/SimsARMAHeteroRenyiReg.Rda
SIMSARMAHETERORENYIRESID=data/SimsARMAHeteroRenyiResid.Rda
SIMSARMAHETEROCUSUM=data/SimsARMAHeteroCUSUM.Rda
SIMSARMAHETEROHS=data/SimsARMAHeteroHS.Rda
SIMSGARCHRENYI=data/SimsGARCHRenyiReg.Rda
SIMSGARCHRENYIRESID=data/SimsGARCHRenyiResid.Rda
SIMSGARCHRENYINOKERN=data/SimsGARCHRenyiRegNoKern.Rda
SIMSGARCHCUSUM=data/SimsGARCHCUSUM.Rda
SIMSGARCHHS=data/SimsGARCHHS.Rda
SIMSGARCHHETERORENYI=data/SimsGARCHHeteroRenyiReg.Rda
SIMSGARCHHETERORENYINOKERN=data/SimsGARCHHeteroRenyiRegNoKern.Rda
SIMSGARCHHETERORENYIRESID=data/SimsGARCHHeteroRenyiResid.Rda
SIMSGARCHHETEROCUSUM=data/SimsGARCHHeteroCUSUM.Rda
SIMSGARCHHETEROHS=data/SimsGARCHHeteroHS.Rda
ALLSIMS=$(SIMSNORMALRENYI) $(SIMSNORMALRENYIRESID) $(SIMSNORMALRENYINOKERN) \
        $(SIMSNORMALCUSUM) $(SIMSNORMALHS) $(SIMSNORMALHETERORENYI) \
        $(SIMSNORMALHETERORENYINOKERN) $(SIMSNORMALHETERORENYIRESID) \
        $(SIMSNORMALHETEROCUSUM) $(SIMSNORMALHETEROHS) $(SIMSARMARENYI) \
        $(SIMSARMARENYIRESID) $(SIMSARMACUSUM) $(SIMSARMAHS) \
        $(SIMSARMAHETERORENYI) $(SIMSARMAHETERORENYIRESID) \
        $(SIMSARMAHETEROCUSUM) $(SIMSARMAHETEROHS) $(SIMSGARCHRENYI) \
        $(SIMSGARCHRENYIRESID) $(SIMSGARCHRENYINOKERN) $(SIMSGARCHCUSUM) \
        $(SIMSGARCHHS) $(SIMSGARCHHETERORENYI) $(SIMSGARCHHETERORENYINOKERN) \
        $(SIMSGARCHHETERORENYIRESID) $(SIMSGARCHHETEROCUSUM) \
        $(SIMSGARCHHETEROHS)

SIMSNORMALRENYIDF=$(subst .Rda,DataFrame.Rda,$(SIMSNORMALRENYI))
SIMSNORMALRENYIRESIDDF=$(subst .Rda,DataFrame.Rda,$(SIMSNORMALRENYIRESID))
SIMSNORMALRENYINOKERNDF=$(subst .Rda,DataFrame.Rda,$(SIMSNORMALRENYINOKERN))
SIMSNORMALCUSUMDF=$(subst .Rda,DataFrame.Rda,$(SIMSNORMALCUSUM))
SIMSNORMALHSDF=$(subst .Rda,DataFrame.Rda,$(SIMSNORMALHS))
SIMSNORMALHETERORENYIDF=$(subst .Rda,DataFrame.Rda,$(SIMSNORMALHETERORENYI))
SIMSNORMALHETERORENYINOKERNDF=$(subst .Rda,DataFrame.Rda, \
                                      $(SIMSNORMALHETERORENYINOKERN))
SIMSNORMALHETERORENYIRESIDDF=$(subst .Rda,DataFrame.Rda, \
                                     $(SIMSNORMALHETERORENYIRESID))
SIMSNORMALHETEROCUSUMDF=$(subst .Rda,DataFrame.Rda,$(SIMSNORMALHETEROCUSUM))
SIMSNORMALHETEROHSDF=$(subst .Rda,DataFrame.Rda,$(SIMSNORMALHETEROHS))
SIMSARMARENYIDF=$(subst .Rda,DataFrame.Rda, $(SIMSARMARENYI))
SIMSARMARENYIRESIDDF=$(subst .Rda,DataFrame.Rda,$(SIMSARMARENYIRESID))
SIMSARMACUSUMDF=$(subst .Rda,DataFrame.Rda,$(SIMSARMACUSUM))
SIMSARMAHSDF=$(subst .Rda,DataFrame.Rda,$(SIMSARMAHS))
SIMSARMAHETERORENYIDF=$(subst .Rda,DataFrame.Rda,$(SIMSARMAHETERORENYI))
SIMSARMAHETERORENYIRESIDDF=$(subst .Rda,DataFrame.Rda, \
                                     $(SIMSARMAHETERORENYIRESID))
SIMSARMAHETEROCUSUMDF=$(subst .Rda,DataFrame.Rda,$(SIMSARMAHETEROCUSUM))
SIMSARMAHETEROHSDF=$(subst .Rda,DataFrame.Rda,$(SIMSARMAHETEROHS))
SIMSGARCHRENYIDF=$(subst .Rda,DataFrame.Rda,$(SIMSGARCHRENYI))
SIMSGARCHRENYIRESIDDF=$(subst .Rda,DataFrame.Rda,$(SIMSGARCHRENYIRESID))
SIMSGARCHRENYINOKERNDF=$(subst .Rda,DataFrame.Rda,$(SIMSGARCHRENYINOKERN))
SIMSGARCHCUSUMDF=$(subst .Rda,DataFrame.Rda,$(SIMSGARCHCUSUM))
SIMSGARCHHSDF=$(subst .Rda,DataFrame.Rda,$(SIMSGARCHHS))
SIMSGARCHHETERORENYIDF=$(subst .Rda,DataFrame.Rda,$(SIMSGARCHHETERORENYI))
SIMSGARCHHETERORENYINOKERNDF=$(subst .Rda,DataFrame.Rda, \
                                     $(SIMSGARCHHETERORENYINOKERN))
SIMSGARCHHETERORENYIRESIDDF=$(subst .Rda,DataFrame.Rda, \
                                     $(SIMSGARCHHETERORENYIRESID))
SIMSGARCHHETEROCUSUMDF=$(subst .Rda,DataFrame.Rda,$(SIMSGARCHHETEROCUSUM))
SIMSGARCHHETEROHSDF=$(subst .Rda,DataFrame.Rda,$(SIMSGARCHHETEROHS))

SIMSNORMALDFPREREQ=$(SIMSNORMALRENYIRESIDDF) $(SIMSNORMALRENYINOKERNDF) \
                   $(SIMSNORMALCUSUMDF) $(SIMSNORMALHSDF)
SIMSNORMALDF=data/SimsNormal.Rda
SIMSNORMALHETERODFPREREQ=$(SIMSNORMALHETERORENYINOKERNDF) \
                         $(SIMSNORMALHETERORENYIRESIDDF) \
                         $(SIMSNORMALHETEROCUSUMDF) $(SIMSNORMALHETEROHSDF)
SIMSNORMALHETERODF=data/SimsNormalHetero.Rda
SIMSARMADFPREREQ=$(SIMSARMARENYIRESIDDF) $(SIMSARMARENYIDF) $(SIMSARMACUSUMDF) \
                 $(SIMSARMAHSDF)
SIMSARMADF=data/SimsARMA.Rda
SIMSARMAHETERODFPREREQ=$(SIMSARMAHETERORENYIRESIDDF) $(SIMSARMAHETERORENYIDF) \
                       $(SIMSARMAHETEROCUSUMDF) $(SIMSARMAHETEROHSDF)
SIMSARMAHETERODF=data/SimsARMAHetero.Rda
SIMSGARCHDFPREREQ=$(SIMSGARCHRENYIRESIDDF) $(SIMSGARCHRENYINOKERNDF) \
                  $(SIMSGARCHCUSUMDF) $(SIMSGARCHHSDF)
SIMSGARCHDF=data/SimsGARCH.Rda
SIMSGARCHHETERODFPREREQ=$(SIMSGARCHHETERORENYIRESIDDF) \
                        $(SIMSGARCHHETERORENYINOKERNDF) \
                        $(SIMSGARCHHETEROCUSUMDF) $(SIMSGARCHHETEROHSDF)
SIMSGARCHHETERODF=data/SimsGARCHHetero.Rda

ALLDISTDF=$(SIMSNORMALDF) $(SIMSNORMALHETERODF) $(SIMSARMADF) \
          $(SIMSARMAHETERODF) $(SIMSGARCHDF) $(SIMSGARCHHETERODF)

ALLSIMSDATAFRAME=$(subst .Rda,DataFrame.Rda,$(ALLSIMS))
POWERPLOTPREFIX=vignettes/power_plot_
POWERPLOTNORMALPREFIX=$(POWERPLOTPREFIX)norm_
POWERPLOTARMAPREFIX=$(POWERPLOTPREFIX)ARMA_
POWERPLOTGARCHPREFIX=$(POWERPLOTPREFIX)GARCH_
POWERPLOTNORMALHETEROPREFIX=$(POWERPLOTPREFIX)norm_hetero_
POWERPLOTARMAHETEROPREFIX=$(POWERPLOTPREFIX)ARMA_hetero_
POWERPLOTGARCHHETEROPREFIX=$(POWERPLOTPREFIX)GARCH_hetero_
POWERPLOTWIDTH=3
POWERPLOTHEIGHT=2
POWERPLOTLEVELLINE=dotted
POWERPLOTS=$(wildcard $(POWERPLOTNORMALPREFIX)*.pdf) \
           $(wildcard $(POWERPLOTNORMALHETEROPREFIX)*.pdf) \
           $(wildcard $(POWERPLOTARMAPREFIX)*.pdf) \
           $(wildcard $(POWERPLOTARMAHETEROPREFIX)*.pdf) \
           $(wildcard $(POWERPLOTGARCHPREFIX)*.pdf) \
           $(wildcard $(POWERPLOTGARCHHETEROPREFIX)*.pdf)

VIGNETTES=doc/CollectedPlots.pdf

.PHONY : all
all : inst/Makefile inst/package $(POWERPLOTS) $(ALLSIMSDATAFRAME) $(VIGNETTES)

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

$(SIMSNORMALRENYI) : data/$(CONTEXTPREFIX)Main.Rda \
                     data/$(SIMDATAPREFIX)NormalXY.Rda \
                     data/$(SIMSTATPREFIX)RenyiTypeReg.Rda
$(SIMSNORMALRENYIRESID) : data/$(CONTEXTPREFIX)Main.Rda \
                          data/$(SIMDATAPREFIX)NormalXY.Rda \
                          data/$(SIMSTATPREFIX)RenyiTypeResid.Rda
$(SIMSNORMALRENYINOKERN) : data/$(CONTEXTPREFIX)Main.Rda \
                           data/$(SIMDATAPREFIX)NormalXY.Rda \
                           data/$(SIMSTATPREFIX)RenyiTypeRegNoKern.Rda
$(SIMSNORMALCUSUM) : data/$(CONTEXTPREFIX)Main.Rda \
                     data/$(SIMDATAPREFIX)NormalXY.Rda \
                     data/$(SIMSTATPREFIX)CUSUM.Rda
$(SIMSNORMALHS) : data/$(CONTEXTPREFIX)Main.Rda \
                  data/$(SIMDATAPREFIX)NormalXY.Rda \
                  data/$(SIMSTATPREFIX)HS.Rda
$(SIMSNORMALHETERORENYI) : data/$(CONTEXTPREFIX)Main.Rda \
                           data/$(SIMDATAPREFIX)NormalXYHetero.Rda \
                           data/$(SIMSTATPREFIX)RenyiTypeReg.Rda
$(SIMSNORMALHETERORENYINOKERN) : data/$(CONTEXTPREFIX)Main.Rda \
                                 data/$(SIMDATAPREFIX)NormalXYHetero.Rda \
                                 data/$(SIMSTATPREFIX)RenyiTypeRegNoKern.Rda
$(SIMSNORMALHETERORENYIRESID) : data/$(CONTEXTPREFIX)Main.Rda \
                                data/$(SIMDATAPREFIX)NormalXYHetero.Rda \
                                data/$(SIMSTATPREFIX)RenyiTypeResid.Rda
$(SIMSNORMALHETEROCUSUM) : data/$(CONTEXTPREFIX)Main.Rda \
                     data/$(SIMDATAPREFIX)NormalXYHetero.Rda \
                     data/$(SIMSTATPREFIX)CUSUM.Rda
$(SIMSNORMALHETEROHS) : data/$(CONTEXTPREFIX)Main.Rda \
                        data/$(SIMDATAPREFIX)NormalXYHetero.Rda \
                        data/$(SIMSTATPREFIX)HS.Rda
$(SIMSARMARENYI) : data/$(CONTEXTPREFIX)Main.Rda \
                   data/$(SIMDATAPREFIX)ARMAXY.Rda \
                   data/$(SIMSTATPREFIX)RenyiTypeReg.Rda
$(SIMSARMARENYIRESID) : data/$(CONTEXTPREFIX)Main.Rda \
                        data/$(SIMDATAPREFIX)ARMAXY.Rda \
                        data/$(SIMSTATPREFIX)RenyiTypeResid.Rda
$(SIMSARMACUSUM) : data/$(CONTEXTPREFIX)Main.Rda \
                   data/$(SIMDATAPREFIX)ARMAXY.Rda \
                   data/$(SIMSTATPREFIX)CUSUM.Rda
$(SIMSARMAHS) : data/$(CONTEXTPREFIX)Main.Rda \
                data/$(SIMDATAPREFIX)ARMAXY.Rda \
                data/$(SIMSTATPREFIX)HS.Rda
$(SIMSARMAHETERORENYI) : data/$(CONTEXTPREFIX)Main.Rda \
                         data/$(SIMDATAPREFIX)ARMAXYHetero.Rda \
                         data/$(SIMSTATPREFIX)RenyiTypeReg.Rda
$(SIMSARMAHETERORENYIRESID) : data/$(CONTEXTPREFIX)Main.Rda \
                              data/$(SIMDATAPREFIX)ARMAXYHetero.Rda \
                              data/$(SIMSTATPREFIX)RenyiTypeResid.Rda
$(SIMSARMAHETEROCUSUM) : data/$(CONTEXTPREFIX)Main.Rda \
                         data/$(SIMDATAPREFIX)ARMAXYHetero.Rda \
                         data/$(SIMSTATPREFIX)CUSUM.Rda
$(SIMSARMAHETEROHS) : data/$(CONTEXTPREFIX)Main.Rda \
                      data/$(SIMDATAPREFIX)ARMAXYHetero.Rda \
                      data/$(SIMSTATPREFIX)HS.Rda
$(SIMSGARCHRENYI) : data/$(CONTEXTPREFIX)Main.Rda \
                    data/$(SIMDATAPREFIX)GARCHXY.Rda \
                    data/$(SIMSTATPREFIX)RenyiTypeReg.Rda
$(SIMSGARCHRENYIRESID) : data/$(CONTEXTPREFIX)Main.Rda \
                         data/$(SIMDATAPREFIX)GARCHXY.Rda \
                         data/$(SIMSTATPREFIX)RenyiTypeResid.Rda
$(SIMSGARCHRENYINOKERN) : data/$(CONTEXTPREFIX)Main.Rda \
                           data/$(SIMDATAPREFIX)GARCHXY.Rda \
                           data/$(SIMSTATPREFIX)RenyiTypeRegNoKern.Rda
$(SIMSGARCHCUSUM) : data/$(CONTEXTPREFIX)Main.Rda \
                    data/$(SIMDATAPREFIX)GARCHXY.Rda \
                    data/$(SIMSTATPREFIX)CUSUM.Rda
$(SIMSGARCHHS) : data/$(CONTEXTPREFIX)Main.Rda \
                 data/$(SIMDATAPREFIX)GARCHXY.Rda \
                 data/$(SIMSTATPREFIX)HS.Rda
$(SIMSGARCHHETERORENYI) : data/$(CONTEXTPREFIX)Main.Rda \
                          data/$(SIMDATAPREFIX)GARCHXYHetero.Rda \
                          data/$(SIMSTATPREFIX)RenyiTypeReg.Rda
$(SIMSGARCHHETERORENYIRESID) : data/$(CONTEXTPREFIX)Main.Rda \
                               data/$(SIMDATAPREFIX)GARCHXYHetero.Rda \
                               data/$(SIMSTATPREFIX)RenyiTypeResid.Rda
$(SIMSGARCHHETERORENYINOKERN) : data/$(CONTEXTPREFIX)Main.Rda \
                                 data/$(SIMDATAPREFIX)GARCHXYHetero.Rda \
                                 data/$(SIMSTATPREFIX)RenyiTypeRegNoKern.Rda
$(SIMSGARCHHETEROCUSUM) : data/$(CONTEXTPREFIX)Main.Rda \
                          data/$(SIMDATAPREFIX)GARCHXYHetero.Rda \
                          data/$(SIMSTATPREFIX)CUSUM.Rda
$(SIMSGARCHHETEROHS) : data/$(CONTEXTPREFIX)Main.Rda \
                       data/$(SIMDATAPREFIX)GARCHXYHetero.Rda \
                       data/$(SIMSTATPREFIX)HS.Rda

$(ALLSIMS) : exec/PowerSimRegression.R R/Utils.R

$(ALLSIMS) :
	make package
	$(RSCRIPT) exec/PowerSimRegression.R -C $(word 1, $^) -S $(word 2, $^) \
		 -T $(word 3, $^) -o $@ -N $(POWERREPLICATIONS) -v \
		 -s $(SIMSEED)$(shell echo $@ $^ | md5sum | grep -Eo "[[:digit:]]{3,6}" | head -n1)

$(SIMSNORMALRENYIDF) : $(SIMSNORMALRENYI) \
                       data/$(SIMSTATPREFIX)RenyiTypeReg.Rda \
                       data/$(CONTEXTPREFIX)Main.Rda
$(SIMSNORMALRENYIRESIDDF) : $(SIMSNORMALRENYIRESID) \
                            data/$(SIMSTATPREFIX)RenyiTypeResid.Rda \
                            data/$(CONTEXTPREFIX)Main.Rda
$(SIMSNORMALRENYINOKERNDF) : $(SIMSNORMALRENYINOKERN) \
                             data/$(SIMSTATPREFIX)RenyiTypeRegNoKern.Rda \
                             data/$(CONTEXTPREFIX)Main.Rda
$(SIMSNORMALCUSUMDF) : $(SIMSNORMALCUSUM) data/$(SIMSTATPREFIX)CUSUM.Rda \
                       data/$(CONTEXTPREFIX)Main.Rda
$(SIMSNORMALHSDF) : $(SIMSNORMALHS) data/$(SIMSTATPREFIX)HS.Rda \
                    data/$(CONTEXTPREFIX)Main.Rda
$(SIMSNORMALHETERORENYIDF) : $(SIMSNORMALHETERORENYI) \
                             data/$(SIMSTATPREFIX)RenyiTypeReg.Rda \
                             data/$(CONTEXTPREFIX)Main.Rda
$(SIMSNORMALHETERORENYINOKERNDF) : $(SIMSNORMALHETERORENYINOKERN) \
                                   data/$(SIMSTATPREFIX)RenyiTypeReg.Rda \
                                   data/$(CONTEXTPREFIX)Main.Rda
$(SIMSNORMALHETERORENYIRESIDDF) : $(SIMSNORMALHETERORENYIRESID) \
                                  data/$(SIMSTATPREFIX)RenyiTypeResid.Rda \
                                  data/$(CONTEXTPREFIX)Main.Rda
$(SIMSNORMALHETEROCUSUMDF) : $(SIMSNORMALHETEROCUSUM) \
                             data/$(SIMSTATPREFIX)CUSUM.Rda \
                             data/$(CONTEXTPREFIX)Main.Rda
$(SIMSNORMALHETEROHSDF) : $(SIMSNORMALHETEROHS) data/$(SIMSTATPREFIX)HS.Rda \
                          data/$(CONTEXTPREFIX)Main.Rda
$(SIMSARMARENYIDF) : $(SIMSARMARENYI) \
                     data/$(SIMSTATPREFIX)RenyiTypeReg.Rda \
                     data/$(CONTEXTPREFIX)Main.Rda
$(SIMSARMARENYIRESIDDF) : $(SIMSARMARENYIRESID) \
                          data/$(SIMSTATPREFIX)RenyiTypeResid.Rda \
                          data/$(CONTEXTPREFIX)Main.Rda
$(SIMSARMACUSUMDF) : $(SIMSARMACUSUM) data/$(SIMSTATPREFIX)CUSUM.Rda \
                     data/$(CONTEXTPREFIX)Main.Rda
$(SIMSARMAHSDF) : $(SIMSARMAHS) data/$(SIMSTATPREFIX)HS.Rda \
                  data/$(CONTEXTPREFIX)Main.Rda
$(SIMSARMAHETERORENYIDF) : $(SIMSARMAHETERORENYI) \
                           data/$(SIMSTATPREFIX)RenyiTypeReg.Rda \
                           data/$(CONTEXTPREFIX)Main.Rda
$(SIMSARMAHETERORENYIRESIDDF) : $(SIMSARMAHETERORENYIRESID) \
                                data/$(SIMSTATPREFIX)RenyiTypeResid.Rda \
                                data/$(CONTEXTPREFIX)Main.Rda
$(SIMSARMAHETEROCUSUMDF) : $(SIMSARMAHETEROCUSUM) \
                           data/$(SIMSTATPREFIX)CUSUM.Rda \
                           data/$(CONTEXTPREFIX)Main.Rda
$(SIMSARMAHETEROHSDF) : $(SIMSARMAHETEROHS) data/$(SIMSTATPREFIX)HS.Rda \
                        data/$(CONTEXTPREFIX)Main.Rda
$(SIMSGARCHRENYIDF) : $(SIMSGARCHRENYI) \
                      data/$(SIMSTATPREFIX)RenyiTypeReg.Rda \
                      data/$(CONTEXTPREFIX)Main.Rda
$(SIMSGARCHRENYIRESIDDF) : $(SIMSGARCHRENYIRESID) \
                           data/$(SIMSTATPREFIX)RenyiTypeResid.Rda \
                           data/$(CONTEXTPREFIX)Main.Rda
$(SIMSGARCHRENYINOKERNDF) : $(SIMSGARCHRENYINOKERN) \
                             data/$(SIMSTATPREFIX)RenyiTypeRegNoKern.Rda \
                             data/$(CONTEXTPREFIX)Main.Rda
$(SIMSGARCHCUSUMDF) : $(SIMSGARCHCUSUM) data/$(SIMSTATPREFIX)CUSUM.Rda \
                      data/$(CONTEXTPREFIX)Main.Rda
$(SIMSGARCHHSDF) : $(SIMSGARCHHS) data/$(SIMSTATPREFIX)HS.Rda \
                   data/$(CONTEXTPREFIX)Main.Rda
$(SIMSGARCHHETERORENYIDF) : $(SIMSGARCHHETERORENYI) \
                            data/$(SIMSTATPREFIX)RenyiTypeReg.Rda \
                            data/$(CONTEXTPREFIX)Main.Rda
$(SIMSGARCHHETERORENYINOKERNDF) : $(SIMSGARCHHETERORENYINOKERN) \
                                  data/$(SIMSTATPREFIX)RenyiTypeReg.Rda \
                                  data/$(CONTEXTPREFIX)Main.Rda
$(SIMSGARCHHETERORENYIRESIDDF) : $(SIMSGARCHHETERORENYIRESID) \
                                 data/$(SIMSTATPREFIX)RenyiTypeResid.Rda \
                                 data/$(CONTEXTPREFIX)Main.Rda
$(SIMSGARCHHETEROCUSUMDF) : $(SIMSGARCHHETEROCUSUM) \
                            data/$(SIMSTATPREFIX)CUSUM.Rda \
                            data/$(CONTEXTPREFIX)Main.Rda
$(SIMSGARCHHETEROHSDF) : $(SIMSGARCHHETEROHS) data/$(SIMSTATPREFIX)HS.Rda \
                         data/$(CONTEXTPREFIX)Main.Rda

$(ALLSIMSDATAFRAME) : exec/Aggregator.R R/Utils.R

$(ALLSIMSDATAFRAME) :
	make package
	$(RSCRIPT) exec/Aggregator.R -i $(word 1, $^) -o $@ -a $(LEVEL) \
		 -T $(word 2, $^) -C $(word 3, $^)

$(SIMSNORMALDF) : $(SIMSNORMALDFPREREQ)
$(SIMSNORMALHETERODF) : $(SIMSNORMALHETERODFPREREQ)
$(SIMSARMADF) : $(SIMSARMADFPREREQ)
$(SIMSARMAHETERODF) : $(SIMSARMAHETERODFPREREQ)
$(SIMSGARCHDF) : $(SIMSGARCHDFPREREQ)
$(SIMSGARCHHETERODF) : $(SIMSGARCHHETERODFPREREQ)
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

$(POWERPLOTNORMALHETEROPREFIX)%.pdf : $(SIMSNORMALHETERODF) \
                                      exec/UnivariatePlotter.R R/Utils.R
	make package
	$(RSCRIPT) exec/UnivariatePlotter.R $< -p $(POWERPLOTNORMALHETEROPREFIX) \
		 --width $(POWERPLOTWIDTH) --height $(POWERPLOTHEIGHT) \
		 -l $(LEVEL) --levellinetype $(POWERPLOTLEVELLINE)

$(POWERPLOTARMAPREFIX)%.pdf : $(SIMSARMADF) exec/UnivariatePlotter.R \
                              R/Utils.R
	make package
	$(RSCRIPT) exec/UnivariatePlotter.R $< -p $(POWERPLOTARMAPREFIX) \
		 --width $(POWERPLOTWIDTH) --height $(POWERPLOTHEIGHT) \
		 -l $(LEVEL) --levellinetype $(POWERPLOTLEVELLINE)

$(POWERPLOTARMAHETEROPREFIX)%.pdf : $(SIMSARMAHETERODF) \
                                    exec/UnivariatePlotter.R R/Utils.R
	make package
	$(RSCRIPT) exec/UnivariatePlotter.R $< -p $(POWERPLOTARMAHETEROPREFIX) \
		 --width $(POWERPLOTWIDTH) --height $(POWERPLOTHEIGHT) \
		 -l $(LEVEL) --levellinetype $(POWERPLOTLEVELLINE)

$(POWERPLOTGARCHPREFIX)%.pdf : $(SIMSGARCHDF) exec/UnivariatePlotter.R \
                               R/Utils.R
	make package
	$(RSCRIPT) exec/UnivariatePlotter.R $< -p $(POWERPLOTGARCHPREFIX) \
		 --width $(POWERPLOTWIDTH) --height $(POWERPLOTHEIGHT) \
		 -l $(LEVEL) --levellinetype $(POWERPLOTLEVELLINE)

$(POWERPLOTGARCHHETEROPREFIX)%.pdf : $(SIMSGARCHHETERODF) \
                                     exec/UnivariatePlotter.R R/Utils.R
	make package
	$(RSCRIPT) exec/UnivariatePlotter.R $< -p $(POWERPLOTGARCHHETEROPREFIX) \
		 --width $(POWERPLOTWIDTH) --height $(POWERPLOTHEIGHT) \
		 -l $(LEVEL) --levellinetype $(POWERPLOTLEVELLINE)

doc/CollectedPlots.pdf : vignettes/CollectedPlots.ltx $(POWERPLOTS)

$(VIGNETTES) :
	$(RSCRIPT) -e "devtools::build_vignettes()"

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
	echo "I'm empty for now" > $(POWERPLOTNORMALHETEROPREFIX)n50.pdf
	echo "I'm empty for now" > $(POWERPLOTARMAPREFIX)n50.pdf
	echo "I'm empty for now" > $(POWERPLOTARMAHETEROPREFIX)n50.pdf
	echo "I'm empty for now" > $(POWERPLOTGARCHPREFIX)n50.pdf
	echo "I'm empty for now" > $(POWERPLOTGARCHHETEROPREFIX)n50.pdf

.PHONY : dependencies
dependencies :
	$(RSCRIPT) exec/GetPackages.R

.PHONY : small
small :
	make POWERREPLICATIONS=$(SMALLREPLICATIONS)
