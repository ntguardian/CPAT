################################################################################
# Makefile
################################################################################
# Curtis Miller
# 2019-04-05
################################################################################

################################################################################
# CONTENTS
################################################################################
# VARIABLES
# EXECUTION-VARS
# SIMULATION-DATA-VARS
# POWER-DF-VARS
# COLLECTION-DF-VARS
# POWER-PLOT-VARS
# VIGNETTE-VARS
# EXAMPLE-VARS
# RECIPES
# CONFIG-RECIPES
# SIMULATION-DATA-DEPENDS
# SIMULATION-DATA-RECIPE
# POWER-DF-DEPENDS
# POWER-DF-RECIPE
# COLLECTION-DF-DEPENDS
# COLLECTION-DF-RECIPE
# POWER-PLOT-RECIPE
# EXAMPLE-RECIPES
# VIGNETTE-DEPENDS-RECIPE
# MISC-RECIPE
################################################################################

################################################################################
# VARIABLES
################################################################################

################################################################################
# EXECUTION-VARS
################################################################################

# Makefile execution variables
SHELL=/bin/sh
RSCRIPT=Rscript

# Testing replication count
SMALLREPLICATIONS=20

# Simulation definition files and constants for simulation
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

################################################################################
# SIMULATION-DATA-VARS
################################################################################

# Files containing simulated test statistics
# Normal
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
SIMSPROPNORMALRENYI=data/SimsPropNormalRenyiReg.Rda
SIMSPROPNORMALRENYIRESID=data/SimsPropNormalRenyiResid.Rda
SIMSPROPNORMALRENYINOKERN=data/SimsPropNormalRenyiRegNoKern.Rda
SIMSPROPNORMALCUSUM=data/SimsPropNormalCUSUM.Rda
SIMSPROPNORMALHS=data/SimsPropNormalHS.Rda
SIMSPROPNORMALHETERORENYI=data/SimsPropNormalHeteroRenyiReg.Rda
SIMSPROPNORMALHETERORENYIRESID=data/SimsPropNormalHeteroRenyiResid.Rda
SIMSPROPNORMALHETERORENYINOKERN=data/SimsPropNormalHeteroRenyiRegNoKern.Rda
SIMSPROPNORMALHETEROCUSUM=data/SimsPropNormalHeteroCUSUM.Rda
SIMSPROPNORMALHETEROHS=data/SimsPropNormalHeteroHS.Rda
# ARMA
SIMSARMARENYI=data/SimsARMARenyiReg.Rda
SIMSARMARENYIRESID=data/SimsARMARenyiResid.Rda
SIMSARMACUSUM=data/SimsARMACUSUM.Rda
SIMSARMAHS=data/SimsARMAHS.Rda
SIMSARMAHETERORENYI=data/SimsARMAHeteroRenyiReg.Rda
SIMSARMAHETERORENYIRESID=data/SimsARMAHeteroRenyiResid.Rda
SIMSARMAHETEROCUSUM=data/SimsARMAHeteroCUSUM.Rda
SIMSARMAHETEROHS=data/SimsARMAHeteroHS.Rda
# GARCH
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
# All
NOTHETEROSIMS=$(SIMSNORMALRENYI) $(SIMSNORMALRENYIRESID) \
              $(SIMSNORMALRENYINOKERN) $(SIMSNORMALCUSUM) $(SIMSNORMALHS) \
              $(SIMSPROPNORMALRENYI) $(SIMSPROPNORMALRENYIRESID) \
              $(SIMSPROPNORMALRENYINOKERN) $(SIMSPROPNORMALCUSUM) \
              $(SIMSPROPNORMALHS) $(SIMSARMARENYI) $(SIMSARMARENYIRESID) \
              $(SIMSARMACUSUM) $(SIMSARMAHS) $(SIMSGARCHRENYI) \
              $(SIMSGARCHRENYIRESID) $(SIMSGARCHRENYINOKERN) $(SIMSGARCHCUSUM) \
              $(SIMSGARCHHS) 
HETEROSIMS=$(SIMSNORMALHETERORENYI) $(SIMSNORMALHETERORENYINOKERN) \
           $(SIMSNORMALHETERORENYIRESID) $(SIMSNORMALHETEROCUSUM) \
           $(SIMSNORMALHETEROHS) $(SIMSPROPNORMALHETERORENYI) \
           $(SIMSPROPNORMALHETERORENYIRESID) \
           $(SIMSPROPNORMALHETERORENYINOKERN) $(SIMSPROPNORMALHETEROCUSUM) \
           $(SIMSPROPNORMALHETEROHS) $(SIMSARMAHETERORENYI) \
           $(SIMSARMAHETERORENYIRESID) $(SIMSARMAHETEROCUSUM) \
           $(SIMSARMAHETEROHS) $(SIMSGARCHHETERORENYI) \
           $(SIMSGARCHHETERORENYINOKERN) $(SIMSGARCHHETERORENYIRESID) \
           $(SIMSGARCHHETEROCUSUM) $(SIMSGARCHHETEROHS)
ALLSIMS=$(NOTHETEROSIMS) $(HETEROSIMS)

################################################################################
# POWER-DF-VARS
################################################################################

# Files containing power estimates computed from simulated test statistics
# Normal
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
SIMSPROPNORMALRENYIDF=$(subst .Rda,DataFrame.Rda,$(SIMSPROPNORMALRENYI))
SIMSPROPNORMALRENYIRESIDDF=$(subst .Rda,DataFrame.Rda, \
                           $(SIMSPROPNORMALRENYIRESID))
SIMSPROPNORMALRENYINOKERNDF=$(subst .Rda,DataFrame.Rda, \
                            $(SIMSPROPNORMALRENYINOKERN))
SIMSPROPNORMALCUSUMDF=$(subst .Rda,DataFrame.Rda, $(SIMSPROPNORMALCUSUM))
SIMSPROPNORMALHSDF=$(subst .Rda,DataFrame.Rda,$(SIMSPROPNORMALHS))
SIMSPROPNORMALHETERORENYIDF=$(subst .Rda,DataFrame.Rda, \
                            $(SIMSPROPNORMALHETERORENYI))
SIMSPROPNORMALHETERORENYIRESIDDF=$(subst .Rda,DataFrame.Rda, \
                                 $(SIMSPROPNORMALHETERORENYIRESID))
SIMSPROPNORMALHETERORENYINOKERNDF=$(subst .Rda,DataFrame.Rda, \
                                  $(SIMSPROPNORMALHETERORENYINOKERN))
SIMSPROPNORMALHETEROCUSUMDF=$(subst .Rda,DataFrame.Rda, \
                            $(SIMSPROPNORMALHETEROCUSUM))
SIMSPROPNORMALHETEROHSDF=$(subst .Rda,DataFrame.Rda, \
                         $(SIMSPROPNORMALHETEROHS))
# ARMA
SIMSARMARENYIDF=$(subst .Rda,DataFrame.Rda, $(SIMSARMARENYI))
SIMSARMARENYIRESIDDF=$(subst .Rda,DataFrame.Rda,$(SIMSARMARENYIRESID))
SIMSARMACUSUMDF=$(subst .Rda,DataFrame.Rda,$(SIMSARMACUSUM))
SIMSARMAHSDF=$(subst .Rda,DataFrame.Rda,$(SIMSARMAHS))
SIMSARMAHETERORENYIDF=$(subst .Rda,DataFrame.Rda,$(SIMSARMAHETERORENYI))
SIMSARMAHETERORENYIRESIDDF=$(subst .Rda,DataFrame.Rda, \
                                     $(SIMSARMAHETERORENYIRESID))
SIMSARMAHETEROCUSUMDF=$(subst .Rda,DataFrame.Rda,$(SIMSARMAHETEROCUSUM))
SIMSARMAHETEROHSDF=$(subst .Rda,DataFrame.Rda,$(SIMSARMAHETEROHS))
# GARCH
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
# All
ALLSIMSDATAFRAME=$(subst .Rda,DataFrame.Rda,$(ALLSIMS))

################################################################################
# COLLECTION-DF-VARS
################################################################################

# Files collecting power statistics for DGP contexts
# Normal
SIMSNORMALDFPREREQ=$(SIMSNORMALRENYIRESIDDF) $(SIMSNORMALRENYINOKERNDF) \
                   $(SIMSNORMALCUSUMDF) $(SIMSNORMALHSDF)
SIMSNORMALDF=data/SimsNormal.Rda
SIMSNORMALHETERODFPREREQ=$(SIMSNORMALHETERORENYINOKERNDF) \
                         $(SIMSNORMALHETERORENYIRESIDDF) \
                         $(SIMSNORMALHETEROCUSUMDF) $(SIMSNORMALHETEROHSDF)
SIMSNORMALHETERODF=data/SimsNormalHetero.Rda
SIMSPROPNORMALDFPREREQ=$(SIMSPROPNORMALRENYIRESIDDF) \
                       $(SIMSPROPNORMALRENYINOKERNDF) $(SIMSPROPNORMALCUSUMDF) \
                       $(SIMSPROPNORMALHSDF)
SIMSPROPNORMALDF=data/SimsPropNormal.Rda
SIMSPROPNORMALHETERODFPREREQ=$(SIMSPROPNORMALHETERORENYIRESIDDF) \
                             $(SIMSPROPNORMALHETERORENYINOKERNDF) \
                             $(SIMSPROPNORMALHETEROCUSUMDF) \
                             $(SIMSPROPNORMALHETEROHSDF)
SIMSPROPNORMALHETERODF=data/SimsPropNormalHetero.Rda
# ARMA
SIMSARMADFPREREQ=$(SIMSARMARENYIRESIDDF) $(SIMSARMARENYIDF) $(SIMSARMACUSUMDF) \
                 $(SIMSARMAHSDF)
SIMSARMADF=data/SimsARMA.Rda
SIMSARMAHETERODFPREREQ=$(SIMSARMAHETERORENYIRESIDDF) $(SIMSARMAHETERORENYIDF) \
                       $(SIMSARMAHETEROCUSUMDF) $(SIMSARMAHETEROHSDF)
SIMSARMAHETERODF=data/SimsARMAHetero.Rda
# GARCH
SIMSGARCHDFPREREQ=$(SIMSGARCHRENYIRESIDDF) $(SIMSGARCHRENYINOKERNDF) \
                  $(SIMSGARCHCUSUMDF) $(SIMSGARCHHSDF)
SIMSGARCHDF=data/SimsGARCH.Rda
SIMSGARCHHETERODFPREREQ=$(SIMSGARCHHETERORENYIRESIDDF) \
                        $(SIMSGARCHHETERORENYINOKERNDF) \
                        $(SIMSGARCHHETEROCUSUMDF) $(SIMSGARCHHETEROHSDF)
SIMSGARCHHETERODF=data/SimsGARCHHetero.Rda

# All DGP power collections in one location
ALLDISTDF=$(SIMSNORMALDF) $(SIMSNORMALHETERODF) $(SIMSPROPNORMALDF) \
          $(SIMSPROPNORMALHETERODF) $(SIMSARMADF) $(SIMSARMAHETERODF) \
          $(SIMSGARCHDF) $(SIMSGARCHHETERODF)

################################################################################
# POWER-PLOT-VARS
################################################################################

# Power plot variables, giving plots names and defining plot dimensions
POWERPLOTPREFIX=vignettes/power_plot_
POWERPLOTNORMALPREFIX=$(POWERPLOTPREFIX)norm_
POWERPLOTPROPNORMALPREFIX=$(POWERPLOTPREFIX)prop_norm_
POWERPLOTARMAPREFIX=$(POWERPLOTPREFIX)ARMA_
POWERPLOTGARCHPREFIX=$(POWERPLOTPREFIX)GARCH_
POWERPLOTNORMALHETEROPREFIX=$(POWERPLOTPREFIX)norm_hetero_
POWERPLOTPROPNORMALHETEROPREFIX=$(POWERPLOTPREFIX)prop_norm_hetero_
POWERPLOTARMAHETEROPREFIX=$(POWERPLOTPREFIX)ARMA_hetero_
POWERPLOTGARCHHETEROPREFIX=$(POWERPLOTPREFIX)GARCH_hetero_
POWERPLOTWIDTH=3
POWERPLOTHEIGHT=2
POWERPLOTLEVELLINE=dotted
POWERPLOTS=$(wildcard $(POWERPLOTNORMALPREFIX)*.pdf) \
           $(wildcard $(POWERPLOTPROPNORMALPREFIX)*.pdf) \
           $(wildcard $(POWERPLOTNORMALHETEROPREFIX)*.pdf) \
           $(wildcard $(POWERPLOTPROPNORMALHETEROPREFIX)*.pdf) \
           $(wildcard $(POWERPLOTARMAPREFIX)*.pdf) \
           $(wildcard $(POWERPLOTARMAHETEROPREFIX)*.pdf) \
           $(wildcard $(POWERPLOTGARCHPREFIX)*.pdf) \
           $(wildcard $(POWERPLOTGARCHHETEROPREFIX)*.pdf)

################################################################################
# EXAMPLE-VARS
################################################################################

# Data example definition files
EXAMPLEPREFIX=RealData
EXAMPLEGERMANM1=data/$(EXAMPLEPREFIX)M1Germany.Rda
ALLEXAMPLES=$(EXAMPLEGERMANM1)

# Processed data examples
EXAMPLEGERMANM1PVALS=data/$(EXAMPLEPREFIX)M1GermanypVals.Rda
ALLEXAMPLESPVALS=$(EXAMPLEGERMAN1PVALS)

################################################################################
# VIGNETTE-VARS
################################################################################

# Vignettes to be created
VIGNETTES=doc/CollectedPlots.pdf

################################################################################
# RECIPES
################################################################################

# Recipes
.PHONY : all
all : inst/Makefile inst/package $(POWERPLOTS) $(ALLSIMSDATAFRAME) \
      $(ALLEXAMPLESPVALS) $(VIGNETTES)

.PRECIOUS : $(ALLSIMS)

################################################################################
# CONFIG-RECIPES
################################################################################

# Configuation and simulation definition file recipes
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

################################################################################
# SIMULATION-DATA-DEPENDS
################################################################################

# Dependencies for simulations
# Normal
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
$(SIMSPROPNORMALRENYI) : data/$(CONTEXTPREFIX)PropChange.Rda \
                         data/$(SIMDATAPREFIX)NormalXY.Rda \
                         data/$(SIMSTATPREFIX)RenyiTypeReg.Rda
$(SIMSPROPNORMALRENYIRESID) : data/$(CONTEXTPREFIX)PropChange.Rda \
                              data/$(SIMDATAPREFIX)NormalXY.Rda \
                              data/$(SIMSTATPREFIX)RenyiTypeResid.Rda
$(SIMSPROPNORMALRENYINOKERN) : data/$(CONTEXTPREFIX)PropChange.Rda \
                               data/$(SIMDATAPREFIX)NormalXY.Rda \
                               data/$(SIMSTATPREFIX)RenyiTypeRegNoKern.Rda
$(SIMSPROPNORMALCUSUM) : data/$(CONTEXTPREFIX)PropChange.Rda \
                         data/$(SIMDATAPREFIX)NormalXY.Rda \
                         data/$(SIMSTATPREFIX)CUSUM.Rda
$(SIMSPROPNORMALHS) : data/$(CONTEXTPREFIX)PropChange.Rda \
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
$(SIMSNORMALHETEROHS) : data/$(CONTEXTPREFIX)PropChange.Rda \
                        data/$(SIMDATAPREFIX)NormalXYHetero.Rda \
                        data/$(SIMSTATPREFIX)HS.Rda
$(SIMSPROPNORMALHETERORENYI) : data/$(CONTEXTPREFIX)PropChange.Rda \
                               data/$(SIMDATAPREFIX)NormalXYHetero.Rda \
                               data/$(SIMSTATPREFIX)RenyiTypeReg.Rda
$(SIMSPROPNORMALHETERORENYINOKERN) : data/$(CONTEXTPREFIX)PropChange.Rda \
                                     data/$(SIMDATAPREFIX)NormalXYHetero.Rda \
                                     data/$(SIMSTATPREFIX)RenyiTypeRegNoKern.Rda
$(SIMSPROPNORMALHETERORENYIRESID) : data/$(CONTEXTPREFIX)PropChange.Rda \
                                    data/$(SIMDATAPREFIX)NormalXYHetero.Rda \
                                    data/$(SIMSTATPREFIX)RenyiTypeResid.Rda
$(SIMSPROPNORMALHETEROCUSUM) : data/$(CONTEXTPREFIX)PropChange.Rda \
                               data/$(SIMDATAPREFIX)NormalXYHetero.Rda \
                               data/$(SIMSTATPREFIX)CUSUM.Rda
$(SIMSPROPNORMALHETEROHS) : data/$(CONTEXTPREFIX)PropChange.Rda \
                            data/$(SIMDATAPREFIX)NormalXYHetero.Rda \
                            data/$(SIMSTATPREFIX)HS.Rda
# ARMA
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
# GARCH
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

################################################################################
# SIMULATION-DATA-RECIPE
################################################################################

# Recipe for simulations
$(ALLSIMS) :
	make package
	$(RSCRIPT) exec/PowerSimRegression.R -C $(word 1, $^) -S $(word 2, $^) \
		 -T $(word 3, $^) -o $@ -N $(POWERREPLICATIONS) -v \
		 -s $(SIMSEED)$(shell echo $@ $^ | md5sum | grep -Eo "[[:digit:]]{3,6}" | head -n1)

################################################################################
# POWER-DF-DEPENDS
################################################################################

# Power estimate collections dependencies
# Normal
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
$(SIMSPROPNORMALRENYIDF) : $(SIMSPROPNORMALRENYI) \
                           data/$(SIMSTATPREFIX)RenyiTypeReg.Rda \
                           data/$(CONTEXTPREFIX)PropChange.Rda
$(SIMSPROPNORMALRENYIRESIDDF) : $(SIMSPROPNORMALRENYIRESID) \
                                data/$(SIMSTATPREFIX)RenyiTypeResid.Rda \
                                data/$(CONTEXTPREFIX)PropChange.Rda
$(SIMSPROPNORMALRENYINOKERNDF) : $(SIMSPROPNORMALRENYINOKERN) \
                                 data/$(SIMSTATPREFIX)RenyiTypeRegNoKern.Rda \
                                 data/$(CONTEXTPREFIX)PropChange.Rda
$(SIMSPROPNORMALCUSUMDF) : $(SIMSPROPNORMALCUSUM) \
                           data/$(SIMSTATPREFIX)CUSUM.Rda \
                           data/$(CONTEXTPREFIX)PropChange.Rda
$(SIMSPROPNORMALHSDF) : $(SIMSPROPNORMALHS) data/$(SIMSTATPREFIX)HS.Rda \
                        data/$(CONTEXTPREFIX)PropChange.Rda
$(SIMSNORMALHETERORENYIDF) : $(SIMSNORMALHETERORENYI) \
                             data/$(SIMSTATPREFIX)RenyiTypeReg.Rda \
                             data/$(CONTEXTPREFIX)Main.Rda
$(SIMSNORMALHETERORENYINOKERNDF) : $(SIMSNORMALHETERORENYINOKERN) \
                                   data/$(SIMSTATPREFIX)RenyiTypeRegNoKern.Rda \
                                   data/$(CONTEXTPREFIX)Main.Rda
$(SIMSNORMALHETERORENYIRESIDDF) : $(SIMSNORMALHETERORENYIRESID) \
                                  data/$(SIMSTATPREFIX)RenyiTypeResid.Rda \
                                  data/$(CONTEXTPREFIX)Main.Rda
$(SIMSNORMALHETEROCUSUMDF) : $(SIMSNORMALHETEROCUSUM) \
                             data/$(SIMSTATPREFIX)CUSUM.Rda \
                             data/$(CONTEXTPREFIX)Main.Rda
$(SIMSNORMALHETEROHSDF) : $(SIMSNORMALHETEROHS) data/$(SIMSTATPREFIX)HS.Rda \
                          data/$(CONTEXTPREFIX)Main.Rda
$(SIMSPROPNORMALHETERORENYIDF) : $(SIMSPROPNORMALHETERORENYI) \
                                 data/$(SIMSTATPREFIX)RenyiTypeReg.Rda \
                                 data/$(CONTEXTPREFIX)PropChange.Rda
$(SIMSPROPNORMALHETERORENYINOKERNDF) : $(SIMSPROPNORMALHETERORENYINOKERN) \
                                       data/$(SIMSTATPREFIX)RenyiTypeRegNoKern.Rda \
                                       data/$(CONTEXTPREFIX)PropChange.Rda
$(SIMSPROPNORMALHETERORENYIRESIDDF) : $(SIMSPROPNORMALHETERORENYIRESID) \
                                      data/$(SIMSTATPREFIX)RenyiTypeResid.Rda \
                                      data/$(CONTEXTPREFIX)PropChange.Rda
$(SIMSPROPNORMALHETEROCUSUMDF) : $(SIMSPROPNORMALHETEROCUSUM) \
                                 data/$(SIMSTATPREFIX)CUSUM.Rda \
                                 data/$(CONTEXTPREFIX)PropChange.Rda
$(SIMSPROPNORMALHETEROHSDF) : $(SIMSPROPNORMALHETEROHS) \
                              data/$(SIMSTATPREFIX)HS.Rda \
                              data/$(CONTEXTPREFIX)PropChange.Rda
# ARMA
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
# GARCH
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
                                  data/$(SIMSTATPREFIX)RenyiTypeRegNoKern.Rda \
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

################################################################################
# POWER-DF-RECIPE
################################################################################

# Recipe for power estimate data frame creation
$(ALLSIMSDATAFRAME) :
	make package
	$(RSCRIPT) exec/Aggregator.R -i $(word 1, $^) -o $@ -a $(LEVEL) \
		 -T $(word 2, $^) -C $(word 3, $^)

################################################################################
# COLLECTION-DF-DEPENDS
################################################################################

# Dependencies for variable collection data frames
$(SIMSNORMALDF) : $(SIMSNORMALDFPREREQ)
$(SIMSPROPNORMALDF) : $(SIMSPROPNORMALDFPREREQ)
$(SIMSNORMALHETERODF) : $(SIMSNORMALHETERODFPREREQ)
$(SIMSPROPNORMALHETERODF) : $(SIMSPROPNORMALHETERODFPREREQ)
$(SIMSARMADF) : $(SIMSARMADFPREREQ)
$(SIMSARMAHETERODF) : $(SIMSARMAHETERODFPREREQ)
$(SIMSGARCHDF) : $(SIMSGARCHDFPREREQ)
$(SIMSGARCHHETERODF) : $(SIMSGARCHHETERODFPREREQ)
$(ALLDISTDF) : exec/Appender.R

################################################################################
# COLLECTION-DF-RECIPE
################################################################################

# Recipe for combining statistics for simulations of DGP into one file
$(ALLDISTDF) :
	make package
	$(RSCRIPT) exec/Appender.R -i $(filter-out exec/Appender.R, $^) -o $@

################################################################################
# POWER-PLOT-RECIPE
################################################################################

# Power plot recipes
# Normal
$(POWERPLOTNORMALPREFIX)%.pdf : $(SIMSNORMALDF) exec/UnivariatePlotter.R \
                                R/Utils.R
	make package
	$(RSCRIPT) exec/UnivariatePlotter.R $< -p $(POWERPLOTNORMALPREFIX) \
		 --width $(POWERPLOTWIDTH) --height $(POWERPLOTHEIGHT) \
		 -l $(LEVEL) --levellinetype $(POWERPLOTLEVELLINE)

$(POWERPLOTPROPNORMALPREFIX)%.pdf : $(SIMSPROPNORMALDF) \
                                    exec/UnivariatePlotter.R R/Utils.R
	make package
	$(RSCRIPT) exec/UnivariatePlotter.R $< -p $(POWERPLOTPROPNORMALPREFIX) \
		 --width $(POWERPLOTWIDTH) --height $(POWERPLOTHEIGHT) \
		 -l $(LEVEL) --levellinetype $(POWERPLOTLEVELLINE)

$(POWERPLOTNORMALHETEROPREFIX)%.pdf : $(SIMSNORMALHETERODF) \
                                      exec/UnivariatePlotter.R R/Utils.R
	make package
	$(RSCRIPT) exec/UnivariatePlotter.R $< -p $(POWERPLOTNORMALHETEROPREFIX) \
		 --width $(POWERPLOTWIDTH) --height $(POWERPLOTHEIGHT) \
		 -l $(LEVEL) --levellinetype $(POWERPLOTLEVELLINE)

$(POWERPLOTPROPNORMALHETEROPREFIX)%.pdf : $(SIMSPROPNORMALHETERODF) \
                                          exec/UnivariatePlotter.R R/Utils.R
	make package
	$(RSCRIPT) exec/UnivariatePlotter.R $< -p \
         $(POWERPLOTPROPNORMALHETEROPREFIX) --width $(POWERPLOTWIDTH) \
         --height $(POWERPLOTHEIGHT) -l $(LEVEL) \
         --levellinetype $(POWERPLOTLEVELLINE)

# ARMA
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

# GARCH
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

################################################################################
# EXAMPLE-RECIPES
################################################################################

data/$(EXAMPLEPREFIX)%.Rda : exec/$(EXAMPLEPREFIX)%.R
	$(RSCRIPT) $< $@

$(EXAMPLEGERMANM1PVALS) : $(EXAMPLEGERMANM1) \
                          data/$(SIMSTATPREFIX)RenyiTypeResid.Rda \
                          data/$(SIMSTATPREFIX)CUSUM.Rda \
                          exec/ExpandingWindowpValComputer.R
	$(RSCRIPT) $(lastword $^) --output $@ \
		 --statistics $(filter-out $< $(lastword $^), $^) --firstright 100 $<

################################################################################
# VIGNETTE-DEPENDS-RECIPE
################################################################################

# Vignette dependencies
doc/CollectedPlots.pdf : vignettes/CollectedPlots.ltx $(POWERPLOTS)

# Vignette recipe
$(VIGNETTES) :
	$(RSCRIPT) -e "devtools::build_vignettes()"

################################################################################
# MISC-RECIPE
################################################################################

# Recipe to create configuration data files
.PHONY : simconfig
simconfig : $(CONTEXTGENERATORS) $(SIMDATAGENERATORS) $(SIMSTATGENERATORS)
	make $(CONTEXTGENERATORS:exec/%.R=data/%.Rda)
	make $(SIMDATAGENERATORS:exec/%.R=data/%.Rda)
	make $(SIMSTATGENERATORS:exec/%.R=data/%.Rda)

# Put Makefile into inst/
inst/Makefile : Makefile
	cp $< $@

# Put package variable (meant to be touched) into inst/
inst/package : package
	cp $< $@

# Recipe to rebuild package if source files have changed
package : R/*.R src/*.cpp
	$(RSCRIPT) exec/RemakePackage.R
	touch package

# Cleaning recipes; mostlyclean for all but simulation data, clean for
# everything
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

# Make empty plots to be recognized as needing updated
.PHONY : init
init :
	make simconfig
	echo "I'm empty for now" > $(POWERPLOTNORMALPREFIX)n50.pdf
	echo "I'm empty for now" > $(POWERPLOTPROPNORMALPREFIX)n50.pdf
	echo "I'm empty for now" > $(POWERPLOTNORMALHETEROPREFIX)n50.pdf
	echo "I'm empty for now" > $(POWERPLOTPROPNORMALHETEROPREFIX)n50.pdf
	echo "I'm empty for now" > $(POWERPLOTARMAPREFIX)n50.pdf
	echo "I'm empty for now" > $(POWERPLOTARMAHETEROPREFIX)n50.pdf
	echo "I'm empty for now" > $(POWERPLOTGARCHPREFIX)n50.pdf
	echo "I'm empty for now" > $(POWERPLOTGARCHHETEROPREFIX)n50.pdf

# Get package dependencies
.PHONY : dependencies
dependencies :
	$(RSCRIPT) exec/GetPackages.R

# Make everything, but small; for testing
.PHONY : small
small :
	make POWERREPLICATIONS=$(SMALLREPLICATIONS)
