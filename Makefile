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
SIMSNORMALRENYI=data-sim/Normal/SimsNormalRenyiReg.Rda
SIMSNORMALRENYIRESID=data-sim/Normal/SimsNormalRenyiResid.Rda
SIMSNORMALRENYINOKERN=data-sim/Normal/SimsNormalRenyiRegNoKern.Rda
SIMSNORMALCUSUM=data-sim/Normal/SimsNormalCUSUM.Rda
SIMSNORMALHS=data-sim/Normal/SimsNormalHS.Rda
SIMSNORMALHETERORENYI=data-sim/Normal/SimsNormalHeteroRenyiReg.Rda
SIMSNORMALHETERORENYINOKERN=data-sim/Normal/SimsNormalHeteroRenyiRegNoKern.Rda
SIMSNORMALHETERORENYIRESID=data-sim/Normal/SimsNormalHeteroRenyiResid.Rda
SIMSNORMALHETEROCUSUM=data-sim/Normal/SimsNormalHeteroCUSUM.Rda
SIMSNORMALHETEROHS=data-sim/Normal/SimsNormalHeteroHS.Rda
SIMSPROPNORMALRENYI=data-sim/PropNormal/SimsPropNormalRenyiReg.Rda
SIMSPROPNORMALRENYIRESID=data-sim/PropNormal/SimsPropNormalRenyiResid.Rda
SIMSPROPNORMALRENYINOKERN=data-sim/PropNormal/SimsPropNormalRenyiRegNoKern.Rda
SIMSPROPNORMALCUSUM=data-sim/PropNormal/SimsPropNormalCUSUM.Rda
SIMSPROPNORMALHS=data-sim/PropNormal/SimsPropNormalHS.Rda
SIMSPROPNORMALHETERORENYI=data-sim/PropNormal/SimsPropNormalHeteroRenyiReg.Rda
SIMSPROPNORMALHETERORENYIRESID=data-sim/PropNormal/SimsPropNormalHeteroRenyiResid.Rda
SIMSPROPNORMALHETERORENYINOKERN=data-sim/PropNormal/SimsPropNormalHeteroRenyiRegNoKern.Rda
SIMSPROPNORMALHETEROCUSUM=data-sim/PropNormal/SimsPropNormalHeteroCUSUM.Rda
SIMSPROPNORMALHETEROHS=data-sim/PropNormal/SimsPropNormalHeteroHS.Rda
SIMS4DNORMALRENYI=data-sim/4DNormal/Sims4DNormalRenyiReg.Rda
SIMS4DNORMALRENYIRESID=data-sim/4DNormal/Sims4DNormalRenyiResid.Rda
SIMS4DNORMALRENYINOKERN=data-sim/4DNormal/Sims4DNormalRenyiRegNoKern.Rda
SIMS4DNORMALCUSUM=data-sim/4DNormal/Sims4DNormalCUSUM.Rda
SIMS4DNORMALHS=data-sim/4DNormal/Sims4DNormalHS.Rda
SIMS4DNORMALHETERORENYI=data-sim/4DNormal/Sims4DNormalHeteroRenyiReg.Rda
SIMS4DNORMALHETERORENYIRESID=data-sim/4DNormal/Sims4DNormalHeteroRenyiResid.Rda
SIMS4DNORMALHETERORENYINOKERN=data-sim/4DNormal/Sims4DNormalHeteroRenyiRegNoKern.Rda
SIMS4DNORMALHETEROCUSUM=data-sim/4DNormal/Sims4DNormalHeteroCUSUM.Rda
SIMS4DNORMALHETEROHS=data-sim/4DNormal/Sims4DNormalHeteroHS.Rda
SIMS6DNORMALRENYI=data-sim/6DNormal/Sims6DNormalRenyiReg.Rda
SIMS6DNORMALRENYIRESID=data-sim/6DNormal/Sims6DNormalRenyiResid.Rda
SIMS6DNORMALRENYINOKERN=data-sim/6DNormal/Sims6DNormalRenyiRegNoKern.Rda
SIMS6DNORMALCUSUM=data-sim/6DNormal/Sims6DNormalCUSUM.Rda
SIMS6DNORMALHS=data-sim/6DNormal/Sims6DNormalHS.Rda
SIMS6DNORMALHETERORENYI=data-sim/6DNormal/Sims6DNormalHeteroRenyiReg.Rda
SIMS6DNORMALHETERORENYIRESID=data-sim/6DNormal/Sims6DNormalHeteroRenyiResid.Rda
SIMS6DNORMALHETERORENYINOKERN=data-sim/6DNormal/Sims6DNormalHeteroRenyiRegNoKern.Rda
SIMS6DNORMALHETEROCUSUM=data-sim/6DNormal/Sims6DNormalHeteroCUSUM.Rda
SIMS6DNORMALHETEROHS=data-sim/6DNormal/Sims6DNormalHeteroHS.Rda
# AR1
SIMSAR1RENYI=data-sim/AR1/SimsAR1RenyiReg.Rda
SIMSAR1RENYIRESID=data-sim/AR1/SimsAR1RenyiResid.Rda
SIMSAR1CUSUM=data-sim/AR1/SimsAR1CUSUM.Rda
SIMSAR1HS=data-sim/AR1/SimsAR1HS.Rda
SIMSAR1HETERORENYI=data-sim/AR1/SimsAR1HeteroRenyiReg.Rda
SIMSAR1HETERORENYIRESID=data-sim/AR1/SimsAR1HeteroRenyiResid.Rda
SIMSAR1HETEROCUSUM=data-sim/AR1/SimsAR1HeteroCUSUM.Rda
SIMSAR1HETEROHS=data-sim/AR1/SimsAR1HeteroHS.Rda
# ARMA
SIMSARMARENYI=data-sim/ARMA/SimsARMARenyiReg.Rda
SIMSARMARENYIRESID=data-sim/ARMA/SimsARMARenyiResid.Rda
SIMSARMACUSUM=data-sim/ARMA/SimsARMACUSUM.Rda
SIMSARMAHS=data-sim/ARMA/SimsARMAHS.Rda
SIMSARMAHETERORENYI=data-sim/ARMA/SimsARMAHeteroRenyiReg.Rda
SIMSARMAHETERORENYIRESID=data-sim/ARMA/SimsARMAHeteroRenyiResid.Rda
SIMSARMAHETEROCUSUM=data-sim/ARMA/SimsARMAHeteroCUSUM.Rda
SIMSARMAHETEROHS=data-sim/ARMA/SimsARMAHeteroHS.Rda
# GARCH
SIMSGARCHRENYI=data-sim/GARCH/SimsGARCHRenyiReg.Rda
SIMSGARCHRENYIRESID=data-sim/GARCH/SimsGARCHRenyiResid.Rda
SIMSGARCHRENYINOKERN=data-sim/GARCH/SimsGARCHRenyiRegNoKern.Rda
SIMSGARCHCUSUM=data-sim/GARCH/SimsGARCHCUSUM.Rda
SIMSGARCHHS=data-sim/GARCH/SimsGARCHHS.Rda
SIMSGARCHHETERORENYI=data-sim/GARCH/SimsGARCHHeteroRenyiReg.Rda
SIMSGARCHHETERORENYINOKERN=data-sim/GARCH/SimsGARCHHeteroRenyiRegNoKern.Rda
SIMSGARCHHETERORENYIRESID=data-sim/GARCH/SimsGARCHHeteroRenyiResid.Rda
SIMSGARCHHETEROCUSUM=data-sim/GARCH/SimsGARCHHeteroCUSUM.Rda
SIMSGARCHHETEROHS=data-sim/GARCH/SimsGARCHHeteroHS.Rda
# EXAMPLE-LIKE
SIMSEXAMPLELIKERENYI=data-sim/EL/SimsELRenyiReg.Rda
SIMSEXAMPLELIKERENYIRESID=data-sim/EL/SimsELRenyiResid.Rda
SIMSEXAMPLELIKERENYINOKERN=data-sim/EL/SimsELRenyiRegNoKern.Rda
SIMSEXAMPLELIKECUSUM=data-sim/EL/SimsELCUSUM.Rda
SIMSEXAMPLELIKEHS=data-sim/EL/SimsELHS.Rda
# All
NOTHETEROSIMS=$(SIMSNORMALRENYI) $(SIMSNORMALRENYIRESID) \
              $(SIMSNORMALRENYINOKERN) $(SIMSNORMALCUSUM) $(SIMSNORMALHS) \
              $(SIMSPROPNORMALRENYI) $(SIMSPROPNORMALRENYIRESID) \
              $(SIMSPROPNORMALRENYINOKERN) $(SIMSPROPNORMALCUSUM) \
              $(SIMSPROPNORMALHS) $(SIMS4DMORNALRENYI) \
              $(SIMS4DNORMALRENYIRESID) $(SIMS4DNORMALRENYINOKERN) \
              $(SIMS4DNORMALCUSUM) $(SIMS4DNORMALHS) $(SIMS6DMORNALRENYI) \
              $(SIMS6DNORMALRENYIRESID) $(SIMS6DNORMALRENYINOKERN) \
              $(SIMS6DNORMALCUSUM) $(SIMS6DNORMALHS) $(SIMSARMARENYI) \
              $(SIMSARMARENYIRESID) $(SIMSARMACUSUM) $(SIMSARMAHS) \
              $(SIMSAR1RENYI) $(SIMSAR1RENYIRESID) $(SIMSAR1CUSUM) \
              $(SIMSAR1HS) $(SIMSGARCHRENYI) $(SIMSGARCHRENYIRESID) \
              $(SIMSGARCHRENYINOKERN) $(SIMSGARCHCUSUM) $(SIMSGARCHHS)
HETEROSIMS=$(SIMSNORMALHETERORENYI) $(SIMSNORMALHETERORENYINOKERN) \
           $(SIMSNORMALHETERORENYIRESID) $(SIMSNORMALHETEROCUSUM) \
           $(SIMSNORMALHETEROHS) $(SIMSPROPNORMALHETERORENYI) \
           $(SIMSPROPNORMALHETERORENYIRESID) \
           $(SIMSPROPNORMALHETERORENYINOKERN) $(SIMSPROPNORMALHETEROCUSUM) \
           $(SIMSPROPNORMALHETEROHS) $(SIMS4DNORMALHETERORENYI) \
           $(SIMS4DNORMALHETERORENYIRESID) \
           $(SIMS4DNORMALHETERORENYINOKERN) $(SIMS4DNORMALHETEROCUSUM) \
           $(SIMS4DNORMALHETEROHS) $(SIMS6DNORMALHETERORENYI) \
           $(SIMS6DNORMALHETERORENYIRESID) \
           $(SIMS6DNORMALHETERORENYINOKERN) $(SIMS6DNORMALHETEROCUSUM) \
           $(SIMS6DNORMALHETEROHS) $(SIMSARMAHETERORENYI) \
           $(SIMSARMAHETERORENYIRESID) $(SIMSARMAHETEROCUSUM) \
           $(SIMSARMAHETEROHS) $(SIMSAR1HETERORENYI) \
           $(SIMSAR1HETERORENYIRESID) $(SIMSAR1HETEROCUSUM) \
           $(SIMSAR1HETEROHS) $(SIMSGARCHHETERORENYI) \
           $(SIMSGARCHHETERORENYINOKERN) $(SIMSGARCHHETERORENYIRESID) \
           $(SIMSGARCHHETEROCUSUM) $(SIMSGARCHHETEROHS)
EXAMPLESIMS=$(SIMSEXAMPLELIKERENYI) $(SIMSEXAMPLELIKERENYIRESID) \
            $(SIMSEXAMPLELIKERENYINOKERN) $(SIMSEXAMPLELIKECUSUM) \
            $(SIMSEXAMPLELIKEHS)
ALLSIMS=$(NOTHETEROSIMS) $(HETEROSIMS) $(EXAMPLESIMS)

################################################################################
# POWER-DF-VARS
################################################################################

# Files containing power estimates computed from simulated test statistics
# Normal
SIMSNORMALRENYIDF=$(subst .Rda,DataFrame.Rda,$(SIMSNORMALRENYI))
1IMSNORMALRENYIRESIDDF=$(subst .Rda,DataFrame.Rda,$(SIMSNORMALRENYIRESID))
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
SIMS4DNORMALRENYIDF=$(subst .Rda,DataFrame.Rda,$(SIMS4DNORMALRENYI))
SIMS4DNORMALRENYIRESIDDF=$(subst .Rda,DataFrame.Rda, \
                           $(SIMS4DNORMALRENYIRESID))
SIMS4DNORMALRENYINOKERNDF=$(subst .Rda,DataFrame.Rda, \
                            $(SIMS4DNORMALRENYINOKERN))
SIMS4DNORMALCUSUMDF=$(subst .Rda,DataFrame.Rda, $(SIMS4DNORMALCUSUM))
SIMS4DNORMALHSDF=$(subst .Rda,DataFrame.Rda,$(SIMS4DNORMALHS))
SIMS4DNORMALHETERORENYIDF=$(subst .Rda,DataFrame.Rda, \
                            $(SIMS4DNORMALHETERORENYI))
SIMS4DNORMALHETERORENYIRESIDDF=$(subst .Rda,DataFrame.Rda, \
                                 $(SIMS4DNORMALHETERORENYIRESID))
SIMS4DNORMALHETERORENYINOKERNDF=$(subst .Rda,DataFrame.Rda, \
                                  $(SIMS4DNORMALHETERORENYINOKERN))
SIMS4DNORMALHETEROCUSUMDF=$(subst .Rda,DataFrame.Rda, \
                            $(SIMS4DNORMALHETEROCUSUM))
SIMS4DNORMALHETEROHSDF=$(subst .Rda,DataFrame.Rda, \
                         $(SIMS4DNORMALHETEROHS))
SIMS6DNORMALRENYIDF=$(subst .Rda,DataFrame.Rda,$(SIMS6DNORMALRENYI))
SIMS6DNORMALRENYIRESIDDF=$(subst .Rda,DataFrame.Rda, \
                           $(SIMS6DNORMALRENYIRESID))
SIMS6DNORMALRENYINOKERNDF=$(subst .Rda,DataFrame.Rda, \
                            $(SIMS6DNORMALRENYINOKERN))
SIMS6DNORMALCUSUMDF=$(subst .Rda,DataFrame.Rda, $(SIMS6DNORMALCUSUM))
SIMS6DNORMALHSDF=$(subst .Rda,DataFrame.Rda,$(SIMS6DNORMALHS))
SIMS6DNORMALHETERORENYIDF=$(subst .Rda,DataFrame.Rda, \
                            $(SIMS6DNORMALHETERORENYI))
SIMS6DNORMALHETERORENYIRESIDDF=$(subst .Rda,DataFrame.Rda, \
                                 $(SIMS6DNORMALHETERORENYIRESID))
SIMS6DNORMALHETERORENYINOKERNDF=$(subst .Rda,DataFrame.Rda, \
                                  $(SIMS6DNORMALHETERORENYINOKERN))
SIMS6DNORMALHETEROCUSUMDF=$(subst .Rda,DataFrame.Rda, \
                            $(SIMS6DNORMALHETEROCUSUM))
SIMS6DNORMALHETEROHSDF=$(subst .Rda,DataFrame.Rda, \
                         $(SIMS6DNORMALHETEROHS))
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
# AR1
SIMSAR1RENYIDF=$(subst .Rda,DataFrame.Rda, $(SIMSAR1RENYI))
SIMSAR1RENYIRESIDDF=$(subst .Rda,DataFrame.Rda,$(SIMSAR1RENYIRESID))
SIMSAR1CUSUMDF=$(subst .Rda,DataFrame.Rda,$(SIMSAR1CUSUM))
SIMSAR1HSDF=$(subst .Rda,DataFrame.Rda,$(SIMSAR1HS))
SIMSAR1HETERORENYIDF=$(subst .Rda,DataFrame.Rda,$(SIMSAR1HETERORENYI))
SIMSAR1HETERORENYIRESIDDF=$(subst .Rda,DataFrame.Rda, \
                                     $(SIMSAR1HETERORENYIRESID))
SIMSAR1HETEROCUSUMDF=$(subst .Rda,DataFrame.Rda,$(SIMSAR1HETEROCUSUM))
SIMSAR1HETEROHSDF=$(subst .Rda,DataFrame.Rda,$(SIMSAR1HETEROHS))
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
# EXAMPLE-LIKE
SIMSEXAMPLELIKERENYIDF=$(subst .Rda,DataFrame.Rda,$(SIMSEXAMPLELIKERENYI))
SIMSEXAMPLELIKERENYIRESIDDF=$(subst .Rda,DataFrame.Rda, \
                                    $(SIMSEXAMPLELIKERENYIRESID))
SIMSEXAMPLELIKERENYINOKERNDF=$(subst .Rda,DataFrame.Rda, \
                                     $(SIMSEXAMPLELIKERENYINOKERN))
SIMSEXAMPLELIKECUSUMDF=$(subst .Rda,DataFrame.Rda,$(SIMSEXAMPLELIKECUSUM))
SIMSEXAMPLELIKEHSDF=$(subst .Rda,DataFrame.Rda,$(SIMSEXAMPLELIKEHS))
# All
ALLSIMSDATAFRAME=$(subst .Rda,DataFrame.Rda,$(ALLSIMS))

################################################################################
# COLLECTION-DF-VARS
################################################################################

# Files collecting power statistics for DGP contexts
# Normal
SIMSNORMALDFPREREQ=$(SIMSNORMALRENYIRESIDDF) $(SIMSNORMALRENYINOKERNDF) \
                   $(SIMSNORMALCUSUMDF) $(SIMSNORMALHSDF)
SIMSNORMALDF=data-sim/Normal/SimsNormal.Rda
SIMSNORMALHETERODFPREREQ=$(SIMSNORMALHETERORENYINOKERNDF) \
                         $(SIMSNORMALHETERORENYIRESIDDF) \
                         $(SIMSNORMALHETEROCUSUMDF) $(SIMSNORMALHETEROHSDF)
SIMSNORMALHETERODF=data-sim/Normal/SimsNormalHetero.Rda
SIMSPROPNORMALDFPREREQ=$(SIMSPROPNORMALRENYIRESIDDF) \
                       $(SIMSPROPNORMALRENYINOKERNDF) $(SIMSPROPNORMALCUSUMDF) \
                       $(SIMSPROPNORMALHSDF)
SIMSPROPNORMALDF=data-sim/PropNormal/SimsPropNormal.Rda
SIMSPROPNORMALHETERODFPREREQ=$(SIMSPROPNORMALHETERORENYIRESIDDF) \
                             $(SIMSPROPNORMALHETERORENYINOKERNDF) \
                             $(SIMSPROPNORMALHETEROCUSUMDF) \
                             $(SIMSPROPNORMALHETEROHSDF)
SIMSPROPNORMALHETERODF=data-sim/PropNormal/SimsPropNormalHetero.Rda
SIMS4DNORMALDFPREREQ=$(SIMS4DNORMALRENYIRESIDDF) \
                       $(SIMS4DNORMALRENYINOKERNDF) $(SIMS4DNORMALCUSUMDF) \
                       $(SIMS4DNORMALHSDF)
SIMS4DNORMALDF=data-sim/4DNormal/Sims4DNormal.Rda
SIMS4DNORMALHETERODFPREREQ=$(SIMS4DNORMALHETERORENYIRESIDDF) \
                             $(SIMS4DNORMALHETERORENYINOKERNDF) \
                             $(SIMS4DNORMALHETEROCUSUMDF) \
                             $(SIMS4DNORMALHETEROHSDF)
SIMS4DNORMALHETERODF=data-sim/4DNormal/Sims4DNormalHetero.Rda
SIMS6DNORMALDFPREREQ=$(SIMS6DNORMALRENYIRESIDDF) \
                       $(SIMS6DNORMALRENYINOKERNDF) $(SIMS6DNORMALCUSUMDF) \
                       $(SIMS6DNORMALHSDF)
SIMS6DNORMALDF=data-sim/6DNormal/Sims6DNormal.Rda
SIMS6DNORMALHETERODFPREREQ=$(SIMS6DNORMALHETERORENYIRESIDDF) \
                             $(SIMS6DNORMALHETERORENYINOKERNDF) \
                             $(SIMS6DNORMALHETEROCUSUMDF) \
                             $(SIMS6DNORMALHETEROHSDF)
SIMS6DNORMALHETERODF=data-sim/6DNormal/Sims6DNormalHetero.Rda
# ARMA
SIMSARMADFPREREQ=$(SIMSARMARENYIRESIDDF) $(SIMSARMARENYIDF) $(SIMSARMACUSUMDF) \
                 $(SIMSARMAHSDF)
SIMSARMADF=data-sim/ARMA/SimsARMA.Rda
SIMSARMAHETERODFPREREQ=$(SIMSARMAHETERORENYIRESIDDF) $(SIMSARMAHETERORENYIDF) \
                       $(SIMSARMAHETEROCUSUMDF) $(SIMSARMAHETEROHSDF)
SIMSARMAHETERODF=data-sim/ARMA/SimsARMAHetero.Rda
# AR1
SIMSAR1DFPREREQ=$(SIMSAR1RENYIRESIDDF) $(SIMSAR1RENYIDF) $(SIMSAR1CUSUMDF) \
                 $(SIMSAR1HSDF)
SIMSAR1DF=data-sim/AR1/SimsAR1.Rda
SIMSAR1HETERODFPREREQ=$(SIMSAR1HETERORENYIRESIDDF) $(SIMSAR1HETERORENYIDF) \
                       $(SIMSAR1HETEROCUSUMDF) $(SIMSAR1HETEROHSDF)
SIMSAR1HETERODF=data-sim/AR1/SimsAR1Hetero.Rda
# GARCH
SIMSGARCHDFPREREQ=$(SIMSGARCHRENYIRESIDDF) $(SIMSGARCHRENYINOKERNDF) \
                  $(SIMSGARCHCUSUMDF) $(SIMSGARCHHSDF)
SIMSGARCHDF=data-sim/GARCH/SimsGARCH.Rda
SIMSGARCHHETERODFPREREQ=$(SIMSGARCHHETERORENYIRESIDDF) \
                        $(SIMSGARCHHETERORENYINOKERNDF) \
                        $(SIMSGARCHHETEROCUSUMDF) $(SIMSGARCHHETEROHSDF)
SIMSGARCHHETERODF=data-sim/GARCH/SimsGARCHHetero.Rda
# EXAMPLE-LIKE
SIMSEXAMPLELIKEDFPREREQ=$(SIMSEXAMPLELIKERENYIRESIDDF) \
                        $(SIMSEXAMPLELIKERENYIDF) \
                        $(SIMSEXAMPLELIKECUSUMDF) $(SIMSEXAMPLELIKEHSDF)
SIMSEXAMPLELIKEDF=data-sim/EL/SimsEL.Rda

# All DGP power collections in one location
ALLDISTDF=$(SIMSNORMALDF) $(SIMSNORMALHETERODF) $(SIMSPROPNORMALDF) \
          $(SIMS4DNORMALDF) $(SIMS6DNORMALDF) $(SIMSPROPNORMALHETERODF) \
          $(SIMS4DNORMALHETERODF) $(SIMS6DNORMALHETERODF) $(SIMSARMADF) \
          $(SIMSARMAHETERODF) $(SIMSAR1DF) $(SIMSAR1HETERODF) $(SIMSGARCHDF) \
          $(SIMSGARCHHETERODF) $(SIMSEXAMPLELIKEDF)

################################################################################
# POWER-PLOT-VARS
################################################################################

# Power plot variables, giving plots names and defining plot dimensions
POWERPLOTPREFIX=vignettes/power_plot_
POWERPLOTNORMALPREFIX=$(POWERPLOTPREFIX)norm_
POWERPLOTPROPNORMALPREFIX=$(POWERPLOTPREFIX)prop_norm_
POWERPLOT4DNORMALPREFIX=$(POWERPLOTPREFIX)4D_norm_
POWERPLOT6DNORMALPREFIX=$(POWERPLOTPREFIX)6D_norm_
POWERPLOTARMAPREFIX=$(POWERPLOTPREFIX)ARMA_
POWERPLOTAR1PREFIX=$(POWERPLOTPREFIX)AR1_
POWERPLOTGARCHPREFIX=$(POWERPLOTPREFIX)GARCH_
POWERPLOTNORMALHETEROPREFIX=$(POWERPLOTPREFIX)norm_hetero_
POWERPLOTPROPNORMALHETEROPREFIX=$(POWERPLOTPREFIX)prop_norm_hetero_
POWERPLOT4DNORMALHETEROPREFIX=$(POWERPLOTPREFIX)4D_norm_hetero_
POWERPLOT6DNORMALHETEROPREFIX=$(POWERPLOTPREFIX)6D_norm_hetero_
POWERPLOTARMAHETEROPREFIX=$(POWERPLOTPREFIX)ARMA_hetero_
POWERPLOTAR1HETEROPREFIX=$(POWERPLOTPREFIX)AR1_hetero_
POWERPLOTGARCHHETEROPREFIX=$(POWERPLOTPREFIX)GARCH_hetero_
POWERPLOTEXAMPLELIKEPREFIX=$(POWERPLOTPREFIX)Example_
POWERPLOTWIDTH=3
POWERPLOTHEIGHT=2
POWERPLOTLEVELLINE=dotted
POWERPLOTS=$(wildcard $(POWERPLOTNORMALPREFIX)*.pdf) \
           $(wildcard $(POWERPLOTPROPNORMALPREFIX)*.pdf) \
           $(wildcard $(POWERPLOT4DNORMALPREFIX)*.pdf) \
           $(wildcard $(POWERPLOT6DNORMALPREFIX)*.pdf) \
           $(wildcard $(POWERPLOTNORMALHETEROPREFIX)*.pdf) \
           $(wildcard $(POWERPLOTPROPNORMALHETEROPREFIX)*.pdf) \
           $(wildcard $(POWERPLOT4DNORMALHETEROPREFIX)*.pdf) \
           $(wildcard $(POWERPLOT6DNORMALHETEROPREFIX)*.pdf) \
           $(wildcard $(POWERPLOTARMAPREFIX)*.pdf) \
           $(wildcard $(POWERPLOTARMAHETEROPREFIX)*.pdf) \
           $(wildcard $(POWERPLOTAR1PREFIX)*.pdf) \
           $(wildcard $(POWERPLOTAR1HETEROPREFIX)*.pdf) \
           $(wildcard $(POWERPLOTGARCHPREFIX)*.pdf) \
           $(wildcard $(POWERPLOTGARCHHETEROPREFIX)*.pdf) \
           $(wildcard $(POWERPLOTEXAMPLELIKEPREFIX)*.pdf)

################################################################################
# EXAMPLE-VARS
################################################################################

# Data example definition files
EXAMPLEPREFIX=RealData
EXAMPLEGERMANM1=data/$(EXAMPLEPREFIX)M1Germany.Rda
EXAMPLECXW=data/$(EXAMPLEPREFIX)CXWDonaldTrump.Rda
EXAMPLEGDPPMI=data/$(EXAMPLEPREFIX)GDPPMIPredict.Rda
EXAMPLEOILDOLLAR=data/$(EXAMPLEPREFIX)OilDollar.Rda
ALLEXAMPLES=$(EXAMPLEGERMANM1) $(EXAMPLECXW) $(EXAMPLEGDPPMI) \
            $(EXAMPLEOILDOLLAR)

# Processed data examples
EXAMPLEGERMANM1PVALS=data/$(EXAMPLEPREFIX)M1GermanypVals.Rda
EXAMPLECXWPVALS=data/$(EXAMPLEPREFIX)CXWDonaldTrumppVals.Rda
EXAMPLEGDPPMIPVALS=data/$(EXAMPLEPREFIX)GDPPMIPredictpVals.Rda
EXAMPLEOILDOLLARPVALS=data/$(EXAMPLEPREFIX)OilDollarpVals.Rda
ALLEXAMPLESPVALS=$(EXAMPLEGERMAN1PVALS) $(EXAMPLECXWPVALS) \
                 $(EXAMPLEGDPPMIPVALS) $(EXAMPLEOILDOLLARPVALS)

# Data example plots
EXAMPLEGERMANM1PLOT=vignettes/$(EXAMPLEPREFIX)M1Germany
EXAMPLECXWPLOT=vignettes/$(EXAMPLEPREFIX)CXW
EXAMPLEGDPPMIPLOT=vignettes/$(EXAMPLEPREFIX)GDPPMI
EXAMPLEOILDOLLARPLOT=vignettes/$(EXAMPLEPREFIX)OilDollar
ALLEXAMPLEPLOTS=$(EXAMPLEGERMANM1PLOT) $(EXAMPLECXWPLOT) $(EXAMPLEGDPPMIPLOT) \
                $(EXAMPLEOILDOLLARPLOT)

################################################################################
# VIGNETTE-VARS
################################################################################

# Vignettes to be created
VIGNETTES=doc/CollectedPlots.pdf doc/ZnTable.pdf \
          doc/CXWDonaldTrumpElection.pdf doc/TypeIErrors.pdf \
          doc/GDPPMIStructuralChange.pdf \
          doc/OilDollarStructuralChange.pdf

################################################################################
# RECIPES
################################################################################

# Recipes
.PHONY : all
all : inst/Makefile inst/package $(POWERPLOTS) $(ALLSIMSDATAFRAME) \
      $(ALLEXAMPLEPLOTS:=.pdf) $(VIGNETTES)

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
$(SIMS4DNORMALRENYI) : data/$(CONTEXTPREFIX)4D.Rda \
                         data/$(SIMDATAPREFIX)NormalXY.Rda \
                         data/$(SIMSTATPREFIX)RenyiTypeReg.Rda
$(SIMS4DNORMALRENYIRESID) : data/$(CONTEXTPREFIX)4D.Rda \
                              data/$(SIMDATAPREFIX)NormalXY.Rda \
                              data/$(SIMSTATPREFIX)RenyiTypeResid.Rda
$(SIMS4DNORMALRENYINOKERN) : data/$(CONTEXTPREFIX)4D.Rda \
                               data/$(SIMDATAPREFIX)NormalXY.Rda \
                               data/$(SIMSTATPREFIX)RenyiTypeRegNoKern.Rda
$(SIMS4DNORMALCUSUM) : data/$(CONTEXTPREFIX)4D.Rda \
                         data/$(SIMDATAPREFIX)NormalXY.Rda \
                         data/$(SIMSTATPREFIX)CUSUM.Rda
$(SIMS4DNORMALHS) : data/$(CONTEXTPREFIX)4D.Rda \
                      data/$(SIMDATAPREFIX)NormalXY.Rda \
                      data/$(SIMSTATPREFIX)HS.Rda
$(SIMS6DNORMALRENYI) : data/$(CONTEXTPREFIX)6D.Rda \
                         data/$(SIMDATAPREFIX)NormalXY.Rda \
                         data/$(SIMSTATPREFIX)RenyiTypeReg.Rda
$(SIMS6DNORMALRENYIRESID) : data/$(CONTEXTPREFIX)6D.Rda \
                              data/$(SIMDATAPREFIX)NormalXY.Rda \
                              data/$(SIMSTATPREFIX)RenyiTypeResid.Rda
$(SIMS6DNORMALRENYINOKERN) : data/$(CONTEXTPREFIX)6D.Rda \
                               data/$(SIMDATAPREFIX)NormalXY.Rda \
                               data/$(SIMSTATPREFIX)RenyiTypeRegNoKern.Rda
$(SIMS6DNORMALCUSUM) : data/$(CONTEXTPREFIX)6D.Rda \
                         data/$(SIMDATAPREFIX)NormalXY.Rda \
                         data/$(SIMSTATPREFIX)CUSUM.Rda
$(SIMS6DNORMALHS) : data/$(CONTEXTPREFIX)6D.Rda \
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
$(SIMS4DNORMALHETERORENYI) : data/$(CONTEXTPREFIX)4D.Rda \
                               data/$(SIMDATAPREFIX)NormalXYHetero.Rda \
                               data/$(SIMSTATPREFIX)RenyiTypeReg.Rda
$(SIMS4DNORMALHETERORENYINOKERN) : data/$(CONTEXTPREFIX)4D.Rda \
                                     data/$(SIMDATAPREFIX)NormalXYHetero.Rda \
                                     data/$(SIMSTATPREFIX)RenyiTypeRegNoKern.Rda
$(SIMS4DNORMALHETERORENYIRESID) : data/$(CONTEXTPREFIX)4D.Rda \
                                    data/$(SIMDATAPREFIX)NormalXYHetero.Rda \
                                    data/$(SIMSTATPREFIX)RenyiTypeResid.Rda
$(SIMS4DNORMALHETEROCUSUM) : data/$(CONTEXTPREFIX)4D.Rda \
                               data/$(SIMDATAPREFIX)NormalXYHetero.Rda \
                               data/$(SIMSTATPREFIX)CUSUM.Rda
$(SIMS4DNORMALHETEROHS) : data/$(CONTEXTPREFIX)4D.Rda \
                            data/$(SIMDATAPREFIX)NormalXYHetero.Rda \
                            data/$(SIMSTATPREFIX)HS.Rda
$(SIMS6DNORMALHETERORENYI) : data/$(CONTEXTPREFIX)6D.Rda \
                               data/$(SIMDATAPREFIX)NormalXYHetero.Rda \
                               data/$(SIMSTATPREFIX)RenyiTypeReg.Rda
$(SIMS6DNORMALHETERORENYINOKERN) : data/$(CONTEXTPREFIX)6D.Rda \
                                     data/$(SIMDATAPREFIX)NormalXYHetero.Rda \
                                     data/$(SIMSTATPREFIX)RenyiTypeRegNoKern.Rda
$(SIMS6DNORMALHETERORENYIRESID) : data/$(CONTEXTPREFIX)6D.Rda \
                                    data/$(SIMDATAPREFIX)NormalXYHetero.Rda \
                                    data/$(SIMSTATPREFIX)RenyiTypeResid.Rda
$(SIMS6DNORMALHETEROCUSUM) : data/$(CONTEXTPREFIX)6D.Rda \
                               data/$(SIMDATAPREFIX)NormalXYHetero.Rda \
                               data/$(SIMSTATPREFIX)CUSUM.Rda
$(SIMS6DNORMALHETEROHS) : data/$(CONTEXTPREFIX)6D.Rda \
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
# AR1
$(SIMSAR1RENYI) : data/$(CONTEXTPREFIX)Main.Rda \
                  data/$(SIMDATAPREFIX)AR1XY.Rda \
                  data/$(SIMSTATPREFIX)RenyiTypeReg.Rda
$(SIMSAR1RENYIRESID) : data/$(CONTEXTPREFIX)Main.Rda \
                       data/$(SIMDATAPREFIX)AR1XY.Rda \
                       data/$(SIMSTATPREFIX)RenyiTypeResid.Rda
$(SIMSAR1CUSUM) : data/$(CONTEXTPREFIX)Main.Rda \
                  data/$(SIMDATAPREFIX)AR1XY.Rda \
                  data/$(SIMSTATPREFIX)CUSUM.Rda
$(SIMSAR1HS) : data/$(CONTEXTPREFIX)Main.Rda \
               data/$(SIMDATAPREFIX)AR1XY.Rda \
               data/$(SIMSTATPREFIX)HS.Rda
$(SIMSAR1HETERORENYI) : data/$(CONTEXTPREFIX)Main.Rda \
                        data/$(SIMDATAPREFIX)AR1XYHetero.Rda \
                        data/$(SIMSTATPREFIX)RenyiTypeReg.Rda
$(SIMSAR1HETERORENYIRESID) : data/$(CONTEXTPREFIX)Main.Rda \
                             data/$(SIMDATAPREFIX)AR1XYHetero.Rda \
                             data/$(SIMSTATPREFIX)RenyiTypeResid.Rda
$(SIMSAR1HETEROCUSUM) : data/$(CONTEXTPREFIX)Main.Rda \
                        data/$(SIMDATAPREFIX)AR1XYHetero.Rda \
                        data/$(SIMSTATPREFIX)CUSUM.Rda
$(SIMSAR1HETEROHS) : data/$(CONTEXTPREFIX)Main.Rda \
                     data/$(SIMDATAPREFIX)AR1XYHetero.Rda \
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
# EXAMPLE-LIKE
$(SIMSEXAMPLELIKERENYI) : data/$(CONTEXTPREFIX)Example.Rda \
                          data/$(SIMDATAPREFIX)Example.Rda \
                          data/$(SIMSTATPREFIX)RenyiTypeRegAlt.Rda
$(SIMSEXAMPLELIKERENYIRESID) : data/$(CONTEXTPREFIX)Example.Rda \
                               data/$(SIMDATAPREFIX)Example.Rda \
                               data/$(SIMSTATPREFIX)RenyiTypeResidAlt.Rda
$(SIMSEXAMPLELIKERENYINOKERN) : data/$(CONTEXTPREFIX)Example.Rda \
                                data/$(SIMDATAPREFIX)Example.Rda \
                                data/$(SIMSTATPREFIX)RenyiTypeRegNoKern.Rda
$(SIMSEXAMPLELIKECUSUM) : data/$(CONTEXTPREFIX)Example.Rda \
                          data/$(SIMDATAPREFIX)Example.Rda \
                          data/$(SIMSTATPREFIX)CUSUM.Rda
$(SIMSEXAMPLELIKEHS) : data/$(CONTEXTPREFIX)Example.Rda \
                       data/$(SIMDATAPREFIX)Example.Rda \
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
		 -s $(SIMSEED)$(shell echo "data/$(notdir $@) $^" | md5sum | grep -Eo "[[:digit:]]{3,6}" | head -n1)

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
$(SIMS4DNORMALRENYIDF) : $(SIMS4DNORMALRENYI) \
                           data/$(SIMSTATPREFIX)RenyiTypeReg.Rda \
                           data/$(CONTEXTPREFIX)4D.Rda
$(SIMS4DNORMALRENYIRESIDDF) : $(SIMS4DNORMALRENYIRESID) \
                                data/$(SIMSTATPREFIX)RenyiTypeResid.Rda \
                                data/$(CONTEXTPREFIX)4D.Rda
$(SIMS4DNORMALRENYINOKERNDF) : $(SIMS4DNORMALRENYINOKERN) \
                                 data/$(SIMSTATPREFIX)RenyiTypeRegNoKern.Rda \
                                 data/$(CONTEXTPREFIX)4D.Rda
$(SIMS4DNORMALCUSUMDF) : $(SIMS4DNORMALCUSUM) \
                           data/$(SIMSTATPREFIX)CUSUM.Rda \
                           data/$(CONTEXTPREFIX)4D.Rda
$(SIMS4DNORMALHSDF) : $(SIMS4DNORMALHS) data/$(SIMSTATPREFIX)HS.Rda \
                        data/$(CONTEXTPREFIX)4D.Rda
$(SIMS6DNORMALRENYIDF) : $(SIMS6DNORMALRENYI) \
                           data/$(SIMSTATPREFIX)RenyiTypeReg.Rda \
                           data/$(CONTEXTPREFIX)6D.Rda
$(SIMS6DNORMALRENYIRESIDDF) : $(SIMS6DNORMALRENYIRESID) \
                                data/$(SIMSTATPREFIX)RenyiTypeResid.Rda \
                                data/$(CONTEXTPREFIX)6D.Rda
$(SIMS6DNORMALRENYINOKERNDF) : $(SIMS6DNORMALRENYINOKERN) \
                                 data/$(SIMSTATPREFIX)RenyiTypeRegNoKern.Rda \
                                 data/$(CONTEXTPREFIX)6D.Rda
$(SIMS6DNORMALCUSUMDF) : $(SIMS6DNORMALCUSUM) \
                           data/$(SIMSTATPREFIX)CUSUM.Rda \
                           data/$(CONTEXTPREFIX)6D.Rda
$(SIMS6DNORMALHSDF) : $(SIMS6DNORMALHS) data/$(SIMSTATPREFIX)HS.Rda \
                        data/$(CONTEXTPREFIX)6D.Rda
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
$(SIMS4DNORMALHETERORENYIDF) : $(SIMS4DNORMALHETERORENYI) \
                                 data/$(SIMSTATPREFIX)RenyiTypeReg.Rda \
                                 data/$(CONTEXTPREFIX)4D.Rda
$(SIMS4DNORMALHETERORENYINOKERNDF) : $(SIMS4DNORMALHETERORENYINOKERN) \
                                       data/$(SIMSTATPREFIX)RenyiTypeRegNoKern.Rda \
                                       data/$(CONTEXTPREFIX)4D.Rda
$(SIMS4DNORMALHETERORENYIRESIDDF) : $(SIMS4DNORMALHETERORENYIRESID) \
                                      data/$(SIMSTATPREFIX)RenyiTypeResid.Rda \
                                      data/$(CONTEXTPREFIX)4D.Rda
$(SIMS4DNORMALHETEROCUSUMDF) : $(SIMS4DNORMALHETEROCUSUM) \
                                 data/$(SIMSTATPREFIX)CUSUM.Rda \
                                 data/$(CONTEXTPREFIX)4D.Rda
$(SIMS4DNORMALHETEROHSDF) : $(SIMS4DNORMALHETEROHS) \
                              data/$(SIMSTATPREFIX)HS.Rda \
                              data/$(CONTEXTPREFIX)4D.Rda
$(SIMS6DNORMALHETERORENYIDF) : $(SIMS6DNORMALHETERORENYI) \
                                 data/$(SIMSTATPREFIX)RenyiTypeReg.Rda \
                                 data/$(CONTEXTPREFIX)6D.Rda
$(SIMS6DNORMALHETERORENYINOKERNDF) : $(SIMS6DNORMALHETERORENYINOKERN) \
                                       data/$(SIMSTATPREFIX)RenyiTypeRegNoKern.Rda \
                                       data/$(CONTEXTPREFIX)6D.Rda
$(SIMS6DNORMALHETERORENYIRESIDDF) : $(SIMS6DNORMALHETERORENYIRESID) \
                                      data/$(SIMSTATPREFIX)RenyiTypeResid.Rda \
                                      data/$(CONTEXTPREFIX)6D.Rda
$(SIMS6DNORMALHETEROCUSUMDF) : $(SIMS6DNORMALHETEROCUSUM) \
                                 data/$(SIMSTATPREFIX)CUSUM.Rda \
                                 data/$(CONTEXTPREFIX)6D.Rda
$(SIMS6DNORMALHETEROHSDF) : $(SIMS6DNORMALHETEROHS) \
                              data/$(SIMSTATPREFIX)HS.Rda \
                              data/$(CONTEXTPREFIX)6D.Rda
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
# AR1
$(SIMSAR1RENYIDF) : $(SIMSAR1RENYI) \
                    data/$(SIMSTATPREFIX)RenyiTypeReg.Rda \
                    data/$(CONTEXTPREFIX)Main.Rda
$(SIMSAR1RENYIRESIDDF) : $(SIMSAR1RENYIRESID) \
                         data/$(SIMSTATPREFIX)RenyiTypeResid.Rda \
                         data/$(CONTEXTPREFIX)Main.Rda
$(SIMSAR1CUSUMDF) : $(SIMSAR1CUSUM) data/$(SIMSTATPREFIX)CUSUM.Rda \
                    data/$(CONTEXTPREFIX)Main.Rda
$(SIMSAR1HSDF) : $(SIMSAR1HS) data/$(SIMSTATPREFIX)HS.Rda \
                 data/$(CONTEXTPREFIX)Main.Rda
$(SIMSAR1HETERORENYIDF) : $(SIMSAR1HETERORENYI) \
                          data/$(SIMSTATPREFIX)RenyiTypeReg.Rda \
                          data/$(CONTEXTPREFIX)Main.Rda
$(SIMSAR1HETERORENYIRESIDDF) : $(SIMSAR1HETERORENYIRESID) \
                               data/$(SIMSTATPREFIX)RenyiTypeResid.Rda \
                               data/$(CONTEXTPREFIX)Main.Rda
$(SIMSAR1HETEROCUSUMDF) : $(SIMSAR1HETEROCUSUM) \
                          data/$(SIMSTATPREFIX)CUSUM.Rda \
                          data/$(CONTEXTPREFIX)Main.Rda
$(SIMSAR1HETEROHSDF) : $(SIMSAR1HETEROHS) data/$(SIMSTATPREFIX)HS.Rda \
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
# EXAMPLE-LIKE
$(SIMSEXAMPLELIKERENYIDF) : $(SIMSEXAMPLELIKERENYI) \
                            data/$(SIMSTATPREFIX)RenyiTypeReg.Rda \
                            data/$(CONTEXTPREFIX)Example.Rda
$(SIMSEXAMPLELIKERENYIRESIDDF) : $(SIMSEXAMPLELIKERENYIRESID) \
                                 data/$(SIMSTATPREFIX)RenyiTypeResid.Rda \
                                 data/$(CONTEXTPREFIX)Example.Rda
$(SIMSEXAMPLELIKERENYINOKERNDF) : $(SIMSEXAMPLELIKERENYINOKERN) \
                                  data/$(SIMSTATPREFIX)RenyiTypeRegNoKern.Rda \
                                  data/$(CONTEXTPREFIX)Example.Rda
$(SIMSEXAMPLELIKECUSUMDF) : $(SIMSEXAMPLELIKECUSUM) \
                            data/$(SIMSTATPREFIX)CUSUM.Rda \
                            data/$(CONTEXTPREFIX)Example.Rda
$(SIMSEXAMPLELIKEHSDF) : $(SIMSEXAMPLELIKEHS) data/$(SIMSTATPREFIX)HS.Rda \
                         data/$(CONTEXTPREFIX)Example.Rda

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
$(SIMS4DNORMALDF) : $(SIMS4DNORMALDFPREREQ)
$(SIMS6DNORMALDF) : $(SIMS6DNORMALDFPREREQ)
$(SIMSNORMALHETERODF) : $(SIMSNORMALHETERODFPREREQ)
$(SIMSPROPNORMALHETERODF) : $(SIMSPROPNORMALHETERODFPREREQ)
$(SIMS4DNORMALHETERODF) : $(SIMS4DNORMALHETERODFPREREQ)
$(SIMS6DNORMALHETERODF) : $(SIMS6DNORMALHETERODFPREREQ)
$(SIMSARMADF) : $(SIMSARMADFPREREQ)
$(SIMSARMAHETERODF) : $(SIMSARMAHETERODFPREREQ)
$(SIMSAR1DF) : $(SIMSAR1DFPREREQ)
$(SIMSAR1HETERODF) : $(SIMSAR1HETERODFPREREQ)
$(SIMSGARCHDF) : $(SIMSGARCHDFPREREQ)
$(SIMSGARCHHETERODF) : $(SIMSGARCHHETERODFPREREQ)
$(SIMSEXAMPLELIKEDF) : $(SIMSEXAMPLELIKEDFPREREQ)
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

$(POWERPLOT4DNORMALPREFIX)%.pdf : $(SIMS4DNORMALDF) \
                                    exec/UnivariatePlotter.R R/Utils.R
	make package
	$(RSCRIPT) exec/UnivariatePlotter.R $< -p $(POWERPLOT4DNORMALPREFIX) \
		 --width $(POWERPLOTWIDTH) --height $(POWERPLOTHEIGHT) \
		 -l $(LEVEL) --levellinetype $(POWERPLOTLEVELLINE)

$(POWERPLOT6DNORMALPREFIX)%.pdf : $(SIMS6DNORMALDF) \
                                    exec/UnivariatePlotter.R R/Utils.R
	make package
	$(RSCRIPT) exec/UnivariatePlotter.R $< -p $(POWERPLOT6DNORMALPREFIX) \
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

$(POWERPLOT4DNORMALHETEROPREFIX)%.pdf : $(SIMS4DNORMALHETERODF) \
                                          exec/UnivariatePlotter.R R/Utils.R
	make package
	$(RSCRIPT) exec/UnivariatePlotter.R $< -p \
         $(POWERPLOT4DNORMALHETEROPREFIX) --width $(POWERPLOTWIDTH) \
         --height $(POWERPLOTHEIGHT) -l $(LEVEL) \
         --levellinetype $(POWERPLOTLEVELLINE)

$(POWERPLOT6DNORMALHETEROPREFIX)%.pdf : $(SIMS6DNORMALHETERODF) \
                                          exec/UnivariatePlotter.R R/Utils.R
	make package
	$(RSCRIPT) exec/UnivariatePlotter.R $< -p \
         $(POWERPLOT6DNORMALHETEROPREFIX) --width $(POWERPLOTWIDTH) \
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

# AR1
$(POWERPLOTAR1PREFIX)%.pdf : $(SIMSAR1DF) exec/UnivariatePlotter.R \
                             R/Utils.R
	make package
	$(RSCRIPT) exec/UnivariatePlotter.R $< -p $(POWERPLOTAR1PREFIX) \
		 --width $(POWERPLOTWIDTH) --height $(POWERPLOTHEIGHT) \
		 -l $(LEVEL) --levellinetype $(POWERPLOTLEVELLINE)

$(POWERPLOTAR1HETEROPREFIX)%.pdf : $(SIMSAR1HETERODF) \
                                   exec/UnivariatePlotter.R R/Utils.R
	make package
	$(RSCRIPT) exec/UnivariatePlotter.R $< -p $(POWERPLOTAR1HETEROPREFIX) \
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

# EXAMPLE-LIKE
$(POWERPLOTEXAMPLELIKEPREFIX)%.pdf : $(SIMSEXAMPLELIKEDF) \
                                     exec/UnivariatePlotter.R \
                                     R/Utils.R
	make package
	$(RSCRIPT) exec/UnivariatePlotter.R $< -p $(POWERPLOTEXAMPLELIKEPREFIX) \
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

$(EXAMPLECXWPVALS) : $(EXAMPLECXW) data/$(SIMSTATPREFIX)RenyiTypeResid.Rda \
                     data/$(SIMSTATPREFIX)CUSUM.Rda \
                     data/$(SIMSTATPREFIX)RenyiTypeReg.Rda \
                     data/$(SIMSTATPREFIX)HS.Rda \
                     exec/ExpandingWindowpValComputer.R
	$(RSCRIPT) $(lastword $^) --output $@ \
		 --statistics $(filter-out $< $(lastword $^), $^) --left 12965 \
		 --firstright 13362 --lastright 13385 $<

$(EXAMPLEOILDOLLARPVALS) : $(EXAMPLEOILDOLLAR) \
                           data/$(SIMSTATPREFIX)RenyiTypeResidAlt.Rda \
                           data/$(SIMSTATPREFIX)CUSUMAlt.Rda \
                           data/$(SIMSTATPREFIX)RenyiTypeRegAlt.Rda \
                           data/$(SIMSTATPREFIX)HS.Rda \
                           exec/ExpandingWindowpValComputer.R
	$(RSCRIPT) $(lastword $^) --output $@ \
		 --statistics $(filter-out $< $(lastword $^), $^) \
		 --firstright 120 --lastright 138 $<

$(EXAMPLEGDPPMIPVALS) : $(EXAMPLEGDPPMI) \
                        data/$(SIMSTATPREFIX)RenyiTypeResidAlt.Rda \
                        data/$(SIMSTATPREFIX)CUSUM.Rda \
                        data/$(SIMSTATPREFIX)RenyiTypeRegAlt.Rda \
                        data/$(SIMSTATPREFIX)HS.Rda \
                        exec/ExpandingWindowpValComputer.R
	$(RSCRIPT) $(lastword $^) --output $@ \
		 --statistics $(filter-out $< $(lastword $^), $^) \
		 --firstright 200 $<

$(EXAMPLEGERMANM1PLOT).pdf : $(EXAMPLEGERMANM1PVALS)
$(EXAMPLECXWPLOT).pdf : $(EXAMPLECXWPVALS)
$(EXAMPLEGDPPMIPLOT).pdf : $(EXAMPLEGDPPMIPVALS)
$(EXAMPLEOILDOLLARPLOT).pdf : $(EXAMPLEOILDOLLARPVALS)
$(ALLEXAMPLEPLOTS:=.pdf) : exec/TimeSeriespValPlotter.R

$(ALLEXAMPLEPLOTS:=.pdf) :
	$(RSCRIPT) $(lastword $^) $< -o $(basename $@) -d -t

################################################################################
# VIGNETTE-DEPENDS-RECIPE
################################################################################

# Vignette dependencies
doc/CollectedPlots.pdf : vignettes/CollectedPlots.ltx $(POWERPLOTS)
doc/ZnTable.pdf : vignettes/ZnTable.ltx
doc/TypeIErrors.pdf : vignettes/TypeIErrors.ltx
doc/CXWDonaldTrumpElection.pdf : vignettes/CXWDonaldTrumpElection.ltx \
                                 $(EXAMPLECXWPLOT).pdf
doc/GDPPMIStructuralChange.pdf: vignettes/GDPPMIStructuralChange.ltx \
                                $(EXAMPLEGDPPMIPLOT).pdf \
                                vignettes/GDPPMIPowerPlot.pdf
doc/OilDollarStructuralChange.pdf : vignettes/OilDollarStructuralChange.ltx \
                                    vignettes/OilDollar.pdf \
                                    $(EXAMPLEOILDOLLARPLOT).pdf \
                                    $(POWERPLOTEXAMPLELIKEPREFIX)n130.pdf
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
	touch $(POWERPLOTNORMALPREFIX)n50.pdf
	touch $(POWERPLOTPROPNORMALPREFIX)n50.pdf
	touch $(POWERPLOT4DNORMALPREFIX)n50.pdf
	touch $(POWERPLOT6DNORMALPREFIX)n50.pdf
	touch $(POWERPLOTNORMALHETEROPREFIX)n50.pdf
	touch $(POWERPLOTPROPNORMALHETEROPREFIX)n50.pdf
	touch $(POWERPLOT4DNORMALHETEROPREFIX)n50.pdf
	touch $(POWERPLOT6DNORMALHETEROPREFIX)n50.pdf
	touch $(POWERPLOTARMAPREFIX)n50.pdf
	touch $(POWERPLOTARMAHETEROPREFIX)n50.pdf
	touch $(POWERPLOTAR1PREFIX)n50.pdf
	touch $(POWERPLOTAR1HETEROPREFIX)n50.pdf
	touch $(POWERPLOTGARCHPREFIX)n50.pdf
	touch $(POWERPLOTGARCHHETEROPREFIX)n50.pdf
	touch $(POWERPLOTEXAMPLELIKEPREFIX)n250.pdf
	make simconfig

# Get package dependencies
.PHONY : dependencies
dependencies :
	$(RSCRIPT) exec/GetPackages.R

# Make everything, but small; for testing
.PHONY : small
small :
	make POWERREPLICATIONS=$(SMALLREPLICATIONS)
