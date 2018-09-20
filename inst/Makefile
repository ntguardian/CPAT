SHELL=/bin/sh
RSCRIPT=Rscript

<<<<<<< HEAD
=======
LEVEL=0.05
>>>>>>> bc3bf4917a89e7c9cf000d20599441013fc48e4b
POWERPLOTPREFIX=inst/plots/PowerPlot
POWERPLOT=$(wildcard $(POWERPLOTPREFIX)_*.pdf)
POWERDAT=data/PowerSimStat95Data.csv
POWERREPLICATIONS=20
POWERSEED=20180910
POWERSIMPREFIX=data/XOUTPowerSimulations_
<<<<<<< HEAD
POWERSIMRDA=$(wildcard $(POWERSIMPREFIX)_*.Rda)
POWERSIMFILEMETA=data/PowerSimulationFileMetadata.csv
POWERSIMTEMPFILEMETAPREFIX=data/PowerSimulationsMetadataDist
POWERSIMTEMPFILEMETA=$(wildcard $(POWERSIMTEMPFILEMETAPREFIX)*.csv)
=======
POWERSIMFILEMETA=data/PowerSimulationFileMetadata.csv
POWERSIMTEMPFILEMETAPREFIX=data/PowerSimulationsMetadataDist
POWERSIMPARAMSUFF=PowerSimulationParameters.R
POWERSIMMETARDASUFF=PowerSimulationParameters.Rda
POWERSIMPARAM=$(wildcard exec/*$(POWERSIMPARAMSUFF))
POWERSIMTEMPFILEMETA:= \
  $(POWERSIMPARAM:exec/%$(POWERSIMPARAMSUFF)=$(POWERSIMTEMPFILEMETAPREFIX)%.csv)
POWERSIMMETARDA:= \
  $(POWERSIMPARAM:exec/%$(POWERSIMPARAMSUFF)=data/%$(POWERSIMMETARDASUFF))
>>>>>>> bc3bf4917a89e7c9cf000d20599441013fc48e4b
POWERSIMSTATMETA=data/PowerSimulationStatsMetadata.csv

LRVPLOTPREFIX=inst/plots/LRVEstPlot
LRVPLOT=$(wildcard $(LRVPLOTPREFIX)*.pdf)
LRVREPLICATIONS=20
LRVSEED=20180912
LRVDAT=data/LRVSimulations.Rda

ZNDAT=data/ZnSimulations.Rda
ZNSIMSEED=20180911
ZNSIMREP=20
ZNCONVPLOTPREFIX=inst/plots/DistConv
ZNCONVPLOT=$(wildcard $(ZNCONVPLOTPREFIX)*.pdf)

CAPMPLOT=inst/plots/BankCAPMChange.pdf
CAPMDAT=data/BankCAPMPValues.Rda
FFFILE=data/FF_factors.csv
BANKFILE=data/Portfolios.csv

.PHONY : all
all : $(POWERPLOT) $(LRVPLOT) $(ZNCONVPLOT) $(CAPMPLOT) inst/Makefile \
      inst/package

<<<<<<< HEAD
$(POWERSIMPREFIX)_garch11_a0.1_b0.7_o0.5_*.Rda : \
                               exec/GARCHPowerSimulationParameters.R \
                               exec/PowerSimulations.R R/ProbabilityFunctions.R
	make package
	$(RSCRIPT) $< -f data/GARCHPowerSimulationParameters.Rda
	$(RSCRIPT) exec/PowerSimulations.R -N $(POWERREPLICATIONS) -s $(POWERSEED) \
		-p $(POWERSIMPREFIX) -i data/GARCHPowerSimulationParameters.Rda \
		-o $(POWERSIMTEMPFILEMETAPREFIX)GARCH.csv

$(POWERSIMPREFIX)_ar1_0.5_*.Rda : exec/AR1PowerSimulationParameters.R \
                                  exec/PowerSimulations.R \
                                  R/ProbabilityFunctions.R
	make package
	$(RSCRIPT) $< -f data/AR1PowerSimulationParameters.Rda
	$(RSCRIPT) exec/PowerSimulations.R -N $(POWERREPLICATIONS) -s $(POWERSEED) \
		-p $(POWERSIMPREFIX) -i data/AR1PowerSimulationParameters.Rda \
		-o $(POWERSIMTEMPFILEMETAPREFIX)AR1.csv

$(POWERSIMPREFIX)_norm_*.Rda : exec/NormPowerSimulationParameters.R \
                               exec/PowerSimulations.R R/ProbabilityFunctions.R
	make package
	$(RSCRIPT) $< -f data/NormPowerSimulationParameters.Rda
	$(RSCRIPT) exec/PowerSimulations.R -N $(POWERREPLICATIONS) -s $(POWERSEED) \
		-p $(POWERSIMPREFIX) -i data/NormPowerSimulationParameters.Rda \
		-o $(POWERSIMTEMPFILEMETAPREFIX)Norm.csv

$(POWERSIMFILEMETA) : $(POWERSIMTEMPFILEMETA) $(POWERSIMRDA)
=======
$(POWERSIMMETARDA) : data/%$(POWERSIMMETARDASUFF) : exec/%$(POWERSIMPARAMSUFF)
	$(RSCRIPT) $< -f $@

$(POWERSIMTEMPFILEMETA) : $(POWERSIMTEMPFILEMETAPREFIX)%.csv : \
  data/%$(POWERSIMMETARDASUFF) exec/PowerSimulations.R R/ProbabilityFunctions.R
	make package
	$(RSCRIPT) exec/PowerSimulations.R -N $(POWERREPLICATIONS) -s $(POWERSEED) \
		-p $(POWERSIMPREFIX) -i $< -o $@

$(POWERSIMFILEMETA) : $(POWERSIMTEMPFILEMETA)
>>>>>>> bc3bf4917a89e7c9cf000d20599441013fc48e4b
	head -1 $< > $@
	for filename in $(POWERSIMTEMPFILEMETA); do sed 1d $$filename >> $@; done

$(POWERDAT) : $(POWERSIMFILEMETA) $(POWERSIMSTATMETA) \
              R/ProbabilityFunctions.R R/SimulationUtils.R \
              exec/PowerSimStatDataGenerator.R
<<<<<<< HEAD
	$(RSCRIPT) exec/PowerSimStatDataGenerator.R -f $(POWERSIMFILEMETA) \
		-s $(POWERSIMSTATMETA) -o $@ -a 0.05
=======
	make package
	$(RSCRIPT) exec/PowerSimStatDataGenerator.R -f $(POWERSIMFILEMETA) \
		-s $(POWERSIMSTATMETA) -o $@ -a $(LEVEL)
>>>>>>> bc3bf4917a89e7c9cf000d20599441013fc48e4b

$(POWERPLOT) : $(POWERDAT) exec/PowerPlot.R R/Plotting.R
	make package
	$(RSCRIPT) exec/PowerPlot.R -f $< -o $(POWERPLOTPREFIX) -v
<<<<<<< HEAD
	mv $(notdir $(basename $@).pdf) $@
	rm *.log *.aux
=======
	mv $(notdir $(POWERPLOTPREFIX))*.pdf $(dir $(POWERPLOTPREFIX))
>>>>>>> bc3bf4917a89e7c9cf000d20599441013fc48e4b

$(LRVDAT) : exec/LRVEstAnalysisParallel.R R/ChangePointTests.R
	make package
	$(RSCRIPT) $< -N $(LRVREPLICATIONS) -s $(LRVSEED) -f $@

$(LRVPLOT) : $(LRVDAT) exec/LRVPlot.R R/Plotting.R
	make package
	$(RSCRIPT) exec/LRVPlot.R -f $< -o $(LRVPLOTPREFIX) -v
<<<<<<< HEAD
	mv $(notdir $(basename $@).pdf) $@
	rm *.log *.aux
=======
	mv $(notdir $(basename $(LRVPLOTPREFIX)))*.pdf $(dir $(LRVPLOTPREFIX))
>>>>>>> bc3bf4917a89e7c9cf000d20599441013fc48e4b

$(ZNDAT) : exec/ZnSimulations.R R/ProbabilityFunctions.R
	make package
	$(RSCRIPT) $< -f $@ -s $(ZNSIMSEED) -r $(ZNSIMREP)

$(ZNCONVPLOT) : $(ZNDAT) exec/DistConvPlot.R R/Plotting.R
	make package
	$(RSCRIPT) exec/DistConvPlot.R -f $< -o $(ZNCONVPLOTPREFIX) -v
<<<<<<< HEAD
	mv $(notdir $(basename $@).pdf) $@
	rm *.log *.aux
=======
	mv $(notdir $(basename $(ZNCONVPLOTPREFIX)))*.pdf $(dir $(ZNCONVPLOTPREFIX))
>>>>>>> bc3bf4917a89e7c9cf000d20599441013fc48e4b

$(CAPMDAT) : exec/BankTestPvalComputeEW.R $(FFFILE) $(BANKFILE) \
             R/ChangePointTests.R R/ProbabilityFunctions.R \
             R/SimulationUtils.R
	make package
	$(RSCRIPT) $< -f $(FFFILE) -b $(BANKFILE) -o $@

$(CAPMPLOT) : $(CAPMDAT) exec/CAPMExamplePlot.R R/Plotting.R
	make package
	$(RSCRIPT) exec/CAPMExamplePlot.R -f $< -o $(basename $@) -v
	mv $(notdir $(basename $@).pdf) $@
<<<<<<< HEAD
	rm *.log *.aux
=======
>>>>>>> bc3bf4917a89e7c9cf000d20599441013fc48e4b

R/ChangePointTests.R : src/ChangePointTests.cpp
	touch $@

inst/Makefile : Makefile
	cp $< $@

inst/package : package
	cp $< $@

package : R/*.R
	$(RSCRIPT) exec/RemakePackage.R
	touch package

.PHONY : clean
clean :
<<<<<<< HEAD
	rm $(POWERPLOT)
	rm $(POWERSIMPREFIX)*
	rm $(POWERSIMFILEMETA)
	rm $(POWERDAT)
	rm $(LRVDAT)
	rm $(LRVPLOT)
	rm $(ZNDAT)
	rm $(ZNCONVPLOT)
	rm $(CAPMDAT)
	rm $(CAPMPLOT)

.PHONY : mostlyclean
mostlyclean :
	rm $(POWERPLOT)
	rm $(LRVPLOT)
	rm $(ZNCONVPLOT)
	rm $(CAPMPLOT)

.PHONY : init
init :
	echo "Empty power plot" > $(POWERPLOTPREFIX)_norm_n50.pdf
=======
	-rm $(POWERPLOT)
	-rm $(POWERSIMPREFIX)*
	-rm $(POWERSIMFILEMETA)
	-rm $(POWERSIMTEMPFILEMETA)
	-rm $(POWERSIMMETARDA)
	-rm $(POWERDAT)
	-rm $(LRVDAT)
	-rm $(LRVPLOT)
	-rm $(ZNDAT)
	-rm $(ZNCONVPLOT)
	-rm $(CAPMDAT)
	-rm $(CAPMPLOT)
	-rm inst/plots/*.tex

.PHONY : mostlyclean
mostlyclean :
	-rm $(POWERPLOT)
	-rm $(LRVPLOT)
	-rm $(ZNCONVPLOT)
	-rm $(CAPMPLOT)
	-rm inst/plots/*.tex

.PHONY : init
init :
	echo "Empty power plot" > $(POWERPLOTPREFIX)_norm_n50_log_c4rt.pdf
>>>>>>> bc3bf4917a89e7c9cf000d20599441013fc48e4b
	echo "Empty LRV plot" > $(LRVPLOTPREFIX)_bartlett_garch_50.pdf
	echo "Empty dist. conv. plot" > $(ZNCONVPLOTPREFIX)_norm_n50_log.pdf

