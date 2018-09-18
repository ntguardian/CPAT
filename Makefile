SHELL=/bin/sh
RSCRIPT=Rscript

POWERPLOTPREFIX=inst/plots/PowerPlot
POWERPLOT=$(wildcard $(POWERPLOTPREFIX)_*.pdf)
POWERDAT=data/PowerSimStat95Data.csv
POWERREPLICATIONS=20
POWERSEED=20180910
POWERSIMPREFIX=data/XOUTPowerSimulations_
POWERSIMFILEMETA=data/PowerSimulationFileMetadata.csv
POWERSIMTEMPFILEMETAPREFIX=data/PowerSimulationsMetadataDist
POWERSIMPARAMSUFF=PowerSimulationParameters.R
POWERSIMMETARDASUFF=PowerSimulationParameters.Rda
POWERSIMPARAM=$(wildcard exec/*$(POWERSIMPARAMSUFF))
POWERSIMTEMPFILEMETA:= \
  $(POWERSIMPARAM:exec/%$(POWERSIMPARAMSUFF)=$(POWERSIMTEMPFILEMETAPREFIX)%.csv)
POWERSIMMETARDA:= \
  $(POWERSIMPARAM:exec/%$(POWERSIMPARAMSUFF)=data/%$(POWERSIMMETARDASUFF))
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

$(POWERSIMMETARDA) : data/%$(POWERSIMMETARDASUFF) : exec/%$(POWERSIMPARAMSUFF)
	$(RSCRIPT) $< -f $@

$(POWERSIMTEMPFILEMETA) : $(POWERSIMTEMPFILEMETAPREFIX)%.csv : \
  data/%$(POWERSIMMETARDASUFF) exec/PowerSimulations.R R/ProbabilityFunctions.R
	make package
	$(RSCRIPT) exec/PowerSimulations.R -N $(POWERREPLICATIONS) -s $(POWERSEED) \
		-p $(POWERSIMPREFIX) -i $< -o $@

$(POWERSIMFILEMETA) : $(POWERSIMTEMPFILEMETA)
	head -1 $< > $@
	for filename in $(POWERSIMTEMPFILEMETA); do sed 1d $$filename >> $@; done

$(POWERDAT) : $(POWERSIMFILEMETA) $(POWERSIMSTATMETA) \
              R/ProbabilityFunctions.R R/SimulationUtils.R \
              exec/PowerSimStatDataGenerator.R
	make package
	$(RSCRIPT) exec/PowerSimStatDataGenerator.R -f $(POWERSIMFILEMETA) \
		-s $(POWERSIMSTATMETA) -o $@ -a 0.05

$(POWERPLOT) : $(POWERDAT) exec/PowerPlot.R R/Plotting.R
	make package
	$(RSCRIPT) exec/PowerPlot.R -f $< -o $(POWERPLOTPREFIX) -v
	mv $(notdir $(POWERPLOTPREFIX))*.pdf $(dir $(POWERPLOTPREFIX))

$(LRVDAT) : exec/LRVEstAnalysisParallel.R R/ChangePointTests.R
	make package
	$(RSCRIPT) $< -N $(LRVREPLICATIONS) -s $(LRVSEED) -f $@

$(LRVPLOT) : $(LRVDAT) exec/LRVPlot.R R/Plotting.R
	make package
	$(RSCRIPT) exec/LRVPlot.R -f $< -o $(LRVPLOTPREFIX) -v
	mv $(notdir $(basename $@))*.pdf $(dir $@)

$(ZNDAT) : exec/ZnSimulations.R R/ProbabilityFunctions.R
	make package
	$(RSCRIPT) $< -f $@ -s $(ZNSIMSEED) -r $(ZNSIMREP)

$(ZNCONVPLOT) : $(ZNDAT) exec/DistConvPlot.R R/Plotting.R
	make package
	$(RSCRIPT) exec/DistConvPlot.R -f $< -o $(ZNCONVPLOTPREFIX) -v
	mv $(notdir $(basename $@).pdf) $@

$(CAPMDAT) : exec/BankTestPvalComputeEW.R $(FFFILE) $(BANKFILE) \
             R/ChangePointTests.R R/ProbabilityFunctions.R \
             R/SimulationUtils.R
	make package
	$(RSCRIPT) $< -f $(FFFILE) -b $(BANKFILE) -o $@

$(CAPMPLOT) : $(CAPMDAT) exec/CAPMExamplePlot.R R/Plotting.R
	make package
	$(RSCRIPT) exec/CAPMExamplePlot.R -f $< -o $(basename $@) -v
	mv $(notdir $(basename $@).pdf) $@

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
	-rm $(POWERPLOT)
	-rm $(POWERSIMPREFIX)*
	-rm $(POWERSIMFILEMETA)
	-rm $(POWERDAT)
	-rm $(LRVDAT)
	-rm $(LRVPLOT)
	-rm $(ZNDAT)
	-rm $(ZNCONVPLOT)
	-rm $(CAPMDAT)
	-rm $(CAPMPLOT)

.PHONY : mostlyclean
mostlyclean :
	rm $(POWERPLOT)
	rm $(LRVPLOT)
	rm $(ZNCONVPLOT)
	rm $(CAPMPLOT)

.PHONY : init
init :
	echo "Empty power plot" > $(POWERPLOTPREFIX)_norm_n50_log_c4rt.pdf
	echo "Empty LRV plot" > $(LRVPLOTPREFIX)_bartlett_garch_50.pdf
	echo "Empty dist. conv. plot" > $(ZNCONVPLOTPREFIX)_norm_n50_log.pdf

