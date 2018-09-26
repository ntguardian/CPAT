SHELL=/bin/sh
RSCRIPT=Rscript

SMALLREPLICATIONS=20

LEVEL=0.05    # Level of test
POWERPLOTPREFIX=inst/plots/PowerPlot    # Prefix to files containing power plots
POWERPLOT=$(wildcard $(POWERPLOTPREFIX)_*.pdf)    # All power plot files
POWERDAT=inst/extdata/PowerSimStat95Data.csv    # Power data file
POWERREPLICATIONS=5000    # Number of replications for each power experiment
POWERSEED=20180910    # Seed for power simulation experiments
POWERSIMPREFIX=data/XOUTPowerSimulations_    # Prefix for simulation files
# Power plot file metadata
POWERSIMFILEMETA=inst/extdata/PowerSimulationFileMetadata.csv
# Prefix for files that hold only some file metadata
POWERSIMTEMPFILEMETAPREFIX=inst/extdata/PowerSimulationsMetadataDist
POWERSIMPARAMSUFF=PowerSimulationParameters.R    # Suffix for parameter R files
POWERSIMMETARDASUFF=PowerSimulationParameters.Rda    # Suffix for param. Rda
POWERSIMPARAM=$(wildcard exec/*$(POWERSIMPARAMSUFF))    # All R param. files
POWERSIMTEMPFILEMETA::= \    # Name of CSV files holding only some file meta.
  $(POWERSIMPARAM:exec/%$(POWERSIMPARAMSUFF)=$(POWERSIMTEMPFILEMETAPREFIX)%.csv)
POWERSIMMETARDA::= \    # Name of all Rda files holding simulation metadata
  $(POWERSIMPARAM:exec/%$(POWERSIMPARAMSUFF)=data/%$(POWERSIMMETARDASUFF))
# Power simulation statistic metadata
POWERSIMSTATMETA=inst/extdata/PowerSimulationStatsMetadata.csv

LRVPLOTPREFIX=inst/plots/LRVEstPlot    # Prefix for LRV simulation plots
LRVPLOT=$(wildcard $(LRVPLOTPREFIX)*.pdf)    # All LRV plot files
LRVREPLICATIONS=10000    # Number of replications for each LRV experiment
LRVSEED=20180912    # Seed for LRV simulations
LRVDAT=data/LRVSimulations.Rda    # File holding LRV simulations

ZNDAT=data/ZnSimulations.Rda    # File holding Rényi-type convergence simulation
ZNSIMSEED=20180911    # Seed for convergence simulations
ZNSIMREP=100000    # Number of replications for each convergence experiment
ZNCONVPLOTPREFIX=inst/plots/DistConv    # Prefix for convergence plot files
ZNCONVPLOT=$(wildcard $(ZNCONVPLOTPREFIX)*.pdf)    # All convergence sim. plots

CAPMPLOT=inst/plots/BankCAPMChange.pdf    # File name for CAPM example plot
CAPMDAT=data/BankCAPMPValues.Rda    # Rda file containing CAPM test results
FFFILE=data/ff.rda    # rda file containing Fama-French 5-factor data
BANKFILE=data/banks.rda    # rda file containing bank portfolio data

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
		-s $(POWERSIMSTATMETA) -o $@ -a $(LEVEL)

$(POWERPLOT) : $(POWERDAT) exec/PowerPlot.R R/Plotting.R
	make package
	-mkdir $(dir POWERPLOTPREFIX)
	$(RSCRIPT) exec/PowerPlot.R -f $< -o $(POWERPLOTPREFIX) -v
	mv $(notdir $(POWERPLOTPREFIX))*.pdf $(dir $(POWERPLOTPREFIX))

$(LRVDAT) : exec/LRVEstAnalysisParallel.R R/ChangePointTests.R
	make package
	$(RSCRIPT) $< -N $(LRVREPLICATIONS) -s $(LRVSEED) -f $@

$(LRVPLOT) : $(LRVDAT) exec/LRVPlot.R R/Plotting.R
	make package
	-mkdir $(dir LRVPLOTPREFIX)
	$(RSCRIPT) exec/LRVPlot.R -f $< -o $(LRVPLOTPREFIX) -v
	mv $(notdir $(basename $(LRVPLOTPREFIX)))*.pdf $(dir $(LRVPLOTPREFIX))

$(ZNDAT) : exec/ZnSimulations.R R/ProbabilityFunctions.R
	make package
	$(RSCRIPT) $< -f $@ -s $(ZNSIMSEED) -r $(ZNSIMREP) -v

$(ZNCONVPLOT) : $(ZNDAT) exec/DistConvPlot.R R/Plotting.R
	make package
	-mkdir $(dir ZNCONVPLOT)
	$(RSCRIPT) exec/DistConvPlot.R -f $< -o $(ZNCONVPLOTPREFIX) -v
	mv $(notdir $(basename $(ZNCONVPLOTPREFIX)))*.pdf $(dir $(ZNCONVPLOTPREFIX))

$(CAPMDAT) : exec/BankTestPvalComputeEW.R $(FFFILE) $(BANKFILE) \
             R/ChangePointTests.R R/ProbabilityFunctions.R \
             R/SimulationUtils.R
	make package
	$(RSCRIPT) $< -o $@

$(CAPMPLOT) : $(CAPMDAT) exec/CAPMExamplePlot.R R/Plotting.R
	make package
	-mkdir $(dir CAPMPLOT)
	$(RSCRIPT) exec/CAPMExamplePlot.R -f $< -o $(basename $@) -v
	mv $(notdir $(basename $@).pdf) $@

R/ChangePointTests.R : src/ChangePointTests.cpp
	touch $@

inst/Makefile : Makefile
	cp $< $@

inst/package : package
	cp $< $@

package : R/*.R src/*.cpp $(FFFILE) $(BANKFILE)
	$(RSCRIPT) exec/RemakePackage.R
	touch package

.PHONY : clean
clean :
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
	make initpower
	make initlrv
	make initconv

.PHONY : initpower
initpower :
	-mkdir $(dir $(POWERPLOTPREFIX))
	-mkdir $(dir $(POWERSIMSTATMETA))
	echo "Empty power plot" > $(POWERPLOTPREFIX)_norm_n50_log_c4rt.pdf

.PHONY : initlrv
initlrv :
	-mkdir $(dir $(LRVPLOTPREFIX))
	echo "Empty LRV plot" > $(LRVPLOTPREFIX)_bartlett_garch_50.pdf

.PHONY : initconv
initconv :
	-mkdir $(dir $(ZNCONVPLOTPREFIX))
	echo "Empty dist. conv. plot" > $(ZNCONVPLOTPREFIX)_norm_n50_log.pdf

.PHONY : small
small :
	make POWERREPLICATIONS=$(SMALLREPLICATIONS) \
		LRVREPLICATIONS=$(SMALLREPLICATIONS) \
		ZNSIMREP=$(SMALLREPLICATIONS)

.PHONY : dependencies
dependencies :
	$(RSCRIPT) exec/GetPackages.R
