SHELL=/bin/sh

POWERPLOTPREFIX=inst/plots/power_plot
POWERPLOT=$(wildcard $(POWERPLOTPREFIX)_*.pdf)
POWERDAT=data/PowerSimStat95Data.csv
POWERREPLICATIONS=5000
POWERSEED=20180910
POWERSIMPREFIX=data/XOUTPowerSimulations_
POWERSIMRDA=$(wildcard $(POWERSIMPREFIX)_*.Rda)
POWERSIMFILEMETA=data/PowerSimulationFileMetadata.csv
POWERSIMTEMPFILEMETAPREFIX=data/PowerSimulationsMetadataDist
POWERSIMTEMPFILEMETA=$(wildcard $(POWERSIMTEMPFILEMETAPREFIX)*.csv)
POWERSIMSTATMETA=data/PowerSimulationStatsMetadata.csv

LRVPLOTPREFIX=inst/plots/LRVEstPlot
LRVPLOT=$(wildcard $(LRVPLOTPREFIX)*.pdf)
LRVREPLICATIONS=10000
LRVSEED=20180912
LRVDAT=data/LRVSimulations.Rda

ZNDAT=data/ZnSimulations.Rda
ZNSIMSEED=20180911
ZNSIMREP=100000
ZNCONVPLOTPREFIX=inst/plots/DistConv
ZNCONVPLOT=$(wildcard $(ZNCONVPLOTPREFIX)*.pdf)

CAPMPLOT=inst/plots/BankCAPMChange.pdf
CAPMDAT=BankCAPMPValues.Rda
FFFILE=data/FF_factors.csv
BANKFILE=data/Portfolios.csv

.PHONY : all
all : $(POWERPLOT) $(LRVPLOT) $(ZNCONVPLOT) $(CAPMPLOT) inst/Makefile

$(POWERSIMPREFIX)_garch11_a0.1_b0.7_o0.5_*.Rda : \
                               exec/GARCHPowerSimulationParameters.R \
                               exec/PowerSimulations.R R/ProbabilityFunctions.R
	Rscript exec/RemakePackage.R
	Rscript $< -f data/GARCHPowerSimulationParameters.Rda
	Rscript exec/PowerSimulations.R -N $(POWERREPLICATIONS) -s $(POWERSEED) \
		-p $(POWERSIMPREFIX) -i data/GARCHPowerSimulationParameters.Rda \
		-o $(POWERSIMTEMPFILEMETAPREFIX)GARCH.csv

$(POWERSIMPREFIX)_ar1_0.5_*.Rda : exec/AR1PowerSimulationParameters.R \
                                  exec/PowerSimulations.R \
                                  R/ProbabilityFunctions.R
	Rscript exec/RemakePackage.R
	Rscript $< -f data/AR1PowerSimulationParameters.Rda
	Rscript exec/PowerSimulations.R -N $(POWERREPLICATIONS) -s $(POWERSEED) \
		-p $(POWERSIMPREFIX) -i data/AR1PowerSimulationParameters.Rda \
		-o $(POWERSIMTEMPFILEMETAPREFIX)AR1.csv

$(POWERSIMPREFIX)_norm_*.Rda : exec/NormPowerSimulationParameters.R \
                               exec/PowerSimulations.R R/ProbabilityFunctions.R
	Rscript exec/RemakePackage.R
	Rscript $< -f data/NormPowerSimulationParameters.Rda
	Rscript exec/PowerSimulations.R -N $(POWERREPLICATIONS) -s $(POWERSEED) \
		-p $(POWERSIMPREFIX) -i data/NormPowerSimulationParameters.Rda \
		-o $(POWERSIMTEMPFILEMETAPREFIX)Norm.csv

$(POWERSIMFILEMETA) : $(POWERSIMTEMPFILEMETA) $(POWERSIMRDA)
	head -1 $< > $@
	for filename in $(POWERSIMTEMPFILEMETA); do sed 1d $$filename >> $@; done

$(POWERDAT) : $(POWERSIMFILEMETA) $(POWERSIMSTATMETA) \
              R/ProbabilityFunctions.R R/SimulationUtils.R \
              exec/PowerSimStatDataGenerator.R
	Rscript exec/PowerSimStatDataGenerator.R -f $(POWERSIMFILEMETA) \
		-s $(POWERSIMSTATMETA) -o $@ -a 0.05

$(POWERPLOT) : $(POWERDAT) exec/PowerPlot.R R/Plotting.R
	Rscript exec/RemakePackage.R
	Rscript exec/PowerPlot.R -f $< -o $(POWERPLOTPREFIX) -v
	mv $(basename $@) $@

$(LRVDAT) : exec/LRVEstAnalysisParallel.R R/ChangePointTests.R
	Rscript exec/RemakePackage.R
	Rscript $< -N $(LRVREPLICATIONS) -s $(LRVSEED) -f $@

$(LRVPLOT) : $(LRVDAT) exec/LRVPlot.R R/Plotting.R
	Rscript exec/RemakePackage.R
	Rscript exec/LRVPlot.R -f $< -o $(LRVPLOTPREFIX) -v
	mv $(basename $@) $@

$(ZNDAT) : exec/ZnSimulations.R R/ProbabilityFunctions.R
	Rscript exec/RemakePackage.R
	Rscript $< -f $@ -s $(ZNSIMSEED) -r $(ZNSIMREP)

$(ZNCONVPLOT) : $(ZNDAT) exec/DistConvPlot.R R/Plotting.R
	Rscript exec/RemakePackage.R
	Rscript exec/DistConvPlot.R -f $< -o $(ZNCONVPLOTPREFIX) -v
	mv $(basename $@) $@

$(CAPMDAT) : exec/BankTestPvalComputeEW.R $(FFFILE) $(BANKFILE) \
             R/ChangePointTests.R R/ProbabilityFunctions.R \
             R/SimulationUtils.R
	Rscript exec/RemakePackage.R
	Rscript $< -f $(FFFIILE) -b $(BANKFILE) -o $@

$(CAPMPLOT) : $(CAPMDAT) exec/CAPMExamplePlot.R R/Plotting.R
	Rscript exec/RemakePackage.R
	Rscript exec/CAPMExamplePlot.R -f $< -o $(basename $@) -v
	mv $(basename $@) $@

R/ChangePointTests.R : src/ChangePointTests.cpp
	touch $@

inst/Makefile : Makefile
	cp $< $@

.PHONY : clean
clean:
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
mostlyclean:
	rm $(POWERPLOT)
	rm $(LRVPLOT)
	rm $(ZNCONVPLOT)
	rm $(CAPMPLOT)
