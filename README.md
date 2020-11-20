# A new test statistic for treelikeness
##### March 6th 2020

This repository contains code to simulate a series of sequence alignments and estimate a series of treelikeness and introgression metrics, including a new test for treelikeness (tree proportion). 

***
### Instructions to reproduce the simulations and analyses
1. Clone the `caitlinch/treelikeness` repo
    * The repo contains 2 folders: `code` and `trees`
        * `code`: contains all code necessary to run or replicate the analysis.
            * Numbered files are a step in the pipeline. To replicate the analysis, run each of the four numbered files sequentially:
                * `1_simulateAlignmentsRunTestStatistics.R`
                * `2_CollateData.R`
                * `3_PlotData.R`
            * The `func` prefix indicates files that contain functions used for simulation and analysis
        * `trees`: contains .txt files of the Newick trees used to fix topology in the simulations.
2. Download necessary software:
    * IQ-Tree (http://www.iqtree.org/)
    * 3SEQ (http://mol.ax/software/3seq/)
    * PhiPack (https://www.maths.otago.ac.nz/~dbryant/software.html)
    * SplitsTree v4 (https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/splitstree/)
3. Prepare scripts
    * Create an output folder (to store alignments and other files that are created during the simulations) and a results folder (to store the test statistic and statistical test results and the plots)
    * Open scripts numbered 1-4 in the `code` folder. Open each file, go to `Step 2` and update the file paths for your own machine
4. Run Part 1
    * Determine how many cores to use and set `num_cores`. `num_cores = 1` means the script will run entirely sequentially. Edit this variable to increae the number of cores and therefore the number of simultaneous analyses
    * By default, the number of parametric bootstrap replicates performed is 199. To change this, edit the `num_reps` variable
    * Run the file `1_simulateAlignmentsRunTestStatistics.R`
        * This will generate all of the simulation data for all three experiments, and apply each test statistic or metric to each alignment
        * This will perform a parametric bootstrap on the experiment 3 simulated alignments to determine whether the test statistics perform well as a statistical test.
5. Run Part 2
    * Run the file `2_CollateData.R`
        * The raw test statistic data has been collated into one csv file per experiment, plus one csv for the p values from the parametric bootstrap performed in experiment 3.
6. Run Part 3
    * Run the file `3_PublicationPlots.R`
        * This will generate a number of plots depicting the results of the experiments

***
### Attribution and citations
If you use this method, please reference:

1. This repository (http://github.com/caitlinch/treelikeness)
2. The programs used within this analysis
    * IQ-Tree
        * If using IQ-Tree 1:
            * Nguyen L-T, Schmidt HA, von Haeseler A, Minh BQ. 2015. IQ-TREE: A fast and effective stochastic algorithm for estimating maximum likelihood phylogenies, _Mol. Biol. Evol._, 32(1):268-274. https://doi.org/10.1093/molbev/msu300
        * If using IQ-Tree 2 and when using likelihood mapping implementation in IQ-Tree:
            * Minh BQ, Schmidt HA, Chernomor O, Schrempf D, Woodhams MD, von Haeseler A, Lanfear R. 2020. IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era, _Mol. Biol. Evol._, 37(5):1530-1534
        * __Additionally__, cite ModelFinder (used for model selection) 
            * Kalyaanamoorthy S, Minh BQ, Wong TKF, von Haeseler A, Jermiin LS. 2017. ModelFinder: Fast model selection for accurate phylogenetic estimates. _Nat. Methods_, 14:587-589. https://doi.org/10.1038/nmeth.4285
    * δ plots
        * When using δ plots:
            * Holland BR , Huber KT, Dress A, Moulton V. 2002. δ Plots: A Tool for Analyzing Phylogenetic Distance Data, _Mol Biol Evol_, 19(12):2051–2059
        * When using the delta.plot function from the R package "ape"
            * Paradis E, Claude J, Strimmer K. 2004. APE: analyses of phylogenetics and evolution in R language. _Bioinformatics_ 20(2):289-290.
    * 3SEQ
        * When using 3SEQ: 
            * Lam HM, Ratmann O, Boni MF. 2018. Improved algorithmic complexity for the 3SEQ recombination detection algorithm.  _Mol Biol Evol_, 35(1):247-251
        * When referring to core parts of the statistic:
            * Boni MF, Posada D, Feldman MW. 2007. An exact nonparametric method for inferring mosaic structure in sequence triplets.  _Genetics_, 176(2):1035-1047
    * PhiPack
        * When using PhiPack:
            * Bruen T, Phillipe H and Bryant D. 2006. A quick and robust statistical test to detect the presence of recombination. _Genetics_ 172(4):2665-2681
    * SplitsTree
        * When using SplitsTree:
            * Huson DH and Bryant D. 2006. Application of phylogenetic networks in evolutionary studies. _Mol Biol Evol_, 23(2):254–267
3. The R packages used for this analysis
    * ape
        * Paradis E, Claude J, Strimmer K. 2004. APE: analyses of phylogenetics and evolution in R language. _Bioinformatics_ 20(2):289-290.
    * ggplot2
        * Wickham H. 2016. ggplot2: Elegant Graphics for Data Analysis. New York: Springer-Verlag.
    * phangorn
        * Schliep KP. 2011. phangorn: phylogenetic analysis in R. _Bioinformatics_ 27(4):592-592.
    * phytools
        * Revell LJ. 2012. phytools: An R package for phylogenetic comparative biology (and other things). _Methods in Ecology and Evolution_ 3:217-223.
    * seqinr
        * Charif D, Lobry JR. 2007. SeqinR 1.0-2: a contributed package to the R project for statistical computing devoted to biological sequences retrieval and analysis. In: Bastolla U, Porto M, Roman HE, Vendruscolo M, editors. Structural approaches to sequence evolution: Molecules, networks, populations. New York: Springer Verlag. p. 207-232.
    * stringR
        * Wickham H. 2019. stringr: Simple, Consistent Wrappers for Common String Operations. R package version 1.4.0. https://CRAN.R-project.org/package=stringr.
    * TreeSim
        * Stadler T. 2017. TreeSim: simulating phylogenetic trees. R package version 2.4. http://CRAN.R-project.org/package=TreeSim 




