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
    * Open scripts numbered 1-4 in the `code` folder. In each file, go to Step 2 and update the file paths for your own machine
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
    * Run the file `3_PlotData.R`
        * This will generate a number of plots depicting the results of the experiments

***
### Attribution and citations
If you use this method, please reference:

1. This repository (http://github.com/caitlinch/treelikeness)
2. The programs used within this analysis
    * IQ-Tree
        * If using IQ-Tree 1:
            * L.-T. Nguyen, H.A. Schmidt, A. von Haeseler, B.Q. Minh (2015) IQ-TREE: A fast and effective stochastic algorithm for estimating maximum likelihood phylogenies, _Mol. Biol. Evol._, 32:268-274. https://doi.org/10.1093/molbev/msu300
        * If using IQ-Tree 2:
            * B.Q. Minh, H.A. Schmidt, O. Chernomor, D. Schrempf, M.D. Woodhams, A. von Haeseler, R. Lanfear (2019) IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era, _Mol. Biol. Evol._, in press:. https://doi.org/10.1093/molbev/msaa015
        * __Additionally__, cite ModelFinder (used for model selection) 
            * S. Kalyaanamoorthy, B.Q. Minh, T.K.F. Wong, A. von Haeseler, L.S. Jermiin (2017) ModelFinder: Fast model selection for accurate phylogenetic estimates. _Nat. Methods_, 14:587-589. https://doi.org/10.1038/nmeth.4285
    * 3SEQ
        * When using 3SEQ: 
            * Lam HM, Ratmann O, Boni MF.  Improved algorithmic complexity for the 3SEQ recombination detection algorithm.  _Mol Biol Evol_, 35(1):247-251, 2018.
        * When referring to core parts of the statistic:
            * Boni MF, Posada D, Feldman MW.  An exact nonparametric method for inferring mosaic structure in sequence triplets.  _Genetics_, 176:1035-1047, 2007.
    * PhiPack
        * When using PhiPack:
            * Bruen, T., Phillipe, H. and Bryant, D. 2006. A quick and robust statistical test to detect the presence of recombination. _Genetics_ 172, 2665--2681
    * SplitsTree
        * When using SplitsTree:
            * D. H. Huson and D. Bryant. Application of phylogenetic networks in evolutionary studies.
_Molecular Biology and Evolution_, 23:254–267, 2006.



