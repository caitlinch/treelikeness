# A test statistic to quantify trelikeness in phylogenetics
#### Caitlin Cherryh, Bui Quang Minh, Rob Lanfear
##### March 6th 2020

This repository contains code to simulate a series of sequence alignments and estimate a series of treelikeness and introgression metrics, including a new test for treelikeness called tree proportion. Tree proportion represents the proportion of information from an alignment that is included in the phylogenetic tree estimated from that alignment. 

Our preprint with more details about our methods and results can be found [here](https://www.biorxiv.org/content/10.1101/2021.02.16.431544v1)

***
### How to calculate the tree proportion
**1. Load the function**

The tree proportion function is in the file `treelikeness/code/func_test_statistic.R`. 

```{r}
# If I cloned the treelikeness repository to the location "~/Repos"
source("~/Repos/treelikeness/code/func_test_statistic.R")
```
**2. Prepare filepaths**

There are two example alignments in the folder `treelikeness/worked_example/`:

* `one_tree_example_alignment.nexus`: 1000 base pairs simulated along a balanced 8 taxon tree (treelike)
* `two_tree_example_alignment.nexus`: 1000 base pairs simulated along a balanced 8 taxon tree containing one divergent introgression event

You will also need to have SplitsTree and IQ-Tree downloaded and installed. 

* On Linux, the path to SplitsTree is `path/to/splitstree4/SplitsTree`
* On MacOS, it will be `path/to/splitstree4/SplitsTree.app/Contents/MacOS/JavaApplicationStub`

```{r}
# Paths to programs (if I store the programs in the location "~/Executables")
splitstree_path <- "~/Executables/SplitsTree.app/Contents/MacOS/JavaApplicationStub"
iqtree_path <- "~/Executables/iqtree-2.0/bin/iqtree"

# Paths to alignments
one_tree_alignment <- "~/Repos/treelikeness/worked_example/one_tree_example_alignment.nexus" 
two_tree_alignment <- "~/Repos/treelikeness/worked_example/two_tree_example_alignment.nexus" 
```

**3. Calculate the tree proportion**

```{r}
tree.proportion(iqpath, splitstree_path, path, network_algorithm = "neighbournet", trimmed = TRUE, 
                tree_path = FALSE, run_IQTREE = FALSE, seq_type)
```

There are a few parameters to the tree proportion function: 

* `iqpath`: path to IQ-Tree program
* `splitstree_path`: path to SplitsTree program
* `path`: path to alignment. Currently takes FASTA (.fasta, .fa, .fna, .ffn, .faa, .frn, .fas), PHYLIP (.phy) and Nexus (.nexus, .nex) files.
* `network_algorithm`: Determines which method is used to calculate a network from the alignment. Can be set to `"split decomposition"` or `"neighbournet"` (default). We used `"neighbornet"` for all our analyses.
* `trimmed`: Whether to remove trivial splits (i.e. terminal branches on a phylogeny). Can be `TRUE` (default) or `FALSE`. We used `TRUE` for all our analyses.
* `tree_path`: Tree proportion compares a tree to a network, to calculate the proportion of splits from the network that are included in the tree. You can either provide a tree or allow the function to estimate a tree from your alignment using IQ-Tree. This argument can be `FALSE` or a file path to a tree (the tree must have the same set of taxa as the alignment). 
    + If you are not providing a tree, set `tree_path = FALSE, run_IQTree = TRUE`. The function will then run IQ-Tree to estimate a tree from the alignment and use that tree to calculate the tree proportion. 
    + If you want to provide a specific tree, set `tree_path = "/path/to/tree.nex", run_IQTREE = FALSE`, and the function will use the specified tree. 
* `run_IQTREE`: Whether to run IQ-Tree to estimate a tree from the alignment or not. Can be `TRUE` (meaning the function runs IQ-Tree on the alignment to estimate a tree) or `FALSE` (meaning IQ-Tree will not be run)
* `seq_type`: can be `"dna"` or `"protein"`

To calculate the tree proportion for the two example alignments:
```{r}
tp_one_tree <- tree.proportion(iqtree_path, splitstree_path, one_tree_alignment, 
                network_algorithm = "neighbournet", trimmed = TRUE, 
                tree_path = FALSE, run_IQTREE = TRUE, seq_type = "dna")
                
tp_two_trees <- tree.proportion(iqtree_path, splitstree_path, two_tree_alignment, 
                network_algorithm = "neighbournet", trimmed = TRUE, 
                tree_path = FALSE, run_IQTREE = TRUE, seq_type = "dna")
```

Calculating the tree proportion will generate a number of files:

* a "XXXXX_blocks.nexus" file (where XXXXX is the name of the alignment) - contains a nexus formatted version of the alignment with a taxa block (needed for SplitsTree)
* a "XXXXX_splits.nex" file (where XXXXX is the name of the alignment) - contains the splits calculated in SplitsTree
* any files generated from an IQ-Tree run

***
### Instructions to reproduce the simulations and analyses
1. Clone the `caitlinch/treelikeness` repo
    * The repo contains 2 folders: `code` and `trees`
        * `code`: contains all code necessary to run or replicate the analysis.
            * Numbered files are a step in the pipeline. To replicate the analysis, run each of the three numbered files sequentially:
                * `1_simulateAlignmentsRunTestStatistics.R`
                * `2_CollateData.R`
                * `3_PublicationsPlots.R`
            * The `func` prefix indicates files that contain functions used for simulation and analysis
        * `trees`: contains .txt files of the Newick trees used in the simulations.
2. Download necessary software:
    * IQ-Tree (http://www.iqtree.org/)
    * 3SEQ (http://mol.ax/software/3seq/)
    * PhiPack (https://www.maths.otago.ac.nz/~dbryant/software.html)
    * SplitsTree v4 (https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/splitstree/)
3. Prepare scripts
    * Create an output folder (to store alignments and other files that are created during the simulations) and a results folder (to store the test statistic and statistical test results and the plots)
    * Open scripts numbered 1-3 in the `code` folder. Open each file, go to `Step 2` and update the file paths for your own machine.
4. Run Part 1
    * Determine how many cores to use and set `num_cores`. `num_cores = 1` means the script will run entirely sequentially. Edit this variable to increase the number of cores and therefore the number of simultaneous analyses
    * By default, the number of parametric bootstrap replicates performed is 199. To change this, edit the `num_reps` variable
    * Run the file `1_simulateAlignmentsRunTestStatistics.R`
        * This will generate all of the simulation data for all experiments, and apply each test statistic or metric to each alignment
        * This will perform a parametric bootstrap on the simulated alignments to calculate the statistical test for tree proportion.
5. Run Part 2
    * Run the file `2_CollateData.R`
        * The raw test statistic data has been collated into one csv file per experiment. The p-value data has also been collated into one csv file per experiment.
6. Run Part 3
    * Run the file `3_PublicationPlots.R`
        * This will generate a number of plots depicting the results of the experiments

***
### Attribution and citations
If you use this method, please reference:

1. This repository (http://github.com/caitlinch/treelikeness)
2. The programs and methods used within this analysis
    * IQ-Tree
        * If using IQ-Tree 1:
            * Nguyen L-T, Schmidt HA, von Haeseler A, Minh BQ. 2015. IQ-TREE: A fast and effective stochastic algorithm for estimating maximum likelihood phylogenies, _Mol. Biol. Evol._, 32(1):268-274.
        * If using IQ-Tree 2:
            * Minh BQ, Schmidt HA, Chernomor O, Schrempf D, Woodhams MD, von Haeseler A, Lanfear R. 2020. IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era, _Mol. Biol. Evol._, 37(5):1530-1534
        * __Additionally__, cite ModelFinder (used for model selection) 
            * Kalyaanamoorthy S, Minh BQ, Wong TKF, von Haeseler A, Jermiin LS. 2017. ModelFinder: Fast model selection for accurate phylogenetic estimates. _Nat. Methods_, 14:587-589.
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
    * Likelihood mapping
        * When using likelihood mapping:
            * Strimmer K, von Haeseler A. 1997. Likelihood-mapping: a simple method to visualize phylogenetic content of a sequence alignment. _Proc Natl Acad Sci U S A_, 94(13):6815-6819.
        * When using the IQTree 1 implementation of likelihood mapping, also cite:
            * Nguyen L-T, Schmidt HA, von Haeseler A, Minh BQ. 2015. IQ-TREE: A fast and effective stochastic algorithm for estimating maximum likelihood phylogenies, _Mol. Biol. Evol._, 32(1):268-274.
        * When using the IQTree 2 implementation of likelihood mapping, also cite:
            * Minh BQ, Schmidt HA, Chernomor O, Schrempf D, Woodhams MD, von Haeseler A, Lanfear R. 2020. IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era, _Mol. Biol. Evol._, 37(5):1530-1534
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




