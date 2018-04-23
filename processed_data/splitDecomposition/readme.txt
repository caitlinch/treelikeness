{\rtf1\ansi\ansicpg1252\cocoartf1561\cocoasubrtf200
{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\paperw11900\paperh16840\margl1440\margr1440\vieww10800\viewh8400\viewkind0
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0

\f0\fs24 \cf0 All alignments from roblanf/BenchmarkAlignments GitHub. Each has had the split decomposition function run on them.\
\
mldist files - output distance matrix from IQ-Tree.\
splitDecomposition files - output from my R code.\
To find out how much information was lost (unable to be included in a split), subtract final matrix in splitDecomposition file from mldist matrix.\
\
splitDecomposition file format:\
[[1]] <- list of splits\
[[1]][[x]] <- information for one split\
[[1]][[x]]$isolation_index <- split weight (weight of that d-split)\
[[1]][[x]]$partition <- the partition of the set of all taxa into two subsets. Each part of the partition (a and b) is monophyletic.\
[[1]][[x]]$matrix <- the matrix for the split (should be 0 where the species for the the row and the species for the column are in the same subset, otherwise the cell will equal the isolation index - the distance between the two subsets)\
\
[[2]] <- matrix obtained by summing the splits matrix for each split.}