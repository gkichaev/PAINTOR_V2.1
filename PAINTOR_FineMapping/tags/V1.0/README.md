# PAINTOR

Probabilistic Annotation INtegraTOR 

## Description


We provide a command line implementation of the PAINTOR framework described in Kichaev et al. (2014, under review). Briefly, PAINTOR is a statistical fine-mapping method that integrates functional genomic data with association strength to prioritize variants for follow-up analysis. The software runs on multiple fine-mapping loci simulatenously and takes as input the following data for each set of SNPs at a locus

1. Summary Association Statistics (Z-scores)
2. Linkage Disequlibirium Matrix (Pairwise Pearson correlations coefficients between each SNP)
3. Functional Annotation Matrix (Binary indicator of annotation memebership (i.e. if entry {i,k} = 1, then SNP i is a member of annotation K)  

#### Key Features

1. Outputs a probability for a SNP to be causal which can subsequently be used to prioritize variants
2. Can model multiple causal variants at any risk locus
3. Leverage functional genomic data as a prior probability to improve prioritization
  - This prior probability is not pre-specified, but rather, learned directly from the data via an Empirical Bayes approach
4. Quantify enrichment of causal variants within functional classes
  - Enables users to unbiasedly select functional annotations from a large pool that are most phenotypically relevant 

## Installation
The software has two dependencies: [1] Eigen v3.2 (matrix library) [2] NLopt v2.4.1 (optimization library) which are packaged with PAINTOR in order to simplify installation. Please see the [Eigen homepage](http://eigen.tuxfamily.org/index.php?title=Main_Page) and [NLopt homepage](http://ab-initio.mit.edu/wiki/index.php/NLopt) for more information. 

