# PAINTOR
Probabilistic Annotation INtegraTOR

## UPDATE 06/8/15
Announcing PAINTOR v2.0! We have expanded PAINTOR functionality to conduct fine-mappng with multi-ethnic cohorts. This latest release inclucdes updates to the input formats and command line flags. Please see the new User Manual for complete details.  

## Description

We provide a command line implementation of the PAINTOR and PAINTOR Trans-ethnic frameworks described in [Kichaev et al. (PLOS Genetics, 2014)](http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1004722) and Kichaev and Pasaniuc (In Review, 2015). Briefly, PAINTOR is a statistical fine-mapping method that integrates functional genomic data with association strength from potentially multiple populations to prioritize variants for follow-up analysis. The software runs on multiple fine-mapping loci and/or populations simulatenously and takes as input the following data for each set of SNPs at a locus


1. Summary Association Statistics (Z-scores)
2. Linkage Disequlibirium Matrix/Matrices (Pairwise Pearson correlations coefficients between each SNP)
3. Functional Annotation Matrix (Binary indicator of annotation memebership (i.e. if entry {i,k} = 1, then SNP i is a member of annotation K)

#### Key Features

1. Outputs a probability for a SNP to be causal which can subsequently be used to prioritize variants
2. Can model multiple causal variants at any risk locus
3. Leverage functional genomic data as a prior probability to improve prioritization
  - This prior probability is not pre-specified, but rather, learned directly from the data via an Empirical Bayes approach
4. Quantify enrichment of causal variants within functional classes
  - Enables users to unbiasedly select from a (potentially) large pool functional annotations that are most phenotypically relevant
5. Model population-specific LD patterns. 

## Installation
The software has two dependencies: [1] Eigen v3.2 (matrix library) [2] NLopt v2.4.2 (optimization library) which are packaged with PAINTOR in order to simplify installation. Please see the [Eigen homepage](http://eigen.tuxfamily.org/index.php?title=Main_Page) and [NLopt homepage](http://ab-initio.mit.edu/wiki/index.php/NLopt) for more information.

Download the latest [version](https://github.com/gkichaev/PAINTOR_FineMapping/releases) of the software into your desired target directory. Then unpack and install the software with the following commands:

`$ tar -xvf PAINTOR_FineMapping-2.0.tar`  
`$ cd PAINTOR_FineMapping`  
`$ bash install.sh`  

This will create an executable "PAINTOR". Sample data is provided with the package. To test that the installation worked properly, type:

`$ ./PAINTOR -input SampleData/input.files -in SampleData/ -out SampleData/ -Zhead Zscore.p1,Zscore.p2,Zscore.p3 -LDname LD.p1,LD.p2,LD.p3 -c 2 -annotations Coding,DHS i`

If everything worked correctly the final log-likelihod should be: `-2129.21235`

For quick start simply type:

`$ ./PAINTOR` 

For detailed information on input files and command line flags see the user manual provided.

