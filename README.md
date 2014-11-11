# PAINTOR
Probabilistic Annotation INtegraTOR

## UPDATE 11/10/14
PAINTOR V1.1 has major updates to improve performance and usablity. We now perform a Cholesky Decomposition of the covariance matrix to increase computational efficiency and stability while simultaneously reducing the overall memory footprint. In addition, we have modified the file input routine to make it more user friendly. See updated User Manual for details.   

## Description

We provide a command line implementation of the PAINTOR framework described in [Kichaev et al. (PLOS Genetics, 2014)](http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1004722). Briefly, PAINTOR is a statistical fine-mapping method that integrates functional genomic data with association strength to prioritize variants for follow-up analysis. The software runs on multiple fine-mapping loci simulatenously and takes as input the following data for each set of SNPs at a locus


1. Summary Association Statistics (Z-scores)
2. Linkage Disequlibirium Matrix (Pairwise Pearson correlations coefficients between each SNP)
3. Functional Annotation Matrix (Binary indicator of annotation memebership (i.e. if entry {i,k} = 1, then SNP i is a member of annotation K)

#### Key Features

1. Outputs a probability for a SNP to be causal which can subsequently be used to prioritize variants
2. Can model multiple causal variants at any risk locus
3. Leverage functional genomic data as a prior probability to improve prioritization
  - This prior probability is not pre-specified, but rather, learned directly from the data via an Empirical Bayes approach
4. Quantify enrichment of causal variants within functional classes
  - Enables users to unbiasedly select from a (potentially) large pool functional annotations that are most phenotypically relevant

## Installation
The software has two dependencies: [1] Eigen v3.2 (matrix library) [2] NLopt v2.4.1 (optimization library) which are packaged with PAINTOR in order to simplify installation. Please see the [Eigen homepage](http://eigen.tuxfamily.org/index.php?title=Main_Page) and [NLopt homepage](http://ab-initio.mit.edu/wiki/index.php/NLopt) for more information.

Download the latest [version](https://github.com/gkichaev/PAINTOR_FineMapping/releases) of the software into your desired target directory. Then unpack and install the software with the following commands:

`$ tar -xvf PAINTOR_FineMapping-1.1.tar`  
`$ cd PAINTOR`  
`$ bash install.sh`  

This will create an executable PAINTOR.exe. Sample data is provided with the package. To test that the installation worked properly, type:

`$ ./PAINTOR -o ./SampleData/ -d ./SampleData/ -input input.files -c 3 -i 1`

If everything worked correctly the final log-likelihod should be: `-977.846493`

For information on input files and command line flags see the user manual provided.
