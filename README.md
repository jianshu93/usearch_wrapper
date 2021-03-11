# usearch_wrapper: wsearch.bash
usearch wrapper for amplicon analysis. Support for 16S and ITS. This is still under active development and will be ready in a few months. Please contact jianshu(jianshu.zhao@gatech.edu)
This shell script has the following dependencies:

1. Falco v0.2.2, a C++ implentation of FastQC, way faster and user friendly than FastQC and more importantly, it is under active development (https://github.com/smithlabcode/falco)
2. Vsearch and Usearch v11 (https://drive5.com/usearch/download.html),both 32 bit and 64 bit works but we strongly suggest 64 bit for large dataset. for example when you have more then 300 samples or even more.
3. Mafft v7.0 for multiple sequence alignment (MSA).
4. FastTreeMP (http://www.microbesonline.org/fasttree/#OpenMP) for building maximum likehood tree and bootsrtraping.

falco, vsearch and mafft can be installed via conda install. FastTreeMP can be downloaded and compiled from source.




Useage: chmod a+x ./wsearch.bash

nohup time ./wsearch.bash -d ./ -o ./output -T 16 -p primer1.fa -S OTU -t NBC -e T &

OTU clustering will be performed using UPARSE algorithm implemented in usearch (Edgar 2016, Nat. Method) and taxonomy classfication will performed either using the Näive bayesian classifier (NBC) (Wang et.al. 2007, Appl. Env. Micro) or sintax algorithm inplemented in vsearch(Edgar 2016, bioRxive, https://doi.org/10.1101/074161).

Exact sequence variance (e.g. ASV) will be generated using the unoise2/3 algorithm (Edgar 2016,bioRxive, https://doi.org/10.1101/081257). We do not recommend this method because ASVs can artificially split bacterial genomes into separate clusters (Scholss 2021, bioRxive, https://doi.org/10.1101/2021.02.26.433139) and thus overestimation of diversity (https://aem.asm.org/content/79/19/5962) due to 16S copy numer while OTU clustering at 97% identity is less easily subjected to this issue.

Options:
-d directory contains raw forward and reverse reads, must end with _R1.fq and _R2.fq, .1.fastq and .2.fastq, .R1.fastq and .R2.fastq, .R1.fq and .R2.fq, .1.fq and .2.fq

-o output directory

-t taxonomy classification method, NBC or sintax, default NBC -S species definition, ASV or OTU or both, by default is OTU only. Both will take more time

-p primers used for ampflication, should be in fasta format with fasta header >forward and >reverse, respectively. 


`>forward`
`GTGARTCATCGARTCTTTG`
`>reverse`
`TCCTCCGCTTATTGATATGC`

-e phylogenetric tree building, T or F, defaulf F

-T threads used for vsearch and usearch, default all the threads available

-b database path for taxonomy classification, default none and download database from usearch website for each classifier, use corresponding database for NBC or Sintax

-i identity for OTU clustering, default 0.97
