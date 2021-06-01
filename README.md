# usearch_wrapper: wsearch.bash
usearch wrapper for amplicon analysis. Support for 16S and ITS. This is still under active development and will be ready in a few months. Please contact jianshu(jianshu.zhao@gatech.edu)
This shell script has the following dependencies:

1. Falco v0.2.4 for reads quality, a C++ implentation of FastQC, way faster and user friendly than FastQC and more importantly, it is under active development (https://github.com/smithlabcode/falco)
2. vsearch v2.17 and Usearch v11 of core sequence analysis (https://drive5.com/usearch/download.html),both 32 bit and 64 bit usearch works but we strongly suggest 64 bit for large dataset. For example, when you have more then 300 samples or even more. Note that usearch 32-bit is not supported on MacOS Catillina or later but only Mojave or before. You will not be able to run usearch on MacOS Catillina or latter. But it should works on all major linux distributions. Let me know when you have large dataset and want to use an 64-bit. I can help with that
3. mothur v1.44.3 or latter, tested on ubuntu 18.0.4 and MacOS Mojave or latter.
4. Clustal Omega v1.2.4 for multiple sequence alignment (MSA).
5. FastTreeMP (http://www.microbesonline.org/fasttree/#OpenMP) for building maximum likehood tree and bootsrtraping.

falco, vsearch and clustal omega and fasttreeMP can be installed via conda install. FastTreeMP can be downloaded and compiled from source (those four dependencies are provided in the latest version and directly called without the need to install).

Now this script supported Linux and MacOS, tested on Ubuntu 18.0.4, RHEL 7 and MacOS Mojave (usearch 32 bit not suported on later versions). For MacOS Mojave, you need to install GNU awk (brew install gawk core-untils) after installing brew (https://brew.sh)



```

### install git first, on MacOS just run 'brew install git' after installing brew here: https://brew.sh
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

brew install gawk
brew install git
brew install libgit2

git clone https://github.com/jianshu93/usearch_wrapper

### On linux please run
sudo apt-get install git
sudo apt-get install libgit2
### You can also install it via conda
conda install git

## you can only run script in this directory because there are hard-coded dependencies in this directory
cd usearch_wrapper
chmod a+x ./wsearch.bash
### primer not removed
nohup time bash ./wsearch.bash -d ./demo_input -o ./output -T 16 -p primer1.fa -S OTU -t SINTAX -e T &
### you want to see each step
bash ./wsearch.bash -d ./demo_input -o ./output -T 16 -p primer1.fa -S ASV -t SINTAX -e T

### if primers have been removed
bash ./wsearch.bash -d ./demo_input -o ./output -T 16 -S ASV -t SINTAX -e T

### use mothur NBC instead of default usearch -nbc (memory limits and only 1 threads is used)
bash ./wsearch.bash -d ./demo_input -o ./output -T 16 -S ASV -t mothur -e T

### if you have you own usearch installed (64 bit for example on MacOS and Linux)
bash ./wsearch.bash -d ./demo_input -o ./output -T 16 -p primer1.fa -S ASV -t SINTAX -e T -u /usr/local/bin/usearch



```


OTU clustering will be performed using UPARSE algorithm implemented in usearch (Edgar 2016, Nat. Method,https://www.nature.com/articles/nmeth.2604) and taxonomy classfication will performed either using the Näive bayesian classifier (NBC) (Wang et.al. 2007, Appl. Env. Micro,https://aem.asm.org/content/73/16/5261.short) or sintax algorithm inplemented in vsearch(Edgar 2016, bioRxive, https://doi.org/10.1101/074161).

Exact sequence variance (e.g. ASV) will be generated using the unoise2/3 algorithm (Edgar 2016,bioRxive, https://doi.org/10.1101/081257). We do not recommend this method because ASVs can artificially split bacterial genomes into separate clusters (Scholss 2021, bioRxive, https://doi.org/10.1101/2021.02.26.433139) and thus overestimation of species diversity (https://aem.asm.org/content/79/19/5962) due to 16S copy numer while OTU clustering at 97% identity is less easily subjected to this issue. Or you can cluster ASV sequences into OTUs after the unoise3 step before creating OTU tables using -usearch_global. This option is not available now but will certainly be available in the near future.

There are a number of ways to correct for 16S copy numbers via searching in the microbial genomic database. However, all those methods perform poorly because of large variation in 16S copies even in close related genomes(More than 99.9% genome similarity) (https://link.springer.com/article/10.1186/s40168-018-0420-9). What we suggest is to use presence-absence OTU table at the same time though this will not take into account the abundance of the same genome species in each sample. I am very interested in developing tools to correct the 16S copy bias considering the fact that we have a much larger complete microbial genome database.

Options:
-d directory contains raw forward and reverse reads, must end with _R1.fq and _R2.fq, .1.fastq and .2.fastq, .R1.fastq and .R2.fastq, .R1.fq and .R2.fq, .1.fq and .2.fq

-o output directory

-t taxonomy classification method, NBC or sintax, default NBC -S species definition, ASV or OTU or both, by default is OTU only. Both will take more time

-p primers used for ampflication, should be in fasta format with fasta header >forward and >reverse, respectively. default is none, without removing primers but do quality control directly

```
>forward 
GTGARTCATCGARTCTTTG
>reverse
TCCTCCGCTTATTGATATGC`
```

-e phylogenetric tree building, T or F, defaulf F

-T threads used for vsearch and usearch, default all the threads available

-b database path for taxonomy classification, default none and download database from usearch website for each classifier, use corresponding database for NBC or Sintax

-i identity for OTU clustering, default 0.97, and I strongly recommend the default one

-u usearch binary path, /usr/local/bin/usearch, for example, make it executable first, by default is the ../dependencies/usearch_linux , a 32 bit version


# More related to taxonomy classification

There are a number of improvments/new implementations of RDP algorithm such as IDTAXA (https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0521-5) and microclass package (https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1583-2) but unfortunately they are all R packages (implemented in Rcpp) and there are no independent execuables available. There are some improvments like using both multinomial distribution models of kmer for training and classification.

# Reference
Edgar, R. C., & Flyvbjerg, H. (2015). Error filtering, pair assembly and error correction for next-generation sequencing reads. Bioinformatics, 1–7. http://doi.org/10.1093/bioinformatics/btv401

Rognes, T., Flouri, T., Ben Nichols, Quince, C., & Mahé, F. (2016). VSEARCH: a versatile open source tool for metagenomics. PeerJ, 4(17), e2584–22. http://doi.org/10.7717/peerj.2584

Sievers, Fabian et al. 2011. “Fast, Scalable Generation of High-Quality Protein Multiple Sequence Alignments Using Clustal Omega.” Molecular Systems Biology 7(1):1–6. 

Edgar, Robert C. 2010. “Search and Clustering Orders of Magnitude Faster Than BLAST.” Bioinformatics 26(19):2460–61.

Edgar, Robert C. 2016. “UNOISE2: Improved Error-Correction for Illumina 16S and ITS Amplicon Sequencing.” bioRxiv 1–21.

Price, M. N., Deha, P. S., & Arkin, A. P. (2010). FastTree 2 – Approximately Maximum-Likelihood Trees for Large Alignments, 5(3), e9490. http://doi.org/10.1371/journal.pone.0009490

Edgar, Robert C. 2016. “SINTAX: a Simple Non-Bayesian Taxonomy Classifier for 16S and ITS Sequences.” bioRxiv 1–20.

Wang, Qiong, George M. Garrity, James M. Tiedje, and James R. Cole. 2007. “Naive Bayesian Classifier for Rapid Assignment of rRNA Sequences Into the New Bacterial Taxonomy.” Applied and environmental microbiology 73(16):5261–67.
