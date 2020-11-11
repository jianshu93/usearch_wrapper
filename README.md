# usearch_wrapper
usearch wrapper for amplicon analysis. Support for 16S and ITS. This is still under active development and will be ready in a few months. Please contact jianshu(jianshu.zhao@gatech.edu) and Qi Qi for details (qiq17@mails.tsinghua.edu.cn)
This shell script has the following dependencies:

1. Falco v0.2.2, a C++ implentation of FastQC, way faster and user friendly than FastQC and more importantly, it is under active development (https://github.com/falcosecurity/falco)
2. Usearch v11,both 32 bit and 64 bit works but we strongly suggest 64 bit for large dataset. for example when you have more then 100 samples or even more.
3. Mafft v7.0 for multiple sequence alignment (MSA).
4. FasttreeMP for building maximum likehood tree and bootsrtraping.
