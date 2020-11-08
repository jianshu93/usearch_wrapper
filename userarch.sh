
#### This is the usearch11 wrapper (v1.0) for fast analysis of amplicon sequencing data. 
#### jianshuzhao@yahoo.com, School of Environment, Tsinghua University and Center for Bioformatics and Computational Biology, Georgia Institute of technology

while getopts ":(in):(tmp):(out):T:(clr):(tax):(tre):h" option
do
	case $option in
		in) input_dir=$OPTARG;;
		tmp) tmp_dir=$OPTARG;;
		out) out_dir=$OPTARG;;
		prm) primer==$OPTARG;;
		T) num_threads=$OPTARG;;
		clr) cluster_method=$OPTARG;;
		tax) taxonomy_assign=$OPTARG;;
		tre) tree=$OPTARG;;
		\?) echo "Invalid option: -$OPTARG" >&2
			exit 1
			;;
		:) echo "Option -$OPTARG requires an argument." >&2
			exit 1
			;;
		h) 
           echo "This is usearch wrapper (v1.0) for fast analysis of amplicon sequencing data, the current
           		 version supports both 16S gene and ITS gene (presuamablly V3-V4 and ITS2)

           		 Usage: usearch.sh -in ~/usearh_dir -tmp ~/tmp -out ~/usearch_out -clr OTU -tax NBC -tre True

           		 Options:
           		 -in input directory for demultiplex reads, with *R1.fastq and *R2.fastq in the director (*must
           		  be the same), you can use usearch -demultiplex command to demultiplex you raw reads according 
           		  to barcodes if you do not have them
           		 -tmp temporary folder for storing intermediate files such as merged reads, clean reads et.al.
           		 -out output directory for final output including OTU table with taxonomy, phylogenetic tree
           		 -pmr primer for amplicon sequencing in the format of fasta, first for forward and second for reverse
           		 -T number of threads to use for usearch, by default is 8
           		 -clr clustering method, by default is tradition OTU clusting method using uclust, you can also
           		  use exact sequence variance by using denoise3 algorithm, defalut is uclust for OTU
           		 -tax taxonomy assignment method, Näive Bayesian Classifier algorithm offered by RDP (NBC) or 
           		 Sintax offered in usearch, by default is NBC
           		 -tre build phylogenetic tree or not, it take much longer when you have more than 10000 sequences,
           		 this script use mafft to align sequences and then build phylogenetic tree using fasttreeMP,
           		 evolutionary model is the defaulst GTR"

		