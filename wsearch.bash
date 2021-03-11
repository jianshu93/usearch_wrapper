#!/bin/bash

threads=$(nproc)
dir=./
output=./output
spe_def="ASV"
identity=0.97
tax="NBC"
tre="F"
primer=primer.fa
echo $threads
while getopts ":d:o:t:S:e:T:p:b:i:h" option
do
	case $option in
		d) dir=$OPTARG;;
        o) output=$OPTARG;;
		t) tax=$OPTARG;;
        S) spe_def=$OPTARG;;
		e) tre=$OPTARG;;
		T) threads=$OPTARG;;
		p) primer=$OPTARG;;
		b) db=$OPTARG;;
        i) identity=$OPTARG;;
		\?) echo "Invalid option: -$OPTARG" >&2
			exit 1
			;;
		:) echo "Option -$OPTARG requires an argument." >&2
		    exit 1
			;;
		h) 
           echo "Amplicon sequencing analysis pipeline by Jianshu Zhao, Georgia Institute of Technology, 
           		jianshu.zhao@gatech.edu
                Lisense MIT

           		usage: usearch.bash -d ~/dir -S OTU -tax nbc -p primer.fa -tre T -T 8 -db ~/rdp_v18.fa -o ./results -id 0.97

				OTU clustering will be performed using UPARSE algorithm implemented in usearch 
				(Edgar 2016, Nat. Method) and taxonomy classfication will performed either 
				using the Näive bayesian classifier (NBC) (Wang et.al. 2007, Appl. Env. Micro) or sintax 
				algorithm inplemented in vsearch(Edgar 2016, bioRxive, https://doi.org/10.1101/074161).

                Exact sequence variance (e.g. ASV) will be generated using the unoise2/3 algorithm 
                (Edgar 2016,bioRxive, https://doi.org/10.1101/081257). We do not recommend this method 
                because ASVs can artificially split bacterial genomes into separate clusters 
                (Scholss 2021, bioRxive, https://doi.org/10.1101/2021.02.26.433139) while OTU clustering 
                at 97% identity is less easily subjected to this issue.

           		Options:
           		-d directory contains raw forward and reverse reads, must end with _R1.fq and _R2.fq, 
				   .1.fastq and .2.fastq, .R1.fastq and .R2.fastq, .R1.fq and .R2.fq, .1.fq and .2.fq
                -o output directory
				-t taxonomy classification method, NBC or sintax, default NBC
                -S species definition, ASV or OTU or both, by default is OTU only. Both will take more time
           		-p primers used for ampflication, should be in fasta format with fasta header >forward and
				   >reverse, respectively
                -e phylogenetric tree building, T or F, defaulf F
           		-T threads used for vsearch and usearch, default all the threads available
				-b database path for taxonomy classification, default none and download database from 
					usearch website for each classifier, use corresponding database for NBC or Sintax
                -i identity for OTU clustering, default 0.97"
           exit 1
           ;;
	esac
done

if ! [ -x "$(command -v vsearch)" ]; then
  echo 'Error: vsearch is not installed. Please install it' >&2
  exit 1
fi

if ! [ -x "$(command -v usearch)" ]; then
  echo 'Error: usearch is not installed. Please install it' >&2
  exit 1
fi

### fastq file name transformation
for f in *.R1.fastq; do
    a="$(echo $f | sed s/.R1.fastq/_R1.fq/)"
    mv "$f" "$a"
done

for f in *.R2.fastq; do
    a="$(echo $f | sed s/.R2.fastq/_R2.fq/)"
    mv "$f" "$a"
done

for f in *.1.fastq; do
    a="$(echo $f | sed s/.1.fastq/_R1.fq/)"
    mv "$f" "$a"
done

for f in *.2.fastq; do
    a="$(echo $f | sed s/.2.fastq/_R2.fq/)"
    mv "$f" "$a"
done

for f in *.R1.fq; do
    a="$(echo $f | sed s/.R1.fq/_R1.fq/)"
    mv "$f" "$a"
done

for f in *.R2.fq; do
    a="$(echo $f | sed s/.R2.fq/_R2.fq/)"
    mv "$f" "$a"
done

for f in *.1.fq; do
    a="$(echo $f | sed s/.1.fq/_R1.fq/)"
    mv "$f" "$a"
done

for f in *.2.fq; do
    a="$(echo $f | sed s/.2.fq/_R2.fq/)"
    mv "$f" "$a"
done

if [ -d "$dir" ] 
then
    echo "directory exists" 
    $(cd $dir)
else
    echo "$dir does not exist, please give a right directory"
    exit 1
fi


if [ -d "$output" ] 
then
    echo "Directory $output exists. Please give a new directory"
    exit 1
else
    echo "making directory $output ..."
    $(mkdir $output)
fi


### reads merging
echo "I am merging pair-end amplicon reads using vsearch with $threads threads"

for F in *_R1.fq; do
	BASE=${F##*/}
	SAMPLE=${BASE%_*}
	$(vsearch --fastq_mergepairs $F --reverse ${SAMPLE}_R2.fq --fastqout $output/${SAMPLE}_merged.fq -relabel ${SAMPLE}. --threads $threads)
done
for F in $output/*_merged.fq; do
    $(cat $F >> $output/all_samples_merged.fq)
    $(rm $F)
done
echo "reads merging done"

echo "I am removing primers"
$(usearch -fastx_subsample $output/all_samples_merged.fq -sample_size 5000 -fastqout $output/all_sub_for_primer_check.fq -threads $threads)
$(usearch -search_oligodb $output/all_sub_for_primer_check.fq -db $primer -strand both -userout $output/primer_hits.txt -userfields query+qlo+qhi+qstrand -threads 1)
variable1=`expr $(awk '$4 == "+"' $output/primer_hits.txt | awk '{print $3}' | sort | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }')`
variable2=`expr $(awk '$4 == "-"' $output/primer_hits.txt | awk '{print $3 - $2}' | sort | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }') + 1`
echo $variable1
$(vsearch --fastq_filter $output/all_samples_merged.fq --fastq_stripleft $variable1 --fastq_stripright $variable2 -fastq_maxee 1 --fastq_qmax 42 --fastq_maxlen 290 --fastq_minlen 220 --fastaout $output/QCd_merged.fa)

echo "primer removing done"


echo "I am dereplicating sequences"
$(vsearch --derep_fulllength $output/QCd_merged.fa -sizeout -relabel Uniq -output $output/uniques_vsearch.fa --threads $threads)
echo "Sequence dereplication done"

if [[ "$spe_def" == "ASV" ]]; then
    echo "I am clustering ASVs"
    $(usearch -unoise3 $output/uniques_vsearch.fa -zotus $output/ASVs.fa -minsize 2 -threads $threads)
    $(usearch -usearch_global $output/QCd_merged.fa -db $output/ASVs.fa -id 0.99 -otutabout $output/ASV_counts.txt -threads $threads -strand both)
    echo "I am done generating ASVs"
    echo "I am doing taxonomy assignment of ASVs"
    if [[ "$tax" == "NBC" ]]; then
	    if [[ "$db" == "" ]]; then
		    $(wget https://www.drive5.com/sintax/rdp_16s_v18.fa.gz)
		    $(gunzip rdp_16s_v18.fa.gz)
            $(usearch -nbc_tax $output/ASVs.fa --db rdp_16s_v18.fa -strand plus --threads $threads -tabbedout $output/asv_tax_rdp.txt)
            $(rm rdp_16s_v18.fa)
	    else
            $(usearch -nbc_tax $output/ASVs.fa --db $db -strand plus --threads $threads -tabbedout $output/asv_tax_rdp.txt)
	    fi
        echo "taxonomy asignment of ASVs using Näive bayesian classifier done"
    else
        if [[ "$db" == "" ]]; then
            $(wget https://www.drive5.com/sintax/silva_16s_v123.fa.gz)
            $(gunzip silva_16s_v123.fa.gz)
            $(vsearch --sintax $output/ASVs.fa --db silva_16s_v123.fa --tabbedout $output/asv_tax_sintax.txt --threads $threads --sintax_cutoff 0.8)
            $(rm silva_16s_v123.fa)
        else
            $(vsearch --sintax $output/ASVs.fa --db $db --tabbedout $output/asv_tax_sintax.txt --threads $threads --sintax_cutoff 0.8)
        fi
        echo "taxonomy asignment of ASVs using sintax done"
    fi
else
    if [[ "$spe_def" == "both" ]]; then
        echo "I am clustering OTUs and generate ASVs"
        $(usearch -cluster_otus $output/uniques_vsearch.fa -otus $output/otus.fa -relabel Otu -threads $threads)
        $(usearch -usearch_global $output/QCd_merged.fa -db $output/otus.fa -id 0.97 -otutabout $output/otu_counts.txt -threads $threads -strand both)
        $(usearch -unoise3 $output/uniques_vsearch.fa -zotus ASVs.fa -minsize 2 -threads $threads)
        $(usearch -usearch_global $output/QCd_merged.fa -db $output/ASVs.fa -id 0.99 -otutabout $output/ASV_counts.txt -threads $threads -strand both)

        echo "I am done clustering OTUs and generating ASVs"
        echo "I am doing taxonomy assignment of OTUs and ASVs"
        if [[ "$tax" == "NBC" ]]; then
            echo "I am done clustering OTUs and generating ASVs"
            echo "I am doing taxonomy assignment of OTUs and also ASVs"
	        if [[ "$db" == "" ]]; then
		        $(wget https://www.drive5.com/sintax/rdp_16s_v18.fa.gz)
		        $(gunzip rdp_16s_v18.fa.gz)
                $(usearch -nbc_tax $output/ASVs.fa --db rdp_16s_v18.fa -strand plus --threads $threads -tabbedout $output/asv_tax_rdp.txt)
                $(usearch -nbc_tax $output/otus.fa --db rdp_16s_v18.fa -strand plus --threads $threads -tabbedout $output/otu_tax_rdp.txt)
                $(rm rdp_16s_v18.fa)
	        else
                $(usearch -nbc_tax $output/ASVs.fa --db $db -strand plus --threads $threads -tabbedout $output/asv_tax_rdp.txt)
                $(usearch -nbc_tax $output/otus.fa --db $db -strand plus --threads $threads -tabbedout $output/otu_tax_rdp.txt)
	        fi
            echo "taxonomy asignment of OTUs and ASVs using Näive bayesian classifier done"
        else
            if [[ "$db" == "" ]]; then
                $(wget https://www.drive5.com/sintax/silva_16s_v123.fa.gz)
                $(gunzip silva_16s_v123.fa.gz)
                $(vsearch --sintax $output/ASVs.fa --db silva_16s_v123.fa --tabbedout $output/asv_tax_sintax.txt --threads $threads --sintax_cutoff 0.8)
                $(vsearch --sintax $output/otus.fa --db silva_16s_v123.fa --tabbedout $output/otu_tax_sintax.txt --threads $threads --sintax_cutoff 0.8)
                $(rm silva_16s_v123.fa)
            else
                $(vsearch --sintax $output/ASVs.fa --db $db --tabbedout $output/asv_tax_sintax.txt --threads $threads --sintax_cutoff 0.8)
                $(vsearch --sintax $output/otus.fa --db $db --tabbedout $output/otu_tax_sintax.txt --threads $threads --sintax_cutoff 0.8)
            fi
            echo "taxonomy asignment of OTUs and ASVs using sintax done"
        fi
    else
        echo "I am clustering OTUs"
        $(usearch -cluster_otus $output/uniques_vsearch.fa -otus $output/otus.fa -relabel Otu -threads $threads)
        $(usearch -usearch_global $output/QCd_merged.fa -db $output/otus.fa -id 0.97 -otutabout $output/otu_counts.txt -threads $threads -strand both)
        echo "OTUs clustering done"
        if [[ "$tax" == "NBC" ]]; then
	        if [[ "$db" == "" ]]; then
		        $(wget https://www.drive5.com/sintax/rdp_16s_v18.fa.gz)
		        $(gunzip rdp_16s_v18.fa.gz)
                $(usearch -nbc_tax $output/otus.fa --db rdp_16s_v18.fa -strand plus --threads $threads -tabbedout $output/otu_tax_rdp.txt)
                $(rm rdp_16s_v18.fa)
	        else
                $(usearch -nbc_tax $output/otus.fa --db $db -strand plus --threads $threads -tabbedout $output/otu_tax_rdp.txt)
	        fi
            echo "taxonomy asignment of OTUs using Näive bayesian classifier done"
        else
            if [ "$db" == "" ]; then
                $(wget https://www.drive5.com/sintax/silva_16s_v123.fa.gz)
                $(gunzip silva_16s_v123.fa.gz)
                $(vsearch --sintax $output/otus.fa --db silva_16s_v123.fa --tabbedout $output/otu_tax_sintax.txt --threads $threads --sintax_cutoff 0.8)
                $(rm silva_16s_v123.fa)
            else
                $(vsearch --sintax $output/ASVs.fa --db $db --tabbedout $output/asv_tax_sintax.txt --threads $threads --sintax_cutoff 0.8)
            fi
            echo "taxonomy asignment of OTUs using sintax done"
        fi
    fi
fi
### multiple sequence alignment and phylogenetic tree building
if [[ "$tre" == "T" ]] ; then
    if ! [ -x "$(command -v mafft)" ]; then
        echo 'Error: mafft is not installed. Please install it' >&2
        exit 1
    else
        echo "I am doing multiple sequence alignment. This may take long, please wait..."
        if [[ "$spe_def" == "ASV" ]] ; then
            $(mafft --thread $threads --reorder --auto $output/ASVs.fa > $output/ASVs_aligned.fa)
            echo "multiple sequence alignment for ASVs done"
        else
            if [ "$spe_def" == "both" ] ; then
                $(mafft --thread $threads --reorder --auto $output/ASVs.fa > $output/ASVs_aligned.fa)
                $(mafft --thread $threads --reorder --auto $output/otus.fa > $output/otus_aligned.fa)
                echo "sequence alignment for ASVs and OTUs done"
            else
                $(mafft --thread $threads --reorder --auto $output/otus.fa > $output/otus_aligned.fa)
                echo "multiple sequence alignment for OTUs done"
            fi
        fi
    fi
    if ! [ -x "$(command -v FastTreeMP)" ]; then
        echo 'Error: FastTreeMP is not installed. Please install it' >&2
        exit 1
    else
        echo "I am building phylogenetic tree using FastTreeMP, please wait..."
        if [[ "$spe_def" == "ASV" ]] ; then
            $(FastTreeMP -out $output/ASVs_aligned.tre -nt -gtr $output/ASVs_aligned.fa)
            echo "phylogenetic tree building for ASVs done"
        else
            if [[ "$spe_def" == "both" ]] ; then
                $(FastTreeMP -out $output/ASVs_aligned.tre -nt -gtr $output/ASVs_aligned.fa)
                $(FastTreeMP -out $output/otus_aligned.tre -nt -gtr $output/otus_aligned.fa)
                echo "phylogenetic tree building for ASVs and OTUs done"
            else
                $(FastTreeMP -out $output/otus_aligned.tre -nt -gtr $output/otus_aligned.fa)
                echo "phylogenetic tree building for OTUs done"
            fi
        fi
    fi
fi

$(rm rdp_16s_v18.fa)
echo "Amplicon sequence analysis done, output files are in $output"