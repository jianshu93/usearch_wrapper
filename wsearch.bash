#!/bin/bash
### Jianshu Zhao (jianshu.zhao@gatech.edu)
### Amplicon analysis pipeline based on open source sofware vsearch and usearch.

threads=$(nproc)
dir=./
output=./output
spe_def="ASV"
identity=0.97
tax="NBC"
tre="F"
usearch_bin=./dependencies/usearch11.0.667_i86linux32
echo $threads
while getopts ":d:o:t:S:e:T:p:b:i:u:h" option
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
        u) usearch_bin=$OPTARG;;
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
                -i identity for OTU clustering, default 0.97
                -u usearch binary path, ~/usearch, for example, make it executable first, by default is the 
                    ../dependencies/usearch_linux , a 32 bit version"
           exit 1
           ;;
	esac
done

echo "$primer will be used for primer removing"
echo $usearch_bin

if [ -d "$dir" ] 
then
    echo "" 
else
    echo "$dir does not exist, please offer a directory that exists"
    exit 1
fi

if [ -d "$output" ] 
then
    echo "Directory $output already exists. Please offer a new directory"
    exit 1
else
    echo "making directory $output ..."
    $(mkdir $output)
fi


### reads merging
echo "I am merging pair-end amplicon reads using vsearch with $threads threads"
dfiles="${dir}/*_R1.fastq"

for F in $dfiles; do
	BASE=${F##*/}
	SAMPLE=${BASE%_*}
	$(./dependencies/vsearch_linux --fastq_mergepairs $F --reverse ${dir}/${SAMPLE}_R2.fastq --fastqout $output/${SAMPLE}_merged.fq -relabel ${SAMPLE}. --threads $threads)
    $(./dependencies/falco_linux -o $output/${SAMPLE}_falco_output $F ${dir}/${SAMPLE}_R2.fastq)
done

ofiles="${output}/*_merged.fq"

for F in $ofiles; do
    $(cat $F >> $output/all_samples_merged.fq)
    $(rm $F)
done
echo "reads merging done"

if [[ ! -z "$primer" ]]; then
    echo "I am removing primers"
    $($usearch_bin -fastx_subsample $output/all_samples_merged.fq -sample_size 5000 -fastqout $output/all_sub_for_primer_check.fq -threads $threads)
    $($usearch_bin -search_oligodb $output/all_sub_for_primer_check.fq -db $primer -strand both -userout $output/primer_hits.txt -userfields query+qlo+qhi+qstrand -threads 1)
    variable1=`expr $(awk '$4 == "+"' $output/primer_hits.txt | awk '{print $3}' | sort | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }')`
    variable2=`expr $(awk '$4 == "-"' $output/primer_hits.txt | awk '{print $3 - $2}' | sort | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }') + 1`
    echo "$variable1 bps will be removed at the beginning of the merged reads"
    echo "$variable2 bps will be removed at the end of the merged reads"
    $(./dependencies/vsearch_linux --fastq_filter $output/all_samples_merged.fq --fastq_stripleft $variable1 --fastq_stripright $variable2 -fastq_maxee 1 --fastq_qmax 42 --fastq_maxlen 290 --fastq_minlen 220 --fastaout $output/QCd_merged.fa)
    echo "primer removing done"
else
    $(./dependencies/vsearch_linux --fastq_filter $output/all_samples_merged.fq -fastq_maxee 1 --fastq_qmax 42 --fastq_maxlen 290 --fastq_minlen 220 --fastaout $output/QCd_merged.fa)
    echo "No primers to remove, reads are filtered"
fi

echo "I am dereplicating sequences"
$(./dependencies/vsearch_linux --derep_fulllength $output/QCd_merged.fa -sizeout -relabel Uniq -output $output/uniques_vsearch.fa --threads $threads)
echo "Sequence dereplication done"

if [[ "$spe_def" == "ASV" ]]; then
    echo "I am clustering ASVs"
    $($usearch_bin -unoise3 $output/uniques_vsearch.fa -zotus $output/ASVs.fa -minsize 2 -threads $threads)
    $(./dependencies/vsearch_linux --usearch_global $output/QCd_merged.fa --db $output/ASVs.fa --id 0.99 --otutabout $output/ASV_counts.txt --threads $threads --strand both)
    echo "I am done generating ASVs"
    echo "I am doing taxonomy assignment of ASVs"
    if [[ "$tax" == "NBC" ]]; then
	    if [[ ! -z "$db" ]]; then
		    $(wget https://www.drive5.com/sintax/rdp_16s_v18.fa.gz)
		    $(gunzip rdp_16s_v18.fa.gz)
            $(usearch_bin -nbc_tax $output/ASVs.fa --db rdp_16s_v18.fa -strand plus --threads $threads -tabbedout $output/asv_tax_rdp.txt)
            $(rm rdp_16s_v18.fa)
	    else
            $(usearch_bin -nbc_tax $output/ASVs.fa --db $db -strand plus --threads $threads -tabbedout $output/asv_tax_rdp.txt)
	    fi
        $(echo -e "#OTU ID\ttaxonomy" > $output/asv_tax_rdp_0.5.txt)
        $($awk 'BEGIN {FS="\t"}; {print $1,$4}' OFS='\t' $output/asv_tax_rdp.txt >> $output/asv_tax_rdp_0.5.txt)
        $($awk 'BEGIN {FS="\t"}; FNR==NR { a[$1]=$0; next } $1 in a { print a[$1], $2}' OFS='\t' $output/ASV_counts.txt $output/asv_tax_rdp_0.5.txt > $output/asv_table_rdp.txt)
        echo "taxonomy asignment of ASVs using Näive bayesian classifier done"
    else
        if [[ "$tax" == "SINTAX" ]]; then
            if [[ ! -z "$db" ]]; then
                $(wget https://www.drive5.com/sintax/silva_16s_v123.fa.gz)
                $(gunzip silva_16s_v123.fa.gz)
                $($usearch_bin --sintax $output/ASVs.fa --db silva_16s_v123.fa --tabbedout $output/asv_tax_sintax.txt --threads $threads --sintax_cutoff 0.8 -strand plus)
                $(rm silva_16s_v123.fa)
            else
                $($usearch_bin --sintax $output/ASVs.fa --db $db --tabbedout $output/asv_tax_sintax.txt --threads $threads --sintax_cutoff 0.8 -strand plus)
            fi
            $(echo -e "#OTU ID\ttaxonomy" > $output/asv_tax_sintax_0.8.txt)
            $($awk 'BEGIN {FS="\t"}; {print $1,$4}' OFS='\t' $output/asv_tax_sintax.txt > $output/asv_tax_sintax_0.8.txt)
            $($awk 'BEGIN {FS="\t"}; FNR==NR { a[$1]=$0; next } $1 in a { print a[$1], $2}' OFS='\t' $output/ASV_counts.txt $output/asv_tax_sintax_0.8.txt > $output/asv_table_sintax.txt)
        else
            echo "Not supported"
        fi
        echo "taxonomy asignment of ASVs using sintax done"
    fi
else
    if [[ "$spe_def" == "both" ]]; then
        echo "I am clustering OTUs and generate ASVs"
        $($usearch_bin -cluster_otus $output/uniques_vsearch.fa -otus $output/otus.fa -relabel Otu -threads $threads)
        $(./dependencies/vsearch_linux --usearch_global $output/QCd_merged.fa --db $output/otus.fa --id 0.97 --otutabout $output/otu_counts.txt --threads $threads --strand both)
        $($usearch_bin -unoise3 $output/uniques_vsearch.fa -zotus ASVs.fa -minsize 2 -threads $threads)
        $(./dependencies/vsearch_linux --usearch_global $output/QCd_merged.fa --db $output/ASVs.fa --id 0.99 --otutabout $output/ASV_counts.txt --threads $threads --strand both)

        echo "I am done clustering OTUs and generating ASVs"
        echo "I am doing taxonomy assignment of OTUs and ASVs"
        if [[ "$tax" == "NBC" ]]; then
	        if [[ ! -z "$db" ]]; then
		        $(wget https://www.drive5.com/sintax/rdp_16s_v18.fa.gz)
		        $(gunzip rdp_16s_v18.fa.gz)
                $($usearch_bin -nbc_tax $output/ASVs.fa --db rdp_16s_v18.fa -strand plus --threads $threads -tabbedout $output/asv_tax_rdp.txt)
                $($usearch_bin -nbc_tax $output/otus.fa --db rdp_16s_v18.fa -strand plus --threads $threads -tabbedout $output/otu_tax_rdp.txt)
                $(rm rdp_16s_v18.fa)
	        else
                $($usearch_bin -nbc_tax $output/ASVs.fa --db $db -strand plus --threads $threads -tabbedout $output/asv_tax_rdp.txt)
                $($usearch_bin -nbc_tax $output/otus.fa --db $db -strand plus --threads $threads -tabbedout $output/otu_tax_rdp.txt)
	        fi
            $(echo -e "#OTU ID\ttaxonomy" > $output/asv_tax_rdp_0.5.txt)
            $(echo -e "#OTU ID\ttaxonomy" > $output/otu_tax_rdp_0.5.txt)
            $($awk 'BEGIN {FS="\t"}; {print $1,$4}' OFS='\t' $output/otu_tax_rdp.txt >> $output/otu_tax_rdp_0.5.txt)
            $($awk 'BEGIN {FS="\t"}; {print $1,$4}' OFS='\t' $output/asv_tax_rdp.txt >> $output/asv_tax_rdp_0.5.txt)
            $($awk 'BEGIN {FS="\t"}; FNR==NR { a[$1]=$0; next } $1 in a { print a[$1], $2}' OFS='\t' $output/ASV_counts.txt $output/asv_tax_rdp_0.5.txt > $output/asv_table_rdp.txt)
            $($awk 'BEGIN {FS="\t"}; FNR==NR { a[$1]=$0; next } $1 in a { print a[$1], $2}' OFS='\t' $output/otu_counts.txt $output/otu_tax_rdp_0.5.txt > $output/otu_table_rdp.txt)
            echo "taxonomy asignment of OTUs and ASVs using Näive bayesian classifier done"
        else
            if [[ "$tax" == "SINTAX" ]]; then
                if [[ ! -z "$db" ]]; then
                    $(wget https://www.drive5.com/sintax/silva_16s_v123.fa.gz)
                    $(gunzip silva_16s_v123.fa.gz)
                    $($usearch_bin --sintax $output/ASVs.fa --db silva_16s_v123.fa --tabbedout $output/asv_tax_sintax.txt --threads $threads --sintax_cutoff 0.8 -strand plus)
                    $($usearch_bin --sintax $output/otus.fa --db silva_16s_v123.fa --tabbedout $output/otu_tax_sintax.txt --threads $threads --sintax_cutoff 0.8 -strand plus)
                    $(rm silva_16s_v123.fa)
                else
                    $($usearch_bin --sintax $output/ASVs.fa --db $db --tabbedout $output/asv_tax_sintax.txt --threads $threads --sintax_cutoff 0.8 -strand plus)
                    $($usearch_bin --sintax $output/otus.fa --db $db --tabbedout $output/otu_tax_sintax.txt --threads $threads --sintax_cutoff 0.8 -strand plus)
                fi
                $(echo -e "#OTU ID\ttaxonomy" > $output/asv_tax_sintax_0.8.txt)
                $(echo -e "#OTU ID\ttaxonomy" > $output/otu_tax_sintax_0.8.txt)
                $($awk 'BEGIN {FS="\t"}; {print $1,$4}' OFS='\t' $output/asv_tax_sintax.txt >> $output/asv_tax_sintax_0.8.txt)
                $($awk 'BEGIN {FS="\t"}; {print $1,$4}' OFS='\t' $output/otu_tax_sintax.txt >> $output/otu_tax_sintax_0.8.txt)
                $($awk 'BEGIN {FS="\t"}; FNR==NR { a[$1]=$0; next } $1 in a { print a[$1], $2}' OFS='\t' $output/ASV_counts.txt $output/asv_tax_sintax_0.8.txt > $output/asv_table_sintax.txt)
                $($awk 'BEGIN {FS="\t"}; FNR==NR { a[$1]=$0; next } $1 in a { print a[$1], $2}' OFS='\t' $output/otu_counts.txt $output/otu_tax_sintax_0.8.txt > $output/otu_table_sintax.txt)
            else
                echo "Taxonomy assignment not supported"
            fi
            echo "Taxonomy asignment of OTUs and ASVs using sintax done"
        fi
    else
        if [[ "$spe_def" == "OTU" ]]; then
            echo "I am clustering OTUs"
            $($usearch_bin -cluster_otus $output/uniques_vsearch.fa -otus $output/otus.fa -relabel Otu -threads $threads)
            $(./dependencies/vsearch_linux --usearch_global $output/QCd_merged.fa --db $output/otus.fa --id 0.97 --otutabout $output/otu_counts.txt --threads $threads --strand both)
            echo "OTUs clustering done"
            echo "I am doing taxonomy assignment of OTUs"
            if [[ "$tax" == "NBC" ]]; then
	            if [[ ! -z "$db" ]]; then
		            $(wget https://www.drive5.com/sintax/rdp_16s_v18.fa.gz)
		            $(gunzip rdp_16s_v18.fa.gz)
                    $($usearch_bin -nbc_tax $output/otus.fa --db rdp_16s_v18.fa -strand plus --threads $threads -tabbedout $output/otu_tax_rdp.txt)
                    $(rm rdp_16s_v18.fa)
	            else
                    $($usearch_bin -nbc_tax $output/otus.fa --db $db -strand plus --threads $threads -tabbedout $output/otu_tax_rdp.txt)
	            fi
                $(echo -e "#OTU ID\ttaxonomy" > $output/otu_tax_rdp_0.5.txt)
                $($awk 'BEGIN {FS="\t"}; {print $1,$4}' OFS='\t' $output/otu_tax_rdp.txt >> $output/otu_tax_rdp_0.5.txt)
                $($awk 'BEGIN {FS="\t"}; FNR==NR { a[$1]=$0; next } $1 in a { print a[$1], $2}' OFS='\t' $output/otu_counts.txt $output/otu_tax_rdp_0.5.txt > $output/otu_table_rdp.txt)
                echo "taxonomy asignment of OTUs using Näive bayesian classifier done"
            else
                if [[ "$tax" == "SINTAX" ]]; then
                    if [[ ! -z "$db" ]]; then
                        $(wget https://www.drive5.com/sintax/silva_16s_v123.fa.gz)
                        $(gunzip silva_16s_v123.fa.gz)
                        $($usearch_bin --sintax $output/otus.fa --db silva_16s_v123.fa --tabbedout $output/otu_tax_sintax.txt --threads $threads --sintax_cutoff 0.8 -strand plus)
                        $(rm silva_16s_v123.fa)
                    else
                        $($usearch_bin --sintax $output/ASVs.fa --db $db --tabbedout $output/asv_tax_sintax.txt --threads $threads --sintax_cutoff 0.8 -strand plus)
                    fi
                    $(echo -e "#OTU ID\ttaxonomy" > $output/otu_tax_sintax_0.8.txt)
                    $($awk 'BEGIN {FS="\t"}; {print $1,$4}' OFS='\t' $output/otu_tax_sintax.txt >> $output/otu_tax_sintax_0.8.txt)
                    $($awk 'BEGIN {FS="\t"}; FNR==NR { a[$1]=$0; next } $1 in a { print a[$1], $2}' OFS='\t' $output/otu_counts.txt $output/otu_tax_sintax_0.8.txt > $output/otu_table_sintax.txt)
                    echo "taxonomy asignment of OTUs using SINTAX classifier done"
                else
                    echo "Taxonomy assignment method not supported"
                fi
            fi
        else
            echo "Only ASV and OTU are supported"
        fi
    fi
fi

### multiple sequence alignment and phylogenetic tree building
if [[ "$tre" == "T" ]] ; then
    echo "I am doing multiple sequence alignment. This may take long, please wait..."
    if [[ "$spe_def" == "ASV" ]] ; then
        $(./dependencies/clustalo-1.2.4-Ubuntu-x86_64 -i $output/ASVs.fa -o $output/ASVs_aligned.fa -t DNA --threads $threads)
        echo "multiple sequence alignment for ASVs done"
    else
        if [ "$spe_def" == "both" ] ; then
            $(./dependencies/clustalo-1.2.4-Ubuntu-x86_64 -i $output/ASVs.fa -o $output/ASVs_aligned.fa -t DNA --threads $threads)
            $(./dependencies/clustalo-1.2.4-Ubuntu-x86_64 -i $output/otus.fa -o $output/otus_aligned.fa -t DNA --threads $threads)
            echo "sequence alignment for ASVs and OTUs done"
        else
            $(./dependencies/clustalo-1.2.4-Ubuntu-x86_64 -i $output/otus.fa -o $output/otus_aligned.fa -t DNA --threads $threads)
            echo "multiple sequence alignment for OTUs done"
        fi
    fi
    echo "I am building phylogenetic tree using FastTreeMP, please wait..."
    if [[ "$spe_def" == "ASV" ]] ; then
        $(./dependencies/FastTreeMP_Linux -out $output/ASVs_aligned.tre -nt -gtr $output/ASVs_aligned.fa)
        echo "phylogenetic tree building for ASVs done"
    else
        if [[ "$spe_def" == "both" ]] ; then
            $(./dependencies/FastTreeMP_Linux -out $output/ASVs_aligned.tre -nt -gtr $output/ASVs_aligned.fa)
            $(./dependencies/FastTreeMP_Linux -out $output/otus_aligned.tre -nt -gtr $output/otus_aligned.fa)
            echo "phylogenetic tree building for ASVs and OTUs done"
        else
            $(./dependencies/FastTreeMP_Linux -out $output/otus_aligned.tre -nt -gtr $output/otus_aligned.fa)
            echo "phylogenetic tree building for OTUs done"
        fi
    fi
fi
echo "Amplicon sequence analysis done, output files are in $output"