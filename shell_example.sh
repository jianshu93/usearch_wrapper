#!/bin/bash

#### author: Jianshu Zhao, jianshu.zhao@gatech.edu

#### SNV calling pipeline for human genome, potentially works for genomes of any other organism.
#### required tools (dependencies): fqtools (https://github.com/alastair-droop/fqtools), samtools (>1.9.0), bwa,minimap2, GTAK v3.7.0, bcftools, picard (picard.jar)
#### It works on *nix system and also Darwin based kernel
#### note that this script requires java 1.8 and jdk 8.0 (oracle or openjdk version). You will also need the GATK .jar file (GenomeAnalysisTK.jar),picard (picard.jar)
#### Liscense: MIT 


### predifined variables for default values
picard_path="./picard.jar"


### input options and help information

while getopts ":a:b:r:M:T:m:d:o:G:e:f:z:i:v:h" option
do
	case $option in
		a) reads1=$OPTARG;;
		b) reads2=$OPTARG;;
		r) ref=$OPTARG;;
		M) memory=$OPTARG;;
		T) n_threads=$OPTARG;;
		m) mapping=$OPTARG;;
		d) temp_dir=$OPTARG;;
		o) out_filename=$OPTARG;;
		G) GATK_path=$OPTARG;;
		e) realign=$OPTARG;;
		f) millsFile=$OPTARG;;
		p) picard_path=$OPTARG;;
		z) gunzip1=$OPTARG;;
		i) index=$OPTARG;;
		v) VERBOSE=true;;
		\?) echo "Invalid option: -$OPTARG" >&2
			exit 1
			;;
		:) echo "Option -$OPTARG requires an argument." >&2
			exit 1
			;;
		h) 
           echo "Single Nucleitide Variant calling pipeline by Jianshu Zhao, Georgia Institute of Technology, 
           		 jianshu.zhao@gatech.edu


           		usage: SNV.sh -a Reads_R1.fastq.gz -b Reads_R2.fastq.gz -r Genome.fasta 
           		-T 4 -M 4g -d ./ -m minimap2 -G ./GenomeAnalysisTK.jar -e 1 
           		-f ./Mills_and_1000G_gold_standard.indels.hg38.vcf -p ./picard.jar -z 1 -i 1 -v 1

           		Options:

           		-a forward reads of a fastq file, must be in gzipped format
           		-b Reverse reads of a fastq file, must be in gzipped format
           		-r genome to map reads to, can be in fasta format
           		-M Memory to use for this pipeline, by default 1g
           		-T number of threads to run this pipeline, by default 1
           		-m reads mapping method,can be bwa or minimap2, by default is bwa
           		-d temporatory directory SNV works in, by default current dirctory
           		-o output file name, by default study.filtered.vcf.gz
           		-G GTAK java file (GenomeAnalysisTK.jar) path
           		-e do realign or not, 0 or 1
           		-f Mills_and_1000G_gold_standard.indels.hg38.vcf file path
           		-p picard jar file path, by default current directory
           		-z zip output file or not 0 or 1
           		-i whether index bam files or not, 0 or 1
           		-v verbose mode, true or false"
           		exit 1
           ;;
	esac
done

### File path checking
if [ "$VERBOSE" == "true" ]; then
    echo "I will start to check files and directories and 
    	create directory if not exist, checking ..."
fi

if test -f "$reads1"; then
    echo "$reads1 exists."
else
	echo "$reads1 does not exists."
	exit 1
fi

if test -f "$reads2"; then
    echo "$reads2 exists."
else
	echo "$reads2 does not exists."
	exit 1
fi

if test -f "$ref"; then
    echo "$ref exists."
else
	echo "$ref does not exists."
fi

if test -f "$GATK_path"; then
    echo "$GATK_path exists."
else
	echo "$ref does not exists."
	exit 1
fi

if [ -d "$temp_dir" ] 
then
    echo "Directory $temp_dir exists." 
else
    echo "Error: Directory $temp_dir does not exists, making $temp_dir ..."
    mkdir $temp_dir
fi

if [ "$VERBOSE" == "true" ]; then
    echo "I have finished checking input files and directories. 
    	I will checking dependent softwares, checking ..."
fi

### dependency checking
if ! [ -x "$(command -v fqtools)" ]; then
  echo 'Warming: fqtools is not installed. Fastq file format checking will not be performed' >&2
fi

if ! [ -x "$(command -v bwa)" ]; then
  echo 'Error: bwa is not installed. It is required for reads mapping' >&2
  exit 1
fi

if ! [ -x "$(command -v minimap2)" ]; then
  echo 'Suggestion: minimap2 is not installed. Its mapping speed is much faster than bwa, 
  		I strongly suggest you install it and use it' >&2
fi

if ! [ -x "$(command -v samtools)" ]; then
  echo 'Error: samtools is not installed. Please install it and then try again' >&2
  exit 1
fi

if ! [ -x "$(command -v bcftools)" ]; then
  echo 'Error: bcftools is not installed. Please install it and then try again' >&2
  exit 1
fi

if ! [ -x "$(command -v java)" ]; then
  echo 'Error: java is not installed. Please install it and then try again' >&2
  exit 1
fi

if [ "$VERBOSE" == "true" ]; then
    echo "I have finished checking software dependencies. 
    		I will validate file format, validating ..."
fi

###### file format checking
## reads fastq file format checking

fq1_test=$(fqtools validate $reads1)
fq2_test=$(fqtools validate $reads2)

if [ "$fq1_test" == "OK" && "$fq2_test" == "OK" ]; then
	echo 'fastq file format test passed.' >&2
else
	echo 'fqtools is not installed or the fastq file format validate fails, please check fastq file format'
fi

### reference fasta format checking



### reads mapping
cd $temp_dir

ref_base=$(basename $ref)

if [ "$VERBOSE" == "true" ]; then
    echo "I will start to map reads to reference, mapping ..."
fi
if [ $mapping == "bwa" ]; then
	$(bwa index $ref)
	$(bwa mem -R '@RG\tID:foo\tSM:bar\tLB:library1' -t $n_threads -o ./mapping2.sam $ref <(gunzip -c $reads1) <(gunzip -c $reads2) >&2)
else
	$(minimap2 -a -t $n_threads -x sr -o mapping2.sam $ref $reads1 $reads1 >&2 )
fi

if [ "$VERBOSE" == "true" ]; then
	echo "I have finished reads mapping. I will start transforming file formats"
fi
### fix mapping index error
$(samtools fixmate -@ $n_threads -O bam ./mapping2.sam ./mapping2.bam >&2)

### sort bam file
$(samtools sort -@ $n_threads -O bam -o ./mapping2.sorted.bam ./mapping2.bam >&2)

### index bam file for realignment

if [ $index == 1 ]; then
	$(samtools index ./mapping2.sorted.bam >&2)
fi


if [ "$VERBOSE" == "true" ]; then
	echo "I have finished file tranformation, I will do realignment using GATK 
		if you ask me to, otherwise I will directly do variant calling"
fi

#### relignment
### index reference (create ref.dist and ref.fai file required by next step)
if [ $realign == 1 ]; then
	echo "Doing realignment"
	rm *.dist
	$(java -jar $picard_path CreateSequenceDictionary -R $ref -O ${ref%.*}.dict >&2)
	$(samtools faidx $(basename $ref) >&2)
	### index relign
	$(java -Xmx2g -jar $GATK_path -T RealignerTargetCreator -R $(basename $ref) -I ./mapping2.sorted.bam -o ./lane.intervals --known $millsFile >&2)
	$(java -Xmx4g -jar $GATK_path -T IndelRealigner -R $(basename $ref) -I ./mapping2.sorted.bam -targetIntervals ./lane.intervals -o ./lane_realigned.bam >&2)
	if [ index == 1 ]; then
		$(samtools index ./lane_realigned.bam >&2)
	fi
	echo "I have finished realignment"
fi
#### calling variants
if [ "$VERBOSE" == "true" ]; then
	echo "I will start to call variant. Calling ... Please be patient."
fi

if [ $realign == 1 ]; then
	$(bcftools mpileup -Ou -f $ref ./lane_realigned.bam | bcftools call -vmO z -o $out_filename.vcf.gz >&2)
else
	$(bcftools mpileup -Ou -f $ref ./mapping2.sorted.bam | bcftools call -vmO z -o $out_filename.vcf.gz >&2)
fi

$(bcftools filter -O z -o $out_filename.filtered.vcf.gz -s LOWQUAL -i'%QUAL>10' $out_filename.vcf.gz >&2)


if [ "$VERBOSE" == "true" ]; then
	echo "I have finished variant calling. I will clean variant 
		calling output file to ses SNPs and insertions/deletions"
fi
 
### create human readable output file

$(gunzip -c $out_filename.filtered.vcf.gz > $out_filename.filted.vcf)
$(grep -vE '^#' $out_filename.filted.vcf | awk '{print $1,$2,$2+length ($5)-length ($4),length ($5)-length ($4)}' | sed 's/chr//' > file.bed)
$(grep -vE '^#' $out_filename.filted.vcf | awk '{print $1,$2,$2+length ($5)-length ($4),length ($5)-length ($4)}' | sed 's/chr//' | awk '{if ($4 != 0){print $1,$2,$3,$4};}' > indel.txt)
$(grep -vE '^#' $out_filename.filted.vcf | awk '{print $1,$2,$2+length ($5)-length ($4),length ($5)-length ($4)}' | sed 's/chr//' | awk '{if ($4 == 0){print $1,$2,$3,$4};}' > snps.txt)
$(rm $out_filename.filted.vcf)


if [ $gunzip1 == 0 ]; then
	gunzip -c $out_filename.filtered.vcf.gz > $out_filename.filted.vcf
	cp $out_filename.filted.vcf $temp_dir
else
	cp $out_filename.filtered.vcf.gz $temp_dir
fi

if [ "$VERBOSE" == "true" ]; then
	echo "done"
fi