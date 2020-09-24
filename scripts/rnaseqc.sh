#!/bin/bash

# this script runs picard preprocessing and RNASeQC

# $1 = bamfile
# $2 = samplename 
# $3 = path to annotation file 
# $4 = path to genome file
# $5 = output file directory

module load java/1.8.0_11
module load picard

java -jar $PICARD_JAR AddOrReplaceReadGroups VALIDATION_STRINGENCY=SILENT INPUT=$1.bam OUTPUT=$1_RG.bam SORT_ORDER=coordinate RGLB=$2 RGPU=$2 RGPL=ILLUMINA RGSM=$2 RGCN=Compass

java -jar $PICARD_JAR MarkDuplicates AS=true M=$1_markdp.txt I=$1_RG.bam O=$1_dd.bam REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT

module load samtools 

samtools index $1_dd.bam
module load rnaseqc

java -jar $RNASEQCPATH/RNA-SeQC_v1.1.8.jar -r $4  -o $5/rnaseqc -s "$2|$1_dd.bam|$2" -t $3

#rm $1_RG.bam $1_markdp.txt $1_dd.bam $1_dd.bam.bai

#touch $4
