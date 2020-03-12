#!/bin/bash

module load java/1.8.0_11
module load picard

java -jar $PICARD_JAR AddOrReplaceReadGroups VALIDATION_STRINGENCY=SILENT INPUT=$1.bam OUTPUT=$1_RG.bam SORT_ORDER=coordinate RGLB=$2 RGPU=$2 RGPL=ILLUMINA RGSM=$2 RGCN=Compass

java -jar $PICARD_JAR MarkDuplicates AS=true M=$1_markdp.txt I=$1_RG.bam O=$1_dd.bam REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT

module load samtools 

samtools index $1_dd.bam
module load rnaseqc

java -jar $RNASEQCPATH/RNA-SeQC_v1.1.8.jar -r /data/Compass/Ref/hg19/TSO170_Res/genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa  -rRNA /data/Compass/Tools/TSO500_pipeline_vg/scripts/rRNA_interval_TSO170.txt -o $3 -s "$2|$1_dd.bam|$2" -t /data/Compass/Ref/hg19/TSO170_Res/genomes/Homo_sapiens/UCSC/hg19/Annotation/gencode.v19.annotation_filtered_20160301.gtf

rm $1_RG.bam $1_markdp.txt $1_dd.bam $1_dd.bam.bai
