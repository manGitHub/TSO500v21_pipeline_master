#!/bin/bash


module load bedtools/2.22.0 samtools/0.1.19

cd $1

sample=$(echo $1 |awk -F"/" '{ print $(NF-1) }')

echo -e "chr\tstart\tend\tgene\tposition\tdepth" > $sample.depth_per_base

samtools view -hF 0x400 -q 30 -L $3  $2 | samtools view -ShF 0x4 - | samtools view -SuF 0x200 - | bedtools coverage -abam - -b $3 -d >> $sample.depth_per_base

perl $4/failed_Exon_Final.pl $sample.depth_per_base 50 ${sample}_${5}.failExons ${sample}_${5}.failGenes


slopBed -i $6 -g $7 -b 50 > ${sample}_Region.bed
samtools view -hF 0x400 -q 30 -L ${sample}_Region.bed $2 | samtools view -ShF 0x4 - | samtools view -SuF 0x200 - | bedtools coverage -abam - -b $6 > ${sample}_${5}.hotspot.depth
rm $sample.depth_per_base ${sample}_Region.bed
