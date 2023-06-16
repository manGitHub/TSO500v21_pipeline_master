#!/bin/bash

# This script runs fail exons, fail genes and hotspot depth for the TSO500 bam files.

# python3 {params.script}/DNA_qc.py {params.dir}/Logs_Intermediates/StitchedRealigned/*/*.bam {params.dir} {params.bed} {params.hotspot} {params.size} {params.script}

# $1 = bam file

# $2 = path to output directory

# $3 = bed file

# $4 = hotspot file

# $5 = genome size

# $6 = path to the failed_Exon_Final.pl script


#module load bedtools/2.22.0 samtools/0.1.19
module load samtools/0.1.19
cd $1

echo -e "chr\tstart\tend\tgene\tposition\tdepth" > $2.depth_per_base

samtools view -hF 0x400 -q 30 -L $5  $4 | samtools view -ShF 0x4 - | samtools view -SuF 0x200 - | /data/Compass/local/software/bedtools/2.22.0/bin/bedtools coverage -abam - -b $5 -d >> $2.depth_per_base

perl $6/failed_Exon_Final.pl $2.depth_per_base 50 $2.failExons $2.failGenes


slopBed -i $7 -g $8 -b 50 > $2_Region.bed
samtools view -hF 0x400 -q 30 -L $2_Region.bed $4 | samtools view -ShF 0x4 - | samtools view -SuF 0x200 - | /data/Compass/local/software/bedtools/2.22.0/bin/bedtools coverage -abam - -b $7 > $2.hotspot.depth
rm $2.depth_per_base $2_Region.bed

