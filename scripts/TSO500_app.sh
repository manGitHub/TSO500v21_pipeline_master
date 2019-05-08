#!/bin/bash

sample=$(echo $1 | sed 's/.*\(FastqFolder\)//'|sed 's/\///g')

path=$(echo $1 |sed 's/\(FastqFolder\).*/\1/g')

cd $path

#for file in $sample/*fastq.gz ;do dir=${file%/*fastq.gz};mkdir  "$dir/$dir"; mv "$file" "$dir/$dir" ;done
for file in $sample/*fastq.gz ; do 
    dir=${file%/*fastq.gz}
    mkdir -p "$dir/$dir"
    mv "$file" "$dir/$dir"
done


module load singularity
cd $2
if [ -d "$3/$sample" ]; then rm -Rf $3/$sample; fi
mkdir -p /$3/$sample
./TruSight_Oncology_500.sh --engine singularity --resourcesFolder $4 --fastqFolder $1 --analysisFolder $3/$sample

