#!/bin/bash

snakemake -rp --nolock --latency-wait 60 -k -j 3000 -s $PIPELINE_HOME/TSO500.snakefile --configfile $YAML -d `pwd` --jobscript $PIPELINE_HOME/scripts/jobscript.sh --jobname {params.rulename}.{jobid} --cluster "sbatch -o logs/{params.rulename}.%j.o -e logs/{params.rulename}.%j.e {params.resources}"
