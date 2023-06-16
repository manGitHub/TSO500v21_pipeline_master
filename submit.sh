#!/bin/bash

snakemake -rp --nolock --latency-wait 120 -k -j 3000 -s $PIPELINE_HOME/TSO500.snakefile --configfile $YAML -d `pwd` --jobname {params.rulename}.{jobid} --cluster "sbatch -o logs/{params.rulename}.%j.o -e logs/{params.rulename}.%j.e {params.resources}"
