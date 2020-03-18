#!/usr/bin/env python3
#!/bin/bash

# coding: utf-8
## this step launches the rnaseqc.sh script to run rnaseqc tool
import glob
import pandas as pd
import csv
import re
import os
import sys
import subprocess

bam = glob.glob(sys.argv[3])
sample = sys.argv[1]
bam = str(bam).strip("['']")
bam = os.path.dirname(bam)
bam = bam + "/" + sample
head, sep, tail = bam.partition('RNA_IntermediateFiles')
qc = head + "rnaseqc"
done = head.rsplit('/',2)
done = done[0] + "/" + sys.argv[2] + "_" + sample + "_RNASeQC.done"
import subprocess

subprocess.run(["/data/Compass/Tools/TSO500_pipeline_vg/scripts/rnaseqc.sh", bam, sample,qc,done])
