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

bam = glob.glob(sys.argv[1])
bam = str(bam).strip("['']")
file = os.path.basename(bam)
sample = str(file).strip("[.bam]")
bam = os.path.dirname(bam)
bamfile = bam + "/" + sample
gtf = sys.argv[2]
genome =sys.argv[3]
dir = sys.argv[4] + "/Results"
rnaseqc = sys.argv[5] + "/rnaseqc.sh" 
subprocess.run([rnaseqc,bamfile,sample,gtf,genome,dir])

