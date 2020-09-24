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

manifest = sys.argv[3]
bam = glob.glob(sys.argv[1])
hotspot = sys.argv[4]
genome = sys.argv[5]
bam = str(bam).strip("['']")
file = os.path.basename(bam)
bam = os.path.dirname(bam)
bamfile = bam + "/" + file
dir = sys.argv[2]
sample = os.path.basename(dir)
#sample = str(pair).strip("[_Pair]")
script = sys.argv[6]
subprocess.run(["/data/Compass/dev/gangalapudiv2/TSO_new/scripts/TSO500_QC.sh",dir,sample,file,bamfile,manifest,script,hotspot,genome])
