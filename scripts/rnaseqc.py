#!/usr/bin/env python3
#!/bin/bash

# coding: utf-8

import glob
import pandas as pd
import csv
import re
import os
import sys
import subprocess

bam = glob.glob(sys.argv[1])
sample = sys.argv[2]
bam = str(bam).strip("['']")
bam = os.path.dirname(bam)
bam = bam + "/" + sample
head, sep, tail = bam.partition('RNA_IntermediateFiles')
#print(bam)
qc = head + "rnaseqc"
#print(qc)
import subprocess
#head = sys.argv[1]
sample = sys.argv[2]

subprocess.run(["/data/Compass/Tools/TSO500_pipeline_vg/scripts/rnaseqc.sh", bam, sample,qc])
