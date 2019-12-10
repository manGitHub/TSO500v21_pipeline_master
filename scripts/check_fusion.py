#! /usr/bin/env python
# coding: utf-8


import glob
import pandas as pd
import csv
import re
import os
import sys
import os.path
import glob








#a = "/data/Compass/dev/gangalapudiv2/ProcessedResults_NexSeq/CP02106_T1R_PS/TruSightTumor170_Analysis_*/RNA_CP02106_T1R_PS/a"
fusion = (sys.argv[1])
dir = os.path.dirname(fusion)
dir = os.path.basename(dir)
sample = dir.strip("RNA_")
#a = "/".join(a.split("/", 2))
fusion = fusion.rsplit('/', 2)
fusion = fusion[0]
fusion = glob.glob(fusion)
fusion= str(fusion).strip("['']")

file = fusion + "/" + dir + "/" + sample + "_Fusions.csv" 

print(file)

if not os.path.exists(file):
    with open(file, 'w'): pass    
