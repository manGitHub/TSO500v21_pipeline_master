#!/usr/bin/env python
import glob
import pandas as pd
import csv
import re
import os
import sys
import os.path

# This script combines the rnaseqc tool output for all the RNA samples in the run
# This script takes three inputs. 1. Name of the output file 2. path to the output directory 3. Space seperated list of the metrics.tsv file for all the RNA samples.
# usage: merge_rnaseqc.py outputfilename.xlsx /path/to/the/results/folder  <list of rnaseqc metrics.tsv files> 

if (sys.argv[3])  == 'None' :
    print "no RNA"
else: 
    path = ','.join(sys.argv[3:])
    my_list = path.split(",")
    mypath= [file for file in my_list if os.path.exists(file)]
    qc= sys.argv[2]
    file = sys.argv[1]
    df = pd.concat([pd.read_csv(f,delimiter='\t') for f in mypath])
    df=df[["Sample","Total Purity Filtered Reads Sequenced","Mapped","Mapping Rate","Mapped Unique","Mapped Unique Rate of Total","Unique Rate of Mapped","Duplication Rate of Mapped","rRNA rate","Exonic Rate","Intronic Rate"]]
    df.to_excel(qc + file,index=False,sheet_name="RNA_QC")
