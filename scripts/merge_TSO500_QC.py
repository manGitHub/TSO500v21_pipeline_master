#!/usr/bin/env python
# coding: utf-8

# This script combines the results from TSO500 app of all the samples.
# sys.argv[1] = excel file name ex:TSO500_QC_{runid}.xlsx {RESULT_DIR}/run_qc/ {TSO500_QC}
# sys.argv[2] = path to output directory  
# sys.argv[3] = Comma seperated list of TSO500 summary files *MetricsOutput.tsv
 
# Usage: merge_TSO500_QC.py TSO500_QC.xlsx /path/to/dir Comma seperated list of files

import glob
import pandas as pd
import csv
import re
import os
import sys
import numpy as np

path = ','.join(sys.argv[3:])
mypath =path.split(",")
qc= sys.argv[2]  
##path to qc report
file = sys.argv[1]
df = pd.concat([pd.read_csv(f,delimiter='\t',skiprows =8,skipfooter=3,engine='python') for f in mypath], axis=1)
df = df.T.drop_duplicates().T   # to remove duplicated columns by column contents
df1 = df.filter(like='Unnamed: 3')
df_header = pd.read_csv("/data/Compass/dev/gangalapudiv2/TSO_new/scripts/QC_template.csv",delimiter='\t',engine='python')
df_final = pd.concat([df_header,df1], axis=1)
df_final.to_excel(qc + file,index=False,header=False,sheet_name="Run_metrics")
