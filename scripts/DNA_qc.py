#! /usr/bin/env python
# coding: utf-8

import glob
import pandas as pd
import csv
import re
import os
import sys


path = ','.join(sys.argv[3:])
mypath =path.split(",")
print mypath
qc= sys.argv[2]  
##path to qc report

file = sys.argv[1]
#print file    #output_file_name

df = pd.concat([pd.read_csv(f,delimiter='\t',skiprows =9,skipfooter=3,engine='python') for f in mypath], axis=1)
df = df.T.drop_duplicates().T   # to remove duplicated columns by column contents

df.to_excel(qc + file,index=False,header=False,sheet_name="Run_metrics")

