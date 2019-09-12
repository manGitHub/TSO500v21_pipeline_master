#! /usr/bin/env python
# coding: utf-8


import glob
import pandas as pd
import csv
import re
import os
import sys

#path = r'/data/Compass/dev/gangalapudiv2/qc_combine_Script/*/*/' # use your path

#all = glob.glob(path + "/RNA_SampleMetricsReport.txt")
path = ','.join(sys.argv[3:])
mypath =path.split(",")
print mypath
qc= sys.argv[2]
file = sys.argv[1]

#print all
#df_2= df.iloc[(df.loc[df[0]=='report field'].index[0]+1):, :].reset_index(drop = True)
df = pd.concat([pd.read_csv(f,delimiter='\t',skiprows=7) for f in mypath])

#df = df.loc[:,~df.columns.duplicated()]
#df = df.T.drop_duplicates().T
#print df

df.to_excel(qc + file,index=False,sheet_name="Run_metrics")
