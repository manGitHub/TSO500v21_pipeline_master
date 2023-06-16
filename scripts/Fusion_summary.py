#! /usr/bin/env python
# coding: utf-8


import glob
import pandas as pd
import csv
import re
import os
import sys
import os.path

#if os.path.isfile("/data/Compass/dev/gangalapudiv2/NextSeq/TSO500_Demux/(sys.argv[3])_RNA/Reports/(sys.argv[3])_RNA.html"):

if (sys.argv[3])  == 'None' :
    print("no RNA")
else: 
    path = ','.join(sys.argv[3:])
    mypath =path.split(",")
#    print(mypath)
    fusionpath= sys.argv[2]
    file = sys.argv[1]
    df = pd.concat([pd.read_csv(f,delimiter='\t',skiprows=9) for f in mypath])
    df = df[df['Summary'].str.contains("Failed at Fusion", na = False)]

#print(df)

    if df.empty:
       print('No fusion failure')
    else:
       df.to_csv(fusionpath + file, header=None, index=None)
       msg = ['TSO500 pipeline 1.0 with TSO500_app 1.3.1 ,TSO170_app 1.0.1.0 completed successfully on rest of the samples']
       newline = ''
       with open(fusionpath + file, 'a') as f:
          writer = csv.writer(f)
          writer.writerow(newline)
          writer.writerow(msg)
