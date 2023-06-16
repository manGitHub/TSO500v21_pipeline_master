#! /usr/bin/env python
# coding: utf-8

import glob
import pandas as pd
import csv
import re
import os
import sys


if (sys.argv[3])  == 'None' :
    print("no DNA")
else: 
    path = ','.join(sys.argv[3:])
    mypath =path.split(",")
    #print(mypath)

    qc= sys.argv[2]

    file = sys.argv[1]

    df = pd.concat([pd.read_csv(f,delimiter='\t',skiprows =6,engine='python') for f in mypath], axis=1)
    df = df.T.drop_duplicates().T
    #print df
    df.to_excel(qc + file,index=False,sheet_name="TMB_MSI")


