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
file = sys.argv[1]

#print all
df = pd.concat([pd.read_csv(f,delimiter='\t',skiprows=7) for f in mypath])


df.to_excel(qc + file,index=False,sheet_name="Run_metrics")
