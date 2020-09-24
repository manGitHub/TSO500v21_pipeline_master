#!/usr/bin/env python
import glob
import pandas as pd
import csv
import re
import os
import sys
import os.path


# This script combines the MSI output from all the DNA samples in the run. 
# This script takes two inputs, 1: Name of the output file 2: Space seperated list of all the *.msi.json files
# Usage:  MSI.py path/to/output/file.xlsx <List of all the json files> 



if (sys.argv[2])  == 'None' :
    print "no DNA"
else:
    items = {}
    path = (sys.argv[1])
    all = ','.join(sys.argv[2:])
    mypath =all.split(",")
#print(mypath)
    for file in mypath:
        name = os.path.basename(file)
        name = name.strip(".msi.json")
        json_file = open(file)
        for line in json_file:
            currentline = line.split(",")
        i = (filter(lambda x:"PercentageUnstableSites" in x,currentline))
        i = [elem.strip('"PercentageUnstableSites":') for elem in i]
        key, value = name,i
        items[key] = value       
    df = pd.DataFrame(items, index=['MSI - PercentageUnstableSites',])
    df = df.replace(r'"PercentageUnstableSites":','')
    df.to_excel(path,sheet_name="MSI_TMB",index=True)
