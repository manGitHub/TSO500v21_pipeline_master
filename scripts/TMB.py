#!/usr/bin/env python
import glob
import pandas as pd
import csv
import re
import os
import sys
import os.path
from openpyxl import load_workbook
import numpy as np


# This script combines the TMB output from all the DNA samples in the run. 
# This script takes two inputs, 1: Name of the output file 2: Space seperated list of all the *.tmb.json files
# Usage:  MSI.py path/to/output/file.xlsx <List of all the json files>



if (sys.argv[2])  == 'None' :
    print("no DNA")
else: 
    items = {}
    path = (sys.argv[1])
    all = ','.join(sys.argv[2:])
    mypath =all.split(",")
    for file in mypath:
        name = os.path.basename(file)
        name = name.strip(".tmb.json")
        json_file = open(file)
        for line in json_file:
            currentline = line.split(",")
        i = (filter(lambda x:"AdjustedNonsynonymousTmbPerMb" in x,currentline))
        i = [elem.strip('"AdjustedNonsynonymousTmbPerMb":') for elem in i]
        key, value = name,i
        items[key] = value       
    df = pd.DataFrame(items, index=['TMB - AdjustedNonsynonymousTmbPerMb',])
    df = df.replace(r'"AdjustedNonsynonymousTmbPerMb":','')
    df = df.astype(float).round(2)
    book = load_workbook(path)
    writer = pd.ExcelWriter(path, engine='openpyxl') 
    writer.book = book
    writer.sheets = {ws.title: ws for ws in book.worksheets}
    for sheetname in writer.sheets:
        df.to_excel(writer,sheet_name=sheetname, startrow=5,index=True)
    writer.save()
