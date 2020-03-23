#! /usr/bin/env python
# coding: utf-8
## this script combined the RNASeQC output for all the samples in a run.

import glob
import pandas as pd
import csv
import re
import os
import sys
import os.path


if (sys.argv[3])  == 'None' :
    print "no RNA"
else: 
    path = ','.join(sys.argv[3:])
    mypath =path.split(",")
    print mypath
    qc= sys.argv[2]
    file = sys.argv[1]

    df = pd.concat([pd.read_csv(f,delimiter='\t') for f in mypath])
    df=df[["Sample","Total Purity Filtered Reads Sequenced","Mapped","Mapping Rate","Mapped Unique","Mapped Unique Rate of Total","Unique Rate of Mapped","Duplication Rate of Mapped","rRNA rate","Exonic Rate","Intronic Rate"]]
    df.to_excel(qc + file,index=False,sheet_name="Run_metrics")

