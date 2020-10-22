#!/usr/bin/env python3
from pprint import pprint as pp
from collections import defaultdict
import glob
import os
import pandas as pd
import numpy as np
from collections import OrderedDict
import re
import sys

analysis = sys.argv[2]

SAMPLESHEET = sys.argv[1]

DATA = defaultdict(list)
PAIR = defaultdict(list)
linefull = []
with open(SAMPLESHEET, 'rt') as S:
    text = S.read()  #reading in the samplesheet file
    lines = text.split('Sample_ID')  #splitting the file at sample_ID which is the first column header, hence creating two elements
    lines2 = lines[1].split('\n')   #splitting the second element which is the table by line
    for L in lines2:
        if L.find('RNA') != -1: # getting list of RNA samples
            lines3 = L.split(',')
            linefull.append(lines3)
            lines3[0] = lines3[0].strip()             
            DATA['RNA'].append(lines3[0])
        if L.find('DNA') != -1:
            lines3 = L.split(',')
            lines3[0] = lines3[0].strip()
            DATA['DNA'].append(lines3[0])

lines3 = []
PAIRDR = dict()
PAIRRD = dict()
PAIRDR2 = dict()
tops = []
with open(SAMPLESHEET, 'rt') as IN:
    text = IN.read()
    text = text.replace('\r', '')
    lines = text.split('[Data]')
    top = text.split('Sample_ID')
    top = top[0].split('\n')
    for i in top:
        tops.append(i.split(','))
    lines2 = lines[1].split('\n')
    header = lines2[1].split(',')
    lines2 = lines2[2:]
    for L in lines2:
        lines3.append(L.split(','))
    top =pd.DataFrame(tops)
    top = top.mask(top.eq('None')).dropna()
    df = pd.DataFrame(lines3, columns = header)
    df.dropna(axis = 0, how = 'all', thresh = 7, inplace = True)
    df['Pair'] = df['Pair'].str.replace('_Pair','')
#    pp(df) this is the complete samplesheet dataframe
    emptyPair =	 df[df.Pair == '']
    if not emptyPair.empty: 
         print("--------------------------------------------")
         print("--------------------------------------------")
         print("The following sample is missing Pairing information. Please fix it, before launching the pipeline.")
         print("\n")
         pp(emptyPair)
         print("--------------------------------------------")
         print("--------------------------------------------")
    df = df[df.Pair != '']
    df0 = df.loc[(df['Sample_Type'] == 'DNA') & (df['Pair'] != 'PENDING') & (df['Pair'] != "NA")]
    df01 = df0.isin(DATA['RNA'])
    otherRNA = df01[df01['Pair'] == False]
    if not otherRNA.empty:
        df02 = df.loc[otherRNA.index.tolist()]        
        oldrna = df02
        df02 =df02.drop('Pair',axis=1)
        df02['Pair'] = df02['Sample_ID'] + "_Pair"
        oldrna = oldrna.rename(columns = {"Sample_ID":"Pair","Pair":"Sample_ID"})
       	oldrna = oldrna.replace({"DNA":"RNA"})       
        oldrna = oldrna[['Sample_ID','Sample_Name','Sample_Plate','Sample_Well','Index_ID','index','index2','Sample_Type','Pair']]
        oldrna['Pair'] = oldrna['Pair'] + "_Pair"
        merge1 = pd.concat([df02,oldrna], axis =0)
        PAIRDR2 = df02.set_index('Sample_ID').T.to_dict() 
#        pp(merge1)This is New DNA and Old RNA
    else:
        merge1 = pd.DataFrame()
    df1 = df.loc[(df['Sample_Type'] == 'DNA') & (df['Pair'] != 'PENDING') & (~df.Sample_ID.isin(PAIRDR2.keys()))] 
    nomatch = df1.loc[(df1['Pair'] == "NA")]
    nomatch = nomatch.drop('Pair',axis=1)
    nomatch['Pair'] = nomatch['Sample_ID'] + "_NoPair"
#    pp(nomatch)  This is DNA only samples
    if not df1.empty:
        newDRpair = df1.loc[(df1['Pair'] != "NA")]
        newDRpair = newDRpair.drop('Pair',axis=1)
        newDRpair['Pair'] = newDRpair['Sample_ID'] + "_Pair"

    PAIRDR = df1.set_index('Sample_ID').T.to_dict()
    df2 = df.loc[df['Sample_Type'] == 'RNA']
    df3 = df2.isin(PAIRDR.keys())
    samernadna = df3[df3['Pair'] == True]
    if not samernadna.empty:
       df7 = df.loc[samernadna.index.tolist()]
       df7['Pair'] = df7['Pair'] + "_Pair"
       merge2 = pd.concat([newDRpair,df7],axis=0)
#       pp(merge2)  DNA RNA sequenced in the same run
    else:
       merge2 = pd.DataFrame()        
    other = df3[df3['Pair'] == False]
#    pp(other)
    if not other.empty:
        df4 = df.loc[other.index.tolist()]
        df6 = df4.loc[df4['Pair'] == "NA"] 
        df6 = df6.drop('Pair',axis=1)
        df6['Pair'] = df6['Sample_ID'] + "_NoPair"
        df4 = df4.loc[(df4['Pair'] != "NA") & (df4['Pair'] != 'PENDING')]
        olddna = df4        
        df4['Pair'] = df4['Pair'] + "_Pair"
        olddna = olddna.rename(columns = {"Sample_ID":"Pair","Pair":"Sample_ID"})
        olddna = olddna.replace({"RNA":"DNA"})
        olddna = olddna[['Sample_ID','Sample_Name','Sample_Plate','Sample_Well','Index_ID','index','index2','Sample_Type','Pair']]
        olddna = olddna.drop('Pair',axis=1)
        olddna['Pair'] = olddna['Sample_ID'] 
        olddna['Sample_ID'] = olddna['Sample_ID'].str.replace('_Pair', '')        
        merge3 = pd.concat([df4,olddna],axis=0)
#        pp(merge3)
#        df5 = df4.loc[(df4['Pair'] != "NA") & (df4['Pair'] != 'PENDING')] #The NA here implies unpaired RNA.
    else:
        merge3 = pd.DataFrame()
        df6 = pd.DataFrame()
head =['Sample_ID','Sample_Name','Sample_Plate','Sample_Well','Index_ID','index','index2', 'Sample_Type', 'Pair_ID']
head = pd.DataFrame(head)
head = head.transpose()
merged = pd.concat([merge1,merge2,merge3,df6,nomatch],axis =0,)
merged.columns = range(merged.shape[1])
mf = pd.concat([top,head,merged],axis =0,) 
#pp(mf)
mf.to_csv(analysis,index=False,sep=',',header=False)

