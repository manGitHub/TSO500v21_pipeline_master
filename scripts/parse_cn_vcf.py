import sys
import glob, os
import pandas as pd
from itertools import groupby


# open input/output files
inFile = sys.argv[1]
geneList = sys.argv[2]

# read in foldChange header/tbl and gene list
with open(inFile) as myfile:
    header = [next(myfile) for x in range(2)]
tbl = pd.read_csv(inFile,delimiter='\t',header=2)
geneList = pd.read_csv(geneList,header=None)

# subset foldChange table
tbl = tbl.loc[tbl.GeneName.isin(geneList[0])]

# get gene start/end pos from manifest
manifest = sys.argv[4]
manifest = pd.read_csv(manifest,delimiter='\t')
manifest = manifest.groupby(['chr','gene']).agg({'start':'min','end':'max'}).reset_index()

# update manifest with new gene coords from QCI
qci_coords = sys.argv[7]
qci_coords = pd.read_csv(qci_coords,delimiter=',',nrows=12)
qci_coords = qci_coords[['Chromosome','Gene','QCI Start','QCI End']]
qci_coords.columns = ['chr','gene','start','end']
manifest = pd.concat([manifest,qci_coords]).drop_duplicates(['chr','gene'],keep='last').sort_values(['chr','start'])

# get vcf header, change FORMAT field FC to CN
vcfFile = sys.argv[3]
with open(vcfFile) as f:
    vcfHeader = [next(f) for x in range(35)]
FORMAT = [heads for heads in vcfHeader if '##FORMAT' in heads]
FORMAT = ['##FORMAT=<ID=CN,Number=1,Type=Float,Description="Copy number variant fold change, calculated as ratio of sample to control">\n']
noFormat = [heads for heads in vcfHeader if '##FORMAT' not in heads]
vcfHeader = noFormat + FORMAT
# read in and empty vcf df
vcfTbl = pd.read_csv(vcfFile,delimiter='\t',header=35)
vcfTbl = vcfTbl[0:0]

# get genes/qualScore/FC from 92 gene subset
vcfTbl['gene'] = tbl.GeneName
vcfTbl.iloc[:,9] = tbl.FC
vcfTbl.QUAL = tbl.qScore
vcfTbl = vcfTbl.round({'QUAL':0})

# merge w manifest, repopulate columns
vcfTbl = pd.merge(vcfTbl,manifest,on='gene',how='left')
vcfTbl.start = vcfTbl.start.astype(int)
vcfTbl.end = vcfTbl.end.astype(int)
vcfTbl['#CHROM'] = vcfTbl.chr
vcfTbl.POS = vcfTbl.start
vcfTbl.ID = '.'
vcfTbl.FILTER = 'PASS'
vcfTbl.INFO = 'END=' + vcfTbl.end.astype(str) + ';ANT=' + vcfTbl.gene.astype(str)
vcfTbl.FORMAT = 'CN'
vcfTbl = vcfTbl.round({vcfTbl.columns[9]:3})

# merge with thresholds, add <dup> <del>
thresh = sys.argv[5]
thresh = pd.read_csv(thresh,delimiter=',')
thresh.columns.values[0] = 'gene'
vcfTbl = pd.merge(vcfTbl,thresh,on='gene',how='left')
vcfTbl['ALT'] = '.'
# del
vcfTbl.loc[vcfTbl.iloc[:,9] < vcfTbl['DEL.Cutoff'], 'ALT'] = '<DEL>'
vcfTbl.loc[vcfTbl.iloc[:,9] < vcfTbl['DEL.Cutoff'], 'INFO'] = 'SVTYPE=CNV;' + vcfTbl['INFO']
vcfTbl.loc[vcfTbl.iloc[:,9] < vcfTbl['DEL.Cutoff'], 'POS'] = vcfTbl['POS'] - 1
# dup
vcfTbl.loc[vcfTbl.iloc[:,9] > vcfTbl['Amp.Cutoff'], 'ALT'] = '<DUP>'
vcfTbl.loc[vcfTbl.iloc[:,9] > vcfTbl['Amp.Cutoff'], 'INFO'] = 'SVTYPE=CNV;' + vcfTbl['INFO']
vcfTbl.loc[vcfTbl.iloc[:,9] > vcfTbl['Amp.Cutoff'], 'POS'] = vcfTbl['POS'] - 1
#remove excess columns
vcfTbl = vcfTbl[vcfTbl.columns[:-8]]

# go through genome and get reference nucleotides
genome = sys.argv[6]
def fasta_iter(fasta_name):
    # first open the file outside
    fh = open(fasta_name)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        headerStr = header.__next__()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.__next__())
        yield(headerStr,seq)

# convert pandas to numpy and get chrs and start pos
vcf_numpy = vcfTbl.to_numpy()
chrs = vcf_numpy[:,0]
starts = vcf_numpy[:,1]-1           # acct for py 0 ind

# fill array REF with nucs, add to df
ref_seq = {}
REFs = []
for header, sequence in fasta_iter(genome):
    ref_seq[header] = sequence
for chro,st in zip(chrs,starts):
    REFs.append((ref_seq[chro][st]).upper())
vcfTbl['REF'] = REFs

outFileVCF = vcfFile.replace('_CopyNumberVariants.vcf','_CopyNumberVariants_92genes.vcf')
outFileVCF = open(outFileVCF,'w')
for heads in vcfHeader:
    outFileVCF.write(heads)
outFileVCF.write(vcfTbl.to_csv(index=False,sep='\t'))
outFileVCF.close()
