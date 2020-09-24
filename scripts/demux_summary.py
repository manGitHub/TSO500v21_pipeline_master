import sys
import glob, os
import pandas as pd

inPath = sys.argv[1] + '/Lane*/Demultiplex_Stats.csv'
outFile = sys.argv[2]
outFile = open(outFile,'w')
sampSheet = sys.argv[3]

# samplesheet
df = pd.read_csv(sampSheet,delimiter=',',skiprows=24)
df = pd.DataFrame(df)
DNA_list = df.loc[df['Sample_Type'] == 'DNA','Sample_ID'].tolist()
RNA_list = df.loc[df['Sample_Type'] == 'RNA','Sample_ID'].tolist()

# append files into one df each for rna/dna
dat = pd.concat(pd.read_csv(dFile, index_col=None, header=0) for dFile in sorted(glob.glob(inPath)))

# filter DNA out of RNA reports ... RNA out of DNA reports
# if no DNA/RNA samples for respective report, print error
if 'DNA_Reports' in inPath:
    dat.loc[~dat.SampleID.isin(DNA_list),'Index']='Undetermined'
    dat.loc[~dat.SampleID.isin(DNA_list),'SampleID']='Undetermined'
elif 'RNA_Reports' in inPath:
    dat.loc[~dat.SampleID.isin(RNA_list),'Index']='Undetermined'
    dat.loc[~dat.SampleID.isin(RNA_list),'SampleID']='Undetermined'

# define col order for lane output
# colOrder = ['Lane', '# Reads', '% of the lane', '# Perfect Index Reads', '% Perfect Index Reads','# One Mismatch Index Reads', '% One Mismatch Index Reads','# of >= Q30 Bases (PF)', '% of >= Q30 Bases (PF)','Mean Quality Score (PF)']
colOrder = ['Lane', '# Reads', '% of the lane', '# Perfect Index Reads', '% Perfect Index Reads','# One Mismatch Index Reads', '% One Mismatch Index Reads','# of >= Q30 Bases (PF)', 'Mean Quality Score (PF)']

# summarize by lane
outFile.write('Lane Summary\n')
laneSum = dat.groupby('Lane').agg({'# Reads':'sum','# Perfect Index Reads':'sum','# One Mismatch Index Reads':'sum','# of >= Q30 Bases (PF)':'sum','Mean Quality Score (PF)':'mean'}).reset_index()
totReads = laneSum['# Reads'].sum()
laneSum['% of the lane'] = (laneSum['# Reads'] / totReads) * 100
laneSum['% Perfect Index Reads'] = (laneSum['# Perfect Index Reads'] / laneSum['# Reads']) * 100
laneSum['% One Mismatch Index Reads'] = (laneSum['# One Mismatch Index Reads'] / laneSum['# Reads']) * 100
# when I uncomment this, use proper column order above
# laneSum['% of >= Q30 Bases (PF)'] = (laneSum['# of >= Q30 Bases (PF)'] / (laneSum['# Reads']*101)) * 100
laneSum = laneSum[colOrder]
outFile.write(laneSum.to_csv(index=False)+'\n')

# define col order for sample output
# colOrder = ['SampleID', 'Index', '# Reads', '% of the lane', '# Perfect Index Reads', '% Perfect Index Reads','# One Mismatch Index Reads', '% One Mismatch Index Reads','# of >= Q30 Bases (PF)', '% of >= Q30 Bases (PF)','Mean Quality Score (PF)']
colOrder = ['SampleID', 'Index', '# Reads', '% of the lane', '# Perfect Index Reads', '% Perfect Index Reads','# One Mismatch Index Reads', '% One Mismatch Index Reads','# of >= Q30 Bases (PF)','Mean Quality Score (PF)']

# summarize by sample
outFile.write('Sample Summary\n')
sampSum = dat.groupby(['SampleID','Index']).agg({'# Reads':'sum','# Perfect Index Reads':'sum','# One Mismatch Index Reads':'sum','# of >= Q30 Bases (PF)':'sum','Mean Quality Score (PF)':'mean'}).reset_index()
totReads = sampSum['# Reads'].sum()
sampSum['% of the lane'] = (sampSum['# Reads'] / totReads) * 100
sampSum['% Perfect Index Reads'] = (sampSum['# Perfect Index Reads'] / sampSum['# Reads']) * 100
sampSum['% One Mismatch Index Reads'] = (sampSum['# One Mismatch Index Reads'] / sampSum['# Reads']) * 100
# when I uncomment this, use proper column order above
# sampSum['% of >= Q30 Bases (PF)'] = (sampSum['# of >= Q30 Bases (PF)'] / (sampSum['# Reads']*101)) * 100
sampSum = sampSum[colOrder]
outFile.write(sampSum.to_csv(index=False)+'\n')

# define col order for disagg output
# colOrder = ['Lane','SampleID', 'Index', '# Reads', '% of the lane', '# Perfect Index Reads', '% Perfect Index Reads','# One Mismatch Index Reads', '% One Mismatch Index Reads','# of >= Q30 Bases (PF)', '% of >= Q30 Bases (PF)','Mean Quality Score (PF)']
colOrder = ['Lane','SampleID', 'Index', '# Reads', '% of the lane', '# Perfect Index Reads', '% Perfect Index Reads','# One Mismatch Index Reads', '% One Mismatch Index Reads','# of >= Q30 Bases (PF)','Mean Quality Score (PF)']

# disagg data
outFile.write('Breakdown by Lane and Barcode\n')
totReads = dat['# Reads'].sum()
dat['% of the lane'] = (dat['# Reads'] / totReads) * 100
dat['% Perfect Index Reads'] = (dat['# Perfect Index Reads'] / dat['# Reads']) * 100
dat['% One Mismatch Index Reads'] = (dat['# One Mismatch Index Reads'] / dat['# Reads']) * 100
# when I uncomment this, use proper column order above
# dat['% of >= Q30 Bases (PF)'] = (dat['# of >= Q30 Bases (PF)'] / (dat['# Reads']*101)) * 100
dat = dat[colOrder]
outFile.write(dat.to_csv(index=False))

outFile.close()
