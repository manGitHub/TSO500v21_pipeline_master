from pprint import pprint as pp
from collections import defaultdict
import glob
import os
import pandas as pd
import numpy as np

PIPELINE = os.environ['PIPELINE_HOME']
RUN_DIR = os.environ['RUN_DIR']
DEMUX_DIR = os.environ['DEMUX_DIR']
configfile: PIPELINE + '/config/cluster.json'
GROUP = config['group']
VERSION = config['version']
APP_DIR = os.environ['APP_DIR']
RESULT_DIR = os.environ['RESULT_DIR']
DATA = defaultdict(list)

F = config['RUNID'].split('_')
flowcell = ''
for f in F:
    if re.match('^\d{4}$', f):
        flowcell = re.sub('^\w', '', F[F.index(f) + 1])

samplesheet = RUN_DIR + '/' + config['RUNID'] + '/SampleSheet.csv'

with open(samplesheet, 'rt') as S:
    text = S.read()
    lines = text.split('Sample_ID')
    lines2 = lines[1].split('\n')
    for L in lines2:
        if L.find('PoolRNA') != -1:
            lines3 = L.split(',')
            lines3[0] = lines3[0].strip()
            DATA['RNA'].append(lines3[0])
        if L.find('PoolDNA') != -1:
            lines3 = L.split(',')
            lines3[0] = lines3[0].strip()
            DATA['DNA'].append(lines3[0])


lines3 = []
PAIRDR = dict()
PAIRRD = dict()
UNPAIRR = []
VCFdone_temp = []
VCFdone = []
with open(samplesheet, 'rt') as IN:
    text = IN.read()
    lines = text.split('[Data]')
    lines2 = lines[1].split('\n')
    header = lines2[1].split(',')
    lines2 = lines2[2:]
    for L in lines2:
        lines3.append(L.split(','))

    df = pd.DataFrame(lines3, columns = header)
    df.dropna(axis = 0, how = 'all', thresh = 10, inplace = True)
    #pp(df)
    df0 = df.loc[(df['Manifest'] == 'PoolDNA') & (df['Matched_RNA'] != 'PENDING') & (df['Matched_RNA'] != 'NA'), ['Matched_RNA']]
    df01 = df0.isin(DATA['RNA'])
    otherRNA = df01[df01['Matched_RNA'] == False]
    if not otherRNA.empty:
        df02 = df.loc[otherRNA.index.tolist(), ['Sample_ID', 'Matched_RNA']]
        PAIRDR2 = df02.set_index('Sample_ID').T.to_dict() 
    
    df1 = df.loc[(df['Manifest'] == 'PoolDNA') & (df['Matched_RNA'] != 'PENDING') & (~df.Sample_ID.isin(PAIRDR2.keys())), ['Sample_ID', 'Matched_RNA']] 
    PAIRDR = df1.set_index('Sample_ID').T.to_dict()
    
    df2 = df.loc[df['Manifest'] == 'PoolRNA', ['Matched_DNA']]
    df3 = df2.isin(PAIRDR.keys())
    other = df3[df3['Matched_DNA'] == False]
    if not other.empty:
        df4 = df.loc[other.index.tolist()]
        df5 = df4.loc[(df4['Matched_DNA'] != 'NA') & (df4['Matched_DNA'] != 'PENDING'), ['Sample_ID', 'Matched_DNA']] #The NA here implies unpaired RNA.
        PAIRRD = df5.set_index('Sample_ID').T.to_dict()
        df6 = df4.loc[df4['Matched_DNA'] == 'NA', ['Sample_ID']]
        UNPAIRR = df6['Sample_ID'].values.tolist()
        
VCFdone_temp = ['{rdir}/{dna}/Results/{dna}_{rna}.{runid}.merge_vcf.done'.format(rdir = RESULT_DIR, rna = PAIRDR[sample]['Matched_RNA'], dna = sample, runid = config['RUNID']) for sample in PAIRDR.keys()]
VCFdone = [t.replace('_NA', '') for t in VCFdone_temp]

DNA_VCF = dict()
if PAIRRD:
    for sample in PAIRRD.keys():
        VCFdone += ['{rdir}/{rna}/Results/{rna}.{runid}.merge_vcf_2runs.done'.format(rdir = RESULT_DIR, rna = sample, runid = config['RUNID'])]
        DNA_VCF[sample] = ['{rdir}/{dna}/Logs_Intermediates/SmallVariantFilter/{dna}/{dna}_SmallVariants.genome.vcf'.format(rdir = RESULT_DIR, dna = PAIRRD[sample]['Matched_DNA'])]

if UNPAIRR:
    VCFdone += ['{rdir}/{rna}/Results/{rna}.{runid}.merge_vcf_RNAonly.done'.format(rdir = RESULT_DIR, rna = sample, runid = config['RUNID']) for sample in UNPAIRR]

RNA_SP = dict()
RNA_FS = dict()
if PAIRDR2:
    for sample in PAIRDR2.keys():
        VCFdone += ['{rdir}/{dna}/Results/{dna}.{runid}.merge_vcf_2runs.done'.format(rdir = RESULT_DIR, dna = sample, runid = config['RUNID'])]
        spath = '{rdir}/{rna}/TruSightTumor170_Analysis_*/RNA_{rna}/{rna}_SpliceVariants.vcf'.format(rdir = RESULT_DIR, rna = PAIRDR2[sample]['Matched_RNA'])
        fuspath = '{rdir}/{rna}/TruSightTumor170_Analysis_*/RNA_{rna}/{rna}_Fusions.csv'.format(rdir = RESULT_DIR, rna = PAIRDR2[sample]['Matched_RNA'])
        RNA_SP[sample] = glob.glob(spath)
        RNA_FS[sample] = glob.glob(fuspath)

#pp(DATA)
#KEY=(DATA.keys())
#pp(key)
runid = config['RUNID']
MAIL=config['mail']
TSO500=config['TSO500_version']
TSO170=config['TSO170_version']
DEMUX_STATS= expand(DEMUX_DIR + '/TSO500_Demux/{runid}_{data}/Reports/{runid}_{data}.html',runid = config['RUNID'],data = DATA.keys())
#DNA_MERGE_QC=expand(RESULT_DIR + '/run_qc/DNA_QC_{runid}.xlsx',runid = config['RUNID'])
#RNA_MERGE_QC=expand(RESULT_DIR + '/run_qc/RNA_QC_{runid}.xlsx',runid = config['RUNID'])
TMB_MSI_MERGE=expand(RESULT_DIR + '/run_qc/TMB_MSI_{runid}.xlsx',runid = config['RUNID'])
QC_STAT=touch(expand(RESULT_DIR + '/run_qc/{data}_QC_{runid}.xlsx',runid = config['RUNID'],data = DATA.keys()))

runqc = RESULT_DIR + '/run_qc'

if not os.path.exists(runqc):
 os.makedirs(runqc)

DNA_QC_PATH= expand(RESULT_DIR + '/{sample}/Results/MetricsReport.tsv',sample=DATA['DNA'])
#pp(DNA_QC_PATH)


#RNA_QC_PATH=expand(RESULT_DIR + '/{sample}/TruSightTumor170_Analysis*/RNA_SampleMetricsReport.txt',sample=DATA['RNA'])
#pp(RNA_QC_PATH)


if 'RNA' in DATA:
 RNA_QC_PATH =expand(RESULT_DIR + '/{sample}/TruSightTumor170_Analysis*/RNA_SampleMetricsReport.txt',sample=DATA['RNA'])
 RNA_summary =expand(RESULT_DIR + '/{sample}/TruSightTumor170_Analysis*/Summary.tsv',sample=DATA['RNA'])
else:
 RNA_QC_PATH = "None"
 RNA_summary = "None"

TMB_MSI= expand(RESULT_DIR + '/{sample}/Results/{sample}_BiomarkerReport.txt',sample=DATA['DNA'])
#pp(TMB_MSI)

samples = []
for sample in DATA.values():
    for s in sample:
        samples.append(s)
#pp(samples)

onstart:
    print('Started workflow')
    shell("echo 'TSO500 pipeline {VERSION} started  on run: {runid}' | mutt -s 'TSO500 Pipeline: {runid}'  {MAIL} ")

onsuccess:
    shell(" {PIPELINE}/scripts/Fusion_summary.py Fusion_error_{runid}.txt {RESULT_DIR}/run_qc/ {RNA_summary} ")
    shell(" {PIPELINE}/scripts/RNA_qc.py RNA_QC_{runid}.xlsx {RESULT_DIR}/run_qc/ {RNA_QC_PATH} ")  
    shell(" {PIPELINE}/scripts/DNA_qc.py DNA_QC_{runid}.xlsx {RESULT_DIR}/run_qc/ {DNA_QC_PATH} ")
    shell(" {PIPELINE}/scripts/TMB_MSI.py TMB_MSI_{runid}.xlsx {RESULT_DIR}/run_qc/ {TMB_MSI} ")
    print('Workflow finished, no error')
    for sample in samples:
        shell("find {RESULT_DIR}/{sample}/ -group $USER -exec chgrp -f {GROUP} {{}} \;")
        shell("find {RESULT_DIR}/{sample}/ \( -type f -user $USER -exec chmod g+rw {{}} \; \) , \( -type d -user $USER -exec chmod g+rwx {{}} \; \)")
    shell("find {RESULT_DIR}/run_qc -group $USER -exec chgrp -f {GROUP} {{}} \;")
    shell("find {RESULT_DIR}/run_qc \( -type f -user $USER -exec chmod g+rw {{}} \; \) , \( -type d -user $USER -exec chmod g+rwx {{}} \; \)")
    shell("find {RUN_DIR}/{runid} -maxdepth 1 -type f -group $USER -name  \"SampleSheet_*.csv\" -exec chgrp -f {GROUP} {{}} \;")
    shell("find {RUN_DIR}/{runid} -maxdepth 1 -type f -user $USER -name  \"SampleSheet_*.csv\" -exec chmod g+rw {{}} \;")
    shell("find {DEMUX_DIR}/TSO500_Demux/{runid}_* -group $USER -exec chgrp -f {GROUP} {{}} \;")
    shell("find {DEMUX_DIR}/TSO500_Demux/{runid}_* \( -type f -user $USER -exec chmod g+rw {{}} \; \) , \( -type d -user $USER -exec chmod g+rwx {{}} \; \)")
    shell("if [ -f {RESULT_DIR}/run_qc/Fusion_error_{runid}.txt ]; then cat {RESULT_DIR}/run_qc/Fusion_error_{runid}.txt ; else echo 'TSO500 pipeline {VERSION} with TSO500_app {TSO500} ,TSO170_app {TSO170} completed successfully  on run' ; fi | mutt -s 'TSO500 Pipeline: {runid}' -a {DEMUX_STATS} {QC_STAT} {TMB_MSI_MERGE} {MAIL} ")
    shell("find .snakemake/ logs {runid}.yaml -group $USER -exec chgrp -f {GROUP} {{}} \;")
    shell("find .snakemake/ logs {runid}.yaml \( -type f -user $USER -exec chmod g+rw {{}} \; \) , \( -type d -user $USER -exec chmod g+rwx {{}} \; \)")

onerror:
    print('An error occured')
    for sample in samples:
        shell("find {RESULT_DIR}/{sample}/ -group $USER -exec chgrp -f {GROUP} {{}} \;")
        shell("find {RESULT_DIR}/{sample}/ \( -type f -user $USER -exec chmod g+rw {{}} \; \) , \( -type d -user $USER -exec chmod g+rwx {{}} \; \)")
    shell("find {RESULT_DIR}/run_qc -group $USER -exec chgrp -f {GROUP} {{}} \;")
    shell("find {RESULT_DIR}/run_qc \( -type f -user $USER -exec chmod g+rw {{}} \; \) , \( -type d -user $USER -exec chmod g+rwx {{}} \; \)")
    shell("find {RUN_DIR}/{runid} -maxdepth 1 -type f -group $USER -name  \"SampleSheet_*.csv\" -exec chgrp -f {GROUP} {{}} \;")
    shell("find {RUN_DIR}/{runid} -maxdepth 1 -type f -user $USER -name  \"SampleSheet_*.csv\" -exec chmod g+rw {{}} \;")
    shell("find {DEMUX_DIR}/TSO500_Demux/{runid}_* -group $USER -exec chgrp -f {GROUP} {{}} \;")
    shell("find {DEMUX_DIR}/TSO500_Demux/{runid}_* \( -type f -user $USER -exec chmod g+rw {{}} \; \) , \( -type d -user $USER -exec chmod g+rwx {{}} \; \)")
    shell("echo 'TSO500 pipeline {VERSION} error occurred  on run: {runid}' | mutt  -s 'TSO500 Pipeline: {runid}' {MAIL} ")
    shell("find .snakemake/ logs {runid}.yaml \( -type f -user $USER -exec chmod g+rw {{}} \; \) , \( -type d -user $USER -exec chmod g+rwx {{}} \; \)")
    shell("find .snakemake/ logs {runid}.yaml -group $USER -exec chgrp -f {GROUP} {{}} \;")


rule all:
    input: 
        expand(DEMUX_DIR + '/TSO500_Demux/{runid}_{data}/Reports/{runid}_{data}.{ext}', runid = config['RUNID'], data = DATA.keys(), ext = ['html', 'csv']),
        expand(DEMUX_DIR + '/TSO500_Demux/{runid}_DNA/{sample}/SampleSheet.csv', runid = config['RUNID'], sample = DATA['DNA']),
        expand(DEMUX_DIR + '/TSO500_Demux/{runid}_RNA/{sample}/SampleSheet.csv', runid = config['RUNID'], sample = DATA['RNA']),
        expand(DEMUX_DIR + '/TSO500_Demux/{runid}_DNA/{sample}.cleanup.done', runid = config['RUNID'], sample = DATA['DNA']),
        expand(DEMUX_DIR + '/TSO500_Demux/{runid}_RNA/{sample}.cleanup.done', runid = config['RUNID'], sample = DATA['RNA']),
        expand(RESULT_DIR + '/{sample}/{runid}_{sample}_DNA.done',runid = config['RUNID'],sample = DATA['DNA']),
        expand(RESULT_DIR + '/{sample}/{runid}_{sample}_RNA.done',runid = config['RUNID'],sample = DATA['RNA']),
        expand(RESULT_DIR + '/{sample}/Results/{sample}_{runid}.failGenes',runid = config['RUNID'],sample = DATA['DNA']),
        expand(RESULT_DIR + '/{sample}/Results/{sample}_{runid}.hotspot.depth',runid = config['RUNID'],sample = DATA['DNA']),
        VCFdone

rule sampleSheet:
    input:
        sampleSheet = RUN_DIR + '/{runid}/SampleSheet.csv'

    output:
        RUN_DIR + '/{runid}/SampleSheet_{data,\w{3}}.csv'
    params:
        rulename = 'all.{runid}.{data}',
        resources = config['default']

    shell:
        '''
        module load python/3.6
        {PIPELINE}/scripts/sampleSheet.py {input.sampleSheet} {output[0]} {wildcards.data}
        '''

rule bcl2fastq:
    input:
        sampleSheet = RUN_DIR + '/{runid}/SampleSheet_{data}.csv'

    output:
        html1 = expand(DEMUX_DIR + '/TSO500_Demux/{{runid}}_{{data,\w\w\w}}/Reports/html/{flowcell}/all/all/all/lane.html', flowcell = flowcell)

    params:
        run = RUN_DIR + '/{runid}',
        bclout = DEMUX_DIR + '/TSO500_Demux/{runid}_{data,\w{3}}',
        rulename = 'bcl2fastq.{runid}.{data}',
        resources = config['bcl2fastq']

    shell:
        '''
        module load bcl2fastq/2.20.0
        bcl2fastq --runfolder-dir {params.run} --output-dir {params.bclout} --interop-dir {params.bclout}/InterOp --sample-sheet {input.sampleSheet} --loading-threads 16 --processing-threads 16 --writing-threads 16 
        '''

rule demuxStats:
    input:
        html1 = expand(DEMUX_DIR + '/TSO500_Demux/{{runid}}_{{data}}/Reports/html/{flowcell}/all/all/all/lane.html', flowcell = flowcell)

    output:
        html = DEMUX_DIR + '/TSO500_Demux/{runid}_{data,\w{3}}/Reports/{runid}_{data}.html',
        csv = DEMUX_DIR + '/TSO500_Demux/{runid}_{data,\w{3}}/Reports/{runid}_{data}.csv',
    params:
        rulename = 'demuxStats.{runid}.{data}',
        resources = config['demuxStats']

    shell:
        '''
        {PIPELINE}/scripts/demux_summary.pl --runid {input.html1} --html {output.html} --csv {output.csv}
        '''

rule sampleFolders:
    input:
        sampleSheet = RUN_DIR + '/{runid}/SampleSheet_{data}.csv',
        html = rules.demuxStats.output.html 

    output:
        DEMUX_DIR + '/TSO500_Demux/{runid}_{data,\w{3}}/{sample}/SampleSheet.csv',
        touch(DEMUX_DIR + '/FastqFolder/{sample}/rsync.{runid}_{data,\w{3}}.done'),                #this file is created so that the wild cards can flow to next rule
    params:
        rulename = 'sampleFolders.{runid}.{data}.{sample}',
        resources = config['sampleFolders']

    shell:
        '''
        module load python/3.6
        {PIPELINE}/scripts/sampleSheet.py {input.sampleSheet} {output[0]} {wildcards.sample}
        for fastq in `find $DEMUX_DIR/TSO500_Demux/{wildcards.runid}_{wildcards.data}/ -maxdepth 1 -name "{wildcards.sample}_*.fastq.gz"`;do ln -sf $fastq $DEMUX_DIR/TSO500_Demux/{wildcards.runid}_{wildcards.data}/{wildcards.sample}/.;done
        rsync -avzL $DEMUX_DIR/TSO500_Demux/{wildcards.runid}_{wildcards.data}/{wildcards.sample} $DEMUX_DIR/FastqFolder
        '''

rule cleanup:
    input:
        rules.sampleFolders.output[1]

    output:
        touch(DEMUX_DIR + '/TSO500_Demux/{runid}_{data,\w{3}}/{sample}.cleanup.done')

    params:
        rulename = 'cleanup.{runid}.{data}.{sample}',
        resources = config['sampleFolders']
       
    shell:
        '''
        find $DEMUX_DIR/TSO500_Demux/{wildcards.runid}_{wildcards.data}/ -maxdepth 1 -name "{wildcards.sample}_*.fastq.gz" -exec rm {{}} \;
        '''

rule launchDNA_App:
    input:
         out = DEMUX_DIR + '/FastqFolder/{sample}/rsync.{runid}_DNA.done'
    output:
#       RESULT_DIR + '/{sample}/Results/MetricsReport.tsv',
        touch(RESULT_DIR + '/{sample}/{runid}_{sample}_DNA.done'),
    params:
        rulename = 'launchDNA_App.{sample}',
        resources = config['app'],
	ref = config['TSO500_reference'],
        DNA_app = config['TSO500_app'],
        out_dir = RESULT_DIR,
	dir = DEMUX_DIR + '/FastqFolder/{sample}'
    shell:
        '''
         {PIPELINE}/scripts/TSO500_app.sh {params.dir} {params.DNA_app} {params.out_dir} {params.ref}
        '''

rule launchRNA_App:
    input:
         out = DEMUX_DIR + '/FastqFolder/{sample}/rsync.{runid}_RNA.done'
    output:
        touch(RESULT_DIR + '/{sample}/{runid}_{sample}_RNA.done')
    params:
        rulename = 'launchRNA_App.{sample}',
        resources = config['app'],
        ref = config['TSO170_reference'],
        RNA_app = config['TSO170_app'],
        out_dir = RESULT_DIR + '/{sample}',
        dir = DEMUX_DIR + '/FastqFolder/{sample}',
        fusion = RESULT_DIR + '/{sample}/TruSightTumor170_Analysis_*/RNA_{sample}/{sample}_Fusions.csv'
    shell:
        '''
         if [ -d "{params.out_dir}" ]; then rm -Rf {params.out_dir}; fi
         {PIPELINE}/scripts/TSO170_app.sh -fastq {params.dir} {params.ref} {params.out_dir}
         {PIPELINE}/scripts/check_fusion.py {params.fusion}         
        '''

rule DNA_QC:
    input:
         out = RESULT_DIR + '/{sample}/{runid}_{sample}_DNA.done'
    output:
        gene = RESULT_DIR + '/{sample}/Results/{sample}_{runid}.failGenes',
	hotspot_depth = RESULT_DIR + '/{sample}/Results/{sample}_{runid}.hotspot.depth',
#        touch(RESULT_DIR + '/run_qc/{runid}.done')        
    params:
        rulename = 'DNA_QC.{sample}',
        bam = RESULT_DIR + '/{sample}/Logs_Intermediates/StitchedReads/{sample}/{sample}.stitched.bam',
        resources = config['DNA_QC'],
        bed = config['TSO500_reference'] + '/TST500C_manifest.bed',
	dir = RESULT_DIR + '/{sample}/Results',
        runid = config['RUNID'],
        script = PIPELINE + '/scripts',
        hotspot = config['hotspot_bed'],
        size = config['genome_size'],
    shell:
        '''

	{PIPELINE}/scripts/TSO500_QC.sh {params.dir} {params.bam} {params.bed}  {params.script} {params.runid} {params.hotspot} {params.size} 
        '''

rule merge_VCF_paired:
    input:
        RESULT_DIR + '/{dna}/{runid}_{dna}_DNA.done',
        RESULT_DIR + '/{rna}/{runid}_{rna}_RNA.done'
    output:
        touch(RESULT_DIR + '/{dna}/Results/{dna}_{rna}.{runid}.merge_vcf.done')
    params:
        rulename = 'merge_VCF_paired.{dna}_{rna}',
        resources = config['merge_VCF']
    shell:
        '''
        module load python/3.7
        module load blat/3.5
        module load seqtk/1.3
        #delete any preexisting fusion_contigs.fa, etc
        cd {RESULT_DIR}/{wildcards.dna}/Results/
        if [ -f fusion_contigs.fa ]; then rm -f fusion_contigs.fa; fi
        if [ -f fusion_ref_pos.bed ]; then rm -f fusion_ref_pos.bed; fi
        if [ -f fusion_contigs.psl ]; then rm -f fusion_contigs.psl; fi    
        {PIPELINE}/scripts/merge_vcfs.py {RESULT_DIR}/{wildcards.dna}/Logs_Intermediates/SmallVariantFilter/{wildcards.dna}/{wildcards.dna}_SmallVariants.genome.vcf {RESULT_DIR}/{wildcards.rna}/TruSightTumor170_Analysis_*/RNA_{wildcards.rna}/{wildcards.rna}_SpliceVariants.vcf {RESULT_DIR}/{wildcards.rna}/TruSightTumor170_Analysis_*/RNA_{wildcards.rna}/{wildcards.rna}_Fusions.csv {RESULT_DIR}/{wildcards.dna}/Results/{wildcards.dna}_{wildcards.rna}.merged.vcf
        gzip -f {RESULT_DIR}/{wildcards.dna}/Results/{wildcards.dna}_{wildcards.rna}.merged.vcf
        ''' 

rule merge_VCF_dna:
    input:
        rules.launchDNA_App.output
        #RESULT_DIR + '/{sample}/{runid}_{sample}_DNA.done'
    output:
        touch(RESULT_DIR + '/{sample}/Results/{sample}.{runid}.merge_vcf.done')
    params:
        rulename = 'merge_VCF_dna.{sample}',
        resources = config['merge_VCF']
    shell:
        '''
        module load python/3.7
        module load blat/3.5
        module load seqtk/1.3
        {PIPELINE}/scripts/merge_vcfs.py {RESULT_DIR}/{wildcards.sample}/Logs_Intermediates/SmallVariantFilter/{wildcards.sample}/{wildcards.sample}_SmallVariants.genome.vcf NULL NULL {RESULT_DIR}/{wildcards.sample}/Results/{wildcards.sample}.merged.vcf
        gzip -f {RESULT_DIR}/{wildcards.sample}/Results/{wildcards.sample}.merged.vcf
        '''

rule merge_VCF_paired_separate_runs_older_dna:
    input:
        #rna_done = RESULT_DIR + '/{rna}/{runid}_{rna}_RNA.done',
        rna_done = rules.launchRNA_App.output,
        dna_vcf = lambda wildcards: DNA_VCF[wildcards.sample]
    output:
        touch(RESULT_DIR + '/{sample}/Results/{sample}.{runid}.merge_vcf_2runs.done')
    params:
        rulename = 'merge_VCF_paired_separate_runs_older_dna.{sample}',
        resources = config['merge_VCF'],
        dna = lambda wildcards: PAIRRD[wildcards.sample]['Matched_DNA']
    shell:
        '''
        module load python/3.7
        module load blat/3.5
        module load seqtk/1.3
        #delete any preexisting fusion_contigs.fa, etc
        cd {RESULT_DIR}/{params.dna}/Results/
        if [ -f fusion_contigs.fa ]; then rm -f fusion_contigs.fa; fi
        if [ -f fusion_ref_pos.bed ]; then rm -f fusion_ref_pos.bed; fi
        if [ -f fusion_contigs.psl ]; then rm -f fusion_contigs.psl; fi    
        {PIPELINE}/scripts/merge_vcfs.py {input.dna_vcf} {RESULT_DIR}/{wildcards.sample}/TruSightTumor170_Analysis_*/RNA_{wildcards.sample}/{wildcards.sample}_SpliceVariants.vcf {RESULT_DIR}/{wildcards.sample}/TruSightTumor170_Analysis_*/RNA_{wildcards.sample}/{wildcards.sample}_Fusions.csv {RESULT_DIR}/{params.dna}/Results/{params.dna}_{wildcards.sample}.merged.vcf
        gzip -f {RESULT_DIR}/{params.dna}/Results/{params.dna}_{wildcards.sample}.merged.vcf
        touch {RESULT_DIR}/{params.dna}/Results/{params.dna}_{wildcards.sample}.merge_vcf.done
        '''

rule merge_VCF_rna:
    input:
        rules.launchRNA_App.output
        #RESULT_DIR + '/{sample}/{runid}_{sample}_RNA.done'
    output:
        touch(RESULT_DIR + '/{sample}/Results/{sample}.{runid}.merge_vcf_RNAonly.done')
    params:
        rulename = 'merge_VCF_rna.{sample}',
        resources = config['merge_VCF']
    shell:
        '''
        module load python/3.7
        module load blat/3.5
        module load seqtk/1.3
        #delete any preexisting fusion_contigs.fa, etc
        cd {RESULT_DIR}/{wildcards.sample}/Results/
        if [ -f fusion_contigs.fa ]; then rm -f fusion_contigs.fa; fi
        if [ -f fusion_ref_pos.bed ]; then rm -f fusion_ref_pos.bed; fi
        if [ -f fusion_contigs.psl ]; then rm -f fusion_contigs.psl; fi    
        {PIPELINE}/scripts/merge_vcfs.py NULL {RESULT_DIR}/{wildcards.sample}/TruSightTumor170_Analysis_*/RNA_{wildcards.sample}/{wildcards.sample}_SpliceVariants.vcf {RESULT_DIR}/{wildcards.sample}/TruSightTumor170_Analysis_*/RNA_{wildcards.sample}/{wildcards.sample}_Fusions.csv {RESULT_DIR}/{wildcards.sample}/Results/{wildcards.sample}.merged.vcf
        gzip -f {RESULT_DIR}/{wildcards.sample}/Results/{wildcards.sample}.merged.vcf
        '''

rule merge_VCF_paired_separate_runs_older_rna:
    input:
        dna_done = rules.launchDNA_App.output,
        rna_sp = lambda wildcards: RNA_SP[wildcards.sample],
        rna_fs = lambda wildcards: RNA_FS[wildcards.sample]
    output:
        touch(RESULT_DIR + '/{sample}/Results/{sample}.{runid}.merge_vcf_2runs.done')
    params:
        rulename = 'merge_VCF_paired_separate_runs_older_rna.{sample}',
        resources = config['merge_VCF'],
        rna = lambda wildcards: PAIRDR2[wildcards.sample]['Matched_RNA']
    shell:
        '''
        module load python/3.7
        module load blat/3.5
        module load seqtk/1.3
        #delete any preexisting fusion_contigs.fa, etc
        cd {RESULT_DIR}/{wildcards.sample}/Results/
        if [ -f fusion_contigs.fa ]; then rm -f fusion_contigs.fa; fi
        if [ -f fusion_ref_pos.bed ]; then rm -f fusion_ref_pos.bed; fi
        if [ -f fusion_contigs.psl ]; then rm -f fusion_contigs.psl; fi    
        {PIPELINE}/scripts/merge_vcfs.py {RESULT_DIR}/{wildcards.sample}/Logs_Intermediates/SmallVariantFilter/{wildcards.sample}/{wildcards.sample}_SmallVariants.genome.vcf {input.rna_sp} {input.rna_fs} {RESULT_DIR}/{wildcards.sample}/Results/{wildcards.sample}_{params.rna}.merged.vcf
        gzip -f {RESULT_DIR}/{wildcards.sample}/Results/{wildcards.sample}_{params.rna}.merged.vcf
        touch {RESULT_DIR}/{wildcards.sample}/Results/{wildcards.sample}_{params.rna}.merge_vcf.done
        '''

