from pprint import pprint as pp
from collections import defaultdict
import glob
import os
import pandas as pd
import numpy as np
from collections import OrderedDict
import re
import sys

PIPELINE = os.environ['PIPELINE_HOME']
RUN_DIR = os.environ['RUN_DIR']
DEMUX_DIR = os.environ['DEMUX_DIR']
configfile: PIPELINE + '/config/cluster.json'
GROUP = config['group']
VERSION = config['version']
APP_DIR = os.environ['APP_DIR']
RESULT_DIR = os.environ['RESULT_DIR']
SAMPLESHEET = os.environ['SAMPLESHEET']
F = config['RUNID'].split('_')
#pp(F)
flowcell = ''
for f in F:
    if re.match('^\d{4}$', f):
        flowcell = re.sub('^\w', '', F[F.index(f) + 1])

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


analysis = RUN_DIR + '/' + config['RUNID'] + '/Analysis_SampleSheet.csv'
demux_done = DEMUX_DIR + '/TSO500_Demux/' + config['RUNID'] + '/' + config['RUNID'] + '_demux.done'


runid = config['RUNID']
MAIL=config['mail']

runqc = RESULT_DIR + '/run_qc'

if not os.path.exists(runqc):
 os.makedirs(runqc)


df = pd.read_csv(analysis,delimiter=',',skiprows=24)
df = pd.DataFrame(df)
df_RNA = df.loc[df['Sample_Type'] == "RNA"]
df_RNA = df_RNA[['Sample_Type','Pair_ID']]
RNA_list = df_RNA['Pair_ID'].tolist()
#pp(RNA_list)
df_DNA = df.loc[df['Sample_Type'] == "DNA"]
DNA_list = df_DNA['Pair_ID'].tolist()
#pp(DNA_list)
df = df.drop_duplicates(subset=['Pair_ID'])
df = df[['Sample_Type','Pair_ID']]
pairs= df['Pair_ID'].tolist()
pp(pairs)


#pp(RNA_list)
#pp(DNA_list)

DATA = dict( [(key,value) for key,value in DATA.items() if len(value)>0])

if not DNA_list:
 pp("There are no DNA samples in this run")
 MSI = "None"
 TMB = "None"
 MSI_TMB = "None"
else:
 MSI =expand(RESULT_DIR + '/{pair}/Logs_Intermediates/Msi/*/*.msi.json',pair=DNA_list)
 TMB =expand(RESULT_DIR + '/{pair}/Logs_Intermediates/Tmb/*/*.tmb.json',pair=DNA_list)
 MSI_TMB =expand(RESULT_DIR + '/run_qc/MSI_TMB_{runid}.xlsx',runid = config['RUNID'])


if not RNA_list:
 RNASeQC = "None"
 RNASeQC_summary = "None"
else:
 RNASeQC =expand(RESULT_DIR + '/{pair}/Results/rnaseqc/metrics.tsv',pair=RNA_list)
 RNASeQC_summary =expand(RESULT_DIR + '/run_qc/RNASeQC_{runid}.xlsx',runid = config['RUNID'])

DATA = dict( [(key,value) for key,value in DATA.items() if len(value)>0])

TSO_WORKFLOW=config['workflow_version']
TSO500_QC =expand(RESULT_DIR + '/{pair}/Results/MetricsOutput.tsv',pair=pairs)
TSO500_QC_Summary =expand(RESULT_DIR + '/run_qc/TSO500_QC_{runid}.xlsx',runid = config['RUNID'])
DEMUX_stats = expand(DEMUX_DIR + '/TSO500_Demux/{runid}/Logs_Intermediates/FastqGeneration/{data}_Reports/{runid}_{data}.csv',runid = config['RUNID'],data = DATA.keys())

TSO_Error = expand(RESULT_DIR + '/{pair}_for_delete',pair=pairs)

onstart:
    print('Started workflow')
    shell("echo 'TSO500 pipeline {VERSION} started  on run: {runid}' | mutt -s 'TSO500 Pipeline: {runid}'  {MAIL} ")

onsuccess:
    shell(" {PIPELINE}/scripts/merge_rnaseqc.py RNASeQC_{runid}.xlsx {RESULT_DIR}/run_qc/ {RNASeQC} ")
    shell(" {PIPELINE}/scripts/merge_TSO500_QC.py TSO500_QC_{runid}.xlsx {RESULT_DIR}/run_qc/ {PIPELINE}/scripts  {TSO500_QC} ")
    shell(" {PIPELINE}/scripts/MSI.py {RESULT_DIR}/run_qc/MSI_TMB_{runid}.xlsx {MSI} ")
    shell(" {PIPELINE}/scripts/TMB.py {RESULT_DIR}/run_qc/MSI_TMB_{runid}.xlsx {TMB} ")
    print('Workflow finished, no error')
    for pair in pairs:
        shell("find {RESULT_DIR}/{pair}/ -group $USER -exec chgrp -f {GROUP} {{}} \;")
        shell("find {RESULT_DIR}/{pair}/ \( -type f -user $USER -exec chmod g+rw {{}} \; \) , \( -type d -user $USER -exec chmod g+rwx {{}} \; \)")
    shell("find {RESULT_DIR}/run_qc -group $USER -exec chgrp -f {GROUP} {{}} \;")
    shell("find {RESULT_DIR}/run_qc \( -type f -user $USER -exec chmod g+rw {{}} \; \) , \( -type d -user $USER -exec chmod g+rwx {{}} \; \)")
    shell("find {RUN_DIR}/{runid} -maxdepth 1 -type f -group $USER -name  \"Analysis_SampleSheet.csv\" -exec chgrp -f {GROUP} {{}} \;")
    shell("find {RUN_DIR}/{runid} -maxdepth 1 -type f -user $USER -name  \"Analysis_SampleSheet.csv\" -exec chmod g+rw {{}} \;")
    shell("find {DEMUX_DIR}/TSO500_Demux/{runid} -group $USER -exec chgrp -f {GROUP} {{}} \;")
    shell("find {DEMUX_DIR}/TSO500_Demux/{runid} \( -type f -user $USER -exec chmod g+rw {{}} \; \) , \( -type d -user $USER -exec chmod g+rwx {{}} \; \)")
#    shell("if [ -f {RESULT_DIR}/run_qc/Fusion_error_{runid}.txt ]; then cat {RESULT_DIR}/run_qc/Fusion_error_{runid}.txt ; else echo 'TSO500 pipeline {VERSION} with TSO500_app {TSO500} ,TSO170_app {TSO170} completed successfully  on run' ; fi | if [[ -f {RESULT_DIR}/run_qc/TMB_MSI_{runid}.xlsx && -f {RESULT_DIR}/run_qc/RNASeQC_{runid}.xlsx ]]; then mutt -s 'TSO500 Pipeline: {runid}' -a {DEMUX_STATS} {QC_STAT} {TMB_MSI_MERGE} {MAIL} {RNASeQC_summary} ; elif [[ -f {RESULT_DIR}/run_qc/TMB_MSI_{runid}.xlsx &&  ! -f {RESULT_DIR}/run_qc/RNASeQC_{runid}.xlsx ]] ; then mutt -s 'TSO500 Pipeline: {runid}' -a {DEMUX_STATS} {QC_STAT} {TMB_MSI_MERGE} {MAIL} ; else mutt -s 'TSO500 Pipeline: {runid}' -a {DEMUX_STATS} {QC_STAT} {MAIL} ; fi ") 
#    shell("if [ -f {RESULT_DIR}/run_qc/RNASeQC_{runid}.xlsx ]; then echo 'TSO500 pipeline {VERSION} completed successfully' | mutt -s 'TSO500 Pipeline: {runid}' -a  {TSO500_QC_Summary} {MAIL} {RNASeQC_summary} ; else echo 'TSO500 pipeline {VERSION} completed successfully' | mutt -s 'TSO500 Pipeline: {runid}' -a {TSO500_QC_Summary} {MAIL} ; fi ")
    shell(" if [[ -f {RESULT_DIR}/run_qc/MSI_TMB_{runid}.xlsx && -f {RESULT_DIR}/run_qc/RNASeQC_{runid}.xlsx ]]; then echo '\nTSO500 pipeline {VERSION} with workflow version: {TSO_WORKFLOW} completed successfully' | mutt -s 'TSO500 Pipeline: {runid}' -a {DEMUX_stats} {TSO500_QC_Summary} {MAIL} {RNASeQC_summary} {MSI_TMB} -i {RESULT_DIR}/run_qc/{runid}_app_complete.txt ; elif [[ -f {RESULT_DIR}/run_qc/MSI_TMB_{runid}.xlsx &&  ! -f {RESULT_DIR}/run_qc/RNASeQC_{runid}.xlsx ]] ; then echo '\nTSO500 pipeline {VERSION} with workflow version: {TSO_WORKFLOW} completed successfully' | mutt -s 'TSO500 Pipeline: {runid}' -a {DEMUX_stats} {TSO500_QC_Summary} {MAIL} {MSI_TMB} -i {RESULT_DIR}/run_qc/{runid}_app_complete.txt ; else echo '\nTSO500 pipeline {VERSION} with workflow version: {TSO_WORKFLOW} completed successfully' | mutt -s 'TSO500 Pipeline: {runid}' -a {DEMUX_stats} {TSO500_QC_Summary} {RNASeQC_summary} {MAIL} -i {RESULT_DIR}/run_qc/{runid}_app_complete.txt ; fi ")
    shell("find .snakemake/ logs {runid}.yaml -group $USER -exec chgrp -f {GROUP} {{}} \;")
    shell("find .snakemake/ logs {runid}.yaml \( -type f -user $USER -exec chmod g+rw {{}} \; \) , \( -type d -user $USER -exec chmod g+rwx {{}} \; \)")
#    shell("for i in `cat {DEMUX_DIR}/TSO500_Demux/{runid}/sample_list` ; do rsync -av {DEMUX_DIR}/TSO500_Demux/{runid}/*demux.done {DEMUX_DIR}/FastqFolder/$i ; done")

onerror:
    print('An error occured')
    for pair in pairs:
        shell("find {RESULT_DIR}/{pair}/ -group $USER -exec chgrp -f {GROUP} {{}} \;")
        shell("find {RESULT_DIR}/{pair}/ \( -type f -user $USER -exec chmod g+rw {{}} \; \) , \( -type d -user $USER -exec chmod g+rwx {{}} \; \)")
    shell("find {RESULT_DIR}/run_qc -group $USER -exec chgrp -f {GROUP} {{}} \;")
    shell("find {RESULT_DIR}/run_qc \( -type f -user $USER -exec chmod g+rw {{}} \; \) , \( -type d -user $USER -exec chmod g+rwx {{}} \; \)") 
    shell("find {RUN_DIR}/{runid} -maxdepth 1 -type f -group $USER -name  \"Analysis_SampleSheet.csv\" -exec chgrp -f {GROUP} {{}} \;")
    shell("find {RUN_DIR}/{runid} -maxdepth 1 -type f -user $USER -name  \"Analysis_SampleSheet.csv\" -exec chmod g+rw {{}} \;")
    shell("find {DEMUX_DIR}/TSO500_Demux/{runid}* -group $USER -exec chgrp -f {GROUP} {{}} \;")
    shell("find {DEMUX_DIR}/TSO500_Demux/{runid}* \( -type f -user $USER -exec chmod g+rw {{}} \; \) , \( -type d -user $USER -exec chmod g+rwx {{}} \; \)")
    shell(" if [ -f {RESULT_DIR}/run_qc/{runid}_app_complete.txt ]; then echo 'TSO500 pipeline {VERSION} error occurred  on run: {runid}' | mutt  -s 'TSO500 Pipeline: {runid}' {MAIL} -i {RESULT_DIR}/run_qc/{runid}_app_complete.txt ; else echo 'TSO500 pipeline {VERSION} error occurred  on run: {runid}' | mutt  -s 'TSO500 Pipeline: {runid}' {MAIL} ; fi ")
    shell("find .snakemake/ logs {runid}.yaml \( -type f -user $USER -exec chmod g+rw {{}} \; \) , \( -type d -user $USER -exec chmod g+rwx {{}} \; \)")
    shell("find .snakemake/ logs {runid}.yaml -group $USER -exec chgrp -f {GROUP} {{}} \;")


runid = config['RUNID']


rule all:
    input:
        expand(DEMUX_DIR + '/TSO500_Demux/{runid}/{runid}_demux.done',runid = config['RUNID']),
        expand(RESULT_DIR + '/{pair}/{runid}_{pair}.done',runid = config['RUNID'],pair=pairs),
        expand(RESULT_DIR + '/{pair}/{runid}_{pair}_DNA_RNA_QC.done',runid = config['RUNID'],pair=pairs),
        expand(DEMUX_DIR + '/TSO500_Demux/{runid}/{runid}_demuxSummary.done',runid = config['RUNID']),
        expand(RESULT_DIR + '/{pair}/{runid}_{pair}_parse_CN.done',runid = config['RUNID'],pair=pairs),
        expand(RESULT_DIR + '/{pair}/{runid}_{pair}_QCI_zip.done',runid = config['RUNID'],pair=pairs),
        expand(RESULT_DIR + '/{pair}/{runid}_{pair}_filter_splice.done',runid = config['RUNID'],pair=pairs)



rule Demultiplexing:
    input:
        run = RUN_DIR + '/{runid}/RunCompletionStatus.xml'
    output:
        demux_done = DEMUX_DIR + '/TSO500_Demux/{runid}/{runid}_demux.done'
    params:
        rulename = 'Demultiplexing',
        resource = config['TSO500_reference'],
        TSO = config['TSO500_script'],
        resources = config['app']
    shell:
        '''
        module load singularity
        cd /data/Compass/Tools/TSO500_App/TSO500_RUO_2.1.0/TSO500_RUO_LocalApp
        if [ -d "{DEMUX_DIR}/TSO500_Demux/{runid}" ]; then rm -rf {DEMUX_DIR}/TSO500_Demux/{runid} ; fi
        mkdir {DEMUX_DIR}/TSO500_Demux/{runid}
        /data/Compass/Tools/TSO500_App/TSO500_RUO_2.1.0/TSO500_RUO_LocalApp/TruSight_Oncology_500_RUO.sh --engine singularity --analysisFolder {DEMUX_DIR}/TSO500_Demux/{runid} --resourcesFolder {params.resource} --runFolder {RUN_DIR}/{runid} --sampleSheet {SAMPLESHEET}  --demultiplexOnly 
        for fastq in `find $DEMUX_DIR/TSO500_Demux/{wildcards.runid}/Logs_Intermediates/FastqGeneration/ -name *fastq.gz`;do path=${{fastq%/*}}; echo ${{path##*/}} ;done | sort -u |sed '/Undetermined$/d' > $DEMUX_DIR/TSO500_Demux/{wildcards.runid}/sample_list
        for i in `cat $DEMUX_DIR/TSO500_Demux/{wildcards.runid}/sample_list` ; do rsync -avz $DEMUX_DIR/TSO500_Demux/{wildcards.runid}/Logs_Intermediates/FastqGeneration/$i $DEMUX_DIR/FastqFolder/ ; done 
        for i in `cat $DEMUX_DIR/TSO500_Demux/{wildcards.runid}/sample_list` ; do rm -rf $DEMUX_DIR/TSO500_Demux/{wildcards.runid}/Logs_Intermediates/FastqGeneration/$i ; done
        touch {output.demux_done}
        touch {RESULT_DIR}/run_qc/{runid}_app_complete.txt
        '''

rule demuxSummary:
    input:
        demux_done = DEMUX_DIR + '/TSO500_Demux/{runid}/{runid}_demux.done',
    output:
        demuxSummary_done =  DEMUX_DIR + '/TSO500_Demux/{runid}/{runid}_demuxSummary.done',
    params:
        rulename = 'demuxSummary',
        resources = config['DNA_QC'],
        dna = DEMUX_DIR + '/TSO500_Demux/{runid}/Logs_Intermediates/FastqGeneration/DNA_Reports',
        rna = DEMUX_DIR + '/TSO500_Demux/{runid}/Logs_Intermediates/FastqGeneration/RNA_Reports',
	dna_out = DEMUX_DIR + '/TSO500_Demux/{runid}/Logs_Intermediates/FastqGeneration/DNA_Reports/{runid}_DNA.csv',
        rna_out = DEMUX_DIR + '/TSO500_Demux/{runid}/Logs_Intermediates/FastqGeneration/RNA_Reports/{runid}_RNA.csv'
    shell:
        '''
        if [ -f {params.dna}/Lane_1/Demultiplex_Stats.csv ]
        then
	python3 {PIPELINE}/scripts/demux_summary.py {params.dna} {params.dna_out} {SAMPLESHEET}
	fi

	if [ -f {params.rna}/Lane_1/Demultiplex_Stats.csv ]
        then
        python3 {PIPELINE}/scripts/demux_summary.py {params.rna} {params.rna_out} {SAMPLESHEET}
	fi
        touch {output.demuxSummary_done}
        '''


rule TSO500_app:
    input:
        {demux_done}
    output:
        appdone = RESULT_DIR + '/{pair}/{runid}_{pair}.done'
    params:
        rulename = 'TSO500_app.{pair}',
        resources = config['app'],
        resource = config['TSO500_reference'],
        fastq = DEMUX_DIR + '/FastqFolder' 
    shell:
        '''
        module load singularity
        cd /data/Compass/Tools/TSO500_App/TSO500_RUO_2.1.0/TSO500_RUO_LocalApp
        if [ ! -f {RESULT_DIR}/{wildcards.pair}/*.done ]; then  echo "processing this sample first time" ; else  mv "{RESULT_DIR}/{wildcards.pair}"  "{RESULT_DIR}/{wildcards.pair}_archive" ; fi
        
#        mkdir {RESULT_DIR}/{wildcards.pair}
         /data/Compass/Tools/TSO500_App/TSO500_RUO_2.1.0/TSO500_RUO_LocalApp/TruSight_Oncology_500_RUO.sh --engine singularity --analysisFolder {RESULT_DIR}/{wildcards.pair} --resourcesFolder {params.resource} --fastqFolder {params.fastq} --samplePairIDs {wildcards.pair} --sampleSheet {RUN_DIR}/{runid}/Analysis_SampleSheet.csv

##  Parsing the MetricsOutput.tsv file to make sure all steps in TSO app completed successfully.

        if [ -f "{RESULT_DIR}/{wildcards.pair}/Results/MetricsOutput.tsv" ]
        then
         grep -i -B1 -A1 "COMPLETED_ALL_STEPS" {RESULT_DIR}/{wildcards.pair}/Results/MetricsOutput.tsv | paste -d " "  - - - >> {RESULT_DIR}/run_qc/{runid}_app_complete.txt
         if grep -i {wildcards.pair} {RESULT_DIR}/run_qc/{runid}_app_complete.txt |grep -q TRUE 
         then
          echo "All the steps in TSO500_app completed successfully" 
         else
          echo "TSO500 app Failed "
         fi
        else
         echo "TSO500 app Failed"
        fi

        if [  -f {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}/*/*_MergedSmallVariants.genome.vcf ]; then  gzip -f {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}/*/*_MergedSmallVariants.genome.vcf ; else echo "This is a RNA only sample" ; fi 
        touch {output.appdone}
        '''    


rule DNA_RNA_QC:
    input:
        RESULT_DIR + '/{pair}/{runid}_{pair}.done'
    output:
       done = RESULT_DIR + '/{pair}/{runid}_{pair}_DNA_RNA_QC.done'
    params:
        rulename = 'DNA_RNA_QC.{pair}',
        resources = config['DNA_QC'],
        resource = config['TSO500_reference'],
        bed = config['manifest'],
        dir = RESULT_DIR + '/{pair}',
        runid = config['RUNID'],
        script = PIPELINE + '/scripts',
        hotspot = config['hotspot_bed'],
        size = config['genome_size'],
        gtf = config['gtf'],
        rna_genome = config['rna'],
    shell:
        '''
        if [ -f {RESULT_DIR}/{wildcards.pair}/Logs_Intermediates/StitchedRealigned/*/*.bam ] 
        then  
           python3 {params.script}/DNA_qc.py {params.dir}/Logs_Intermediates/StitchedRealigned/*/*.bam {params.dir} {params.bed} {params.hotspot} {params.size} {params.script}
        else
           echo "no DNA in this pair"
        fi
        if [ -f  {RESULT_DIR}/{wildcards.pair}/Logs_Intermediates/RnaAlignment/*/*.bam ]
        then
           python3 {params.script}/rnaseqc.py {params.dir}/Logs_Intermediates/RnaAlignment/*/*.bam {params.gtf} {params.rna_genome} {params.dir} {params.script}
           rm {params.dir}/Logs_Intermediates/RnaAlignment/*/*_RG.bam
           rm {params.dir}/Logs_Intermediates/RnaAlignment/*/*_dd.bam*
        else
           echo "no RNA in this pair"
        fi
        touch {output.done}
        '''  


rule parse_CN:
        input:
            RESULT_DIR + '/{pair}/{runid}_{pair}.done'
        output:
           done = RESULT_DIR + '/{pair}/{runid}_{pair}_parse_CN.done'
        params:
            rulename = 'parse_CN.{pair}',
            resources = config['DNA_QC'],
            resource = config['TSO500_reference'],
            dir = RESULT_DIR + '/{pair}',
            script = PIPELINE + '/scripts/parse_cn_vcf.py',
            cnv_genes = config['cnv_genes'],
            craft_manifest = config['craft_manifest'],
            thresholds = config['thresholds'],
            genome = config['dna_genome'],
            qci_coords = config['qci_coords']

        shell:
            '''
            if [ -f {params.dir}/Results/*/*/*_CopyNumberVariants.vcf ]
            then
               python3 {params.script} {params.dir}/Logs_Intermediates/CnvCaller/*/*_foldChange.tsv {params.cnv_genes} {params.dir}/Results/*/*/*_CopyNumberVariants.vcf {params.craft_manifest} {params.thresholds} {params.genome} {params.qci_coords}
            else
               echo "no copy number results found"
            fi
            touch {output.done}
            '''

rule Filter_Splice:
    input:
       RESULT_DIR + '/{pair}/{runid}_{pair}_parse_CN.done'
    output:
       done = RESULT_DIR + '/{pair}/{runid}_{pair}_filter_splice.done'
    params:
        rulename = 'Filter_Splice.{pair}',
        resources = config['DNA_QC'],
        dir = RESULT_DIR + '/{pair}',
        runid = config['RUNID'],
        script = PIPELINE + '/scripts'
    shell:
        '''
        #module load bedtools/2.29.2
        if [ -f {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}/*/*_SpliceVariants.vcf ]
        then
         /data/Compass/local/software/bedtools/2.29.2/bin/intersectBed -a {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}/*/*_SpliceVariants.vcf -b {params.script}/TSO500_snv.bed -wa -header | if awk -v FS='\t' -v OFS='\t' '/1/ {{$6 = "."}} 1' |grep  "PASS" > {params.dir}/Results/{wildcards.pair}/temp.vcf ; then cat {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}/*/*_SpliceVariants.vcf | grep -e  "^#" |cat - {params.dir}/Results/{wildcards.pair}/temp.vcf > {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}/{wildcards.pair}_filtered_splice.vcf ; rm {params.dir}/Results/{wildcards.pair}/temp.vcf ; else rm {params.dir}/Results/{wildcards.pair}/temp.vcf ;  fi 
        else
          echo "no RNA sample in this run"
        fi
        touch {output.done}
        '''

rule QCI_input:
    input:
       RESULT_DIR + '/{pair}/{runid}_{pair}_filter_splice.done'
    output:
       done = RESULT_DIR + '/{pair}/{runid}_{pair}_QCI_zip.done'
    params:
        rulename = 'QCI_input.{pair}',
        resources = config['DNA_QC'],
        dir = RESULT_DIR + '/{pair}',
        runid = config['RUNID'],
        script = PIPELINE + '/scripts'

    shell:
        '''

##  As we handle DNA-RNA pairs, DNAonly sample and RNAonly samples. We check for three files here and create zip files accordingly. *AllFusions.csv *_MergedSmallVariants.genome.vcf.gz *_filtered_splice.vcf.

        if [ -f {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}/*/*_AllFusions.csv ] && [ -f {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}/*/*_MergedSmallVariants.genome.vcf.gz ] && [ -f {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}/{wildcards.pair}_filtered_splice.vcf ]
        then
          zip -j  {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}.zip {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}/{wildcards.pair}_CombinedVariantOutput.tsv {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}/*/*_AllFusions.csv {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}/{wildcards.pair}_filtered_splice.vcf {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}/*/*_MergedSmallVariants.genome.vcf.gz {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}/*/*_CopyNumberVariants_92genes.vcf
        elif [ -f {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}/*/*_AllFusions.csv ] && [ -f {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}/*/*_MergedSmallVariants.genome.vcf.gz ] && [ ! -f {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}/{wildcards.pair}_filtered_splice.vcf ]
        then
          zip -j  {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}.zip {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}/{wildcards.pair}_CombinedVariantOutput.tsv {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}/*/*_AllFusions.csv  {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}/*/*_MergedSmallVariants.genome.vcf.gz {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}/*/*_CopyNumberVariants_92genes.vcf
        elif [ -f {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}/*/*_AllFusions.csv ] && [ -f {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}/{wildcards.pair}_filtered_splice.vcf ] && [ ! -f {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}/*/*_MergedSmallVariants.genome.vcf.gz ]
        then
          zip -j  {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}.zip {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}/{wildcards.pair}_CombinedVariantOutput.tsv {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}/*/*_AllFusions.csv {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}/{wildcards.pair}_filtered_splice.vcf
        elif [ -f {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}/*/*_AllFusions.csv ] && [ ! -f {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}/{wildcards.pair}_filtered_splice.vcf ] && [ ! -f {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}/*/*_MergedSmallVariants.genome.vcf.gz ]
        then
          zip -j  {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}.zip {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}/{wildcards.pair}_CombinedVariantOutput.tsv {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}/*/*_AllFusions.csv 
        elif [ -f {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}/*/*_MergedSmallVariants.genome.vcf.gz ] && [ ! -f {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}/*/*_AllFusions.csv ]
        then
          zip -j  {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}.zip {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}/{wildcards.pair}_CombinedVariantOutput.tsv {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}/*/*_MergedSmallVariants.genome.vcf.gz {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}/*/*_CopyNumberVariants_92genes.vcf
        fi
        if [ -f {RESULT_DIR}/{wildcards.pair}/Results/{wildcards.pair}.zip ] ; then  touch {output.done} ; fi
        '''
