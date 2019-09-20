from pprint import pprint as pp
from collections import defaultdict
import glob
import os

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
#print('RUNID:', F)
flowcell = ''
for f in F:
    if re.match('^\d{4}$', f):
        flowcell = re.sub('^\w', '', F[F.index(f) + 1])

#print(flowcell)

samplesheet = RUN_DIR + '/' + config['RUNID'] + '/SampleSheet.csv'
#print('sample_sheet path', samplesheet)

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

DNA_QC_PATH= expand(RESULT_DIR + '/{sample}/Results/MetricsReport.tsv',sample=DATA['DNA'])
#pp(DNA_QC_PATH)


RNA_QC_PATH=expand(RESULT_DIR + '/{sample}/TruSightTumor170_Analysis*/RNA_SampleMetricsReport.txt',sample=DATA['RNA'])
#pp(RNA_QC_PATH)

TMB_MSI= expand(RESULT_DIR + '/{sample}/Results/{sample}_BiomarkerReport.txt',sample=DATA['DNA'])
#pp(TMB_MSI)


onstart:
    print('Started workflow')
    shell("echo 'TSO500 pipeline {VERSION} started  on run: {runid}' | mutt -s 'TSO500 Pipeline: {runid}'  {MAIL} ")
onsuccess:
    shell(" {PIPELINE}/scripts/RNA_QC.sh RNA_QC_{runid}.xlsx {RESULT_DIR}/run_qc/ {RNA_QC_PATH} ")  
    shell(" {PIPELINE}/scripts/DNA_qc.py DNA_QC_{runid}.xlsx {RESULT_DIR}/run_qc/ {DNA_QC_PATH} ")
    shell(" {PIPELINE}/scripts/TMB_MSI.py TMB_MSI_{runid}.xlsx {RESULT_DIR}/run_qc/ {TMB_MSI} ")
    print('Workflow finished, no error')
    shell("echo 'TSO500 pipeline {VERSION} with TSO500_app {TSO500} ,TSO170_app {TSO170} completed successfully  on run: {runid}' | mutt -s 'TSO500 Pipeline: {runid}' -a {DEMUX_STATS} {QC_STAT} {TMB_MSI_MERGE} {MAIL} ")

onerror:
    print('An error occured')
    shell("echo 'TSO500 pipeline {VERSION} error occurred  on run: {runid}' | mutt  -s 'TSO500 Pipeline: {runid}' {MAIL} ")

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
        expand(RESULT_DIR + '/{sample}/Results/{sample}_{runid}.hotspot.depth',runid = config['RUNID'],sample = DATA['DNA'])

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
#check if bcl2fastq will complain if folder exists
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
        for fastq in `find $DEMUX_DIR/TSO500_Demux/{wildcards.runid}_{wildcards.data}/ -maxdepth 1 -name "{wildcards.sample}_*.fastq.gz"`;do ln -s $fastq $DEMUX_DIR/TSO500_Demux/{wildcards.runid}_{wildcards.data}/{wildcards.sample}/.;done
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
       touch(RESULT_DIR + '/{sample}/{runid}_{sample}_DNA.done')
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

rule LaunchRNA_App:
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
        dir = DEMUX_DIR + '/FastqFolder/{sample}'
    shell:
        '''
         if [ -d "{params.out_dir}" ]; then rm -Rf {params.out_dir}; fi
         {PIPELINE}/scripts/TSO170_app.sh -fastq {params.dir} {params.ref} {params.out_dir}

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







