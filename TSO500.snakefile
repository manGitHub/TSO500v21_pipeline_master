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
            DATA['RNA'].append(lines3[0])
        if L.find('PoolDNA') != -1:
            lines3 = L.split(',')
            DATA['DNA'].append(lines3[0])


#pp(DATA)
#pp(DATA.keys())

onstart:
    print('Started workflow')

onsuccess:
    print('Workflow finished, no error')
    #shell("chmod ...")

onerror:
    print('An error occured')
    #shell("mail -s "an error occurred" youremail@provider.com < {log}")


rule all:
    input: 
        expand(DEMUX_DIR + '/TSO500_Demux/{runid}_{data}/Reports/{runid}_{data}.{ext}', runid = config['RUNID'], data = DATA.keys(), ext = ['html', 'csv']),
        expand(DEMUX_DIR + '/TSO500_Demux/{runid}_DNA/{sample}/SampleSheet.csv', runid = config['RUNID'], sample = DATA['DNA']),
        expand(DEMUX_DIR + '/TSO500_Demux/{runid}_RNA/{sample}/SampleSheet.csv', runid = config['RUNID'], sample = DATA['RNA'])

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
        csv = DEMUX_DIR + '/TSO500_Demux/{runid}_{data,\w{3}}/Reports/{runid}_{data}.csv'

    params:
        rulename = 'demuxStats.{runid}.{data}',
        resources = config['demuxStats']

    shell:
        '''
        {PIPELINE}/scripts/demux_summary.pl --runid {input.html1} --html {output.html} --csv {output.csv}
        #mutt command
        '''

rule sampleFolders:
    input:
        sampleSheet = RUN_DIR + '/{runid}/SampleSheet_{data}.csv',
        html = rules.demuxStats.output.html 

    output:
        DEMUX_DIR + '/TSO500_Demux/{runid}_{data,\w{3}}/{sample}/SampleSheet.csv',
        touch(DEMUX_DIR + '/FastqFolder/{sample}/rsync.{runid}_{data,\w{3}}.done')
 
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
"""
rule launchApps:
    input:
        sampleSheet = '/data/Compass/DATA/NextSeq/{sample}'

    output:
        ' 
"""
