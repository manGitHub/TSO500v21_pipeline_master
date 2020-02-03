#! /usr/bin/env python3

from sys import argv
from pprint import pprint as pp
import os
import subprocess
import pandas as pd

def getTophit(hits, chrm, pos):
    if not hits.empty:
        hits['dist1'] = abs(hits['T_start'] - int(pos))
        hits['dist2'] = abs(hits['T_end'] - int(pos))
        pp(hits)
        hits_sub = hits.loc[(hits['T_name'] == chrm) & ((hits['dist1'] <= 5) | (hits['dist2'] <= 5))]   
    return hits_sub

def getStrand(pos, result, partner):
    #Handle empty data frames from no hits. 
    hits = pd.read_csv('fusion_contigs.psl', sep='\t', skiprows=5, names=('match', 'mis-match', 'rep.match', 'Ns', 'Q_gap_count', 'Q_gap_bases', 'T_gap_count', 'T_gap_bases', 'strand', 'Q_name', 'Q_size', 'Q_start', 'Q_end', 'T_name', 'T_size', 'T_start', 'T_end', 'block_count', 'blockSizes', 'qStarts', 'tStarts'))
    #pp(hits)
    hits1 = hits.loc[hits['Q_name'] == result]
    #pp(hits1)
    pos1 = pos.split(':')
    pp(pos1)
    strand = None
    tophits = getTophit(hits1, pos1[0], pos1[1])
    if not tophits.empty:
        tophit = tophits.loc[tophits['match'].idxmax()]
        pp(tophit)
        #strand = None
        if (partner == 'geneA') & (tophit['strand'] == '+'):
            if (pos1[0] == tophit['T_name']) & (abs(int(pos1[1]) - tophit['T_end']) < 5):
                strand = tophit['strand']
            else:
                print('Alignment location does not match position %s' % pos)
        elif (partner == 'geneA') & (tophit['strand'] == '-'):
            if (pos1[0] == tophit['T_name']) & (abs(int(pos1[1]) - tophit['T_start']) < 5):
                strand = tophit['strand']
            else:
                print('Top hit alignment location does not match position %s' % pos)
        elif (partner == 'geneB') & (tophit['strand'] == '+'):
            if (pos1[0] == tophit['T_name']) & (abs(int(pos1[1]) - tophit['T_start']) < 5):
                strand = tophit['strand']
        elif (partner == 'geneB') & (tophit['strand'] == '-'):
            if (pos1[0] == tophit['T_name']) & (abs(int(pos1[1]) - tophit['T_end']) < 5):
                strand = tophit['strand']
        else:
            print('Error: strand in the blat result is not + or -')  
    else:
        print("No hits returned for %s" % pos)
    return strand
    

def runBlat(fname, result):
    process = subprocess.run(['blat', '-stepSize=5', '-minScore=20', '-minIdentity=0', '/data/Compass/dev/ref/TSO170_resources/blat_3.5/twoBit/hg19.2bit', fname, result], check=True)
    if process.returncode == 0:
        print('Blat completed successfully on %s' % fname)


def revcomp(ref):
    rc = {
    'A':'T',
    'T':'A',
    'G':'C',
    'C':'G'
    }
    ref_revcom = "".join(rc[base] for base in reversed(ref))
    #return rc[ref]
    return(ref_revcom)


def writefasta(seq, name):
    with open('fusion_contigs.fa', 'a+') as Q:
        Q.write(">{}\n{}\n".format(name, seq))


def writebed(pos1):
    pos = pos1.split(':')
    with open('fusion_ref_pos.bed', 'a+') as BED:
        BED.write("{}\t{}\t{}\n".format(pos[0], int(pos[1]) - 1, pos[1]))  

def runseqtk(bed):
    process = subprocess.run(['seqtk', 'subseq', '/data/Compass/Ref/hg19/TSO170_Res/genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa', bed], check=True, capture_output=True)
    if process.returncode == 0:
        print('seqtk completed successfully on %s' % bed)
        return process.stdout

if not (argv[1] == "NULL"):
    VCF1 = open(argv[1])
    all = VCF1.read()
    VCF = all.split('\n')
    index = [i for i, content in enumerate(VCF) if content.startswith('##FORMAT')][-1]
    #pp(index)
else:
    VCF = []

if not (argv[2] == "NULL") & (argv[3] == "NULL"):
    for SP in open(argv[2]):
        if not (argv[1] == "NULL"):
            if SP.startswith('##QUAL') or SP.startswith('##ALT') or SP.startswith('##INFO') or SP.startswith('##FILTER'):
                #pp(SP)
                VCF.insert(index+1, SP.rstrip())
                index+=1
        else:
            if SP.startswith('#'):
                VCF.append(SP.rstrip())
        if not SP.startswith('#') and SP.find('PASS') != -1:
            SPfields1 = SP.split('ANT=')
            SPfields2 = SPfields1[1].split('|')
            if (SPfields2[3] == 'EGFR') or (SPfields2[3] == 'MET') or (SPfields2[3] == 'AR'):
                #pp(SP)
                SPfields = SP.split('\t')
                #SPfields[4] = SPfields[4].replace('DEL', 'CNV')
                SPfields[5] = '.'
                #SPfields[7] = SPfields[7].replace('DEL', 'CNV')
                #SPfields[8] = 'GT:CN'
                #SPfields[9] = './.:1'
                SPfin = '\t'.join(SPfields)
                VCF.append(SPfin.rstrip())
                #VCF.append(SP.rstrip())
   
    VCF = list(filter(None, VCF))
    flag = 0
    for FUS in open(argv[3]):
        if FUS.rstrip().endswith('True'):
            flag = 1
            fields = FUS.split(',')
            idf = '{}-{}'.format(fields[1], fields[2])
            fragA = fields[16][0:int(fields[17])]
            fnameA = idf + '_seqA'
            writefasta(fragA, fnameA)
            fragB = fields[16][int(fields[17]):]
            fnameB = idf + '_seqB'
            writefasta(fragB, fnameB)
            fragArc = revcomp(fragB)
            fnameArc = idf + '_rc_seqA'
            writefasta(fragArc, fnameArc)
            fragBrc = revcomp(fragA)
            fnameBrc = idf + '_rc_seqB'
            writefasta(fragBrc, fnameBrc)
            writebed(fields[3])
        
    if flag == 1:
        index1 = [i for i, content in enumerate(VCF) if content.startswith('##INFO')][-1]
        infotext = '##INFO=<ID=READ_COUNT,Number=.,Type=String,Description="Number of reads aligned to the Target">'
        VCF.insert(index1+1, infotext)
        runBlat('fusion_contigs.fa', 'fusion_contigs.psl')
        refBases = runseqtk('fusion_ref_pos.bed').decode('utf-8')
    #refBases1 = refBases.decode('utf-8')
        REF = refBases.split('\n')
        pp(REF)
 
    for FUS in open(argv[3]):
        if FUS.rstrip().endswith('True'):
            fields = FUS.split(',')
            print('Processing fusion record:')
            pp(fields)
        
            idf = '{}-{}'.format(fields[1], fields[2])
       
            fragA = fields[16][0:int(fields[17])]
            print(fragA)
            fnameA = idf + '_seqA'
            fragB = fields[16][int(fields[17]):]
            print(fragB)
            fnameB = idf + '_seqB'
            fragA_strand = getStrand(fields[3], fnameA, 'geneA')
            fragB_strand = getStrand(fields[4], fnameB, 'geneB')
            print(fragA_strand, fragB_strand)
            if (fragA_strand == None) & (fragB_strand == None):
                fragArc = revcomp(fragB)
                print(fragArc)
                fnameArc = idf + '_rc_seqA'
                fragA_strand = getStrand(fields[3], fnameArc, 'geneA')
                fragBrc = revcomp(fragA)
                print(fragBrc)
                fnameBrc = idf + '_rc_seqB'
                fragB_strand = getStrand(fields[4], fnameBrc, 'geneB')
                fragA = fragArc
                fragB = fragBrc
            print(fragA_strand, fragB_strand)

            chrom, pos = fields[3].split(':')

            indexr = [i for i, content in enumerate(REF) if content.startswith('>' + fields[3])][-1]
            #pp(indexr)
            ref = REF[indexr+1]  

            if (fragA_strand == '-') & (fragB_strand == '+'):
                alt = '[{}[{}'.format(fields[4], ref)
            elif (fragA_strand == '+') & (fragB_strand == '+'):
                alt = '{}[{}['.format(ref, fields[4])
            elif (fragA_strand == '+') & (fragB_strand == '-'):
                alt = '{}]{}]'.format(ref, fields[4])
            elif (fragA_strand == '-') & (fragB_strand == '-'):
                alt = ']{}]{}'.format(fields[4], ref)
            else:
                alt = '.'
        
            reads = int(fields[13]) + int(fields[14])

            record = "{}\t{}\t{}\t{}\t{}\t.\tPASS\tSVTYPE=Fusion;READ_COUNT={}\tGT:GQ\t./.:.".format(chrom, pos, idf, ref, alt, reads)
            print(record) 
            VCF.append(record)

with open (argv[4], 'wt') as OUT:
    for line in VCF:
        OUT.write("%s\n" % line)


