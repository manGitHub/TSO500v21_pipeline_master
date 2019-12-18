#! /usr/bin/env python3
from sys import argv
import re

final = []
with open(argv[1], 'rt') as IN, open(argv[2], 'wt') as OUT:
    text = IN.read()
    lines = text.split('\n')
    if argv[3] == 'RNA':
        for L in lines:
            if L.find('Read1UMI') == -1 and L.find('Read2UMI') == -1 and L.find('PoolDNA') == -1:
                if L.find('PoolRNA') != -1:
                    L = L.replace(' ', '')
                final.append(L)
    elif argv[3] == 'DNA':
        for L in lines:
            if L.find('PoolRNA') == -1:
                if L.find('PoolDNA') != -1:
                    L = L.replace(' ', '')
                final.append(L)

    else:
        index = [i for i, word in enumerate(lines) if word.startswith('Sample_ID')][0] 
        final = lines[0:index+1]
        for L in lines[index:]:
            F = L.split(',')
            m = re.match("^" + argv[3] + "$", F[0])
            if m != None:
                final.append(L)

    OUT.write('\n'.join(final))
                


#print('\n'.join(final))
