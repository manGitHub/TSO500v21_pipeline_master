#! /usr/bin/env python3
from sys import argv

final = []
with open(argv[1], 'rt') as IN, open(argv[2], 'wt') as OUT:
    text = IN.read()
    lines = text.split('\n')
    if argv[3] == 'RNA':
        for L in lines:
            if L.find('Read1UMI') == -1 and L.find('Read2UMI') == -1 and L.find('PoolDNA') == -1:
                final.append(L)
    elif argv[3] == 'DNA':
        for L in lines:
            if L.find('PoolRNA') == -1:
                final.append(L)

    else:
#        index = lines.index('Sample_ID,Sample_Name,Sample_Plate,Sample_Well,Index_ID,index,I7_Index_ID,index2,I5_Index_ID,Manifest')
        index = [i for i, word in enumerate(lines) if word.startswith('Sample_ID')][0] 
#        print(index)
        final = lines[0:index+1]
        for L in lines[index:]:
            if L.find(argv[3]) != -1:
                final.append(L)

    OUT.write('\n'.join(final))
                


#print('\n'.join(final))
