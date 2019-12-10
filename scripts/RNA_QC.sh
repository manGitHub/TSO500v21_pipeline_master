#!/bin/bash


if [[ -f /data/Compass/dev/gangalapudiv2/NextSeq/TSO500_Demux/$3_RNA/Reports/$3_RNA.html ]]

then
     python /data/Compass/Tools/TSO500_pipeline_vg/scripts/RNA_qc.py $1 $2 $3 $4
else
     echo "no RNA"
fi

