#!/bin/bash


if [[ -f "$3" ]]

then
     python /data/Compass/Tools/TSO500_pipeline_vg/scripts/RNA_qc.py $1 $2 $3
else
     echo "no files"
fi

