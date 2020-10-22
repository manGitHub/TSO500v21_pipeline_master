#! /bin/bash

function usage() {
    echo "USAGE: $0 --rundir [optional: /path/to/run/dir] --runid [RUNID] --demuxdir [optional: /path/to/demultiplexing/dir] --resultsdir [optional: /path/to/App/results/dir] --pipeline [optional: /path/to/pipeline/dir] --samplesheet [optional: /path/to/custom/samplesheet] --dryrun "
}
function fail() {
    echo "$@"
    exit 1
}

if [ $# -eq 0 ]; then
    usage
    exit 0
fi

RUN_DIR=/data/Compass/NextSeq_raw
DEMUX_DIR=/data/Compass/DATA/NextSeq
RESULT_DIR=/data/Compass/Analysis/ProcessedResults_NexSeq/TSO500_Results
PIPELINE_HOME=$(dirname "$0")  #/data/Compass/Tools/TSO500_pipeline
APP_DIR=/data/Compass/Tools
#SAMPLESHEET=$RUN_DIR/$runid/SampleSheet.csv

while [ "$1" != "" ]; do
    case $1 in
        --rundir )	shift
			RUN_DIR=$1
			;;
        --runid )	shift
			runid=$1
                        SAMPLESHEET=$RUN_DIR/$runid/SampleSheet.csv
			;;
        --dryrun )	dryrun=1
			;;
        --demuxdir )	shift
			DEMUX_DIR=$1
			;;
        --resultsdir )	shift
			RESULT_DIR=$1
			;;
        --pipeline )	shift
			PIPELINE_HOME=$1
 			;;
        --samplesheet ) shift
                        SAMPLESHEET=$1
                        ;;
        -h | --help )	usage
			exit
    esac
    shift
done

export RUN_DIR
export DEMUX_DIR
export RESULT_DIR
export PIPELINE_HOME
export APP_DIR
export DATE=`date +'%m%d%Y_%H%M%S'`
export SAMPLESHEET

echo SAMPLESHEET: $SAMPLESHEET
echo RUN_DIR: $RUN_DIR
echo DEMUX_DIR: $DEMUX_DIR
echo RESULT_DIR: $RESULT_DIR
echo PIPELINE_HOME: $PIPELINE_HOME
echo APP_DIR: $APP_DIR
echo DATETIME: $DATE

YAML=${runid}.yaml
echo $YAML
echo "RUNID:" > $YAML
echo "    '$runid'" >> $YAML
cat $YAML

export YAML 

if [ ! -d "$RUN_DIR/$runid" ]; then
    fail "Could not find RUN $runid at the input path $RUN_DIR/"
fi

module load snakemake/5.4.0 &> /dev/null || fail "Could not load module snakemake/5.4.0"

if [ "$dryrun" == '1' ];then
    #Dryrun
    echo "Dryrun"
    snakemake -nrp --nolock -k -j 3000 -s $PIPELINE_HOME/TSO500.snakefile --configfile $YAML -d `pwd` 
else
    echo "Executing TSO500 pipeline on RUN: $runid"
    $PIPELINE_HOME/scripts/analysis_samplesheet.py $SAMPLESHEET $RUN_DIR/$runid/Analysis_SampleSheet.csv
    if [ ! -d "logs" ]; then
        mkdir logs
        chgrp -f Compass logs
        chmod g+rwx logs
    fi
    sbatch -e pipeline.%j.%x.e -o pipeline.%j.%x.o --job-name=TSO500.$runid.$DATE --mem=1G --partition=ccr,norm --time=20:00:00 --cpus-per-task=1 $PIPELINE_HOME/submit.sh	
fi


