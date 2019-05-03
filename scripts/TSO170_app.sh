#! /bin/bash
#> USAGE
#>     170.sh [-fastq] [run_path|fastq_path] resources output
#> DESCRIPTION
#>     Wrapper for the TruSight Tumor 170 pipeline.
#> OPTIONS
#>     -fastq    Starting with a fastq folder instead of BCLs

function usage() {
    grep '^#>' $0 | sed -re 's/^#> ?//'
}

function fail() {
    echo "$@"
    exit 1
}

################################################################################
#                                    Setup                                     #
################################################################################

script="$(basename ${BASH_SOURCE[0]})"
path="$(readlink -f $(dirname ${BASH_SOURCE[0]}))"
#image="${path}/TruSight_Tumor_170.sif"
#image="${path}/TruSight_Tumor_170_Local_App_1.0.1.0.sif"
image="/data/Compass/Tools/TSO170_App_Singularity/TruSight_Tumor_170_Local_App_1.0.1.0.sif"

opts=""
if [[ "${1:-none}" == "-fastq" ]]; then
    opts="-fastq"
    shift
fi

data=${1:-none}
resources=${2:-none}
output=${3:-none}

[[ -e "${image}" ]] || fail "Could not find singularity image ${image}"

[[ "${data}" == "none" ]] && { usage; exit 1; }
[[ -d "${data}" ]] || fail "Could not find data folder ${data}"

[[ "${resources}" == "none" ]] && { usage; exit 1; }
[[ -d "${resources}" ]] || fail "Could not find resources folder ${resources}"

[[ "${output}" == "none" ]] && { usage; exit 1; }
if [[ -d "${output}" ]]; then
    fail "The output folder ${output} already exists. Exiting"
else
    mkdir "${output}"
fi

module load singularity &> /dev/null || fail "Could not load singularity module"

################################################################################
#                               Run the pipeline                               #
################################################################################
singularity run -B "${data}":/data -B "${resources}":/genomes -B "${output}":/analysis "${image}" ${opts}

