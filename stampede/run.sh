#!/bin/bash

#SBATCH -J lc
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p normal
#SBATCH -t 24:00:00
#SBATCH -A iPlant-Collabs

module load tacc-singularity
module load launcher

set -u

function ADVANCED_USAGE() {
	echo "Usage: pipeline.py [-h] -i INPUT_DIR -w WORK_DIR"
	echo "		-b BARCODE_LENGTH -m MAPPING_FILE [-p PAIRED_FILE]"
	echo
	echo "Required arguments:"
    echo "  -i INPUT_DIR            path to the input directory"
    echo "  -w WORK_DIR         path to the output directory" 
    echo "  -b BARCODE_LENGTH       length of the barcodes (int)"
    echo "  -m MAPPING_FILE         path to file containing information about the barcodes. Must be in the QIIME mapping file format"
	echo
    echo "Optional arguments:"
	echo "	-h				show this help message and exit"
	echo "	-p PAIRED_FILE			path to the input file's pair"
	echo
	exit 1
}


function USAGE() {
    echo "Usage: pipeline.py [-h] -i INPUT_DIR -w WORK_DIR"
    echo "		-b BARCODE_LENGTH -m MAPPING_FILE [-p PAIRED_FILE]"
	echo "Required arguments:"
    echo "	-b BARCODE_LENGTH"
    echo "	-m MAPPING_FILE"
    echo "  -i INPUT_DIR"
    echo "  -w WORK_DIR"
    echo
    echo "Options:"
    echo "	-h"
	echo "	-p PAIRED_FILE"
	echo
    exit 1
}
INPUT_DIR=""
WORK_DIR=""
PAIRED_ENDS=""
BARCODE_LENGTH=0
MAPPING_FILE=""
IMG="demultiplexer.img"


[[ $# -eq 0 ]] && USAGE 1

while getopts :i:w:b:m:p:h OPT; do
  case $OPT in
    h)
      ADVANCED_USAGE
      ;;
    i)
      INPUT_DIR="$OPTARG"
      ;;
    w)
      WORK_DIR="$OPTARG"
	  ;;
	b)
	  BARCODE_LENGTH=$OPTARG
	  ;;
	m)
	  MAPPING_FILE="$OPTARG"
	  ;;
	p)
	  PAIRED_ENDS="$OPTARG"
	  ;;
    :)
      echo "Error: Option -$OPTARG requires an argument."
      exit 1
      ;;
    \?)
      echo "Error: Invalid option: -${OPTARG:-""}"
      exit 1
  esac
done

if [[ $BARCODE_LENGTH -eq 0 ]]; then
	echo "BARCODE_LENGTH is required"
	exit 1
fi

if [[ $MAPPING_FILE -eq "" ]]; then
	echo "MAPPING_FILE is required"
	exit 1
fi

if [[ $WORK_DIR -eq "" ]]; then
    echo "WORK_DIR is required"
    exit 1
fi

if [[ $INPUT_DIR -eq "" ]]; then
    echo "INPUT_DIR is required"
    exit 1
fi

if [[ ! -f "$IMG" ]]; then
    echo "Missing IMG \"$IMG\""
    exit 1
fi

PARAMRUN="$TACC_LAUNCHER_DIR/paramrun"
export LAUNCHER_PLUGIN_DIR="$TACC_LAUNCHER_DIR/plugins"
export LAUNCHER_WORKDIR="$PWD"
export LAUNCHER_RMI="SLURM"
export LAUNCHER_SCHED="interleaved"

#
# Detect if INPUT is a regular file or directory, expand to list of files
#
INPUT_FILES=$(mktemp)
if [[ -f "$INPUT_DIR" ]]; then
    echo "$INPUT_DIR" > "$INPUT_FILES"
elif [[ -d "$INPUT_DIR" ]]; then
    find "$INPUT_DIR" -type f > "$INPUT_FILES"
else
    echo "-i \"$INPUT_DIR\" is neither file nor directory"
    exit 1
fi

NUM_INPUT=$(wc -l "$INPUT_FILES" | awk '{print $1}')
if [[ $NUM_INPUT -lt 1 ]]; then
    echo "There are no files to process."
    exit 1
fi

echo "I will process NUM_INPUT \"$NUM_INPUT\" files"
cat -n "$INPUT_FILES"

#
# Here is how to use LAUNCHER for parallelization
#
PARAM="$$.param"
cat /dev/null > "$PARAM"
while read -r FILE; do
    if [[ $PAIRED_ENDS -eq "" ]]; then
        echo "singularity run $IMG -i $FILE -w $WORK_DIR -m $MAPPING_FILE -b $BARCODE_LENGTH" >> "$PARAM"
    else
        echo "singularity run $IMG -i $FILE -w $WORK_DIR -m $MAPPING_FILE -b $BARCODE_LENGTH -p $PAIRED_ENDS" >> "$PARAM"
    fi
done < "$INPUT_FILES"

echo "Starting Launcher"
if [[ $NUM_INPUT -lt 16 ]]; then
    LAUNCHER_PPN=$NUM_INPUT
else
    LAUNCHER_PPN=16
fi

LAUNCHER_JOB_FILE="$PARAM"
export LAUNCHER_JOB_FILE
export LAUNCHER_PPN
$PARAMRUN
echo "Ended LAUNCHER"
unset LAUNCHER_PPN

rm "$PARAM"

echo "Done."

