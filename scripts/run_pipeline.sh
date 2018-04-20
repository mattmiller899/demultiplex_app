#!/bin/bash

set -u

function ADVANCED_USAGE() {
	echo "Usage: pipeline.py [-h] [-i INPUT_DIR] [-w WORK_DIR] [-c CORE_COUNT]"
	echo "		-b BARCODE_LENGTH -m MAPPING_FILE [-p]"
	echo
	echo "Required arguments:"
	echo "	-b BARCODE_LENGTH		length of the barcodes (int)"
	echo "	-m MAPPING_FILE			path to file containing information about the barcodes. Must be in the QIIME mapping file format"
	echo
	echo "Optional arguments:"
	echo "	-h						show this help message and exit"
	echo "	-i INPUT_DIR			path to the input directory (default = current directory)"
	echo "	-w WORK_DIR				path to the output directory (default = current directory)"
	echo "	-c CORE_COUNT			number of cores to use (int) (default = 1)"
	echo "	-p						Indicates that the input files are paired end reads"
	echo
	exit 1
}


function USAGE() {
    echo "Usage: pipeline.py [-h] [-i INPUT_DIR] [-w WORK_DIR] [-c CORE_COUNT]"
    echo "		-b BARCODE_LENGTH -m MAPPING_FILE [-p]"
	echo "Required arguments:"
    echo "	-b BARCODE_LENGTH"
    echo "	-m MAPPING_FILE"
    echo
    echo "Options:"
    echo "	-h"
    echo "	-i INPUT_DIR"
	echo "	-w WORK_DIR"
	echo "	-c CORE_COUNT"
	echo "	-p"
	echo
    exit 1
}
echo "TEST"
INPUT_DIR="."
WORK_DIR="."
CORE_COUNT=1
PAIRED_ENDS=false
BARCODE_LENGTH=0
MAPPING_FILE=""

[[ $# -eq 0 ]] && USAGE 1

while getopts :i:w:c:b:m:ph OPT; do
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
    c)
      CORE_COUNT=$OPTARG
      ;;
	b)
	  BARCODE_LENGTH=$OPTARG
	  ;;
	m)
	  MAPPING_FILE="$OPTARG"
	  ;;
	p)
	  PAIRED_ENDS=true
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

if [[ $BARCODE_LENGTH -eq 0 ]] || [[ $MAPPING_FILE -eq "" ]]; then
	USAGE
fi

if [[ PAIRED_ENDS ]]; then
	exec /miniconda3/bin/python pipeline.py -i $INPUT_DIR -w $WORK_DIR -c $CORE_COUNT -b $BARCODE_LENGTH -m $MAPPING_FILE -p
else
	exec /miniconda3/bin/python pipeline.py -i $INPUT_DIR -w $WORK_DIR -c $CORE_COUNT -b $BARCODE_LENGTH -m $MAPPING_FILE
fi
