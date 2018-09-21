# imicrobe-demultiplexer

A demultiplexing and index-removing pipeline for single and paired-end data.

##Introduction

There are three ways to run this pipeline.

    + Clone the repository and run with Python 2.7
    + As a Singularity Container.
    + As a Cyverse/iMicrobe app.

## Clone and Run with Python 2.7
The only requirement to run the pipeline as a Python 2.7+ package is a Python 2.7+ interpreter and `Git`. It is not required but is highly recommended to install the pipeline in a virtual environment. You will need qiime installed to your Python installation. Clone the package:

```
$ git clone git@github.com:mattmiller899/demultiplex_app.git
$ cd demultiplex_app/scripts
$ python pipeline.py \
    -i /path/to/input/directory/ \
    -w /path/to/output/directory/ \
    -m /path/to/mapping/file.txt \
    -b 12 \
    -p /path/to/paired/ends/directory \
    -d /path/to/index/file.txt
```

## Build and Run as a Singularity container

`Singularity`, `Git`, and `make` must be installed to build the pipeline as a Singularity container. In addition, `sudo` privilege is required. If you don't ahve `sudo` privilege, I would recommend using VirtualBox and creating a Linux VM.

Build the pipeline container:

```
$ git clone git@github.com:mattmiller899/demultiplex_app.git
$ cd demultiplex_app/singularity
$ make img
```
This may take 15 minutes or more. The Singularity container will be built in the `demultiplex_app/singularity` directory.

Run the pipeline:
```
$ singularity run singularity/demultiplexer.img \
    -i /path/to/input/directory/ \
    -w /path/to/output/directory/ \
    -m /path/to/mapping/file.txt \
    -b 12 \
    -p /path/to/paired/ends/directory \
    -d /path/to/index/file.txt
```

## Argument Descriptions

```
-i INPUT_PATH (required): path to input directory (or single file)
-w WORK_DIR (required): path to directory that output will be written
-m MAPPING_FILE (required): path to file containing information about the barcodes. Must be in the QIIME mapping file format
-b BARCODE_LENGTH (required): length of barcodes/indices
-p PAIRED_PATH (optional): path to paired-ends directory (or single file)
-d INDEX_PATH (optional): path to index file. Including this argument indicates that reads in the INPUT_PATH and PAIRED_PATH do not have barcodes
```
