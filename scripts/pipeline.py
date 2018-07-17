"""
Design principles:
  1. reading this script should not require a high level of Python experience
  2. pipeline steps should correspond directly to functions
  3. there should be a single function that represents the pipeline
  4. for development the pipeline should be 'restartable'
"""
import argparse
import glob
import gzip
import itertools
import logging
import os
import re
import shutil
import sys
import qiime
import scandir

from pipeline_util import create_output_dir, get_forward_fastq_files, get_associated_reverse_fastq_fp, get_associated_barcodes_fp, run_cmd, PipelineException, get_associated_barcodes_unpaired_fp
from fasta_qual_to_fastq import fasta_qual_to_fastq


def main():
    logging.basicConfig(level=logging.INFO)
    args = get_args()

    Pipeline(**args.__dict__).run(input_file=args.input_file)
    return 0

def get_args():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('-i', '--input-file', required=True,
                            help='path to the input file')
    arg_parser.add_argument('-w', '--work-dir', required=True,
                            help='path to the output directory')

    arg_parser.add_argument('-b', '--barcode-length', required=True, type=int,
                            help='length of the barcodes')

    arg_parser.add_argument('-m', '--mapping-file', required=True,
                            help='path to file containing information about the barcodes. Must be in the QIIME mapping '
                                 'file format')

    arg_parser.add_argument('-p', '--paired-ends', default='', 
                            help='path to paired end directory')

    '''
    arg_parser.add_argument('--uchime-ref-db-fp', default='/16SrDNA/pr2/pr2_gb203_version_4.5.fasta',
                            help='database for vsearch --uchime_ref')


    arg_parser.add_argument('--pear-min-overlap', required=True, type=int,
                            help='-v/--min-overlap for pear')
    arg_parser.add_argument('--pear-max-assembly-length', required=True, type=int,
                            help='-m/--max-assembly-length for pear')
    arg_parser.add_argument('--pear-min-assembly-length', required=True, type=int,
                            help='-m/--min-assembly-length for pear')

    arg_parser.add_argument('--chimera-detection-usearch', action='store_true', default=False)

    arg_parser.add_argument('--vsearch-filter-maxee', required=True, type=int,
                            help='fastq_maxee for vsearch')
    arg_parser.add_argument('--vsearch-filter-trunclen', required=True, type=int,
                            help='fastq_trunclen for vsearch')

    arg_parser.add_argument('--vsearch-derep-minuniquesize', required=True, type=int,
                            help='minimum unique size for vsearch -derep_fulllength')
    '''
    args = arg_parser.parse_args()
    return args


class Pipeline:
    def __init__(self,
            input_file,
            work_dir,
            mapping_file,
            barcode_length,
            paired_ends,
            **kwargs  # allows some command line arguments to be ignored
            ):
        
        self.input_file = input_file
        self.work_dir = work_dir

        self.mapping_file = mapping_file
        self.barcode_length = barcode_length
        self.paired_ends_path = paired_ends
        self.paired_ends = False
        self.paired_ends_dir = False
        if self.paired_ends_path != '':
            self.paired_ends = True
            if os.path.isdir(self.paired_ends_path):
                self.paired_ends_dir = True


    def run(self, input_file):
        output_dir_list = list()
        output_dir_list.append(self.step_01_remove_barcodes(input_file=input_file))
        output_dir_list.append(self.step_02_split_libraries(input_dir=output_dir_list[-1]))
        output_dir_list.append(self.step_03_demultiplex(input_dir=output_dir_list[-1]))
        if self.paired_ends is True:
            output_dir_list.append(self.step_04_make_paired_end_files(input_dir=output_dir_list[-1]))

        return output_dir_list


    def initialize_step(self):
        function_name = sys._getframe(1).f_code.co_name
        log = logging.getLogger(name=function_name)
        #Make step output_dir
        output_dir = create_output_dir(output_dir_name=function_name, parent_dir=self.work_dir)
        #Make specific file output_dir
        name, ext = os.path.splitext(self.input_file)
        fileout_dir = create_output_dir(output_dir_name=os.path.basename(name), parent_dir=output_dir)
        return log, fileout_dir


    def complete_step(self, log, output_dir):
        output_dir_list = sorted(os.listdir(output_dir))
        if len(output_dir_list) == 0:
            raise PipelineException('ERROR: no output files in directory "{}"'.format(output_dir))
        else:
            log.info('output files:\n\t%s', '\n\t'.join(os.listdir(output_dir)))
            # apply FastQC to all .fastq files
            subfolders = [f.path for f in scandir.scandir(output_dir) if f.is_dir()]
            if len(subfolders) == 0 or (len(subfolders) == 1 and os.path.basename(subfolders[0]) == "fastqc_results"):
                fastq_file_list = [
                    os.path.join(output_dir, output_file)
                    for output_file
                    in output_dir_list
                    if re.search(pattern=r'\.fastq(\.gz)?$', string=output_file)
                ]

                if len(fastq_file_list) == 0:
                    log.info('no .fastq files found in "{}"'.format(output_dir))
                else:
                    fastqc_output_dir = os.path.join(output_dir, 'fastqc_results')
                    if not os.path.isdir(fastqc_output_dir):
                        os.mkdir(fastqc_output_dir)
                    for file_fp in fastq_file_list:
                        run_cmd(
                            [
                                'fastqc',
                                '--threads', '1',
                                '--outdir', fastqc_output_dir,
                                file_fp
                            ],
                            log_file=os.path.join(fastqc_output_dir, 'log')
                        )
            else:
                fastqc_output_dir = os.path.join(output_dir, 'fastqc_results')
                if not os.path.isdir(fastqc_output_dir):
                    os.mkdir(fastqc_output_dir)
                for folder in subfolders:
                    folder_basename = os.path.basename(folder)
                    folder_dir_list = sorted(os.listdir(folder))
                    if folder_basename == "fastqc_results":
                        continue
                    folder_output_dir = os.path.join(fastqc_output_dir, folder_basename)
                    os.mkdir(folder_output_dir)
                    fastq_file_list = [
                        os.path.join(folder, output_file)
                        for output_file
                        in folder_dir_list
                        if re.search(pattern=r'\.fastq(\.gz)?$', string=output_file)
                    ]
                    if len(fastq_file_list) == 0:
                        log.info('no .fastq files found in "{}"'.format(output_dir))
                    else:
                        for file_fp in fastq_file_list:
                            run_cmd(
                                [
                                    'fastqc',
                                    '--threads', '1',
                                    '--outdir', folder_output_dir,
                                    file_fp
                                ],
                                log_file=os.path.join(fastqc_output_dir, 'log')
                            )


    def step_01_remove_barcodes(self, input_file):
        log, output_dir = self.initialize_step()
        if len(os.listdir(output_dir)) > 0:
            log.info('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            log.info('Going to remove barcodes')
            if self.paired_ends is True:
                if self.paired_ends_dir is True:
                    paired_end_file = get_associated_reverse_fastq_fp(forward_fp=input_file, reverse_input_dir=self.paired_ends_path)
                else:
                    paired_end_file = self.paired_ends_path
                log.info('removing barcodes from forward reads "%s"', input_file)
                log.info('removing barcodes from reverse reads "%s"', paired_end_file)
                run_cmd([
                        'python', '/miniconda/bin/extract_barcodes.py',
                        '-f', input_file,
                        '-r', paired_end_file,
                        '-c', 'barcode_paired_end',
                        '-m', str(self.mapping_file),
                        '-l', str(self.barcode_length),
                        '-L', str(self.barcode_length),
                        '-o', str(output_dir)
                    ],
                    log_file=os.path.join(output_dir, 'log')
                )
                forward_fastq_basename = os.path.basename(input_file)
                tmp = re.split('_([0R])1', forward_fastq_basename)
                file_name = tmp[0] + re.split('.fastq', tmp[2])[0]
                os.rename(os.path.join(output_dir, 'reads1.fastq'), os.path.join(output_dir, file_name + '_debarcoded_R1.fastq'))
                os.rename(os.path.join(output_dir, 'reads2.fastq'), os.path.join(output_dir, file_name + '_debarcoded_R2.fastq'))
                os.rename(os.path.join(output_dir, 'barcodes.fastq'), os.path.join(output_dir, file_name + '_barcodes.fastq'))
            else:
                log.info('removing barcodes from "%s"', input_file)
                run_cmd([
                    'python', '/miniconda/bin/extract_barcodes.py',
                    '-f', input_file,
                    '-c', 'barcode_single_end',
                    '-m', str(self.mapping_file),
                    '-l', str(self.barcode_length),
                    '-o', str(output_dir)
                   ],
                    log_file=os.path.join(output_dir, 'log')
                )
                file_basename = os.path.basename(input_file)
                file_name = re.split('.fastq', file_basename)[0]
                os.rename(os.path.join(output_dir, 'reads.fastq'),
                          os.path.join(output_dir, file_name + '_debarcoded.fastq'))
                os.rename(os.path.join(output_dir, 'barcodes.fastq'),
                          os.path.join(output_dir, file_name + '_barcodes.fastq'))
        self.complete_step(log, output_dir)
        return output_dir


    def step_02_split_libraries(self, input_dir):
        log, output_dir = self.initialize_step()
        if len(os.listdir(output_dir)) > 0:
            log.info('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            log.info('Splitting library based on barcodes')
            if self.paired_ends is True:
                for forward_fastq_fp in get_forward_fastq_files(input_dir=input_dir):
                    reverse_fastq_fp = get_associated_reverse_fastq_fp(forward_fp=forward_fastq_fp, reverse_input_dir=input_dir)
                    barcodes_fp = get_associated_barcodes_fp(forward_fastq_fp)
                    log.info('Splitting libraries of "%s" and "%s" with "%s"', forward_fastq_fp, reverse_fastq_fp, barcodes_fp)
                    #TODO Make argument for -q and --max_barcode_errors and rev_comp (1-step vs 2-step PCR?)?
                    run_cmd([
                            'python', '/miniconda/bin/split_libraries_fastq.py',
                            '-o', str(output_dir),
                            '-b', str(barcodes_fp) + ',' + str(barcodes_fp),
                            '-i', str(forward_fastq_fp) + ',' + str(reverse_fastq_fp),
                            '-m', str(self.mapping_file),
                            '--barcode_type', str(self.barcode_length),
                            '-q', '0',
                            '--max_barcode_errors', '0',
                            '--rev_comp_barcode',
                            '--phred_offset=33',
                            '--store_demultiplexed_fastq'
                        ],
                        log_file=os.path.join(output_dir, 'log')
                    )
                    forward_fastq_basename = os.path.basename(forward_fastq_fp)
                    file_name = re.split('_([0R])1', forward_fastq_basename)[0]
                    os.rename(os.path.join(output_dir, 'seqs.fastq'),
                              os.path.join(output_dir, file_name + '_seqs.fastq'))
                    os.rename(os.path.join(output_dir, 'histograms.txt'),
                              os.path.join(output_dir, file_name + '_histograms.txt'))
                    os.rename(os.path.join(output_dir, 'split_library_log.txt'),
                              os.path.join(output_dir, file_name + '_split_library_log.txt'))
                    os.remove(os.path.join(output_dir, 'seqs.fna'))
            else:
                input_files_glob = os.path.join(input_dir, '*.fastq*')
                for input_file in glob.glob(input_files_glob):
                    log.info("Input file: %s", input_file)
                    if "barcodes.fastq" not in input_file:
                        barcodes_fp = get_associated_barcodes_unpaired_fp(input_file)
                        log.info('Splitting libraries of "%s" with "%s"', input_file, barcodes_fp)
                        # TODO Make argument for -q and --max_barcode_errors and rev_comp (1-step vs 2-step PCR?)?
                        run_cmd([
                            'python', '/miniconda/bin/split_libraries_fastq.py',
                            '-o', str(output_dir),
                            '-b', str(barcodes_fp),
                            '-i', str(input_file),
                            '-m', str(self.mapping_file),
                            '--barcode_type', str(self.barcode_length),
                            '-q', '0',
                            '--max_barcode_errors', '0',
                            '--rev_comp_barcode',
                            '--phred_offset=33',
                            '--store_demultiplexed_fastq'
                            ],
                            log_file=os.path.join(output_dir, 'log')
                        )
                        file_basename = os.path.basename(input_file)
                        file_name = re.split('.fastq', file_basename)[0]
                        os.rename(os.path.join(output_dir, 'seqs.fastq'),
                                  os.path.join(output_dir, file_name + '_seqs.fastq'))
                        os.rename(os.path.join(output_dir, 'histograms.txt'),
                                  os.path.join(output_dir, file_name + '_histograms.txt'))
                        os.rename(os.path.join(output_dir, 'split_library_log.txt'),
                                  os.path.join(output_dir, file_name + '_split_library_log.txt'))
                        os.remove(os.path.join(output_dir, 'seqs.fna'))
        self.complete_step(log, output_dir)
        return output_dir

    def step_03_demultiplex(self, input_dir):
        log, output_dir = self.initialize_step()
        if len(os.listdir(output_dir)) > 0:
            log.info('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            log.info('Splitting seqs files based on sampleID')

            input_files_glob = os.path.join(input_dir, '*_seqs.fastq')
            log.info('input file glob: "%s"', input_files_glob)
            for split_fastq_fp in glob.glob(input_files_glob):
                log.info('Splitting seq file "%s"', split_fastq_fp)
                split_fastq_basename = os.path.basename(split_fastq_fp)
                file_name = re.split('_seqs', split_fastq_basename)[0]
                run_cmd([
                        'python', '/miniconda/bin/split_sequence_file_on_sample_ids.py',
                        '-i', split_fastq_fp,
                        '-o', str(output_dir),
                        '--file_type', 'fastq'
                    ],
                    log_file=os.path.join(output_dir, 'log')
                )
                for sample_file in glob.glob(os.path.join(output_dir, '*.fastq')):
                    sample_file_basename = os.path.basename(sample_file)
                    os.rename(sample_file, os.path.join(output_dir, file_name + '_' + sample_file_basename))

        self.complete_step(log, output_dir)
        return output_dir


    def step_04_make_paired_end_files(self, input_dir):
        log, output_dir = self.initialize_step()
        if len(os.listdir(output_dir)) > 0:
            log.info('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            log.info('Splitting sample files back into paired end files')
            input_files_glob = os.path.join(input_dir, '*.fastq')
            for input_file in glob.glob(input_files_glob):
                log.info('Making paired end files with "%s"', input_file)
                input_file_basename = os.path.basename(input_file)
                input_file_no_fastq = input_file_basename.split('.fastq')[0]
                out1 = open(os.path.join(output_dir, input_file_no_fastq + '_R1.fastq'), 'w')
                out2 = open(os.path.join(output_dir, input_file_no_fastq + '_R2.fastq'), 'w')
                with open(input_file, 'r') as f:
                    count = 0
                    write_file = 1
                    for l in f:
                        if count % 4 == 0:
                            write_file = l.split(' ')[2].split(':')[0]
                            write_line = l.split(' ')
                            if write_file == '1':
                                out1.write('@' + ' '.join(write_line[1:]))
                            elif write_file == '2':
                                out2.write('@' + ' '.join(write_line[1:]))
                            else:
                                log.info('Bad number for paired end files "%s" in "%s', write_file, l)
                        elif write_file == '1':
                            out1.write(l)
                        elif write_file == '2':
                            out2.write(l)
                        count += 1
                out1.close()
                out2.close()

        self.complete_step(log, output_dir)
        return output_dir
    '''

    def step_05_copy_and_compress(self, input_dir):
        log, output_dir = self.initialize_step()
        if len(os.listdir(output_dir)) > 0:
            log.info('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            log.debug('output_dir: %s', output_dir)
            # Check for fasta and qual files
            fasta_file_glob = os.path.join(input_dir, '*fasta*')
            qual_file_glob = os.path.join(input_dir, '*.qual*')
            fasta_files = sorted(glob.glob(fasta_file_glob))
            qual_files = sorted(glob.glob(qual_file_glob))
            fastas_no_slashes = [fasta.replace('\\', ' ').replace('/', ' ').split()[-1] for fasta in fasta_files]
            quals_no_slashes = [qual.replace('\\', ' ').replace('/', ' ').split()[-1] for qual in qual_files]
            fastas = [fasta.split(".")[0] for fasta in fastas_no_slashes]
            quals = [qual.split('.')[0] for qual in quals_no_slashes]
            for fasta in fastas:
                for qual in quals:
                    if fasta == qual:
                        tmpfasta = os.path.join(input_dir, fasta + ".fasta")
                        tmpqual = os.path.join(input_dir, qual + ".qual")
                        tmpfastq = os.path.join(input_dir, fasta + ".fastq")
                        log.info("Creating %s", tmpfastq)
                        fasta_qual_to_fastq(tmpfasta, tmpqual, tmpfastq)
                        log.info('%s created', tmpfastq)
                        quals.remove(qual)
                        break

            subfolders = [f.path for f in scandir.scandir(input_dir) if f.is_dir()]
            for folder in subfolders:
                folder_basename = os.path.basename(folder)
                if folder_basename == "fastqc_results":
                    continue
                input_file_glob = os.path.join(folder, '*.fastq*')
                log.debug('input_file_glob: %s', input_file_glob)
                input_fp_list = sorted(glob.glob(input_file_glob))
                log.info('input files: %s', input_fp_list)

                if len(input_fp_list) == 0:
                    raise PipelineException('found no fastq files in directory "{}"'.format(folder))

                out_dir = create_output_dir(output_dir_name=folder_basename, parent_dir=output_dir)
                for input_fp in input_fp_list:
                    destination_fp = os.path.join(out_dir, os.path.basename(input_fp))
                    if input_fp.endswith('.gz'):
                        with open(input_fp, 'rb') as f, open(destination_fp, 'wb') as g:
                            shutil.copyfileobj(fsrc=f, fdst=g)
                    else:
                        destination_fp = destination_fp + '.gz'
                        with open(input_fp, 'rt') as f, gzip.open(destination_fp, 'wt') as g:
                            shutil.copyfileobj(fsrc=f, fdst=g)

        self.complete_step(log, output_dir)
        return output_dir


    def step_06_remove_primers(self, input_dir):
        log, output_dir = self.initialize_step()
        if len(os.listdir(output_dir)) > 0:
            log.info('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            if self.paired_ends is True:
                log.info('using cutadapt "%s"', self.cutadapt_executable_fp)
                subfolders = [f.path for f in scandir.scandir(input_dir) if f.is_dir()]
                for folder in subfolders:
                    folder_basename = os.path.basename(folder)
                    if folder_basename == "fastqc_results":
                        continue
                    out_dir = create_output_dir(output_dir_name=folder_basename, parent_dir=output_dir)
                    for forward_fastq_fp in get_forward_fastq_files(input_dir=folder):
                        log.info('removing forward primers from file "%s"', forward_fastq_fp)
                        forward_fastq_basename = os.path.basename(forward_fastq_fp)

                        reverse_fastq_fp = get_associated_reverse_fastq_fp(forward_fp=forward_fastq_fp)
                        log.info('removing reverse primers from file "%s"', reverse_fastq_fp)

                        trimmed_forward_fastq_fp = os.path.join(
                            out_dir,
                            re.sub(
                                string=forward_fastq_basename,
                                pattern='_([0R])1',
                                repl=lambda m: '_trimmed_{{name}}.{}1'.format(m.group(1))))
                        trimmed_reverse_fastq_fp = os.path.join(
                            out_dir,
                            re.sub(
                                string=forward_fastq_basename,
                                pattern='_([0R])1',
                                repl=lambda m: '_trimmed_{{name}}.{}2'.format(m.group(1))))
                        run_cmd([
                            self.cutadapt_executable_fp,
                            '-a', self.forward_primer,
                            '-A', self.reverse_primer,
                            '-o', trimmed_forward_fastq_fp,
                            '-p', trimmed_reverse_fastq_fp,
                            '-m', str(self.cutadapt_min_length),
                            forward_fastq_fp,
                            reverse_fastq_fp
                            ],
                            log_file=os.path.join(out_dir, 'log')
                        )
            #Nanopore porechop
            elif self.nanopore is True:
                log.info('Using Porechop %s', self.porechop_executable_fp)
                subfolders = [f.path for f in scandir.scandir(input_dir) if f.is_dir()]
                for folder in subfolders:
                    folder_basename = os.path.basename(folder)
                    if folder_basename == "fastqc_results":
                        continue
                    out_dir = create_output_dir(output_dir_name=folder_basename, parent_dir=output_dir)
                    input_files_glob = os.path.join(input_dir, '*.fastq*')
                    for input_file in glob.glob(input_files_glob):
                        log.info("Removing primers from %s", input_file)

            #Single reads
            else:
                pass
        self.complete_step(log, output_dir)
        return output_dir


    def step_03_merge_forward_reverse_reads_with_vsearch(self, input_dir):
        log, output_dir = self.initialize_step()
        if len(os.listdir(output_dir)) > 0:
            log.info('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            log.info('vsearch executable: "%s"', self.vsearch_executable_fp)

            for forward_fastq_fp in get_forward_fastq_files(input_dir=input_dir):
                reverse_fastq_fp = get_associated_reverse_fastq_fp(forward_fp=forward_fastq_fp)

                joined_fastq_basename = re.sub(
                    string=os.path.basename(forward_fastq_fp),
                    pattern=r'_([0R]1)',
                    repl=lambda m: '_merged'.format(m.group(1))
                )
                joined_fastq_fp = os.path.join(output_dir, joined_fastq_basename)[:-3]
                log.info('writing joined paired-end reads to "%s"', joined_fastq_fp)

                notmerged_fwd_fastq_basename = re.sub(
                    string=os.path.basename(forward_fastq_fp),
                    pattern=r'_([0R]1)',
                    repl=lambda m: '_notmerged_fwd'.format(m.group(1))
                )
                notmerged_fwd_fastq_fp = os.path.join(output_dir, notmerged_fwd_fastq_basename)[:-3]

                notmerged_rev_fastq_basename = re.sub(
                    string=os.path.basename(forward_fastq_fp),
                    pattern=r'_([0R]1)',
                    repl=lambda m: '_notmerged_rev'.format(m.group(1))
                )
                notmerged_rev_fastq_fp = os.path.join(output_dir, notmerged_rev_fastq_basename)[:-3]

                run_cmd([
                        self.vsearch_executable_fp,
                        '--fastq_mergepairs', forward_fastq_fp,
                        '--reverse', reverse_fastq_fp,
                        '--fastqout', joined_fastq_fp,
                        '--fastqout_notmerged_fwd', notmerged_fwd_fastq_fp,
                        '--fastqout_notmerged_rev', notmerged_rev_fastq_fp,
                        '--fastq_minovlen', str(self.pear_min_overlap),
                        '--fastq_maxlen', str(self.pear_max_assembly_length),
                        '--fastq_minlen', str(self.pear_min_assembly_length),
                        '--threads', str(self.core_count)
                    ],
                    log_file = os.path.join(output_dir, 'log')
                )

                es(glob.glob(os.path.join(output_dir, '*.fastq')))

        self.complete_step(log, output_dir)
        return output_dir

    def step_03_merge_forward_reverse_reads_with_pear(self, input_dir):
        log, output_dir = self.initialize_step()
        if len(os.listdir(output_dir)) > 0:
            log.info('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            log.info('PEAR executable: "%s"', self.pear_executable_fp)

            for compressed_forward_fastq_fp in get_forward_fastq_files(input_dir=input_dir):
                compressed_reverse_fastq_fp = get_associated_reverse_fastq_fp(forward_fp=compressed_forward_fastq_fp)

                forward_fastq_fp, reverse_fastq_fp = unes(
                    compressed_forward_fastq_fp,
                    compressed_reverse_fastq_fp,
                    target_dir=output_dir
                )

                joined_fastq_basename = re.sub(
                    string=os.path.basename(forward_fastq_fp),
                    pattern=r'_([0R]1)',
                    repl=lambda m: '_merged'.format(m.group(1)))[:-6]

                joined_fastq_fp_prefix = os.path.join(output_dir, joined_fastq_basename)
                log.info('joining paired ends from "%s" and "%s"', forward_fastq_fp, reverse_fastq_fp)
                log.info('writing joined paired-end reads to "%s"', joined_fastq_fp_prefix)
                run_cmd([
                        self.pear_executable_fp,
                        '-f', forward_fastq_fp,
                        '-r', reverse_fastq_fp,
                        '-o', joined_fastq_fp_prefix,
                        '--min-overlap', str(self.pear_min_overlap),
                        '--max-assembly-length', str(self.pear_max_assembly_length),
                        '--min-assembly-length', str(self.pear_min_assembly_length),
                        '-j', str(self.core_count)
                    ],
                    log_file = os.path.join(output_dir, 'log')
                )

                # delete the uncompressed input files
                os.remove(forward_fastq_fp)
                os.remove(reverse_fastq_fp)

                es(glob.glob(joined_fastq_fp_prefix + '.*.fastq'))

        self.complete_step(log, output_dir)
        return output_dir

    def step_04_qc_reads_with_vsearch(self, input_dir):
        log, output_dir = self.initialize_step()
        if len(os.listdir(output_dir)) > 0:
            log.info('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            input_files_glob = os.path.join(input_dir, '*.assembled.fastq.gz')
            log.info('input file glob: "%s"', input_files_glob)
            for assembled_fastq_fp in glob.glob(input_files_glob):
                input_file_basename = os.path.basename(assembled_fastq_fp)
                output_file_basename = re.sub(
                    string=input_file_basename,
                    pattern='\.fastq\.gz',
                    repl='.ee{}trunc{}.fastq.gz'.format(self.vsearch_filter_maxee, self.vsearch_filter_trunclen)[:-3]
                )
                output_fastq_fp = os.path.join(output_dir, output_file_basename)

                log.info('vsearch executable: "%s"', self.vsearch_executable_fp)
                log.info('filtering "%s"', assembled_fastq_fp)
                run_cmd([
                        self.vsearch_executable_fp,
                        '-fastq_filter', assembled_fastq_fp,
                        '-fastqout', output_fastq_fp,
                        '-fastq_maxee', str(self.vsearch_filter_maxee),
                        '-fastq_trunclen', str(self.vsearch_filter_trunclen),
                        '-threads', str(self.core_count)
                    ],
                    log_file = os.path.join(output_dir, 'log')
                )

            es(glob.glob(os.path.join(output_dir, '*.assembled.*.fastq')))

        self.complete_step(log, output_dir)
        return output_dir

    def step_05_combine_runs(self, input_dir):
        log, output_dir = self.initialize_step()
        if len(os.listdir(output_dir)) > 0:
            log.info('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            log.info('input directory listing:\n\t%s', '\n\t'.join(os.listdir(input_dir)))
            input_files_glob = os.path.join(input_dir, '*.assembled.*.fastq.gz')
            log.info('input file glob: "%s"', input_files_glob)
            input_fp_list = sorted(glob.glob(input_files_glob))
            log.info('combining files:\n\t%s', '\n\t'.join(input_fp_list))

            output_file_name = get_combined_file_name(input_fp_list=input_fp_list)

            log.info('combined file: "%s"', output_file_name)
            output_fp = os.path.join(output_dir, output_file_name)
            with gzip.open(output_fp, 'wt') as output_file:
                for input_fp in input_fp_list:
                    with gzip.open(input_fp, 'rt') as input_file:
                        shutil.copyfileobj(fsrc=input_file, fdst=output_file)

        self.complete_step(log, output_dir)
        return output_dir

    def step_06_dereplicate_sort_remove_low_abundance_reads(self, input_dir):
        log, output_dir = self.initialize_step()
        if len(os.listdir(output_dir)) > 0:
            log.info('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            log.info('input directory listing:\n\t%s', '\n\t'.join(os.listdir(input_dir)))
            input_files_glob = os.path.join(input_dir, '*.assembled.*.fastq.gz')
            log.info('input file glob: "%s"', input_files_glob)
            input_fp_list = sorted(glob.glob(input_files_glob))

            for input_fp in input_fp_list:
                output_fp = os.path.join(
                    output_dir,
                    re.sub(
                        string=os.path.basename(input_fp),
                        pattern='\.fastq\.gz$',
                        repl='.derepmin{}.fasta'.format(self.vsearch_derep_minuniquesize)))

                uc_fp = os.path.join(
                    output_dir,
                    re.sub(
                        string=os.path.basename(input_fp),
                        pattern='\.fastq\.gz$',
                        repl='.derepmin{}.txt'.format(self.vsearch_derep_minuniquesize)))

                run_cmd([
                        self.vsearch_executable_fp,
                        '-derep_fulllength', input_fp,
                        '-output', output_fp,
                        '-uc', uc_fp,
                        '-sizeout',
                        '-minuniquesize', str(self.vsearch_derep_minuniquesize),
                        '-threads', str(self.core_count)
                    ],
                    log_file = os.path.join(output_dir, 'log')
                )

            es(glob.glob(os.path.join(output_dir, '*.fasta')))

        self.complete_step(log, output_dir)
        return output_dir

    def step_07_cluster_97_percent(self, input_dir):
        log, output_dir = self.initialize_step()
        if len(os.listdir(output_dir)) > 0:
            log.info('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            # can usearch read gzipped files? no
            input_files_glob = os.path.join(input_dir, '*.fasta.gz')
            input_fp_list = sorted(glob.glob(input_files_glob))

            for compressed_input_fp in input_fp_list:

                input_fp, *_ = unes(compressed_input_fp, target_dir=output_dir)
                log.debug('input_fp: "%s"', input_fp)
                otu_output_fp = os.path.join(
                    output_dir,
                    re.sub(
                        string=os.path.basename(input_fp),
                        pattern='\.fasta$',
                        repl='.rad3.fasta'
                    )
                )

                uparse_output_fp = os.path.join(
                    output_dir,
                    re.sub(
                        string=os.path.basename(input_fp),
                        pattern='\.fasta$',
                        repl='.rad3.txt'
                    )
                )

                run_cmd([
                        self.usearch_executable_fp,
                        '-cluster_otus', input_fp,
                        '-otus', otu_output_fp,
                        '-relabel', 'OTU_',
                        # '-sizeout',
                        '-uparseout', uparse_output_fp
                    ],
                    log_file = os.path.join(output_dir, 'log')
                )

                os.remove(input_fp)

        self.complete_step(log, output_dir)
        return output_dir

    def step_08_reference_based_chimera_detection(self, input_dir):
        log, output_dir = self.initialize_step()
        if len([entry for entry in os.scandir(output_dir) if not entry.name.startswith('.')]) > 0:
            log.info('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            input_fps = glob.glob(os.path.join(input_dir, '*.fasta'))
            for input_fp in input_fps:
                uchimeout_fp = os.path.join(
                    output_dir,
                    re.sub(
                        string=os.path.basename(input_fp),
                        pattern='\.fasta$',
                        repl='.uchime.txt'))

                notmatched_fp = os.path.join(
                    output_dir,
                    re.sub(
                        string=os.path.basename(input_fp),
                        pattern='\.fasta$',
                        repl='.uchime.fasta'))

                log.info('starting chimera detection on file "%s"', input_fp)
                if self.chimera_detection_usearch:
                    run_cmd([
                            self.usearch_executable_fp,
                            '-uchime2_ref', input_fp,
                            '-db', self.uchime_ref_db_fp,
                            '-uchimeout', uchimeout_fp,
                            '-mode', 'balanced',
                            '-strand', 'plus',
                            '-notmatched', notmatched_fp,
                            '-threads', str(self.core_count)
                        ],
                        log_file=os.path.join(output_dir, 'log')
                    )

                else:
                    run_cmd([
                            self.vsearch_executable_fp,
                            '-uchime_ref', input_fp,
                            '-db', self.uchime_ref_db_fp,
                            '-uchimeout', uchimeout_fp,
                            '-mode', 'balanced',
                            '-strand', 'plus',
                            '-nonchimeras', notmatched_fp,
                            '-threads', str(self.core_count)
                        ],
                        log_file = os.path.join(output_dir, 'log')
                    )


        self.complete_step(log, output_dir)
        return output_dir

    def step_09_create_otu_table(self, input_dir):
        log, output_dir = self.initialize_step()
        if len([entry for entry in os.scandir(output_dir) if not entry.name.startswith('.')]) > 0:
            log.info('output directory "%s" is not empty, this step will be skipped', output_dir)
        else:
            otus_fp, *_ = glob.glob(os.path.join(input_dir, '*rad3.uchime.fasta'))
            input_fps = glob.glob(os.path.join(self.work_dir, 'step_03*', '*.assembled.fastq.gz'))
            for input_fp in input_fps:
                fasta_fp = os.path.join(
                    output_dir,
                    re.sub(
                        string=os.path.basename(input_fp),
                        pattern='\.fastq\.gz',
                        repl='.fasta'))

                log.info('convert fastq file\n\t%s\nto fasta file\n\t%s', input_fp, fasta_fp)
                run_cmd([
                    self.vsearch_executable_fp,
                    '--fastq_filter', input_fp,
                    '--fastaout', fasta_fp
                    ],
                    log_file=os.path.join(output_dir, 'log')
                )

                otu_table_fp = os.path.join(
                    output_dir,
                    re.sub(
                        string=os.path.basename(input_fp),
                        pattern='\.assembled\.fastq\.gz$',
                        repl='.uchime.otutab.txt'
                    )
                )
                otu_table_biom_fp = os.path.join(
                    output_dir,
                    re.sub(
                        string=os.path.basename(input_fp),
                        pattern='\.assembled\.fastq\.gz$',
                        repl='.uchime.otutab.json'
                    )
                )

                run_cmd([
                        self.vsearch_executable_fp,
                        '--usearch_global', fasta_fp,
                        '--db', otus_fp,
                        '--id', '0.97',
                        '--biomout', otu_table_biom_fp,
                        '--otutabout', otu_table_fp
                    ],
                    log_file = os.path.join(output_dir, 'log')
                )

                # os.remove(fasta_fp)

        self.complete_step(log, output_dir)
        return output_dir

'''

def get_combined_file_name(input_fp_list):
    if len(input_fp_list) == 0:
        raise PipelineException('get_combined_file_name called with empty input')

    def sorted_unique_elements(elements):
        return sorted(set(elements))

    return '_'.join(  # 'Mock_Run3_Run4_V4.fastq.gz'
        itertools.chain.from_iterable(  # ['Mock', 'Run3', 'Run4', 'V4.fastq.gz']
            map(  # [{'Mock'}, {'Run3', 'Run4'}, {'V4.fastq.gz'}]
                sorted_unique_elements,
                zip(  # [('Mock', 'Mock'), ('Run4', 'Run3'), ('V4.fastq.gz', 'V4.fastq.gz')]
                    *[  # [('Mock', 'Run3', 'V4.fastq.gz'), ('Mock', 'Run4', 'V4.fastq.gz')]
                        os.path.basename(fp).split('_')  # ['Mock', 'Run3', 'V4.fastq.gz']
                        for fp  # '/some/data/Mock_Run3_V4.fastq.gz'
                        in input_fp_list  # ['/input/data/Mock_Run3_V4.fastq.gz', '/input_data/Mock_Run4_V4.fastq.gz']
                    ]))))


if __name__ == '__main__':
    main()
