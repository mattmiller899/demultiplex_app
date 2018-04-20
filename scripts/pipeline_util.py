import glob
import gzip
import logging
from operator import attrgetter
import os
import re
import shutil
import subprocess
import traceback


class PipelineException(BaseException):
    pass


def get_sorted_file_list(dir_path):
    return tuple(
        sorted(
            [
                entry
                for entry
                in os.scandir(dir_path)
                if entry.is_file()
            ],
        key=attrgetter('name')
        )
    )


def create_output_dir(output_dir_name, parent_dir=None, input_dir=None):
    if parent_dir is not None and input_dir is None:
        pass
    elif input_dir is not None and parent_dir is None:
        parent_dir, _ = os.path.split(input_dir)
    else:
        raise ValueError('exactly one of parent_dir and input_dir must be None')

    log = logging.getLogger(name=__name__)
    output_dir = os.path.join(parent_dir, output_dir_name)
    if os.path.exists(output_dir):
        log.warning('directory "%s" already exists', output_dir)
    else:
        os.mkdir(output_dir)
    return output_dir


def get_forward_fastq_files(input_dir):
    log = logging.getLogger(name=__name__)
    input_glob = os.path.join(input_dir, '*_[R0]1*.fastq*')
    log.info('searching for forward read files with glob "%s"', input_glob)
    forward_fastq_files = glob.glob(input_glob)
    if len(forward_fastq_files) == 0:
        raise PipelineException('found no forward reads from glob "{}"'.format(input_glob))
    return forward_fastq_files


def get_associated_reverse_fastq_fp(forward_fp):
    forward_input_dir, forward_basename = os.path.split(forward_fp)
    reverse_fastq_basename = re.sub(
        string=forward_basename,
        pattern=r'_([0R])1',
        repl=lambda m: '_{}2'.format(m.group(1)))
    reverse_fastq_fp = os.path.join(forward_input_dir, reverse_fastq_basename)
    return reverse_fastq_fp


def get_associated_barcodes_fp(forward_fp):
    forward_input_dir, forward_basename = os.path.split(forward_fp)
    barcodes_fastq_basename = re.sub(
        string=forward_basename,
        pattern='_debarcoded_R1',
        repl='_barcodes')
    barcodes_fastq_fp = os.path.join(forward_input_dir, barcodes_fastq_basename)
    return barcodes_fastq_fp


def get_associated_barcodes_unpaired_fp(forward_fp):
    forward_input_dir, forward_basename = os.path.split(forward_fp)
    barcodes_fastq_basename = re.sub(
        string=forward_basename,
        pattern='_debarcoded',
        repl='_barcodes')
    barcodes_fastq_fp = os.path.join(forward_input_dir, barcodes_fastq_basename)
    return barcodes_fastq_fp


def run_cmd(cmd_line_list, log_file, **kwargs):
    log = logging.getLogger(name=__name__)
    try:
        with open(log_file, 'at') as log_file:
            if cmd_line_list[0] == 'complete':
                tmp = cmd_line_list[-1]
                tmp = ' '.join(tmp)
                cmd_line_list[-1] = tmp
                cmd_line_list.remove('complete')
            cmd_line_str = ' '.join((str(x) for x in cmd_line_list))
            log.info('executing "%s"', cmd_line_str)
            log_file.write('executing "{}"'.format(cmd_line_str))
            output = subprocess.call(
                cmd_line_list,
                stdout=log_file,
                stderr=subprocess.STDOUT,
                **kwargs)
            log.info(output)
        return output
    except subprocess.CalledProcessError as c:
        logging.exception(c)
        print(c.message)
        print(c.cmd)
        print(c.output)
        raise c
    except Exception as e:
        logging.exception(e)
        print('blarg!')
        print(e)
        traceback.print_exc()
        raise e
