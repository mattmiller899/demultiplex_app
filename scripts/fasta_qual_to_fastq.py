import argparse
import itertools


def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return itertools.zip_longest(*args, fillvalue=fillvalue)


def main():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('--fasta')
    arg_parser.add_argument('--qual')
    arg_parser.add_argument('--fastq')

    args = arg_parser.parse_args()

    fasta_qual_to_fastq(**args.__dict__)


def fasta_qual_to_fastq(fasta, qual, fastq):
    with open(fasta, 'rt') as fasta_file, open(qual, 'rt') as qual_file, open(fastq, 'wt') as fastq_file:
        for (fasta_hdr, fasta_seq), (qual_hdr, qual_seq) in zip(grouper(fasta_file, 2), grouper(qual_file, 2)):
            fasta_hdr = fasta_hdr.strip()
            fasta_seq = fasta_seq.strip()
            qual_hdr = qual_hdr.strip()
            qual_seq = qual_seq.strip()
            if fasta_hdr == qual_hdr and len(fasta_seq) == len(qual_seq):
                fastq_file.write('@{}\n'.format(fasta_hdr[1:]))
                fastq_file.write(fasta_seq.strip() + '\n')
                fastq_file.write('+\n')
                fastq_file.write(qual_seq.strip() + '\n')
            else:
                print('problem with')
                print('  FASTA header: {}'.format(fasta_hdr))
                print('  FASTA seq   : {}'.format(fasta_seq))
                print('  QUAL header : {}'.format(qual_hdr))
                print('  QUAL seq    : {}'.format(qual_seq))
                break


if __name__ == '__main__':
    main()