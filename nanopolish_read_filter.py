#!/usr/bin/env python3

"""
This script takes two arguments:
  1) a Nanopolish-extracted fasta file of reads
  2) a fastq.gz file of acceptable reads
It outputs a fasta file of reads in the Nanopolish format, including only the reads which are
also in the acceptable reads fastq.
"""

import sys
import gzip


def main():
    nanopolish_fasta_filename = sys.argv[1]
    fastq_filename = sys.argv[2]
    fastq_read_names = load_fastq_read_names(fastq_filename)
    with open(nanopolish_fasta_filename, 'rt') as nanopolish_fasta_file:
        for line in nanopolish_fasta_file:
            read_name_line = line.rstrip()
            sequence_line = next(nanopolish_fasta_file).rstrip()
            read_name = read_name_line[1:].split('_Basecall_')[0]
            if read_name in fastq_read_names:
                print(read_name_line)
                print(sequence_line)


def load_fastq_read_names(fastq_filename):
    read_names = set()
    with gzip.open(fastq_filename, 'rt') as fastq_file:
        for line in fastq_file:
            read_names.add(line.strip()[1:].split()[0])
            _ = next(fastq_file)
            _ = next(fastq_file)
            _ = next(fastq_file)
    return read_names


if __name__ == '__main__':
    main()
