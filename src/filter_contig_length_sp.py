#!/usr/local/bioconda/bin/python
#Script taken and modified from Honey-bee project
import sys
from Bio import SeqIO
import re

cov_pattern = re.compile("cov_([0-9.]+)")

min_length, min_coverage, fasta_file_path = sys.argv[1:]

with open(fasta_file_path.replace('fasta', 'fltr.fasta'), 'w') as filtered_fasta:
        with open(fasta_file_path, 'rU') as input_fasta:
                def filtered_contigs_generator(minl, minc):
                        for contig in SeqIO.parse(input_fasta, 'fasta'):
                                result = cov_pattern.search(contig.name)
                                if result:
                                        if len(contig) >= minl and float(result.group(1)) >= minc:
                                                yield contig
                SeqIO.write(filtered_contigs_generator(int(min_length), float(min_coverage)), filtered_fasta, 'fasta')
