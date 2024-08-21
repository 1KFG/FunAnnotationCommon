#!/usr/bin/env python3

from Bio import SeqIO
import sys
import argparse
from natsort import natsorted

parser = argparse.ArgumentParser(description="Clean genome fasta file")
parser.add_argument("-l","--length", help="minimum length of the sequence",
                    type=int, default=1000)
args = parser.parse_args()
sequences = {}
for record in SeqIO.parse(sys.stdin, "fasta"):
    error = 0
    if len(record) > args.length:
        record.description = ""
        characters = {}
        for nuc in record:
            nuc = nuc.upper()
            if nuc not in characters:
                characters[nuc] = 1
            else:
                characters[nuc] += 1
        if len(characters) >= 4:
            sequences[record.id] = record

seqs = []

for seqname in natsorted(sequences.keys()):
    seqs.append(sequences[seqname])

SeqIO.write(seqs, sys.stdout, "fasta")
