import sys
from blast_alignment import BlastAlignment

def parse(blast_file, outfmt=6):
    with open(blast_file, 'r') as f:
        for line in f:
            yield BlastAlignment(line, outfmt)
