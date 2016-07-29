from sam_alignment import SamAlignment

def parse_sam(sam_file):
    with open(sam_file, 'r') as f:
        for line in f:
            yield SamAlignment(line)
