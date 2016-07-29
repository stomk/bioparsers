from cigar import Cigar

class SamAlignment(object):
    def __init__(self, qname='', flag=0, rname='', pos=0, mapq=0, cigar='',
                 rnext='', pnext=0, tlen=0, seq='', qual=''):
        self.qname = qname
        self.flag  = flag
        self.rname = rname
        self.pos   = pos
        self.mapq  = mapq
        self.cigar = cigar
        self.rnext = rnext
        self.pnext = pnext
        self.tlen  = tlen
        self.seq   = seq
        self.qual  = qual

    def set_mandatory_fields(self, array):
        self.qname = array[0]
        self.flag  = array[1]
        self.rname = array[2]
        self.pos   = array[3]
        self.mapq  = array[4]
        self.cigar = array[5]
        self.rnext = array[6]
        self.pnext = array[7]
        self.tlen  = array[8]
        self.seq   = array[9]
        self.qual  = array[10]

    def set_optional_fields(self, array):
        optional_fields = {}
        for field in array:
            key, typ, val = field.split(':', 2)
            optional_fields[key] = val
        self.optional_fields = optional_fields


    ## Class methods

    @classmethod
    def create_from_array(cls, array):
        aln = cls()
        aln.set_mandatory_fields(array[:11])
        aln.set_optional_fields(array[11:])
        return aln

    @classmethod
    def create_from_line(cls, line):
        return cls.create_from_array(line.strip().split())



def parse_sam(sam_file):
    with open(sam_file, 'r') as f:
        for line in f:
            yield SamAlignment.create_from_line(line)
