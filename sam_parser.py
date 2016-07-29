from alignment import Alignment
from cigar import Cigar

class SamAlignment(Alignment):
    def __init__(self, line):
        super(self.__class__, self).__init__()
        self.line = line
        fields = line.strip().split()
        self.set_mandatory_fields(fields[:11])
        self.set_optional_fields(fields[11:])
        self.set_cigar()
        self.set_attributes()

    def set_mandatory_fields(self, array):
        fields = {}
        fields['qname'] = array[0]
        fields['flag']  = array[1]
        fields['rname'] = array[2]
        fields['pos']   = array[3]
        fields['mapq']  = array[4]
        fields['cigar'] = array[5]
        fields['rnext'] = array[6]
        fields['pnext'] = array[7]
        fields['tlen']  = array[8]
        fields['seq']   = array[9]
        fields['qual']  = array[10]
        self.mandatory_fields = fields

    def set_optional_fields(self, array):
        optional_fields = {}
        for field in array:
            key, typ, val = field.split(':', 2)
            optional_fields[key] = val
        self.optional_fields = optional_fields

    def set_cigar(self):
        self.cigar = Cigar(self.mandatory_fields['cigar'])

    def set_attributes(self):
        m = self.mandatory_fields
        c = self.cigar

        self.q_name = m['qname']
        self.q_bgn  = c.aln_bgn()
        self.q_end  = self.q_bgn + c.num_query_bases_in_aln() - 1
        self.q_len  = c.num_query_bases()

        self.t_name = m['rname']
        self.t_bgn  = int(m['pos'])
        self.t_end  = self.t_bgn + c.num_target_bases_in_aln() - 1

        self.aln_len = c.aln_len()
        self.aln_identity = c.aln_identity()



def parse_sam(sam_file):
    with open(sam_file, 'r') as f:
        for line in f:
            yield SamAlignment(line)
