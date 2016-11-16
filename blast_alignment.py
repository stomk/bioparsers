class BlastAlignment:
    def __init__(self, line, outfmt):
        self.raw = line
        self.set_fields(outfmt)
        
    def set_fields(self, outfmt=6):
        fields = self.raw.strip().split()
        if outfmt == 6:            
            self.q_name   = fields[0]
            self.t_name   = fields[1]
            self.identity = float(fields[2])
            self.length   = int(fields[3])
            self.mismatch = int(fields[4])
            self.gap_open = int(fields[5])
            self.q_bgn    = int(fields[6])
            self.q_end    = int(fields[7])
            self.t_bgn    = int(fields[8])
            self.t_end    = int(fields[9])
            self.evalue   = float(fields[10])
            self.bitscore = float(fields[11])
        else:
            raise ValueError('Currently only outfmt=6 is supported')
