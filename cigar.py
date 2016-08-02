class Cigar(object):
    symbols = ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X']

    def __init__(self, cigar_str=''):
        self.cigar_str = cigar_str
        self.cigar_arr = self.parse_cigar(cigar_str)

    def __eq__(self, other):
        return self.cigar_str == other.cigar_str

    def parse_cigar(self, cigar_s):
        cigar_arr = []
        length = ''
        for c in cigar_s:
            if c in self.__class__.symbols:
                cigar_arr.append((c, int(length)))
                length = ''
            else:
                length += c
        return cigar_arr

    def num_total_bases(self):
        return sum(map(lambda x: x[1], self.cigar_arr))

    def num_match(self):
        return sum(map(lambda x: x[1], filter(lambda x: x[0] == '=', self.cigar_arr)))

    def num_mismatch(self):
        return sum(map(lambda x: x[1], filter(lambda x: x[0] == 'X', self.cigar_arr)))

    def num_insertion(self):
        return sum(map(lambda x: x[1], filter(lambda x: x[0] == 'I', self.cigar_arr)))

    def num_deletion(self):
        return sum(map(lambda x: x[1], filter(lambda x: x[0] == 'D', self.cigar_arr)))

    def num_left_softclip(self):
        return self.cigar_arr[0][1] if self.cigar_arr[0][0] == 'S' else 0

    def num_right_softclip(self):
        return self.cigar_arr[-1][1] if self.cigar_arr[-1][0] == 'S' else 0

    def num_softclip(self):
        return self.num_left_softclip(), self.num_right_softclip()

    def num_query_bases(self):
        return self.num_left_softclip() + self.num_query_bases_in_aln() + self.num_right_softclip()

    def num_query_bases_in_aln(self):
        return self.num_match() + self.num_mismatch() + self.num_insertion()

    def num_target_bases_in_aln(self):
        return self.num_match() + self.num_mismatch() + self.num_deletion()

    def reverse(self):
        reversed_cigar_arr = list(reversed(self.cigar_arr))
        return self.__class__.create_from_array(reversed_cigar_arr)

    def strip_softclip(self):
        array = self.cigar_arr
        if array[0][0] == 'S':
            array = array[1:]
        if array[-1][0] == 'S':
            array = array[:-1]
        return self.__class__.create_from_array(array)

    def aln_bgn(self):
        return self.num_left_softclip() + 1

    def aln_end(self):
        return self.num_total_bases() - self.num_right_softclip()

    def aln_len(self):
        return self.num_total_bases() - self.num_left_softclip() - self.num_right_softclip()

    def aln_identity(self):
        return round(1.0 * self.num_match() / self.aln_len(), 4)

    def slice_left(self, slice_len):
        sliced_cigar = []
        total_len = 0
        for block in self.cigar_arr:
            s, l = block
            if total_len + l >= slice_len:
                sliced_cigar.append((s, slice_len - total_len))
                break
            else:
                sliced_cigar.append(block)
                total_len += l
        return self.__class__.create_from_array(sliced_cigar)

    def slice_right(self, slice_len):
        return self.reverse().slice_left(slice_len).reverse()

    def to_legacy_match(self):
        converted_arr = []
        for symbol, length in self.cigar_arr:
            if symbol in ['=', 'X']:
                converted_arr.append(('M', length))
            else:
                converted_arr.append((symbol, length))

        merged_arr = []
        prev_symbol, sum_length = converted_arr[0]
        for symbol, length in converted_arr[1:]:
            if symbol == prev_symbol:
                sum_length += length
            else:
                merged_arr.append((prev_symbol, sum_length))
                prev_symbol = symbol
                sum_length = length

        merged_arr.append((prev_symbol, sum_length))
        return self.__class__.create_from_array(merged_arr)


    ## Class methods

    @classmethod
    def create_from_array(cls, cigar_arr):
        cigar_str = cls.arr_to_str(cigar_arr)
        return cls(cigar_str)

    @classmethod
    def arr_to_str(cls, cigar_arr):
        return ''.join([str(c[1]) + c[0] for c in cigar_arr])

    @classmethod
    def convert_to_legacy_match(cls, cigar_str):
        return cls(cigar_str).to_legacy_match().cigar_str

