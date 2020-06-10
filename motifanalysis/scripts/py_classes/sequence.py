class Sequence:

    # acceptable DNA base pairs
    ValidChars = 'TCAGtcag' + 'YMRWSK' + '[]'
    # construction a complement table to use with str.translate to generate complement sequence
    ComplementTable = str.maketrans('TCAGtcag[]', 'AGTCagtc][')

    @classmethod
    def invalid_char(cls, seqstring):
        return {chr for chr in seqstring if chr not in cls.ValidChars}

    def __init__(self, seqstring):
        invalid = self.invalid_char(seqstring)
        if invalid:
            raise Exception(type(self).__name__ +
                            'Sequence contains one or more invalid characters:'
                            + str(invalid))
        self.__seq = seqstring

    def get_sequence(self):
        return self.__seq

    def complement(self, seq):
        return seq.translate(self.ComplementTable)

    def reverse_complement(self, seq):
        return self.complement(seq)[::-1]
