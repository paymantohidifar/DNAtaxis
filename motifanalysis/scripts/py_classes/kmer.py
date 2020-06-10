from .sequence import Sequence
import re as re
from statistics import median


class Kmer(Sequence):

    @staticmethod
    def extract_positions(matchiter):
        return [match.start() for match in matchiter]

    def __init__(self, seq):
        super().__init__(seq)
        self.count = 0
        self.positions = tuple()
        self.spacings = tuple()
        self.meandist = None

    def get_counts(self):
        return self.count

    def get_positions(self):
        return self.positions

    def get_spacings(self):
        return self.spacings

    def get_mean_dist(self):
        return self.meandist

    def set_positions(self, pos):
        self.positions = pos

    def set_counts(self, motifcount):
        self.count = motifcount

    def set_spacings(self, sp):
        self.spacings = sp

    def set_mean_distance(self, md):
        self.meandist = md

    def search(self, chrseqgen, rc=True, cluster=False):
        seq = self.get_sequence()
        pat = re.compile('|'.join([seq, self.reverse_complement(seq)]), re.I) if rc else re.compile(seq, re.I)
        pos_gen = self.find_positions_generator(chrseqgen, pat)
        allpos = (-len(seq),)
        for elt in pos_gen:
            allpos += tuple(pos for pos in elt if (pos - allpos[-1]) >= len(seq))
        self.set_positions(allpos[1:])
        self.set_counts(len(allpos[1:]))

        if cluster:
            allspacing = tuple()
            for n in range(0, self.count-1):
                allspacing += (self.positions[n+1] - self.positions[n] - len(seq),)
            self.set_spacings(allspacing)
            self.set_mean_distance(median(allspacing))

    def find_positions_generator(self, chrseqgen, pat):
        overlap = len(self.get_sequence()) - 1
        chrseq = next(chrseqgen, '')
        if not chrseq:
            return None
        buf = chrseq + ''
        postracker = 0
        while len(buf) != overlap:
            poslst = self.next_positions(buf, pat)
            if poslst:
                yield tuple(pos + postracker for pos in poslst)
            postracker += len(buf) - overlap
            buf = chrseq[-overlap:]
            chrseq = next(chrseqgen, '')
            buf += chrseq

    def next_positions(self, buffer, pat):
        return self.extract_positions(pat.finditer(buffer))