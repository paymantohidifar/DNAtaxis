from .sequence import Sequence
from .kmer import Kmer
import re as re
from itertools import product
from zipfile import ZipFile
from numpy import mod


class Chromosome(Sequence):

    def __init__(self, filename, kmerseq, kmersize=None):
        super().__init__(kmerseq)
        self.filename = filename
        self.src = None
        self.info = {'acc': [], 'strain': []}
        self.kmersize = kmersize

    def set_chr_info(self, header):
        info, gnomtyp = header.split(',') if ',' in header else (header, '')
        acc, ignore, strain = info.partition(' ')
        self.info['acc'].append(acc.strip('>'))
        self.info['strain'].append(strain)

    def get_chr_info(self):
        return self.info

    def process_motifs(self):
        kmerseqlst = [self.get_sequence()] if self.get_sequence() else self.create_kmer_seqs(self.kmersize)
        len_kmerseqlst = len(kmerseqlst)
        counter = 1
        for kmerseq in kmerseqlst:
            kmerobj = Kmer(kmerseq)
            chrseq_gen = self.chromosome_seq_generator()
            kmerobj.search(chrseq_gen, cluster=True)
            kmerprop = dict()
            kmerprop.setdefault('kmer seq', '/'.join([kmerseq, self.reverse_complement(kmerseq)]))
            kmerprop.setdefault('count', kmerobj.get_counts())
            kmerprop.setdefault('mean distance', kmerobj.get_mean_dist())
            if mod(counter, round(len_kmerseqlst/10)) == 0:
                print('{} % completed ...'.format(round(counter/len_kmerseqlst*100)))
            counter += 1
            yield kmerprop

    def chromosome_seq_generator(self):
        self.extract_zip_file()
        buffer = self.read_next_chunk()
        while buffer:
            yield self.scan_buffer(buffer)
            buffer = self.read_next_chunk()

    def read_next_chunk(self):
        """Reads small chunks of chromosome to
        avoid overwhelming the memory"""
        lncnt, seq = 0, ''
        while lncnt < 50:
            line = self.src.readline().decode('UTF-8')
            if not line:
                return seq
            seq += line
            lncnt += 1
        return seq

    def scan_buffer(self, seqstring):
        pat = re.compile(r'^>.+?$', re.I | re.M)
        matchobj = pat.search(seqstring, 0)
        if not matchobj:
            return seqstring.replace('\n', '')
        seq = seqstring[0:matchobj.start()]
        nxtmatchobj = pat.search(seqstring, matchobj.end()+1)
        while nxtmatchobj:
            if self.parse_header(matchobj.group()):
                start = matchobj.end() + 1
                end = nxtmatchobj.start()
                seq += seqstring[start:end]
            matchobj = nxtmatchobj
            nxtmatchobj = pat.search(seqstring, matchobj.end()+1)
        if self.parse_header(matchobj.group()):
            seq += seqstring[matchobj.end()+1:]

        # check sequence contained allowed nucleotides
        seq = seq.replace('\n', '')
        invalid = self.invalid_char(seq)
        if invalid:
            raise Exception('Sequence contains one or more invalid characters:' + str(invalid))
        else:
            return seq

    def parse_header(self, header):
        keywrds = ['chromosome', 'complete genome', 'complete sequence']
        if any([elt in header.lower() for elt in keywrds]) and 'plasmid' not in header.lower():
            self.set_chr_info(header)
            return header
        else:
            return None

    def create_kmer_seqs(self, size):
        kmers = [''.join(seq) for seq in product('ACGT', repeat=size)]
        kmerscopy = kmers
        for kmer in kmerscopy:
            rc = self.reverse_complement(kmer)
            if rc in kmers:
                kmers.remove(rc)
        return kmers

    def extract_zip_file(self):
        with ZipFile(self.filename) as zipfile:
            fnafilename = [filename for filename in zipfile.namelist()
                           if filename.endswith('.fna')]
            self.src = zipfile.open(fnafilename[0], 'r')

    def process_chromosome(self):
        chrseq_gen = self.chromosome_seq_generator()
        chrseq = next(chrseq_gen, '')
        gc_content = 0
        chr_len = 0
        while chrseq:
            gc_content += chrseq.upper().count('C') + chrseq.upper().count('G')
            chr_len += len(chrseq)
            chrseq = next(chrseq_gen, '')
        if not chr_len:
            raise Exception(ValueError, 'Parser did not find the right genome sequence!')
        return [chr_len, round(gc_content/chr_len*100, 2)]
