from zipfile import ZipFile
import re as re
from itertools import product
from statistics import median
import pandas as pd
import scipy.stats
from numpy import mod
import pylab
import matplotlib.pyplot as plt


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
        for kmerseq in kmerseqlst:
            kmerobj = Kmer(kmerseq)
            chrseq_gen = self.chromosome_seq_generator()
            kmerobj.search(chrseq_gen, cluster=True)
            kmerprop = dict()
            kmerprop.setdefault('kmer seq', '/'.join([kmerseq, self.reverse_complement(kmerseq)]))
            kmerprop.setdefault('count', kmerobj.get_counts())
            kmerprop.setdefault('mean distance', kmerobj.get_mean_dist())
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


def run_analysis(names, seq, size):
    chr_info = []
    kmer_info = []
    for strainname, taxid in names.items():
        chrobj = Chromosome(refdir + taxid + '.zip', seq) if seq \
            else Chromosome(refdir + taxid + '.zip', '', size)
        chr_info.append((strainname, chrobj.process_chromosome()))
        kmer_info.append((strainname, chrobj.process_motifs()))
    return create_data_frames(chr_info, kmer_info)


def create_data_frames(chrinfo, kmerinfo):
    # construct data frame for chromosome sequence properties
    chrinfodict = {}
    for chr in chrinfo:
        chrinfodict.setdefault(chr[0], pd.Series(chr[1], index=['genome length', 'GC content']))
    chrinfo_df = pd.DataFrame(chrinfodict)

    # construct data frames for count and average
    # distance data of each kmer for each chromosome
    cdict = {}
    ddict = {}
    for chr in kmerinfo:
        print('Processing: ', chr[0], 50*'-', sep='\n')
        idxlst = []
        cntlst = []
        distlst = []
        for kmer in chr[1]:
            idxlst.append(kmer['kmer seq'])
            cntlst.append(kmer['count'])
            distlst.append(kmer['mean distance'])
        cdict.setdefault(chr[0], pd.Series(cntlst, index=idxlst))
        ddict.setdefault(chr[0], pd.Series(distlst, index=idxlst))
    count_df = pd.DataFrame(cdict)
    dist_df = pd.DataFrame(ddict)
    return chrinfo_df, count_df, dist_df


def write_df_to_csv(df, filname):
    # compression_opts = dict('method':'zip', 'archive_name':resultsdir+filname+'.csv')
    df.to_csv(resultsdir+filname+'.csv', index=True, header=True)


def generate_corr_coeff_df(covar_df, response):
    y = response.to_numpy()
    idxs = covar_df.index
    corr_data = []
    for idx in idxs:
        x = covar_df.loc[idx].to_numpy()
        lnreg = scipy.stats.linregress(x, y)
        corr_data.append([round(elt, 2) for elt in lnreg])
    corr_df = pd.DataFrame(data=corr_data,
                           columns=['Slope', 'Intercept', 'Pearson\'s r', 'Pearson\'s pvalue', 'Reg. std error'],
                           index=idxs)
    corr_df.sort_values(by=['Pearson\'s r'], inplace=True, ascending=False)
    return corr_df


def main():
    chrnames = {    # plasmid DNA sequences are excluded.
                    # 'strain name':            'tax id'
                    'bacillus subtilis 168':    '224308',
                    'bacillus megaterium':      '545693',
                    'bacillus pumilus':         '1408',
                    'bacillus licheniformis':   '1402',
                    'bacillus halodurans':      '272558',
                    'geobacillus therm':        '634956',
                    'vibrio splendidus':        '29497',
                    'escherichia coli':         '83333',
                    'zymomonas mobilis':        '542',
                    'salmonella enterica':      '99287',
                    'lambda phage':             '10710'
                }
    che_response = pd.Series(data=[1378.22, 1590.78, 1338.22, 708.44, 469.44, 465.33,
                                   560.56, 606.89, 612.00, 550.00, 352.33],
                             index=[name for name in chrnames.keys()])
    kmerseq = ''
    kmersize = 5
    # pdata, cdata, ddata = run_analysis(chrnames, kmerseq, kmersize)
    # write_df_to_csv(pdata, 'pdata_{}'.format(kmerseq if kmerseq else kmersize))
    # write_df_to_csv(cdata, 'cdata_{}'.format(kmerseq if kmerseq else kmersize))
    # write_df_to_csv(ddata, 'ddata_{}'.format(kmerseq if kmerseq else kmersize))

    pdata = pd.read_csv(resultsdir + 'pdata_5.csv', header=0, index_col=0)
    cdata = pd.read_csv(resultsdir + 'cdata_5.csv', header=0, index_col=0)
    ddata = pd.read_csv(resultsdir + 'ddata_5.csv', header=0, index_col=0)

    chrlen = pdata.loc['genome length']
    norm_cdata = cdata / chrlen * 1000
    count_corr_df = generate_corr_coeff_df(norm_cdata, che_response)
    write_df_to_csv(count_corr_df, 'ccorr_{}'.format(kmerseq if kmerseq else kmersize))

    inverted_ddata = 1. / ddata
    dist_corr_df = generate_corr_coeff_df(inverted_ddata, che_response)
    write_df_to_csv(dist_corr_df, 'dcorr_{}'.format(kmerseq if kmerseq else kmersize))

    # update salmonella response data
    # modify run_analysis function using list instead of pd.Series
    # fit distance data to a distribution and calculate mean values

    # scipy.stats.probplot(che_response.to_numpy(), dist="norm", plot=pylab)
    # pylab.show()

    # scipy.stats.probplot(norm_cdata.iloc[100].to_numpy(), dist="norm", plot=pylab)
    # pylab.show()


if __name__ == '__main__':
    refdir = '../reference_data/'
    resultsdir= '../results/'
    main()
