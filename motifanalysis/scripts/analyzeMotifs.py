from py_classes.chromosome import Chromosome
import pandas as pd
import scipy.stats
from numpy import mod
import pylab
import matplotlib.pyplot as plt


def analyze_kmers(names, seq, size):
    kmer_info = []
    for name, taxid in names.items():
        chrobj = Chromosome(refdir + taxid + '.zip', seq) if seq \
            else Chromosome(refdir + taxid + '.zip', '', size)
        kmer_info.append((name, chrobj.process_motifs()))
    return generate_kmers_df(kmer_info)


def analyze_chromosome(names, seq, size):
    chr_info = []
    for name, taxid in names.items():
        chrobj = Chromosome(refdir + taxid + '.zip', seq) if seq \
            else Chromosome(refdir + taxid + '.zip', '', size)
        chr_info.append((name, chrobj.process_chromosome()))
    return generate_chromosome_df(chr_info)


# def generate_kmers_df(kmerinfo):
#     """constructs data frames for count and average
#     distance data of all kmers for each chromosome"""
#     cdict = {}
#     ddict = {}
#     for chr_name, kmers_gen in kmerinfo:
#         print('Processing: ', chr_name, 50 * '-', sep='\n')
#         idxlst = []
#         cntlst = []
#         distlst = []
#         for kmer in kmers_gen:
#             idxlst.append(kmer['kmer seq'])
#             cntlst.append(kmer['count'])
#             distlst.append(kmer['mean distance'])
#         cdict.setdefault(chr_name, pd.Series(cntlst, index=idxlst))
#         ddict.setdefault(chr_name, pd.Series(distlst, index=idxlst))
#     count_df = pd.DataFrame(cdict)
#     dist_df = pd.DataFrame(ddict)
#     return count_df, dist_df


def generate_kmers_df(kmerinfo):
    """constructs data frames for count and average
    distance data of all kmers for each chromosome"""
    idxlst = []
    cntlst = []
    distlst = []
    chr_name, kmers_gen = kmerinfo[0]
    print('Processing: ', chr_name, 50 * '-', sep='\n')
    for kmer in kmers_gen:
        idxlst.append(kmer['kmer seq'])
        cntlst.append(kmer['count'])
        distlst.append(kmer['mean distance'])
    cdict = {chr_name: cntlst}
    ddict = {chr_name: distlst}

    for chr_name, kmers_gen in kmerinfo[1:]:
        print('Processing: ', chr_name, 50 * '-', sep='\n')
        cntlst = []
        distlst = []
        for kmer in kmers_gen:
            cntlst.append(kmer['count'])
            distlst.append(kmer['mean distance'])
        cdict.setdefault(chr_name, cntlst)
        ddict.setdefault(chr_name, distlst)
    count_df = pd.DataFrame(cdict, index=idxlst)
    dist_df = pd.DataFrame(ddict, index=idxlst)
    return count_df, dist_df


def generate_chromosome_df(chr_info):
    """construct data frame for chromosome sequence info"""
    chrinfodict = {}
    for chr in chr_info:
        chrinfodict.setdefault(chr[0], pd.Series(chr[1], index=['genome length', 'GC content']))
    chrinfo_df = pd.DataFrame(chrinfodict)
    return chrinfo_df


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
    # Inputs by user
    kmerseq = 'CACAA'
    kmersize = None

    # Analysis
    chrdata = analyze_chromosome(chrnames, kmerseq, kmersize)
    cdata, ddata = analyze_kmers(chrnames, kmerseq, kmersize)

    # Write data as csv-formatted tables in 'results' subdirectory
    for df, df_name in [(chrdata, 'chr'), (cdata, 'count'), (ddata, 'distance')]:
        write_df_to_csv(df, df_name + '_table_{}'.format(kmerseq if kmerseq else kmersize))

    # pdata = pd.read_csv(resultsdir + 'pdata_5.csv', header=0, index_col=0)
    # cdata = pd.read_csv(resultsdir + 'cdata_5.csv', header=0, index_col=0)
    # ddata = pd.read_csv(resultsdir + 'ddata_5.csv', header=0, index_col=0)
    #
    # chrlen = pdata.loc['genome length']
    # norm_cdata = cdata / chrlen * 1000
    # count_corr_df = generate_corr_coeff_df(norm_cdata, che_response)
    # write_df_to_csv(count_corr_df, 'ccorr_{}'.format(kmerseq if kmerseq else kmersize))
    #
    # inverted_ddata = 1. / ddata
    # dist_corr_df = generate_corr_coeff_df(inverted_ddata, che_response)
    # write_df_to_csv(dist_corr_df, 'dcorr_{}'.format(kmerseq if kmerseq else kmersize))

    # update salmonella response data
    # modify run_analysis function using list instead of pd.Series
    # fit distance data to a distribution and calculate mean values

    ## check normality of kmer data
    # scipy.stats.probplot(che_response.to_numpy(), dist="norm", plot=pylab)
    # pylab.show()

    # scipy.stats.probplot(norm_cdata.iloc[100].to_numpy(), dist="norm", plot=pylab)
    # pylab.show()


if __name__ == '__main__':
    refdir = '../reference_data/'
    resultsdir= '../results/'
    main()
