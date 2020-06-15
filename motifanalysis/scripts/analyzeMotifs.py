from py_classes.chromosome import Chromosome
import numpy as np
from itertools import product
import pandas as pd
import scipy.stats
import pylab
import matplotlib.pyplot as plt


def analyze_chromosome(names, seq, size):
    chr_info = []
    for name, taxid in names.items():
        chrobj = Chromosome(refdir + taxid + '.zip', seq) if seq \
            else Chromosome(refdir + taxid + '.zip', '', size)
        chr_info.append((name, chrobj.process_chromosome()))
    return generate_chromosome_df(chr_info)


def analyze_kmers(names, seq, size):
    kmer_info = []
    for name, taxid in names.items():
        chrobj = Chromosome(refdir + taxid + '.zip', seq) if seq \
            else Chromosome(refdir + taxid + '.zip', '', size)
        kmer_info.append((name, chrobj.process_motifs()))
    return generate_kmers_df(kmer_info)


def generate_chromosome_df(chr_info):
    """construct data frame for chromosome sequence info"""
    chrinfodict = {}
    for chr in chr_info:
        chrinfodict.setdefault(chr[0], pd.Series(chr[1], index=['genome length', 'GC content']))
    chrinfo_df = pd.DataFrame(chrinfodict)
    return chrinfo_df


def generate_kmers_df(kmerinfo):
    """constructs data frames for count and average
    distance data of all kmers for each chromosome"""
    idxlst = []
    cntlst = []
    distlst = []
    chr_name, kmers_gen = kmerinfo[0]
    print('Processing: ', chr_name)
    for kmer in kmers_gen:
        idxlst.append(kmer['kmer seq'])
        cntlst.append(kmer['count'])
        distlst.append(kmer['mean distance'])
    cdict = {chr_name: cntlst}
    ddict = {chr_name: distlst}
    print(25 * '-')

    for chr_name, kmers_gen in kmerinfo[1:]:
        print('Processing: ', chr_name)
        cntlst = []
        distlst = []
        for kmer in kmers_gen:
            cntlst.append(kmer['count'])
            distlst.append(kmer['mean distance'])
        cdict.setdefault(chr_name, cntlst)
        ddict.setdefault(chr_name, distlst)
        print(25 * '-')
    count_df = pd.DataFrame(cdict, index=idxlst)
    dist_df = pd.DataFrame(ddict, index=idxlst)
    return count_df, dist_df


def write_df_to_csv(df, filname):
    df.to_csv(resultsdir+filname+'.csv', index=True, header=True)


def generate_corr_df(covar_df, response):
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


def make_corr_trend_plot(corr_df, kmers, kmersize, plottype):
    selected_labels = [itm for itm in list(corr_df.index) if itm.split('/')[0] in kmers or
                       itm.split('/')[1] in kmers]
    selected_xrange = [n for n, elt in enumerate(list(corr_df.index)) if elt in selected_labels]
    colormap = [tuple(int(elt) for elt in seq) for seq in product('01', repeat=3)]
    fig, ax = plt.subplots()
    ax.plot(range(0, len(corr_df.index)), corr_df.get('Pearson\'s r'),
            linewidth=0, marker='o', markersize=5, color='grey', label='All {} mers'.format(kmersize))
    ax.set_title('Pearson\'s correlation coefficient for all {} mers'.format(kmersize))
    ax.set_xlabel('Kmer number')
    ax.set_ylabel('Pearson\'s correlation coefficient (r)')
    ax.set_ylim([-1, 1])
    for n, elt in enumerate(selected_labels):
        ax.plot(selected_xrange[n], corr_df.loc[elt, 'Pearson\'s r'],
                linewidth=0, marker='o', markersize=5, color=colormap[n], label=elt)
    ax.legend(facecolor='white')
    figname = '{}_corr_{}.eps'.format(plottype, kmersize)
    plt.savefig(figdir + figname, dpi=None, facecolor='w', edgecolor='w',
                format='eps', pad_inches=0.1)
    plt.show()


def make_corr_plot(kmers, df, corr_df, response, plottype):
    y = response.to_numpy()
    selected_labels = [itm for itm in list(df.index) if itm.split('/')[0] in kmers or
                       itm.split('/')[1] in kmers]
    for n, label in enumerate(selected_labels):
        xpoints = df.loc[label]
        xrange = np.linspace(min(xpoints), max(xpoints), 100)
        intercept = corr_df.loc[label, 'Intercept']
        slope = corr_df.loc[label, 'Slope']
        pearsonr = corr_df.loc[label, 'Pearson\'s r']
        line = f'Reg. line: r={pearsonr:.2f}'
        fig, ax = plt.subplots()
        ax.plot(xpoints, y, linewidth=0, marker='o', label='Data points', color=(0, 0, 1))
        ax.plot(xrange, intercept + slope * xrange, linewidth=2, label=line, color=(1, 0, 0))
        ax.set_title('Correlation plot for \'{}\''.format(label))
        ax.set_xlabel('Normalized frequency (per 1000 bp)') if plottype == 'count' \
            else ax.set_xlabel('Inverted distance')
        ax.set_ylabel('Chemotaxis response')
        ax.legend(facecolor='white')
        figname = '{}_corr_{}.eps'.format(plottype, kmers[n])
        plt.savefig(figdir + figname, dpi=None, facecolor='w', edgecolor='w', orientation='portrait',
                    papertype=None,
                    format='eps', transparent=False, bbox_inches=None, pad_inches=0.1, metadata=None)
    plt.show()


def main():
    chrnames = {    # plasmid DNA sequences are excluded.
                    # 'strain name':            'taxonomy id'
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
    # Inputs by the user
    kmerseq = 'CACAA'
    kmersize = None
    selected_kmers = [kmerseq] if kmersize else ['CACAA', 'CATAA', 'CCCAA']

    # Run kmer analysis
    # chrdata = analyze_chromosome(chrnames, kmerseq, kmersize)
    # cdata, ddata = analyze_kmers(chrnames, kmerseq, kmersize)

    # Write kmer analysis results as csv-formatted tables in 'results' subdirectory
    # for df, df_name in [(chrdata, 'chr'), (cdata, 'count'), (ddata, 'distance')]:
    #     write_df_to_csv(df, df_name + '_table_{}'.format(kmerseq if kmerseq else kmersize))


    # temp
    chrdata = pd.read_csv(resultsdir + 'chr_table_5.csv', header=0, index_col=0)
    cdata = pd.read_csv(resultsdir + 'count_table_5.csv', header=0, index_col=0)
    ddata = pd.read_csv(resultsdir + 'distance_table_5.csv', header=0, index_col=0)
    print(chrdata)
    # temp

    # Perform correlation analysis for count data
    chrlen = chrdata.loc['genome length']
    norm_cdata = cdata / chrlen * 1000
    count_corr_df = generate_corr_df(norm_cdata, che_response)
    write_df_to_csv(count_corr_df, 'ccorr_{}'.format(kmerseq if kmerseq else kmersize))
    if kmersize:
        make_corr_trend_plot(count_corr_df, selected_kmers, kmersize=5, plottype='count')
    make_corr_plot(selected_kmers, norm_cdata, count_corr_df, che_response, plottype='count')

    # Perform correlation analysis for distance data
    inverted_ddata = 1. / ddata
    dist_corr_df = generate_corr_df(inverted_ddata, che_response)
    write_df_to_csv(dist_corr_df, 'dcorr_{}'.format(kmerseq if kmerseq else kmersize))
    if kmersize:
        make_corr_trend_plot(dist_corr_df, selected_kmers, kmersize=5, plottype='distance')
    make_corr_plot(selected_kmers, inverted_ddata, dist_corr_df, che_response, plottype='distance')

    # update salmonella response data


if __name__ == '__main__':
    refdir = '../reference_data/'
    resultsdir= '../results/'
    figdir = '../results/figures/'
    main()
