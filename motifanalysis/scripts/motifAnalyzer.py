import numpy as np
import itertools
import matplotlib.pyplot as plt
from Bio.Blast import NCBIWWW

''' modules'''

def mismatch(seq, p):
    """
    Generates sequences with p mismatches in input sequence.

    :param seq:
    :param p:
    :return:
    """

    nucleotides = 'ATCG'
    for locs in itertools.combinations(range(len(seq)), p):
        new_seq = [[nucs] for nucs in seq]
        for loc in locs:
            orig_nuc = seq[loc]
            new_seq[loc] = [l for l in nucleotides if l != orig_nuc]
        for i in itertools.product(*new_seq):
            yield ''.join(i)

def complement(s):

    """
    This function returns the complement of the sequence.

    :param s:
    :return:
    """

    basecomplement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'D': 'D', 'R': 'R', 'S': 'S', 'M': 'M', 'W': 'W',
                      'H': 'H', 'B': 'B', 'N': 'N', 'Y': 'Y', 'K': 'K', 'V': 'V'}
    letters = list(s)
    letters = [basecomplement[base] for base in letters]
    return ''.join(letters)


def reverse_complement(motiflist):

    """
    It reverse-complements sequence.

    :param motiflist:
    :return:
    """

    new_list = []
    for s in motiflist:
        new_list.append(complement(s[::-1]))
    return new_list


def find_hits(gs, mot):

    """
    Returns locations of input sequence in genome sequence (gs)

    :param gs:
    :param mot:
    :return:
    """

    loc_list = []
    found = gs.find(mot)
    loc = found + 1
    loc_list.append(loc)
    while found > -1:
        found = gs.find(mot, found + 1)
        if found == -1:
            break
        loc = found + 1
        loc_list.append(loc)
    return loc_list


def merge(listA, listB):

    """
    Merges and sorts locations of motif sequences in both strands of genome sequence.

    :param listA:
    :param listB:
    :return:
    """

    tot = listA + listB
    tot.sort(key=float)
    return tot


def find_cluster(locs_list, gl, ml, wl):

    """
    Returns total number of motifs within sequence window with user-specified length of wl
    Note: Overlap betweem each adajacent window is wl/2

    :param locs_list:
    :param gl:  genome sequence length
    :param ml:  motif length
    :param wl:  sequence window length
    :return:    total number of motifs within sequence window
    """

    cluster = []
    win_min = 0
    win_max = wl
    jump_len = 3 * wl // 4
    region_counter = 0
    while win_max <= gl:
        counter = 0
        for loc in locs_list:
            if loc >= win_min and loc < win_max - ml:
                counter += 1
            elif loc >= (win_max - ml):
                break
        cluster.append([counter, win_min])
        win_min += jump_len
        win_max += jump_len
        region_counter += 1
    return cluster, region_counter


def find_inert_frag(inert_frag_size, inert_frag_num, motif_list, gs):

    """
    Returns location of inert DNA fragments with zero number of motifs

    :param inert_frag_size:
    :param inert_frag_num:
    :param motif_list:
    :param gs:
    :return:
    """

    gl = len(gs)
    locs_list = []
    counter = 0
    win_min = 0
    win_max = inert_frag_size
    while win_max <= gl and counter < inert_frag_num:
        is_inert = True
        frag_seq = gs[win_win:win_max]
        for motif in motif_list:
            if motif in frag_seq:
                is_inert = False
                break
        if is_inert:
            locs_List.append(win_min)
            counter += 1
        win_min += inert_frag_size
        win_max += inert_frag_size
        return locs_list


def find_enriched_region(cluster_data, gs, enriched_frag_size, num_motif, undesired_motif_list):

    """
    Returns starting position of regions containing motifs with min number of threshold value

    :param cluster_data:
    :param gs:
    :param enriched_frag_size:
    :param num_motif:
    :param undesired_motif_list:
    :return:
    """

    enriched_locs = []
    for data in cluster_data:
        no_undesired_motif = True
        if data[0] == num_motif:
            for undesired_motif in undesired_motif_list:
                if undesired_motif in gs[data[1]:data[1] + enriched_frag_size]:
                    no_undesired_motif = False
                    break
            if no_undesired_motif:
                enriched_locs.append(data[1])
    return enriched_locs


def find_spacing_hist(all_locs, ml):
    """
    Returns distribution of spacing between adjacent motifs with maximum spacing of 500 bp.

    :param all_locs:
    :param ml:
    :return:
    """

    bins = range(0, 60, 1)
    spacing = []
    for i in range(0, len(all_locs) - 1):
        spacing.append(all_locs[i + 1] - all_locs[i] - ml)
    hist, bin_edges = np.histogram(spacing, bins, density=False)
    cum_hist = np.cumsum(hist)
    return hist, cum_hist, bins


''' module __main__'''

help(mismatch)

''' Chromosomes'''
chromosome = 'B.subtilis168'
# chromosome = 'B.pumilus'
# chromosome = 'B.megateriumQM'
# chromosome = 'B.licheniformisATCC14580'
# chromosome = 'B.haloduransC125'
# chromosome = 'G.thermoglucosidasiusC56'
# chromosome = 'P.aeruginosaPAO1'
# chromosome  = 'E.coliK12MG1655'
# chromosome  = 'Z.mobilis'
# chromosome = 'V.splendidus12B01'
# chromosome = 'V.splendidus13B01'
# chromosome = 'L.phage'


''' Print chromosome's information'''
print(chromosome)
fileIn = open('../genomes/' + chromosome, 'r')
genomeSeq = fileIn.readline()
fileIn.close()
genomeLen = len(genomeSeq)
print('Genome length: ', genomeLen)
print('-----------------------------')

''' Print motif's information'''
motifList = ['CACAA', 'CATAA']
motifLen = len(motifList[0])
print('Motifs length: ', motifLen)
motifRCList = reverse_complement(motifList)
print('Motifs: ', motifList)
print('Reverse complement of Motifs: ', motifRCList)

'''Return the locations of motifs in both strands of chromosome'''
locs = []
for motif in motifList:
    locs += find_hits(genomeSeq, motif)

locsRC = []
for motifRC in motifRCList:
    locsRC += find_hits(genomeSeq, motifRC)

locsAll = merge(locs, locsRC)

''' Print results to the console '''
print('Number of motifs in (+) strand: ', len(locs))
print('Number of motifs in (-) strand: ', len(locsRC))
print('Number of motifs in (+/-) strands: ', len(locsAll))


''' Show spacing distributions '''
# spacingHist, cumSpacing, bins = find_spacing_hist(locsAll, motifLen)
# #plt.plot(bins[:-1],spacingHist/float(len(locsAll)))
# plt.figure('Spacing distribution')
# plt.plot(bins[:-1],spacingHist)
# plt.figure('Cumulative distribution of spacing')
# plt.plot(bins[:-1],cumSpacing)
# plt.show()

''' Returns number of motifs in every sub-sequence '''
regionWindow = 500
clusterData, regionNum = find_cluster(locsAll, genomeLen, motifLen, regionWindow)
bins = range(1, 20)
counts = []
for data in clusterData:
    counts.append(data[0])
regionHist, bin_edges = np.histogram(counts, bins, density=False)
print(bin_edges)
print(bins)
print(regionHist)
# print(regionNum)
# plt.bar(bins[:-1], regionHist, align='center')
# plt.show()

''' Return enriched fragments with certain number of motifs in each fragments '''
# undesiredMotifList = ['GACAA','CGCAA','CATAA','CCCAA','CACCA','CTCAA','CACGA']
# undesiredMotifList = ['CATAA']
# undesiredMotifRevComList = []
# for undesiredMotif in undesiredMotifList:
#     undesiredMotifRevComList.append(reverse_complement(undesiredMotif))

# undesiredMotifTotalList = undesiredMotifList + undesiredMotifRevComList

undesiredMotifTotalList = []
enrichedLocs = find_enriched_region(clusterData, genomeSeq, regionWindow, 10, undesiredMotifTotalList)

print(enrichedLocs)
print(len(enrichedLocs))
for seqLoc in enrichedLocs:
    print('>', seqLoc)
    print(genomeSeq[seqLoc:seqLoc + regionWindow])

# efsFileOut = open('enrichedFragmentSeqs9CACAA.txt','w')
# for enrichedLoc in enrichedLocs:
# 	efsFileOut.write('>'+chromosome+'  '+str(enrichedLoc)+' : '+str(enrichedLoc+sizeFrag)+'\n')
# 	efsFileOut.write(genomeSeq[enrichedLoc:enrichedLoc+sizeFrag]+'\n')
# efsFileOut.close()
#
#
# ## Return the inert fragments (zero motifs) sequences
# numInertFrags = 10
# sizeInertFrags = 400
# numMismatch = 1
# #motifVariantList   = list(mismatch(motif,numMismatch))
# motifVariantList = ['GACAA','CGCAA','CATAA','CCCAA','CACCA','CTCAA','CACGA']
# motifVariantList.append(motif)
# motifVariantRevComList = []
# for motifVariant in motifVariantList:
# 	motifVariantRevComList.append(reverse_complement(motifVariant))
#
# motifVariantTotalList = motifVariantList + motifVariantRevComList
#
# print motifVariantList
# print motifVariantRevComList
# print motifVariantTotalList
#
#
# inertFragLocs = find_inert_frag(sizeInertFrags, numInertFrags, motifVariantTotalList, genomeSeq)
# print inertFragLocs
#
# ifsFileOut = open('inertFragSeqs.txt','w')
# for inertFragLoc in inertFragLocs:
# 	ifsFileOut.write('>'+chromosome+'  '+str(inertFragLoc)+' : '+str(inertFragLoc+sizeInertFrags)+'\n')
# 	ifsFileOut.write(genomeSeq[inertFragLoc:inertFragLoc+sizeInertFrags]+'\n')
# ifsFileOut.close()
