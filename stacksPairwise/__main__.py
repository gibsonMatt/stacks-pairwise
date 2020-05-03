# /usr/bin/python3

import sys
from Bio import SeqIO
from os import path
import time
import argparse
import logging as log
import numpy as np 
import pandas as pd 
import warnings
from scipy.spatial.distance import squareform

def makePairs(popmap):
    pairs = []
    for x in popmap:
        for y in popmap:
            pairs.append([x,y])
    return(pairs)

def getLocus(name):
    return(name.split("_")[1])

def getName(description):
    return(description.split()[1].split(';')[0][1:len(description.split(';')[0])])

def getChrAndPos(description):
    chrr = description.split(" ")[2][:-1]
    startPos = description.split(" ")[3][:-1]
    return(chrr, startPos)

def getNs(length):
    n = ""
    for i in range(0, length):
        n += 'N'
    return(n)

def pairwise(records, pairs, names):
    seqlength = len(records[0].seq)
    by_indv = {}
    for record in records:
        name = getName(record.description)
        if (name not in list(by_indv.keys())):
            by_indv[name] = [str(record.seq)]
        else:
            by_indv[name].append(str(record.seq))
    ests = []
    sites_total = []
    for pair in pairs:
        p1 = pair[0]
        p2 = pair[1]

        try:
            p1_seq_0 = by_indv[p1][0]
            p1_seq_1 = by_indv[p1][1]
        except KeyError:
            p1_seq_0 = getNs(seqlength)
            p1_seq_1 = getNs(seqlength)
        try:
            p2_seq_0 = by_indv[p2][0]
            p2_seq_1 = by_indv[p2][1]
        except KeyError:
            p2_seq_0 = getNs(seqlength)
            p2_seq_1 = getNs(seqlength)

        diffs, sites = divergence2(p1_seq_0, p1_seq_1, p2_seq_0, p2_seq_1)

        ests.append(diffs)
        sites_total.append(sites)
    return(ests, sites_total)

def divergence(a,b,c,d):

    totalSitesWithData = 0
    div = 0
    for i in range(0, len(a)):
        aa = a[i]
        bb = b[i]
        cc = c[i]
        dd = d[i]

        if (aa != 'N') and (bb != 'N') and (cc != 'N') and (dd != "N"):
            totalSitesWithData += 1
            if (aa == bb) and (cc == dd) and (aa != cc) and (aa != dd):
                #two different homos
                div += 2
            elif (aa != bb) and (cc != dd) and (aa != cc) and (aa != dd) and (bb != cc) and (bb != dd):
                #two different hets
                div += 2
            elif (aa != bb) and (cc != dd) and (aa == cc) and (aa != dd) and (bb != cc) and (bb != dd):
                #Two hets but a and c same
                div += 1
            elif (aa != bb) and (cc != dd) and (aa != cc) and (aa == dd) and (bb != cc) and (bb != dd):
                #Two hets but a and d the same
                div += 1
            elif (aa != bb) and (cc != dd) and (aa != cc) and (aa != dd) and (bb == cc) and (bb != dd):
                #two hets but b and c the same
                div += 1
            elif (aa != bb) and (cc != dd) and (aa != cc) and (aa != dd) and (bb != cc) and (bb == dd):
                #two hets but b and d the same
                div += 1
            elif (aa != bb) and (cc == dd) and (cc in [aa,bb]):
                #indv 1 is het and indv 2 is homo, s
                div += 1
            elif (aa == bb) and (cc != dd) and (aa in [cc,dd]):
                #indv 2 is het and indv 1 is homo
                div += 1
            elif (aa != bb) and (cc == dd) and (cc not in [aa,bb]):
                #indv 1 is het and indv 2 is homo, s
                div += 2
            elif (aa == bb) and (cc != dd) and (aa not in [cc,dd]):
                #indv 2 is het and indv 1 is homo
                div += 2
    if totalSitesWithData == 0:
        return('NA')
    else:
        return(str(div) + '/' + str(totalSitesWithData))

def d1(x,y):
    d = 0
    n = 0
    for i in range(0, len(x)):
        if x[i] != 'N' and y[i] != 'N':
            n += 1
            if x[i] != y[i]:
                d += 1
    
    return(n, d)


def divergence2(a,b,c,d):
    n_ac, d_ac = d1(a,c)
    n_ad, d_ad = d1(a,d)
    n_bc, d_bc = d1(b,c)
    n_bd, d_bd = d1(b,d)

    num =  d_ac + d_ad + d_bc + d_bd
    den = n_ac + n_ad + n_bc + n_bd
    avg_diff = float(num/4)
    avg_sites = float(den/4)
    return(avg_diff, avg_sites)

def writeDistanceMatrix(df1, df2, df3):
    df1_a = df1.drop(['Chromosome', 'Position'], axis = 1)
    a = df1_a.mean(axis = 0)
    a_indices = a.index.values

    unique = []
    for aa in a_indices:
        aa = aa.split('_')
        if aa[0] not in unique:
            unique.append(aa[0])
        if aa[1] not in unique:
            unique.append(aa[1])

    distance_df1 = np.zeros((len(unique), len(unique)))
    distance_df1 = pd.DataFrame(distance_df1, index = unique, columns = unique)

    for i, aa in enumerate(a_indices):
        aa = aa.split('_')
        val = a[i]

        distance_df1.at[aa[0], aa[1]] = val
    

    df2_a = df2.drop(['Chromosome', 'Position'], axis = 1)
    a = df2_a.sum(axis = 0)
    a_indices = a.index.values

    unique = []
    for aa in a_indices:
        aa = aa.split('_')
        if aa[0] not in unique:
            unique.append(aa[0])
        if aa[1] not in unique:
            unique.append(aa[1])

    distance_df2 = np.zeros((len(unique), len(unique)))
    distance_df2 = pd.DataFrame(distance_df2, index = unique, columns = unique)

    for i, aa in enumerate(a_indices):
        aa = aa.split('_')
        val = a[i]

        distance_df2.at[aa[0], aa[1]] = val

    df3_a = df3.drop(['Chromosome', 'Position'], axis = 1)
    a = df3_a.sum(axis = 0)
    a_indices = a.index.values

    unique = []
    for aa in a_indices:
        aa = aa.split('_')
        if aa[0] not in unique:
            unique.append(aa[0])
        if aa[1] not in unique:
            unique.append(aa[1])

    distance_df3 = np.zeros((len(unique), len(unique)))
    distance_df3 = pd.DataFrame(distance_df3, index = unique, columns = unique)

    for i, aa in enumerate(a_indices):
        aa = aa.split('_')
        val = a[i]

        distance_df3.at[aa[0], aa[1]] = val

    distance_df4 = np.divide(distance_df2, distance_df3)

    return(distance_df1, distance_df2, distance_df3, distance_df4)
    

def main(*args):
    start = time.time()
    parser = argparse.ArgumentParser(
        description="Calculate pairwise divergence (pairwise pi) from Stacks `samples.fa` output fle"
    )

    parser.add_argument(
        "-names",
        "--names",
        metavar="",
        help="Names of samples to analyze. Either a text file or comma seperated list.",
        default=None,
    )

    parser.add_argument(
        "samples",
        metavar="samples",
        help="Path to `samples.fa` file (from Stacks output)"
    )

    parser.add_argument("-o", "--outputdir", metavar="", help="Output directory/prefix", default="./stackspairwise")

    args = parser.parse_args()

    # Setup ###################
    log.basicConfig(level=log.DEBUG)
    logger = log.getLogger()
    ##########################



    pops = args.names
    seqs = args.samples
    output = args.outputdir
    
    #Read in names file
    popmap = {}
    try:
        if (path.exists(pops)):
            names = open(pops, 'r')
           
            for i, line in enumerate(names):
                l = line.replace('\n','').split()
                popmap[l[0]] = 'A'
            names.close()
        else:
            names = pops.split(',')
            if len(names) < 2:
                sys.exit("Error in sample names inputs")
                
            else:
                for n in names:
                    popmap[n] = 'A'
    except:
        sys.exit("Error in sample names inputs")
    


    individuals = list(popmap.keys())
    pairs = makePairs(individuals)

    #Make pairs strings
    pairs2 = []
    for p in pairs:
        pairs2.append(p[0] + '_' + p[1])


    #Calculate pairwise
    loci = []
    names = []
    chromosomes = []
    positions = []
    with open(seqs, 'r') as handle:
        focal = []
        prev = ''
        i = 0
        x = 0
        prev_chr = ''
        for record in SeqIO.parse(handle, "fasta"):

            locus = getLocus(record.name)
            if i == 0:
                focal.append(record)
            else:
                if locus == prev:
                    focal.append(record)
                else:
                    chrr, pos = getChrAndPos(record.description)
                    if chrr != prev_chr:
                        print("Reading chromosome/scaffold " + chrr)
                        
                    prev_chr = chrr

                    if (len(focal)/2) > 1:
                        diffs, sites = pairwise(focal, pairs, individuals)
                        diffs = [float(x) for x in diffs]
                        sites = [float(x) for x in sites]
                        if x == 0:
                            diffs_np = np.array(diffs, dtype=np.float64)
                            sites_np = np.array(sites, dtype=np.float64)
                        else:
                            diffs_np = np.vstack([diffs_np, diffs])
                            sites_np = np.vstack([sites_np, sites])
                        
                    else:
                        if x == 0:
                            diffs_np = np.array([np.nan]*len(pairs2), dtype=np.float64)
                            sites_np = np.array([np.nan]*len(pairs2), dtype=np.float64)
                        else:
                            diffs_np = np.vstack([diffs_np, [np.nan]*len(pairs2)])
                            sites_np = np.vstack([sites_np, [np.nan]*len(pairs2)])
                    focal = [record]
                    loci.append(getLocus(record.name))
                    chromosomes.append(chrr)
                    positions.append(pos)
                    x += 1
                
            prev = locus
            i += 1

    shape = diffs_np.shape
    print("Processed " + str(shape[0]) + " loci across " + str(shape[1]) + " pairwise comparisons")

    warnings.filterwarnings('ignore')

    #Write output files
    final_df = np.divide(diffs_np, sites_np)
    np.savetxt('text.csv', final_df)
    pairs3 = [x for x in pairs2]
    df1 = pd.DataFrame(data=final_df, index = loci, columns = pairs2, dtype='Float64')
    df1['Chromosome'] = chromosomes
    df1['Position'] = positions
    pairs2.insert(0, "Position")
    pairs2.insert(0, "Chromosome")
    df1 = df1[pairs2]
    df1.index.names = ['Locus']
    df1.to_csv(output + '.estimates.csv')

    df2 = pd.DataFrame(data=diffs_np, index = loci, columns = pairs3)
    df2['Chromosome'] = chromosomes
    df2['Position'] = positions
    df2 = df2[pairs2]
    df2.index.names = ['Locus']
    df2.to_csv(output + '.diffs.csv')

    df3 = pd.DataFrame(data=sites_np, index = loci, columns = pairs3)
    df3['Chromosome'] = chromosomes
    df3['Position'] = positions
    df3 = df3[pairs2]
    df3.index.names = ['Locus']
    df3.to_csv(output + '.sites.csv')

    dist1, dist2, dist3, dist4 = writeDistanceMatrix(df1, df2, df3)

    matrix_output = open(output + ".summary.txt", 'w')

    matrix_output.write("Average per-locus estimates:\n")
    matrix_output.write(dist1.to_string() + '\n\n')

    matrix_output.write("Sum of per-site distances:\n")
    matrix_output.write(dist2.to_string() + '\n\n')

    matrix_output.write("Total number of sites with data:\n")
    matrix_output.write(dist3.to_string() + '\n\n')

    matrix_output.write("Genome-wide estimates:\n")
    matrix_output.write(dist4.to_string() + '\n\n')
    matrix_output.close()

if __name__ == "__main__":
    main(*sys.argv)
