##########
#Matthew JS Gibson
#Indiana University, Bloomington, IN, USA
#2019

#Calculate pairwise divergence across RAD loci generated in Stacks.
##########


#input1 = popmap
#input2 = samples.fa
#input3 = output name


import sys
from Bio import SeqIO

def makePairs(individuals):
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
        ests.append(divergence(p1_seq_0, p1_seq_1, p2_seq_0, p2_seq_1))
    return(ests)

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



pops = sys.argv[1]
seqs = sys.argv[2]
output = sys.argv[3]


pops = open(pops, 'r')

popmap = {}
for i, line in enumerate(pops):
    l = line.replace('\n','').split()
    popmap[l[0]] = l[1]
pops.close()

individuals = list(popmap.keys())
pairs = makePairs(individuals)


output = open(output, 'w')

pairs2 = []
for p in pairs:
    pairs2.append(p[0] + '_' + p[1])

output.write("LocusID\tChr\tStartPos\t" + '\t'.join(pairs2) + '\n')
loci = []
with open(seqs, 'r') as handle:
    focal = []
    prev = ''
    i = 0
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
                    print(chrr)
                prev_chr = chrr

                if (len(focal)/2) > 1:
                    row = pairwise(focal, pairs, individuals)
                    
                    output.write(getLocus(record.name) + '\t' + chrr + '\t' + pos + '\t' + '\t'.join(row) + '\n')
                else:
                    output.write(getLocus(record.name) + '\t' + '\t' + chrr + '\t' + pos + '\n')




                focal = [record]
        prev = locus
        i += 1
output.close()