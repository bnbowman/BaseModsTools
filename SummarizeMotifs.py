#! /usr/bin/env python

import sys

from collections import defaultdict

PRECISION = 6

genome = sys.argv[1]
fns = sys.argv[2:]

if genome.lower().startswith("cagg") or genome.lower().startswith("c_agg"):
    motifsDict = {"m6A": ["GATC", "CATG", "GANTC", "CATCC", "CGGAG", "CCCASG", "TTCGAA", "TCCGGAHD"],
                  "m4C": ["CACCC", "GGWCC", "RGATCY", "GCCNNNNNGGC", "CGGCCG", "TGGCCA"]}
elif genome.lower().startswith("rpal") or genome.lower().startswith("r_pal"):
    motifsDict = {"m6A": ["GTYGGAG", "CTGCAG", "GANTC"],
                  "m4C": ["GCCCGCC"]}
elif genome.lower().startswith("daceto") or genome.lower().startswith("d_aceto"):
    motifsDict = {"m6A": ["GACCCA", "CGGAAG", "CTGGAA", "CRGAGG", "GACGA", "GCCATC", "GGGCGA", "GAAYNNNNNGGG"],
                  "m4C": ["CCCNNNNNRTTC", "GTAC", "CGCG", "CTAG", "CCACCC", "ARGCCC"]}
elif genome.lower().startswith("mjann") or genome.lower().startswith("m_jann"):
    motifsDict = {"m6A": ["CSATC", "CCANNNNNNNTTG", "CAANNNNNNNTGG", "YACNNNNNTGG", "CCANNNNNGTR", "GAYNNNNNGTAA", "TTACNNNNNRTC", "TAGNNNNNNTGC"],
                  "m4C": ["GGNCCD", "GTACR"]}
elif genome.lower().startswith("saurus") or genome.lower().startswith("s_aureus"):
    motifsDict = {"m6A": ["AGGNNNNNGAT", "ATCNNNNNCCT", "ACANNNNNNRTGG", "CCAYNNNNNNTGT"],
                  "m4C": []}
elif genome.lower().startswith("bsub") or genome.lower().startswith("b_sub"):
    motifsDict = {"m6A": ["CNCANNNNNNNRTGT", "ACAYNNNNNNNTGNG"],
                  "m4C": []}
elif genome.lower().startswith("ecoli") or genome.lower().startswith("e_coli"):
    motifsDict = {"m6A": ["AACNNNNNNGTGC", "GCACNNNNNNGTT", "GATC"],
                  "m4C": []}
else:
    raise ValueError("Unrecognized genome - no known motifs for '{0}'".format(genome))

def CleanAndSplitLine( string ):
    return [s[1:-1] for s in string.strip().split(',')]

## Example
## "motifString","centerPos","modificationType","fraction","nDetected","nGenome","groupTag","partnerMotifString","meanScore","meanIpdRatio","meanCoverage","objectiveScore"
## "GATC","1","m6A","0.99701947","38134","38248","GATC","GATC","65.82299","5.777036","35.60458","2503359.8"

print "Genome,Movie,TrainingSet,TrainingRep,Coverage,Modification,ModificationType,FracPositive"
for fn in fns:
    fileParts  = fn.split('.')
    movie      = fileParts[0]
    #shortMovie = "_".join(movie.split('_')[:3])
    shortMovie = movie
    coverage  = int(fn.split('.')[-3][1:])

    if len(fileParts) == 6:
        barcode     = fileParts[1].split("-")[0]
        trainingSet = fileParts[2]
        trainingRep = "5k"
    #if len(fileParts) == 6:
    #    trainingSet = fileParts[1]
    #    trainingRep = fileParts[2]
    else:
        trainingSet = "N/A"
        trainingRep = "N/A"

    pos   = defaultdict(int)
    total = defaultdict(int)
    with open(fn) as handle:
        handle.next()
        for line in handle:
            motif, _, mod, _, mPos, mTotal = CleanAndSplitLine( line )[:6]
            if mod not in ["m6A", "m4C"]:
                continue
            if motif not in motifsDict[mod]:
                continue
            pos[mod]   += int(mPos)
            total[mod] += int(mTotal)
            print "{0},{1}_{2},{3},{4},{5},{6},{7},{8}".format(genome, shortMovie, barcode, trainingSet, trainingRep, coverage, motif, mod, round(int(mPos) / float(mTotal), PRECISION))
            #print "{0},{1},{2},{3},{4},{5},{6},{7}".format(genome, shortMovie, trainingSet, trainingRep, coverage, motif, mod, round(int(mPos) / float(mTotal), PRECISION))
    for mod, modTotal in total.iteritems():
        modPos = pos[mod]
        print "{0},{1}_{2},{3},{4},{5},WeightedAvg,{6},{7}".format(genome, shortMovie, barcode, trainingSet, trainingRep, coverage, mod, round(modPos / float(modTotal), PRECISION))
        #print "{0},{1},{2},{3},{4},WeightedAvg,{5},{6}".format(genome, shortMovie, trainingSet, trainingRep, coverage, mod, round(modPos / float(modTotal), PRECISION))
