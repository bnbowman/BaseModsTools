#! /usr/bin/env python

import sys
import re
import string

from collections import defaultdict, Counter
from itertools import groupby

from pbcore.io import *

PRECISION = 6

COMPLEMENT = string.maketrans("ACGTMRWSYKVHDBN-", "TGCAKYWSRMBDHVN-")

refFn       = sys.argv[1]
dataFn      = sys.argv[2]
motifStr    = sys.argv[3]
posStr      = sys.argv[4]
excludedStr = sys.argv[5] if len(sys.argv) > 5 else None

#print "RefFn", refFn
#print "DataFn", dataFn
#print "motifStr", motifStr
#print "posStr", posStr
#print "excluded", excludedStr

def _decodeIpd( ip ):
    if ip <= 64:
        return ip
    elif ip > 64 and ip < 128:
        return 64 + (ip-64) * 2
    elif ip >= 128 and ip < 192:
        return 192 + (ip-128) * 4
    else:
        return 448 + (ip-192) * 8

def _reverseComplement( string ):
    return string.upper().translate(COMPLEMENT)[::-1]

def _rindex(alist, value):
    try:
        return len(alist) - alist[-1::-1].index(value) - 1
    except:
        raise ValueError

def _string_to_re( string ):
    # Base  Meaning  Compl.
    #  M   A or C        K
    #  R   A or G        Y
    #  W   A or T        W
    #  S   C or G        S
    #  Y   C or T        R
    #  K   G or T        M
    #  V   A or C or G   B
    #  H   A or C or T   D
    #  D   A or G or T   H
    #  B   C or G or T   V
    groups = []
    for base, group in groupby( string ):
        ct = len(list(group))
        if base in "ACGT":
            groups.append( base * ct )
        else:
            # Convert an ambiguous base to it's components
            baseStr = base
            if base == "N":
                baseStr = "[ACGT]"
            elif base == "M":
                baseStr = "[AC]"
            elif base == "R":
                baseStr = "[AG]"
            elif base == "W":
                baseStr = "[AT]"
            elif base == "S":
                baseStr = "[CG]"
            elif base == "Y":
                baseStr = "[CT]"
            elif base == "K":
                baseStr = "[GT]"
            elif base == "V":
                baseStr = "[ACG]"
            elif base == "H":
                baseStr = "[ACT]"
            elif base == "D":
                baseStr = "[AGT]"
            elif base == "B":
                baseStr = "[CGT]"
            # Convert the ambiguous base components to an RE
            if ct == 1:
                groups.append( baseStr )
            else: 
                groups.append( baseStr + "{" + str(ct) + "}" )
    return "".join(groups)

def _qStartToPass( qStarts, qStart ):
    if qStart == 0:
        nPass = "First"
    elif qStarts.index(qStart) == 1:
        nPass = "Second"
    elif qStarts.index(qStart) >= 2:
        nPass = "Third+"
    else:
        nPass = "Unknown"
    return nPass

class Motif( object ):
    """
    Represent a possible methylation motif in a genome
    """
    def __init__( self, arg ):
        self._string = arg.upper()
        self._regularExpression = _string_to_re( self.String )

    @property
    def String(self):
        return self._string

    @property
    def regularExpression(self):
        return self._regularExpression
    
    @property
    def definedPositions(self):
        return [p+1 for p in range(len(self.String)) if self.String[p] != "N"]

    def Find(self, sequence):
        length = len(sequence)
        rc = _reverseComplement(sequence)
        fwd = [(m.start(), m.end(), 0) for m in re.finditer(self.regularExpression, sequence)]
        rev = [(length-m.end(), length-m.start(), 1) for m in re.finditer(self.regularExpression, rc)]
        return fwd + rev

    def __len__(self):
        return len(self.String)


class MaskedReference( object ):
    """
    Represent a possibly masked reference
    """

    def __init__( self, rec ):
        self._rec = rec
        self._mask = ["0"] * len(rec.sequence)
        self._unmaskedRegions = []
        self._regionsNeedUpdate = True

    @property
    def id(self):
        return self._rec.id

    @property
    def header(self):
        return self._rec.header
    
    @property
    def sequence(self):
        return None

    @property
    def rawSequence(self):
        return self._rec.sequence

    @property
    def mask(self):
        return self._mask
    
    @property
    def maskString(self):
        return "".join(self._mask)
    
    @property
    def UnmaskedRegions(self):
        if self._regionsNeedUpdate:
            self._updateRegions()
        return self._unmaskedRegions

    @property
    def UnmaskedLength(self):
        if self._regionsNeedUpdate:
            self._updateRegions()
        return sum(e-s for s,e in self._unmaskedRegions)

    def __len__(self):
        return len(self.mask)

    def _updateRegions(self):
        self._unmaskedRegions = [(m.start(), m.end()) for m in re.finditer("0+", self.maskString)]
        self._regionsNeedUpdate = False

    def MaskRepeat(self, motif):
        for start, end, _ in motif.Find( rec.rawSequence ):
            self._regionsNeedUpdate = True
            for p in range(start, end):
                self.mask[p] = "1"

    def Find(self, motif):
        if self._regionsNeedUpdate:
            self._updateRegions()
        unmasked = []
        for mS, mE, strand in motif.Find( rec.rawSequence ):
            for rS, rE in self.UnmaskedRegions:
                if rS >= mE:
                    break
                if mS >= rS and mE <= rE:
                    unmasked.append( (mS, mE, strand) )
                    #if strand == 1:
                    #    print mS, mE, strand, rec.rawSequence[mS:mE]
                    break
        return unmasked

def ReadGffFile( fn ):
    data = {}
    with open(fn) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            #modType = parts[2]
            modStrand = 0 if parts[6] == "+" else 1
            modPos = int(parts[3]) - 1                                          # Adjust because GFF is 1-indexed
            modQv = int(parts[-1].split('=')[-1]) if "Qv=" in parts[-1] else 0  # Default QV is Zero
            data[(modPos, modStrand)] = modQv
    return data

def FilterCalls( calls, maskedRefs ):
    newCalls = {}
    for ref in maskedRefs:
        for rS, rE in ref.UnmaskedRegions:
            for p in range(rS, rE):
                # For each unmasked position P, copy any call on either strand
                try:
                    c = calls[(p, 0)]
                    newCalls[(p, 0)] = c
                except:
                    pass
                try:
                    c = calls[(p, 1)]
                    newCalls[(p, 1)] = c
                except:
                    pass
    return newCalls

# Read input files
refRecs   = [MaskedReference(rec) for rec in list(FastaReader(refFn))]
calls     = ReadGffFile( dataFn )
motifs    = [Motif(m) for m in motifStr.split(',')]
positions = [int(p)-1 for p in posStr.split(',')]
excluded  = [Motif(m) for m in excludedStr.split(',')] if excludedStr else []

#print "RefFn", refRecs
#print "DataFn", len(calls), sorted(calls.keys())[:10]
#print "motifStr", motifs
#print "posStr", positions
#print "excluded", excluded

# First we mask any exlcuded regions and tabulate how much is left
for m in excluded:
    for rec in refRecs:
        rec.MaskRepeat( m )
totalUnmasked = 2*sum(r.UnmaskedLength for r in refRecs)

# Second we filter down to only the calls in unmasked regions
filteredCalls = FilterCalls(calls, refRecs)

# Third, we identify our TP and FNs from the filtered call list
posCalls = set()
nFail = 0
for m, p in zip(motifs, positions):
    currPos, currFail = 0, 0
    for rec in refRecs:
        for mB, mE, mS in sorted(rec.Find( m )):
            if mS == 0:
                try:
                    c = filteredCalls[(mB+p, mS)]
                    posCalls.add((mB+p, mS))
                    currPos += 1
                    #print mS, mB, p, mB+p, c, refRecs[0].rawSequence[mB:mE]
                except:
                    currFail += 1
            if mS == 1:
                try:
                    c = filteredCalls[(mE-p-1, mS)]
                    posCalls.add((mE-p-1, mS))
                    currPos += 1
                    #print mS, mE, (-1*p)-1, mE-p-1, c, refRecs[0].rawSequence[mB:mE]
                except:
                    currFail += 1
    print m.String, p, currPos, currFail, currPos / float(currPos + currFail)
    nFail += currFail

# Fourth, we identify all of our true-negative calls
negCalls = []
for key in filteredCalls.iterkeys():
    if key in posCalls:
        continue
    else:
        negCalls.append( key )

# Fifth, we convert our pos/neg calls into Counters, for rapid threshold testing
posCounts = Counter()
for key in posCalls:
    posCounts[filteredCalls[key]] += 1
## Treat FNs as QV=0 hits
posCounts[-1] += nFail

negCounts = Counter()
for key in negCalls:
    negCounts[filteredCalls[key]] += 1
## Treat uncalled TN bases as QV=-1, to differentiate them from called bases
negCounts[-1] += totalUnmasked - len(filteredCalls) - nFail

qvRange = sorted(list(set(posCounts.keys() + negCounts.keys())))

print
print totalUnmasked
print len(calls)
print len(filteredCalls)
print len(posCalls)
print len(negCalls)
print nFail
print 
print posCounts
print negCounts
print
print qvRange

print "MinQv,TP,FN,FP,TN,Sensitivity,Specificity,Accuracy"
for minQv in qvRange:
    TP = sum(v for k,v in posCounts.iteritems() if k >= minQv)
    FN = sum(v for k,v in posCounts.iteritems() if k < minQv)
    FP = sum(v for k,v in negCounts.iteritems() if k >= minQv)
    TN = sum(v for k,v in negCounts.iteritems() if k < minQv)
    sens = round(TP / float(TP + FN), PRECISION)
    spec = round(TN / float(TN + FP), PRECISION)
    acc  = round((TP + TN) / float(TP + FN + FP + TN), PRECISION)
    print "{0},{1},{2},{3},{4},{5},{6},{7}".format(minQv, TP, FN, FP, TN, sens, spec, acc)
