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

def _reverseComplement( string ):
    return string.upper().translate(COMPLEMENT)[::-1]

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
        self._mask = [0] * len(rec.sequence)
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
        return ''.join("N" if m==1 else b for b,m in zip(self.rawSequence, self._mask))

    @property
    def rawSequence(self):
        return self._rec.sequence

    @property
    def mask(self):
        return self._mask
    
    @property
    def maskString(self):
        return "".join(str(m) for m in self._mask)
    
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
                self.mask[p] = 1

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
                    break
        return unmasked

def ReadGffFile( fn ):
    baseToModType = {'A': 'm6A', 'C': 'm4C', 'G': 'Unknown', 'T': 'Unknown'}

    data = {}
    with open(fn) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            modType = parts[2]
            modStrand = 0 if parts[6] == "+" else 1
            modNotes = parts[-1].split(';')
            modCtx = modNotes[1].split('=')[-1]
            modBase = modCtx[20]
            modType = baseToModType[modBase] if modType == "modified_base" else modType
            modPos = int(parts[3]) - 1                                             # Adjust because GFF is 1-indexed
            modQv = int(modNotes[-1].split('=')[-1]) if "Qv=" in parts[-1] else 0  # Default QV is Zero
            data[(modPos, modStrand)] = (modType, modQv)
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

# First we mask any exlcuded regions and tabulate how much is left
for m in excluded:
    for rec in refRecs:
        rec.MaskRepeat( m )

# Second we filter down to only the calls in unmasked regions
filteredCalls = FilterCalls(calls, refRecs)

# Third, we identify our TP and FNs from the filtered call list
posCalls = defaultdict(set)
nFails   = defaultdict(int)
for m, p in zip(motifs, positions):
    if m.String[p] == "A":
        mType = "m6A"
    elif m.String[p] == "C":
        mType = "m4C"
    else:
        mType = "Unknown"

    currPos, currFail = 0, 0
    for rec in refRecs:
        for mB, mE, mS in sorted(rec.Find( m )):
            if mS == 0:
                try:
                    c = filteredCalls[(mB+p, mS)]
                    posCalls[mType].add((mB+p, mS))
                    currPos += 1
                except:
                    currFail += 1
            if mS == 1:
                try:
                    c = filteredCalls[(mE-p-1, mS)]
                    posCalls[mType].add((mE-p-1, mS))
                    currPos += 1
                except:
                    currFail += 1
    print m.String, p, mType, currPos, currFail, currPos / float(currPos + currFail)
    nFails[mType] += currFail

# Fourth, we identify all of our false-positive calls
negCalls = defaultdict(list)
for key, call in filteredCalls.iteritems():
    isPositive = False
    for mType, callSet in posCalls.iteritems():
        if key in posCalls:
            isPositive = True
            break
    if not isPositive:
        mType, qv = call
        negCalls[mType].append( key )

posCounts = defaultdict(Counter)
for mType, callSet in posCalls.iteritems():
    for key in callSet:
        _, qv = filteredCalls[key]
        posCounts[mType][qv] += 1
    ## Treat  as QV=-1 hits
    posCounts[mType][-1] += nFails[mType]

nBases = Counter(refRecs[0].sequence)
negCounts = defaultdict(Counter)
for mType, callList in negCalls.iteritems():
    for key in callList:
        _, qv = filteredCalls[key]
        negCounts[mType][qv] += 1
    ## Treat uncalled bases as QV=-1 calls
    if mType == "m6A":
        negCounts[mType][-1] += nBases['A'] - len(posCalls[mType]) - len(callList)
    elif mType == "m4C":
        negCounts[mType][-1] += nBases['C'] - len(posCalls[mType]) - len(callList)
    # NOTE: Figure out how to handle Unknown bases here
    #else:
    #    negCounts[mType][-1] += n

print "ModType,MinQv,TP,FN,FP,TN,Sensitivity,Specificity,Accuracy"
for mType in posCounts.iterkeys():
    qvRange = sorted(list(set(posCounts[mType].keys() + negCounts[mType].keys())))
    for minQv in qvRange:
        TP = sum(v for k,v in posCounts[mType].iteritems() if k >= minQv)
        FN = sum(v for k,v in posCounts[mType].iteritems() if k < minQv)
        FP = sum(v for k,v in negCounts[mType].iteritems() if k >= minQv)
        TN = sum(v for k,v in negCounts[mType].iteritems() if k < minQv)
        sens = round(TP / float(TP + FN), PRECISION)
        spec = round(TN / float(TN + FP), PRECISION)
        acc  = round((TP + TN) / float(TP + FN + FP + TN), PRECISION)
        print "{0},{1},{2},{3},{4},{5},{6},{7},{8}".format(mType, minQv, TP, FN, FP, TN, sens, spec, acc)
