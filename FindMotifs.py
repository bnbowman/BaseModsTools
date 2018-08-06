#! /usr/bin/env python

import sys
import re
import string

from collections import defaultdict
from itertools import groupby

from pbcore.io import *

MAX_COUNT = 100000
MAX_PER_READ = 100

COMPLEMENT = string.maketrans("ACGTMRWSYKVHDBN-", "TGCAKYWSRMBDHVN-")

refFn = sys.argv[1]
dataFn = sys.argv[2]
motifStr = sys.argv[3] if len(sys.argv) >= 4 else None
targetMotif = sys.argv[4] if len(sys.argv) >= 5 else None

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
    for base, group in groupby( string):
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
        rc = _reverseComplement(sequence)
        fwd = [(m.start(), m.end(), 0) for m in re.finditer(self.regularExpression, sequence)]
        rev = [(m.start(), m.end(), 1) for m in re.finditer(self.regularExpression, rc)]
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

    def __len__(self):
        return len(self.mask)

    def _updateRegions(self):
        self._unmaskedRegions = [(m.start(), m.end()) for m in re.finditer("0+", self.maskString)]
        self._regionsNeedUpdate = False

    def MaskRepeat(self, motif):
        for start, end, strand in motif.Find( rec.rawSequence ):
            self._regionsNeedUpdate = True
            if strand == 1:
                s = start
                start = len(self) - end
                end = len(self) - s
            for p in range(start, end):
                self.mask[p] = "1"

    def Find(self, motif):
        if self._regionsNeedUpdate:
            self._updateRegions()
        unmasked = []
        for mS, mE, strand in motif.Find( rec.rawSequence ):
            if strand == 1:
                continue
            for rS, rE in self.UnmaskedRegions:
                if rS >= mE:
                    break
                if mS >= rS and mE <= rE:
                    unmasked.append( (mS, mE, strand) )
                    break
        return unmasked

def ReadSubreadStarts( dataFn ):
    starts = defaultdict(list)
    for record in openDataSet(dataFn):
        starts[record.holeNumber].append( record.qStart )
    
    for hn, qS in starts.iteritems():
        starts[hn] = sorted(starts[hn])

    return starts

def PrintUnmaskedIpds( dataFn, refs, maxCount=MAX_COUNT, maxPerRead=MAX_PER_READ ):
    print "Pass,IPD"
    count = 0
    currCount = 0

    qStarts = ReadSubreadStarts( dataFn )
    passCounts = defaultdict(int)

    for record in openDataSet(dataFn):
        nPass = _qStartToPass( qStarts[record.holeNumber], record.qStart )
        if passCounts[nPass] > MAX_COUNT:
            continue

        currCount = 0
        regions = refRecs[record.tId].UnmaskedRegions
        ovls = [r for r in regions if r[0] < record.tEnd and r[1] > record.tStart]
        for start, end in ovls:
            clip = record.clippedTo(start, end)
            for ip in clip.IPD(aligned=True, orientation="native"):
                if ip == 65535:
                    continue
                if nPass == "First" and currCount > MAX_PER_READ:
                    continue
                print "{0},{1}".format(nPass, ip)
                currCount += 1
        passCounts[nPass] += currCount

def PrintMotifIpds( dataFn, refs, motif, maxMotifs=MAX_COUNT ):
    mRegions = [ref.Find(motif) for ref in refs]
    definedPos = motif.definedPositions

    print "Pass,MotifStart,Pos,Base,IPD"
    count = 0

    qStarts = ReadSubreadStarts( dataFn )
    passCounts = defaultdict(int)

    for record in openDataSet(dataFn):
        nPass = _qStartToPass( qStarts[record.holeNumber], record.qStart )
        if passCounts[nPass] > MAX_COUNT:
            continue
        if not record.isForwardStrand:
            continue

        ovls = [r for r in mRegions[record.tId] if r[0] < record.tEnd and r[1] > record.tStart]
        for start, end, strand in ovls:
            clip = record.clippedTo(start, end)
            refBases = clip.reference(aligned=True, orientation="native")
            ipds = clip.IPD(aligned=True, orientation="native")

            p = 0
            for b, ip in zip(refBases, ipds):
                if b == "-":
                    continue
                else:
                    p += 1
                if ip == 65535:
                    continue
                if p in definedPos:
                    print "{0},{1},{2},{3},{4}".format(nPass, start, p, b, ip)
            count += 1

def PrintMotifIpdSums( dataFn, refs, motif, maxMotifs=MAX_COUNT, window=8):
    paddedM = Motif("N" * window + motif.String + "N" * window)
    mRegions = [ref.Find(paddedM) for ref in refs]
    definedPos = motif.definedPositions
    
    qStarts = ReadSubreadStarts( dataFn )

    windowStrs = [str(w+1) for w in range(window)]
    print "MotifStart,nPass,{0},MotifOnly,{1}".format(",".join("n" + w for w in windowStrs) , ",".join("p" + w for w in windowStrs))
    count = 0

    motifStartIdx = window + 1
    motifEndIdx = window + len(motif) + 1
    for record in openDataSet(dataFn):
        #if count > MAX_COUNT:
        #    break
        if not record.isForwardStrand:
            continue

        if record.qStart == 0:
            nPass = "First"
        elif qStarts[record.holeNumber].index(record.qStart) == 1:
            nPass = "Second"
        elif qStarts[record.holeNumber].index(record.qStart) >= 2:
            nPass = "Third+"
        else:
            nPass = "Unknown"

        ovls = [r for r in mRegions[record.tId] if r[0] < record.tEnd and r[1] > record.tStart]
        for start, end, strand in ovls:
            clip = record.clippedTo(start, end)
            refBases = clip.reference(aligned=True, orientation="native")
    
            ipds = [ip if ip != 65535 else 0 for ip in clip.IPD(aligned=True, orientation="native")]
            ipds = [min(_decodeIpd(ip), 952) for ip in ipds]

            p = 0
            pos = []
            for b, ip in zip(refBases, ipds):
                if b == "-":
                    pos.append( p )
                    continue
                else:
                    p += 1
                pos.append( p )

            # If there are alignment issues such that we can't call the begin/end
            #  of the motif, then we skip this read/motif pair
            try:
                motifStart = pos.index(motifStartIdx)
                motifEnd = _rindex(pos, motifEndIdx)
            except ValueError:
                continue

            sums = []
            # Left-side sums
            for i in range(1, window+1):
                try:
                    idx = pos.index(i)
                    sums.append( sum(ipds[idx:motifEnd+1]) )
                except:
                    if len(sums):
                        sums.append( sums[-1] )
                    else:
                        sums.append( 0 )
            # Motif-only sum 
            #print motifStart, motifEnd
            sums.append( sum(ipds[motifStart:motifEnd+1]) )
            # Right-side sums
            for i in range(window + len(motif) + 1, max(pos) + 1):
                try:
                    idx = _rindex(pos, i)
                    sums.append( sum(ipds[motifStart:idx+1]) )
                except:
                    sums.append( sums[-1] )

            print "{0},{1},{2}".format(start + window, nPass, ",".join(str(s) for s in sums))

# Read input files
refRecs = [MaskedReference(rec) for rec in FastaReader(refFn)]
maskedMs = [Motif(m) for m in motifStr.split(',')] if motifStr else []
targetM = Motif(_reverseComplement(targetMotif)) if targetMotif else None

# Mask repeat regions
for motif in maskedMs:
    for rec in refRecs:
        rec.MaskRepeat( motif )

if targetM:
    PrintMotifIpds( dataFn, refRecs, targetM )
    #PrintMotifIpdSums( dataFn, refRecs, targetM )
else:
    PrintUnmaskedIpds( dataFn, refRecs )
