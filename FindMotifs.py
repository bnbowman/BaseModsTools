#! /usr/bin/env python

import sys

from collections import defaultdict

from pbcore.io import *

from resources.utils import createMaskedReferences
from resources.dna import reverseComplement
from resources.regex import reFromString
from resources.Motif import Motif

START_POS = 0
END_POS   = 20000
MAX_COUNT = 100000
MAX_PER_READ = 100

refFn = sys.argv[1]
dataFn = sys.argv[2]
motifStr = sys.argv[3] if len(sys.argv) >= 4 else None
targetMotif = sys.argv[4] if len(sys.argv) >= 5 else None

def _rindex(alist, value):
    try:
        return len(alist) - alist[-1::-1].index(value) - 1
    except:
        raise ValueError

def _qStartToPass( qStarts, qStart ):
    if qStart == 0:
        nPass = "First"
    elif qStarts[0] == 0 and qStarts.index(qStart) == 1:
        nPass = "Second"
    elif qStarts.index(qStart) >= 2:
        nPass = "RC"
    else:
        nPass = "Unknown"
    return nPass

def ReadSubreadStarts( dataFn ):
    starts = defaultdict(list)
    for record in openDataSet(dataFn):
        starts[record.holeNumber].append( record.qStart )
    
    for hn, qS in starts.iteritems():
        starts[hn] = sorted(starts[hn])

    return starts

def PrintUnmaskedIpds( dataFn, refs, maxPosition=END_POS ):
    qStarts = ReadSubreadStarts( dataFn )
    print "Pass,Pos,Base,IPD"
    for record in openDataSet(dataFn):
        if record.tStart > maxPosition:
            continue
        if not record.isForwardStrand:
            continue
        nPass = _qStartToPass( qStarts[record.holeNumber], record.qStart )
        if nPass == "Unknown":
            continue

        regions = refRecs[record.tId].UnmaskedRegions
        ovls = [r for r in regions if r[0] < record.tEnd and r[1] > record.tStart]
        for start, end in ovls:
            clip = record.clippedTo(start, end)
            aStart = max(start, record.tStart)
            readBases = clip.read(aligned=True, orientation="native")
            tplBases = clip.reference(aligned=True, orientation="native")
            ipds = clip.IPD(aligned=True, orientation="native")

            p = aStart
            for r, t, ip in zip(readBases, tplBases, ipds):
                # Count forward our position if it's not an insertion
                if t == "-":
                    continue
                else:
                    p += 1
                # Skip any mismatches or deletions
                if r != t:
                    continue
                if p > maxPosition:
                    break
                print "{0},{1},{2},{3}".format(nPass, p, r, ip)

def PrintMotifIpds( dataFn, refs, motif, maxCount=MAX_COUNT ):
    mRegions = [ref.Find(motif) for ref in refs]
    definedPos = motif.definedPositions

    qStarts = ReadSubreadStarts( dataFn )
    passCounts = defaultdict(int)

    print "Pass,MotifStart,Pos,Base,IPD"
    for record in openDataSet(dataFn):
        nPass = _qStartToPass( qStarts[record.holeNumber], record.qStart )
        if nPass == "Unknown":
            continue
        if passCounts[nPass] > maxCount:
            continue
        if not record.isForwardStrand:
            continue

        ovls = [r for r in mRegions[record.tId] if r[0] >= record.tStart and r[1] <= record.tEnd]
        for start, end, strand in ovls:
            clip = record.clippedTo(start, end)
            readBases = clip.read(aligned=True, orientation="native")
            tplBases = clip.reference(aligned=True, orientation="native")
            ipds = clip.IPD(aligned=True, orientation="native")

            p = 0
            for r, t, ip in zip(readBases, tplBases, ipds):
                # Count forward our position if it's not an insertion
                if t == "-":
                    continue
                else:
                    p += 1
                # Skip any mismatches or deletions
                if r != t:
                    continue
                if p in definedPos:
                    print "{0},{1},{2},{3},{4}".format(nPass, start, p, r, ip)

refRecs = createMaskedReferences( refFn, motifStr )
targetM = Motif(reverseComplement(targetMotif)) if targetMotif else None
if targetM:
    PrintMotifIpds( dataFn, refRecs, targetM )
else:
    PrintUnmaskedIpds( dataFn, refRecs )
