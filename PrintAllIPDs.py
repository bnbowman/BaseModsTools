#! /usr/bin/env python

import sys

from pbcore.io import openDataSet

from collections import defaultdict

MAX_POS = 100000

FWD = defaultdict(list)
REV = defaultdict(list)

fns      = sys.argv[1:]

for fn in fns:
    ds = openDataSet(fn)
    for rec in ds:
        if rec.tStart >= MAX_POS:
            continue

        ipd   = rec.IPD(aligned=True, orientation="native")
        trans = rec.transcript(orientation="native", style="gusfield")
        pos   = rec.referencePositions(aligned=True, orientation="native")

        if rec.isForwardStrand:
            for (ip, t, p) in zip(ipd, trans, pos):
                if t != "M" or p > MAX_POS:
                    continue
                FWD[p].append(ip)
        else:
            for (ip, t, p) in zip(ipd, trans, pos):
                if t != "M" or p > MAX_POS:
                    continue
                REV[p].append(ip)

pos = list(set(FWD.keys() + REV.keys()))

print "Pos,Strand,IPDs"
for p in pos:
    print '{},FWD,"{}"'.format(p, ",".join(str(ip) for ip in sorted(FWD[p])))
    print '{},REV,"{}"'.format(p, ",".join(str(ip) for ip in sorted(REV[p])))
