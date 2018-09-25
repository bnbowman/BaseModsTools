#! /usr/bin/env python

import sys

from collections import defaultdict

# 99-perctile cutoffs range from ~270-330, so this applies
#  an appoximate
MAX_IPD = 300

csvStyle = sys.argv[1]
csvFn = sys.argv[2]

def SummarizeMotifIpds( csvFn, maxIpd=MAX_IPD ):
    data = defaultdict(list)
    with open(csvFn) as handle:
        handle.next()
        for line in handle:
            parts = line.strip().split(',')
            if "G" in parts:
                ipd = int(parts[-1])
                if ipd > maxIpd:
                    continue
                key = (parts[0], parts[1], parts[2])
                data[key].append(ipd)

    print "Pass,MotifPos,Pos,Count,IPD"
    for (nPass, mPos, pos), ipds in data.iteritems():
        ct = len(ipds)
        avg = sum(ipds) / ct
        print "{0},{1},{2},{3},{4}".format(nPass, mPos, pos, ct, avg)

def SummarizeNonMotifIpds( csvFn, maxIpd=MAX_IPD ):
    data = defaultdict(list)
    with open(csvFn) as handle:
        handle.next()
        for line in handle:
            parts = line.strip().split(',')
            if "G" in parts:
                ipd = int(parts[-1])
                if ipd > maxIpd:
                    continue
                key = (parts[0], parts[1])
                data[key].append(ipd)

    print "Pass,Pos,Count,IPD"
    for (nPass, pos), ipds in data.iteritems():
        ct = len(ipds)
        avg = sum(ipds) / ct
        print "{0},{1},{2},{3}".format(nPass, pos, ct, avg)

if csvStyle.lower() == "motif":
    SummarizeMotifIpds( csvFn )
else:
    SummarizeNonMotifIpds( csvFn )
