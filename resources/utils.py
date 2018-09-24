#! /usr/bin/env python

from pbcore.io import FastaReader

from MaskedReference import MaskedReference
from Motif import Motif

def createMaskedReferences( refFn, motifStr ):
    # Read input files
    refRecs = [MaskedReference(rec) for rec in FastaReader(refFn)]
    maskedMs = [Motif(m) for m in motifStr.split(',')] if motifStr else []

    # Mask repeat regions
    for motif in maskedMs:
        for rec in refRecs:
            rec.MaskRepeat( motif )

    return refRecs
