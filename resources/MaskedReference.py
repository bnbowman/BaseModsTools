#! /usr/bin/env python

import re

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
        for start, end, _ in motif.Find( self.rawSequence ):
            self._regionsNeedUpdate = True
            for p in range(start, end):
                self.mask[p] = 1

    def Find(self, motif):
        if self._regionsNeedUpdate:
            self._updateRegions()
        unmasked = []
        for mS, mE, strand in motif.Find( self.rawSequence ):
            for rS, rE in self.UnmaskedRegions:
                if rS >= mE:
                    break
                if mS >= rS and mE <= rE:
                    unmasked.append( (mS, mE, strand) )
                    break
        return unmasked
