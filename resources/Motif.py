#! /usr/bin/env python

import re

from dna import reverseComplement
from regex import reFromString

class Motif( object ):
    """
    Represent a possible methylation motif in a genome
    """
    def __init__( self, arg ):
        self._string = arg.upper()
        self._regularExpression = reFromString( self.String )

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
        rc = reverseComplement(sequence)
        fwd = [(m.start(), m.end(), 0) for m in re.finditer(self.regularExpression, sequence)]
        rev = [(length-m.end(), length-m.start(), 1) for m in re.finditer(self.regularExpression, rc)]
        return fwd + rev

    def __len__(self):
        return len(self.String)
