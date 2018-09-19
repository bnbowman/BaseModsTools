#! /usr/bin/env python

import string

COMPLEMENT = string.maketrans("ACGTMRWSYKVHDBN-", "TGCAKYWSRMBDHVN-")

def reverseComplement( string ):
    return string.upper().translate(COMPLEMENT)[::-1]

