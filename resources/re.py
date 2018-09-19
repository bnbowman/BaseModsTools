#! /usr/bin/env python

def refromString( string ):
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
