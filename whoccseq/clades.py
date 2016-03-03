# -*- Python -*-
# license
# license.

"""
Detects clade by aa sequence.
"""

import os, sys, json
import logging; module_logger = logging.getLogger(__name__)

# ----------------------------------------------------------------------

def clades_of(virus_type, lineage, sequence, shift):
    r = None
    if sequence is not None and shift is not None:
        if virus_type == "B":
            r = clades_B(lineage=lineage, sequence=sequence, shift=shift)
        elif virus_type == "A(H1N1)":
            r = clades_H1(lineage=lineage, sequence=sequence, shift=shift)
        elif virus_type == "A(H3N2)":
            r = clades_H3(lineage=lineage, sequence=sequence, shift=shift)
    return r

# ----------------------------------------------------------------------

def clades_B(lineage, sequence, shift):
    r = None
    if lineage == "YAMAGATA":
        r = clades_B_yamagata(sequence=sequence, shift=shift)
    elif lineage == "VICTORIA":
        r = clades_B_victoria(sequence=sequence, shift=shift)
    elif lineage:
        module_logger.warning('Unsupported B lineage {}'.format(lineage))
    return r

# ----------------------------------------------------------------------

def clades_B_yamagata(sequence, shift):
    # 165N -> Y2, 165Y -> Y3
    r = None
    pos = 164 - shift
    if len(sequence) > pos and pos > 0:
        if sequence[pos] == "N":
            r = ["Y2"]
        elif sequence[pos] == "Y":
            r = ["Y3"]
    return r

# ----------------------------------------------------------------------

def clades_B_victoria(sequence, shift):
    r = None
    return r

# ----------------------------------------------------------------------

def clades_H1(lineage, sequence, shift):
    r = []
    # 84N+162N+216T - 6B.1, 152T+173I+501E - 6B.2
    seq_len = len(sequence)
    pos84 = 83 - shift
    if seq_len > pos84 and pos84 > 0 and sequence[pos84] == "N":
        pos162 = 161 - shift
        pos216 = 215 - shift
        if seq_len > pos216 and sequence[pos162] == "N" and sequence[pos216] == "T":
            r.append("6B1")
        # else:
        #     r.append("6B.2")

    pos152 = 151 - shift
    pos173 = 172 - shift
    pos501 = 500 - shift
    if seq_len > pos501 and sequence[pos152] == "T" and sequence[pos173] == "I" and sequence[pos501] == "E":
        r.append("6B2")
    return r or None

# ----------------------------------------------------------------------

def clades_H3(lineage, sequence, shift):
    r = []
    seq_len = len(sequence)

    # 158N, 159F -> 3C3, 159Y -> 3c2a, 159S -> 3c3a, 62K+83R+261Q -> 3C3b.
    pos158 = 157 - shift
    pos159 = 158 - shift
    if seq_len > pos159 and pos158 > 0 and sequence[pos158] == "N":
        if sequence[pos159] == "F":
            r.append("3C3")
        elif sequence[pos159] == "Y":
            r.append("3C2a")
        elif sequence[pos159] == "S":
            r.append("3C3a")

    pos62 = 61 - shift
    pos83 = 82 - shift
    pos261 = 260 - shift
    if seq_len > pos261 and sequence[pos159] == "F" and sequence[pos83] == "R" and sequence[pos261] == "Q": # and sequence[pos62] == "K"
        r.append("3C3b")
    # if seq_len > pos261 and sequence[pos62] == "K" and sequence[pos83] == "R" and sequence[pos261] == "Q":
    #     r.append("3C3b?")

    # 160S -> gly, 160T -> gly, 160x -> no gly
    pos = 159 - shift
    if seq_len > pos and pos > 0 and sequence[pos] in "ST":
        r.append("gly")
    else:
        r.append("no-gly")
    return r or None

# ----------------------------------------------------------------------

# ======================================================================
### Local Variables:
### eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
### End:
