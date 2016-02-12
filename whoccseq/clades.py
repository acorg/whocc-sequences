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
    return None

# ----------------------------------------------------------------------

def clades_H3(lineage, sequence, shift):
    return None

# ----------------------------------------------------------------------

# ======================================================================
### Local Variables:
### eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
### End:
