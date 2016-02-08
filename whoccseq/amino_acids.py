# -*- Python -*-
# license
# license.

import sys, os, re, csv, pprint
import logging; module_logger = logging.getLogger(__name__)

# ----------------------------------------------------------------------

CODON_TO_PROTEIN = {
    'UGC': 'C', 'GTA': 'V', 'GTG': 'V', 'CCT': 'P', 'CUG': 'L', 'AGG': 'R', 'CTT': 'L', 'CUU': 'L',
    'CTG': 'L', 'GCU': 'A', 'CCG': 'P', 'AUG': 'M', 'GGC': 'G', 'UUA': 'L', 'GAG': 'E', 'UGG': 'W',
    'UUU': 'F', 'UUG': 'L', 'ACU': 'T', 'TTA': 'L', 'AAT': 'N', 'CGU': 'R', 'CCA': 'P', 'GCC': 'A',
    'GCG': 'A', 'TTG': 'L', 'CAT': 'H', 'AAC': 'N', 'GCA': 'A', 'GAU': 'D', 'UAU': 'Y', 'CAC': 'H',
    'AUA': 'I', 'GUC': 'V', 'TCG': 'S', 'GGG': 'G', 'AGC': 'S', 'CTA': 'L', 'GCT': 'A', 'CCC': 'P',
    'ACC': 'T', 'GAT': 'D', 'TCC': 'S', 'UAC': 'Y', 'CAU': 'H', 'UCG': 'S', 'CAA': 'Q', 'UCC': 'S',
    'AGU': 'S', 'TTT': 'F', 'ACA': 'T', 'ACG': 'T', 'CGC': 'R', 'TGT': 'C', 'CAG': 'Q', 'GUA': 'V',
    'GGU': 'G', 'AAG': 'K', 'AGA': 'R', 'ATA': 'I', 'TAT': 'Y', 'UCU': 'S', 'TCA': 'S', 'GAA': 'E',
    'AGT': 'S', 'TCT': 'S', 'ACT': 'T', 'CGA': 'R', 'GGT': 'G', 'TGC': 'C', 'UGU': 'C', 'CUC': 'L',
    'GAC': 'D', 'UUC': 'F', 'GTC': 'V', 'ATT': 'I', 'TAC': 'Y', 'CUA': 'L', 'TTC': 'F', 'GTT': 'V',
    'UCA': 'S', 'AUC': 'I', 'GGA': 'G', 'GUG': 'V', 'GUU': 'V', 'AUU': 'I', 'CGT': 'R', 'CCU': 'P',
    'ATG': 'M', 'AAA': 'K', 'TGG': 'W', 'CGG': 'R', 'AAU': 'N', 'CTC': 'L', 'ATC': 'I',
    # stops
    'TAA': '*', 'UAA': '*', 'TAG': '*', 'UAG': '*', 'TGA': '*', 'UGA': '*', 'TAR': '*', 'TRA': '*', 'UAR': '*', 'URA': '*',
    }

def translate_sequence_to_amino_acid(sequence, min_offset=0, max_offset=2):

    def translate(sequence):
        r = "".join(CODON_TO_PROTEIN.get(codon, "X") for codon in (sequence[i:i+3] for i in range(0, len(sequence), 3)))
        stop = r.find("*")
        if stop >= 0:
            r = r[:stop]
        return r

    best = ""
    for offset in range(min_offset, max_offset + 1):
        trans = translate(sequence[offset:])
        if len(trans) > len(best):
            best = trans
    return best

# ----------------------------------------------------------------------

def align_to_start(sequences, subtype=None):
    """Aligns aa sequences by trimming initial letters if necessary.
    sequences is a dict name -> sequence.
    returns tuple: (aligned sequences as dict, dict of not aligned sequences)
    passed sequences dict is not modified.
    """
    aligner = SequencesAligner()
    aligner.align(sequences, subtype)
    aligner.trim_to_start()
    return (aligner.good, aligner.not_aligned)

# ----------------------------------------------------------------------

class SequenceIsTooShort (Exception): pass

# ----------------------------------------------------------------------

def align(sequence_aa):
    """Returns dict: {"subtype": <detected-flu-subtype>, "lineage": <B, H1 lineage>, "gene": <detected gene>, "shift": <alignment shift>}
    Raises: SequenceIsTooShort
    """
    global MINIMUM_SEQUENCE_AA_LENGTH
    if len(sequence_aa) < MINIMUM_SEQUENCE_AA_LENGTH:
        raise SequenceIsTooShort("sequence is too short: {}".format(len(sequence_aa)))
    ali = aligner()
    r = ali.match(sequence_aa)
    return r

# ----------------------------------------------------------------------

# http://signalpeptide.com/
ALIGNMENT_RAW_DATA = [
    ("MKTIIALCYILCLVFA",  "A(H3N2)", None, "signalpeptide", "HA"),
    ("MKTIIALSHIFCLVLG",  "A(H3N2)", None, "signalpeptide", "HA"),
    ("MKTIIALSYIFCLAFA",  "A(H3N2)", None, "signalpeptide", "HA"),
    ("MKTIIALSYIFCLAFG",  "A(H3N2)", None, "signalpeptide", "HA"),
    ("MKTIIALSYIFCLALG",  "A(H3N2)", None, "signalpeptide", "HA"),
    ("MKTIIALSYIFCLVFA",  "A(H3N2)", None, "signalpeptide", "HA"),
    ("MKTIIALSYIFCLVLG",  "A(H3N2)", None, "signalpeptide", "HA"),
    ("MKTIIALSYIFCQVFA",  "A(H3N2)", None, "signalpeptide", "HA"),
    ("MKTIIALSYIFCQVLA",  "A(H3N2)", None, "signalpeptide", "HA"),
    ("MKTIIALSYILCLVFA",  "A(H3N2)", None, "signalpeptide", "HA"),
    ("MKTIIALSYISCLVFA",  "A(H3N2)", None, "signalpeptide", "HA"),
    ("MKTIIVLSCFFCLAFS",  "A(H3N2)", None, "signalpeptide", "HA"),
    ("MKTLIALSYIFCLVLG",  "A(H3N2)", None, "signalpeptide", "HA"),
    ("MKTTTILILLTHWVHS",  "A(H3N2)", None, "signalpeptide", "HA"),
    ("ATLCLGHHAV",        "A(H3N2)", None, 10,              "HA"),
    ("TNATELVQ",          "A(H3N2)", None, 36,              "HA"),
    ("VERSKAYSN",         "A(H3N2)", None, 87,              "HA"),

    ("MKVKLLVLLCTFTATYA", "A(H1N1)", None,   "signalpeptide", "HA"),
    ("MKVKLLVLLCTFSATYA", "A(H1N1)", "seas", "signalpeptide", "HA"),
    ("MKAILVVLLYTFATANA", "A(H1N1)", "pdm",  "signalpeptide", "HA"),

    ("MKAIIVLLMVVTSNA", "B", None, "signalpeptide", "HA"),               # http://repository.kulib.kyoto-u.ac.jp/dspace/bitstream/2433/49327/1/8_1.pdf
    ("MVVTSNA",         "B", None, "signalpeptide", "HA"),
]

MINIMUM_SEQUENCE_AA_LENGTH = 200          # actually H3 3C3b clade requires 261Q

class Aligner:

    def __init__(self):
        global ALIGNMENT_RAW_DATA
        self.data = {}          # sequence-to-match -> {"virus_type":, "lineage":, "shift":}
        # check for duplicates
        for match_sequence, virus_type, lineage, shift_data, gene  in ALIGNMENT_RAW_DATA:
            if match_sequence in self.data:
                raise ValueError("match_sequence duplicate: {}".format(match_sequence))
            if shift_data == "signalpeptide":
                shift = - len(match_sequence)
            elif isinstance(shift_data, int) and shift_data > 0:   # infix
                shift = shift_data
            else:
                raise ValueError("Unrecognized shift data {!r} for {!r}".format(shift_data, match_sequence))
            self.data[match_sequence] = {"virus_type": virus_type, "shift": shift, "gene": gene}
            if lineage:
                self.data[match_sequence]["lineage"] = lineage

    def match(self, sequence):
        for matcher, matching_data in self.data.items():
            offset = sequence.find(matcher)
            if offset >= 0:
                r = {k: v for k, v in matching_data.items()}
                r["shift"] -= offset
                break
        else:
            r = None
        return r

# ----------------------------------------------------------------------

sAligner = None

def aligner():
    global sAligner
    if not sAligner:
        sAligner = Aligner()
    return sAligner

# ----------------------------------------------------------------------

class SequencesAligner:

    SIGNAL_PEPTIDE = {
        "H3": [                               # all from http://signalpeptide.com/
            "MKTIIALCYILCLVFA","MKTIIALSHIFCLVLG","MKTIIALSYIFCLAFA","MKTIIALSYIFCLAFG","MKTIIALSYIFCLALG",
            "MKTIIALSYIFCLVFA","MKTIIALSYIFCLVLG","MKTIIALSYIFCQVFA","MKTIIALSYIFCQVLA","MKTIIALSYILCLVFA",
            "MKTIIALSYISCLVFA","MKTIIVLSCFFCLAFS","MKTLIALSYIFCLVLG","MKTTTILILLTHWVHS"
            ],
        "H1": [
            "MKVKLLVLLCTFTATYA","MKVKLLVLLCTFSATYA", # H1seas, http://signalpeptide.com/
            "MKAILVVLLYTFATANA",              # H1pdm
            ],
        "B": [
            "MKAIIVLLMVVTSNA",                # http://repository.kulib.kyoto-u.ac.jp/dspace/bitstream/2433/49327/1/8_1.pdf
            "MVVTSNA",                        # http://signalpeptide.com/
            ],
    }

    SIGNAL_PEPTIDE_LEN = {"H3": 16, "H1": 17, "B": 7}

    GOOD_INFIXES = {
        "H3": [
            {"infix": "ATLCLGHHAV", "shift": 10},
            {"infix": "TNATELV", "shift": 36},   # http://www.rcsb.org/
            {"infix": "VERSKAYSN", "shift": 87},
            ],
        "H1": [
            {"infix": "DTLCIGYHA", "shift": 0},
            {"infix": "DTICIGYHANN", "shift": 0},
            {"infix": "DTICMGYHANN", "shift": 0},
            {"infix": "GYHANNSTDTV", "shift": 5},
            {"infix": "GYHANNSADTV", "shift": 5},
            ],
        "B": [
            {"infix": "RICT", "shift": 1, "check": [{"offset": 40, "seq": "FANL"}]},
            {"infix": "SICT", "shift": 1, "check": [{"offset": 40, "seq": "FANL"}]},
            {"infix": "FANLKG", "shift": 40, "check": [{"offset": 58, "seq": "NCTDLDVAL"}]},
            ],
    }

    EXCLUDE = {
        "H3": [],
        "H1": [
            {"seq": "MNPNQKIITIGSVCMTI", "type": "neuraminidase"},   # http://sbkb.org/
            ],
        "B":  [
            {"seq": "MLPSTIQ", "type": "neuraminidase"},
            {"seq": "MADNMTT", "type": "ns1"},
            {"seq": "MANNNMT", "type": "ns1"},
            ],
    }

    def __init__(self):
        self.good = {}
        self.not_aligned = {}
        self.min_len = 100                # to extract short

    def align(self, sequences, subtype):
        if not subtype:
            from . import fasta as fasta_m
            subtype = fasta_m.detect_subtype(sequences)
        signal_peptide = self.SIGNAL_PEPTIDE[subtype]
        signal_peptide_len = self.SIGNAL_PEPTIDE_LEN[subtype]
        self.good_start = set()
        self.to_align = sequences
        module_logger.info('{} sequences to align'.format(len(self.to_align)))
        if self.to_align and isinstance(next(iter(self.to_align.values())), list):
            # merged
            self.extract_short_merged()
            module_logger.info('{} sequences to align upon removing short ones'.format(len(self.to_align)))
            # self._check_good()
            self.extract_by_signal_peptide_merged(signal_peptide)
            module_logger.info('{} extracted by signal peptide, {} left to align'.format(len(self.good), len(self.to_align)))
            # self._check_good()
            self.extract_by_good_start_merged()
            module_logger.info('{} aligned after using good start, {} left to align'.format(len(self.good), len(self.to_align)))
            # self._check_good()
            self.extract_by_mismatched_signal_peptide_at_beginning_merged(signal_peptide, signal_peptide_len)
            module_logger.info('{} aligned after using mismatched_signal_peptide_at_beginning, {} left to align'.format(len(self.good), len(self.to_align)))
            # self._check_good()
            self.extract_excluded_merged(self.EXCLUDE[subtype])
            module_logger.info('{} aligned after excluding sequences for other proteins, {} left to align'.format(len(self.good), len(self.to_align)))
            # self._check_good()
            self.extract_by_good_infixes_merged(self.GOOD_INFIXES[subtype])
            module_logger.info('{} sequences aligned after using good infixes, {} left to align'.format(len(self.good), len(self.to_align)))
            # self._check_good()
        else:
            # regular
            self.extract_short()
            module_logger.info('{} sequences to align upon removing short ones'.format(len(self.to_align)))
            self.extract_by_signal_peptide(signal_peptide)
            module_logger.info('{} sequences extracted by signal peptide, {} left to align'.format(len(self.good), len(self.to_align)))
            self.extract_by_good_start()
            module_logger.info('{} sequences aligned after using good start, {} left to align'.format(len(self.good), len(self.to_align)))
            self.extract_by_mismatched_signal_peptide_at_beginning(signal_peptide, signal_peptide_len)
            module_logger.info('{} sequences aligned after using mismatched_signal_peptide_at_beginning, {} left to align'.format(len(self.good), len(self.to_align)))
            self.extract_excluded(self.EXCLUDE[subtype])
            module_logger.info('{} aligned after excluding sequences for other proteins, {} left to align'.format(len(self.good), len(self.to_align)))
            self.extract_by_good_infixes(self.GOOD_INFIXES[subtype])
            module_logger.info('{} sequences aligned after using good infixes, {} left to align'.format(len(self.good), len(self.to_align)))
        module_logger.info('\n{:5d} good\n{:5d} unrecognized\n{}\n\nUnrecognized:\n{}'.format(len(self.good), len(self.to_align), "\n".join("{:5d} {}".format(len(v), k) for k, v in self.not_aligned.items()), pprint.pformat(self.to_align)))

    def trim_to_start(self):
        def trim(data):
            if data["shift"] < 0:
                return data["sequence"][-data["shift"]:]
            elif data["shift"] > 0:
                return ("X" * data["shift"]) + data["sequence"]
            else:
                return data["sequence"]
        self.good = {name: trim(data) for name, data in self.good.items()}

    def extract_short(self):
        short_target = self.not_aligned.setdefault("short", {})
        def short(name_sequence):
            if len(name_sequence[1]) < self.min_len:
                short_target[name_sequence[0]] = name_sequence[1]
                name_sequence = (None, None)
            return name_sequence

        self.to_align = dict(e for e in map(short, self.to_align.items()) if e[0] is not None and e[1] is not None)

    def extract_short_merged(self):
        short_target = self.not_aligned.setdefault("short", {})
        def short(name_sequence):
            seq_long = list(filter(lambda s: len(s) >= self.min_len, name_sequence[1]))
            if seq_long:
                name_sequence = (name_sequence[0], seq_long)
            else:                     # all are short
                short_target[name_sequence[0]] = name_sequence[1]
                name_sequence = (None, None)
            return name_sequence

        self.to_align = dict(e for e in map(short, self.to_align.items()) if e[0] is not None and e[1] is not None)

    def extract_by_signal_peptide(self, signal_peptide, start_len=16):
        def good(name, sequence):
            offset = -1
            for sp in signal_peptide:
                offset = sequence.find(sp)
                if offset >= 0:
                    len_sp = len(sp)
                    self.good[name] = {"signal_peptide_len": len_sp, "shift": - offset - len_sp, "sequence": sequence}
                    self.good_start.add(sequence[offset + len_sp : offset + len_sp + start_len])
                    break
            return offset >= 0
        self.to_align = {name: sequence for name, sequence in self.to_align.items() if not good(name, sequence)}

    def extract_by_signal_peptide_merged(self, signal_peptide, start_len=16):

        def good(name_sequence):
            offset = -1
            for sp in signal_peptide:
                for seq in name_sequence[1]:
                    offset = seq.find(sp)
                    if offset >= 0:
                        len_sp = len(sp)
                        self.good[name_sequence[0]] = {"signal_peptide_len": len_sp, "shift": - offset - len_sp, "sequence": seq}
                        self.good_start.add(seq[offset + len_sp : offset + len_sp + start_len])
                        break
                if offset >= 0:
                    break
            return (None, None) if offset >= 0 else name_sequence

        self.to_align = dict(e for e in map(good, self.to_align.items()) if e[0] is not None and e[1] is not None)

    def extract_by_good_start(self):
        def good(name, sequence):
            offset = -1
            for gs in self.good_start:
                o = sequence.find(gs)
                if o >= 0 and (offset < 0 or o < offset):
                    offset = o
            if offset >= 0:
                self.good[name] = {"signal_peptide_len": 0, "shift": - offset, "sequence": sequence}
            return offset >= 0
        self.to_align = {name: sequence for name, sequence in self.to_align.items() if not good(name, sequence)}

    def extract_by_good_start_merged(self):
        def good(name_sequence):
            offset = -1
            for gs in self.good_start:
                for seq in name_sequence[1]:
                    o = seq.find(gs)
                    if o >= 0 and (offset < 0 or o < offset):
                        offset = o
                if offset >= 0:
                    self.good[name_sequence[0]] = {"signal_peptide_len": 0, "shift": - offset, "sequence": seq}
                    break
            return (None, None) if offset >= 0 else name_sequence

        self.to_align = dict(e for e in map(good, self.to_align.items()) if e[0] is not None and e[1] is not None)

    def extract_by_mismatched_signal_peptide_at_beginning(self, signal_peptide, signal_peptide_len):
        def good(name, sequence):
            sp_len, diff_sp = self._differences(signal_peptide, sequence)
            gs_len, diff_gs = self._differences(self.good_start, sequence[sp_len:])
            if diff_gs < 2 and diff_sp < 8:
                self.good[name] = {"signal_peptide_len": sp_len, "shift": - sp_len, "sequence": sequence}
                r = True
            else:
                r = False
            return r
        self.to_align = {name: sequence for name, sequence in self.to_align.items() if not good(name, sequence)}

    def extract_by_mismatched_signal_peptide_at_beginning_merged(self, signal_peptide, signal_peptide_len):
        def good(name_sequence):
            for seq in name_sequence[1]:
                sp_len, diff_sp = self._differences(signal_peptide, seq)
                gs_len, diff_gs = self._differences(self.good_start, seq[sp_len:])
                if diff_gs < 2 and diff_sp < 8:
                    self.good[name_sequence[0]] = {"signal_peptide_len": sp_len, "shift": - sp_len, "sequence": seq}
                    name_sequence = (None, None)
                    break
            return name_sequence

        self.to_align = dict(e for e in map(good, self.to_align.items()) if e[0] is not None and e[1] is not None)

    def extract_excluded(self, exlusions):
        def exclude(name, sequence):
            for p in exlusions:
                if sequence[:len(p["seq"])] == p["seq"]:
                    self.not_aligned.setdefault(p["type"], {})[name] = sequence
                    return True
            return False
        self.to_align = {name: sequence for name, sequence in self.to_align.items() if not exclude(name, sequence)}

    def extract_excluded_merged(self, exlusions):

        def exclude(name_sequence):

            exclude_p = []

            def test_p(seq, p):
                r = seq[:len(p["seq"])] == p["seq"]
                if r:
                    exclude_p.append(p)
                return r

            def exclude_seq(seq):
                return any(test_p(seq, p) for p in exlusions)

            seq = list(filter(lambda s: not exclude_seq(s), name_sequence[1]))
            if seq:
                name_sequence = (name_sequence[0], seq)
            else:
                # if len(exclude_p) > 1:
                #     module_logger.info('EXCLUDE {} {}'.format(name_sequence[0], exclude_p))
                if exclude_p:
                    self.not_aligned.setdefault(exclude_p[0]["type"], {})[name_sequence[0]] = name_sequence[1]
                name_sequence = (None, None)   # all variants excluded
            return name_sequence

        self.to_align = dict(e for e in map(exclude, self.to_align.items()) if e[0] is not None and e[1] is not None)

    def extract_by_good_infixes(self, good_infixes):
        def check(sequence, offset, gi):
            valid = True
            if gi.get("check"):
                for c in gi["check"]:
                    b = offset - gi["shift"] + c["offset"]
                    e = b + len(c["seq"])
                    if sequence[b:e] != c["seq"]:
                        valid = False
                        break
            return valid
        def good(name, sequence):
            r = False
            for gi in good_infixes:
                offset = sequence.find(gi["infix"])
                if offset >= 0 and check(sequence, offset, gi):
                    self.good[name] = {"signal_peptide_len": 0, "shift": - offset + gi["shift"], "sequence": sequence}
                    r = True
                    break
            return r
        self.to_align = {name: sequence for name, sequence in self.to_align.items() if not good(name, sequence)}

    def extract_by_good_infixes_merged(self, good_infixes):

        def check(sequence, offset, gi):
            valid = True
            if gi.get("check"):
                for c in gi["check"]:
                    b = offset - gi["shift"] + c["offset"]
                    e = b + len(c["seq"])
                    if sequence[b:e] != c["seq"]:
                        valid = False
                        break
            return valid

        def good(name_sequence):
            for gi in good_infixes:
                offset = -1
                for seq in name_sequence[1]:
                    offset = seq.find(gi["infix"])
                    if offset >= 0 and check(seq, offset, gi):
                        self.good[name_sequence[0]] = {"signal_peptide_len": 0, "shift": - offset + gi["shift"], "sequence": seq}
                        name_sequence = (None, None)
                        break
                if offset >= 0:
                    break
            return name_sequence

        self.to_align = dict(e for e in map(good, self.to_align.items()) if e[0] is not None and e[1] is not None)

    def _differences(self, prefixes, seq):
        """returns tuple (prefix_len, difference)"""
        def diff(prefix):
            return sum(0 if prefix[p] == seq[p] else 1 for p in range(len(prefix)))
        def sorting_key(e):
            return e[1]
        return sorted(((len(prefix), diff(prefix)) for prefix in prefixes), key=sorting_key)[0]

    def _check_good(self):
        if None in self.good:
            raise RuntimeError("self.good has key None")

# ======================================================================
### Local Variables:
### eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
### End:
