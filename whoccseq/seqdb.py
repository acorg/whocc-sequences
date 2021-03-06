# -*- Python -*-
# license
# license.

"""
Class to access sequence database
"""

import os, sys, json, re, bisect, collections, operator
import logging; module_logger = logging.getLogger(__name__)
from . import open_file, acmacs, amino_acids, utility
from .utility import timeit

# TODO
# export fasta with encoded names
# match HI by cdcid and by name
#      If there are multiple sequences per name/passage, choose the last one (Colin's request of Wed, 4 Sep 2013 20:04:22 +0200)
# convert to aa
# detect clades

# ----------------------------------------------------------------------

class Exclude (Exception): pass

# ----------------------------------------------------------------------

class Seq:

    sReTrailingXs = re.compile(r"X+$")

    def __init__(self, db_entry, seq):
        self.db_entry = db_entry
        self.seq = seq

    def name(self):
        return self.db_entry["N"]

    def hi_name(self):
        return self.seq.get("h")

    def name_hi(self):
        return self.hi_name() or (self.name() + " " + self.passage())

    def seq_id(self):
        return self.name() + "__" + self.passage()

    def virus_type(self):
        return self.db_entry["v"]

    def lineage(self):
        return self.db_entry.get("l")

    def passage(self):
        p = self.seq.get("p")
        if p:
            p = p[0]
        else:
            p = ""
        return p

    def shift(self):
        return self.seq.get("s")

    def nuc_shift(self):
        return self.seq.get("t")

    def aa(self):
        return self.seq.get("a")

    def nuc(self):
        return self.seq["n"]

    def date(self):
        d = self.db_entry.get("d") or None
        if d:
            d = d[-1]
        return d

    def aligned(self, amino_acid :bool, truncate :int = None):
        if amino_acid:
            return self.aa_aligned(truncate=truncate)
        else:
            return self.nuc_aligned(truncate=truncate)

    def nuc_aligned(self, truncate=None):
        shift = self.nuc_shift()
        if shift is None:
            raise RuntimeError("Not aligned")
        s = self._shift(self.nuc(), shift, fill="-")
        return self._truncate(s, truncate)

    def aa_aligned(self, truncate=None):
        shift = self.shift()
        if shift is None:
            raise RuntimeError("Not aligned")
        s = self._shift(self.aa(), shift, fill="X")
        if True:
            # find the longest part not having *, replace other parts with X
            parts = sorted(enumerate(s.split("*")), key=lambda e: len(e[1]), reverse=True)
            if len(parts) == 2 and parts[1][1] == "":
                parts = parts[:1]
            s = "X".join(e[1] for e in sorted((p if i == 0 else (p[0], "X" * len(p[1])) for i, p in enumerate(parts)), key=operator.itemgetter(0)))
            if truncate is None:
                # strip trailing Xs
                s = self.sReTrailingXs.sub("", s)
        return self._truncate(s, truncate)

    @classmethod
    def _shift(cls, seq, shift, fill):
        if shift < 0:
            seq = seq[-shift:]
        elif shift > 0:
            seq = (fill * shift) + seq
        return seq

    @classmethod
    def _truncate(cls, seq, truncate):
        if truncate is not None:
            if len(seq) < truncate:
                seq = seq + "-" * (truncate - len(seq))
            elif len(seq) > truncate:
                seq = seq[:truncate]
        return seq

    def hamming_distance(self, other :"Seq", amino_acid :bool, truncate :int = None):
        s1 = self.aligned(amino_acid=amino_acid, truncate=truncate)
        s2 = other.aligned(amino_acid=amino_acid, truncate=truncate)
        l = min(len(s1), len(s2))
        return sum(1 for pos in range(l) if s1[pos] != s2[pos])

    def hamming_distance_many(self, other :list, amino_acid :bool):
        """Returns list of distances (int), order is the same as other"""
        most_common_length = most_common_length(other, amino_acid=amino_acid)
        return [self.hamming_distance(e, amino_acid=amino_acid, truncate=most_common_length) for e in other]

    def gene(self):
        return self.seq.get("g", "HA")

    def lab_matches(self, lab):
        return bool(lab) and lab.upper() in self.seq["l"]

    def labs(self):
        return list(self.seq["l"])

    def lab_ids(self, lab):
        return self.seq["l"][lab]

    def aa_len(self):
        return len(self.aa()) - max(- self.shift(), 0)

    def aa_at_pos0(self, pos):
        """Returns aa at pos, pos is counted starting with 0"""
        shift = self.shift()
        if pos >= shift:
            return self.aa()[pos - shift]
        else:
            return "X"

    def aa_at_pos(self, pos):
        """Returns aa at pos, pos is counted starting with 1"""
        return self.aa_at_pos0(pos - 1)

    def align(self):
        SeqDB.align(sequence=self.aa(), entry_passage=self.seq, data=None, db_entry=self.db_entry, verbose=True)

    def clades(self):
        return self.seq.get("c")

    def set_clades(self, clades):
        self.seq["c"] = clades

# ----------------------------------------------------------------------

def most_common_length(sequences, amino_acid=True):
    """Returns most common length for the passed list of Seq"""
    len_stat = collections.Counter(len(e.aligned(amino_acid=amino_acid)) for e in sequences)
    return len_stat.most_common(1)[0][0]

# ----------------------------------------------------------------------

# self.data is list of dicts sorted by "N":
# {
#     "N": <name>,
#     "s": [
#         {
#             ["p": <list of passages (str)>,]
#             "n": <sequence-nucleotides>,
#             "a": <sequence-amino-acids>,
#             "s": <shift (int) for aa sequence>,
#             "t": <shift (int) for nuc sequence>,
#             "l": {<lab> :[<lab_id>]},
#             "g": <gene: HA|NA>,
#             "h": <hi-name>,
#             "c": <list of clades (str)>,
#         },
#     ],
#     "v": <virus_type>,
#     "l": "VICTORIA",  # VICTORIA, YAMAGATA, 2009PDM, SEASONAL
#     "d": [<date>]
# }

class SeqDB:

    def __init__(self, path_to_db, try_to_load=True):
        self.path_to_db = path_to_db
        if try_to_load:
            self.load()
        else:
            self.data = []

        # --------------------------------------------------

    def add_sequences(self, data):
        """Adds sequences to the database, data is a list of dict {"sequence":, "name":, "date":, "passage":, "lab":, "virus_type":, "gene":, "lab_id":}
        name, passage, date will be normalized before adding.
        """
        self.normalize(data)
        num_added = 0
        for entry in data:
            if self._add_sequence(entry):
                num_added += 1
        module_logger.info('{} new sequences added'.format(num_added))

        # --------------------------------------------------

    def load(self):
        if os.path.isfile(self.path_to_db):
            with timeit('Reading {}'.format(os.path.realpath(self.path_to_db))):
                data = open_file.read_json(self.path_to_db)
                if data.get("  version") != "sequence-database-v2":
                    raise RuntimeError("Unrecognized sequence database version: {}".format(data.get("  version")))
        else:
            data = {}
        self.data = data.get("data", [])

    def save(self):
        with timeit("Writing {}".format(self.path_to_db)):
            open_file.write_json(self.path_to_db, {"  version": "sequence-database-v2", "data": self.data}, indent=1, sort_keys=True)

        # --------------------------------------------------

    def iterate_sequences(self):
        """Yields Seq"""
        for db_entry in self.data:
            for seq in db_entry["s"]:
                yield Seq(db_entry, seq)

    def iterate_sequences_with(self, labs :list = None, virus_types :list = None, lineage :str = None, gene :str = "HA", aligned=False, start_date :str = None, end_date :str = None, with_hi_name=False, name_matcher :str = None):
        """Yields Seq"""
        if not labs:
            labs = None
        elif isinstance(labs, str):
            labs = {self.normalize_lab(labs)}
        elif isinstance(labs, (list, tuple, set)):
            labs = set(self.normalize_lab(l) for l in labs)
        else:
            raise ValueError("Unrecognized labs: {}".format(labs))

        if not virus_types:
            virus_types = None
        elif isinstance(virus_types, str):
            virus_types = [self.normalize_virus_type(virus_types)]
        elif isinstance(virus_types, (list, tuple, set)):
            virus_types = [self.normalize_virus_type(vt) for vt in virus_types]
        else:
            raise ValueError("Unrecognized virus_types: {}".format(virus_types))

        if not lineage:
            lineage = None
        elif isinstance(lineage, str):
            lineage = self.normalize_lineage(lineage)
        else:
            raise ValueError("Unrecognized lineage: {!r}".format(lineage))

        if not name_matcher:
            name_matcher = None
        elif isinstance(name_matcher, str):
            name_matcher = re.compile(name_matcher, re.I)
        else:
            raise ValueError("Unrecognized name_matcher: {!r}".format(name_matcher))

        start_date = self._fix_date(start_date)
        end_date = self._fix_date(end_date)

        for db_entry in self.data:
            for seq in db_entry["s"]:
                if ((not labs or (labs & set(seq["l"])))
                     and (not virus_types or db_entry["v"] in virus_types)
                     and (not lineage or db_entry["l"] == lineage)
                     and gene == seq.get("g", "HA")
                     and (not aligned or seq.get("s") is not None)
                     and (not start_date or (db_entry.get("d") and db_entry["d"][-1] >= start_date))
                     and (not end_date or (db_entry.get("d") and db_entry["d"][-1] < end_date))
                     and (not with_hi_name or seq.get("h"))
                     and (not name_matcher or name_matcher.search(seq.get("h") or "{} {}".format(db_entry["N"], seq.get("p") and seq["p"][0])))
                     ):
                    yield Seq(db_entry, seq)

    def iterate_sequences_aligned_with_virus_type(self, virus_type):
        """Yields Seq"""
        virus_type, _ = utility.fix_virus_type_lineage(virus_type)
        for db_entry in self.data:
            if db_entry["v"] == virus_type:
                for seq in db_entry["s"]:
                    if seq.get("s") is not None:
                        yield Seq(db_entry, seq)

    def iterate_sequences_with_name(self, names):
        """Yields Seq"""
        for raw_name in ([names] if isinstance(names, str) else names):
            name, passage = self._split_name(raw_name)
            db_entry = self.find_by_name(name)
            if db_entry:
                for seq in db_entry["s"]:
                    if not passage or passage in seq.get("p", []):
                        yield Seq(db_entry, seq)
            else:
                module_logger.error('{!r} not found in the database'.format(name))

    def iterate_sequences_not_aligned(self):
        """Yields Seq"""
        for db_entry in self.data:
            for seq in db_entry["s"]:
                if "s" not in seq:
                    yield Seq(db_entry, seq)

        # --------------------------------------------------

    def report(self):
        print("SeqDB names:", len(self.data))

        strange_names = [e["N"] for e in self.data if e["N"].count("/") < 3]
        if strange_names:
            print("Strange names {}:\n\t{}".format(len(strange_names), "\n\t".join(strange_names)))

        multi_seq_per_passage = [ee for ee in (self._multiple_sequences_per_passage(e) for e in self.data) if ee]
        if multi_seq_per_passage:
            # print("Multiple sequences per passage {}:\n\t{}".format(len(multi_seq_per_passage), "\n\t".join("{!r} {!r}".format(n, p) for n, p in multi_seq_per_passage)))
            print("Multiple sequences per passage {}".format(len(multi_seq_per_passage)))

        not_translated_to_aa = sum(1 for ee in (all(bool(e2.get("a")) for e2 in e["s"]) for e in self.data) if not ee)
        if not_translated_to_aa:
            print("Not translated to amino-acids {}".format(not_translated_to_aa))

        aligned = {}                      # by virus_type
        not_aligned = {}
        for e1 in self.data:
            for e2 in e1["s"]:
                if e2.get("g", "HA") == "HA":
                    if e2.get("s") is not None:
                        aligned.setdefault(e1["v"], 0)
                        aligned[e1["v"]] += 1
                    else:
                        not_aligned.setdefault(e1["v"], 0)
                        not_aligned[e1["v"]] += 1
        aligned[" All"] = sum(aligned.values())
        not_aligned[" All"] = sum(not_aligned.values())
        total = {k: (aligned[k] + not_aligned.get(k, 0)) for k in aligned}
        ks = sorted(aligned)
        print("HA aligned:\n  {}\nHA not aligned\n  {}".format("\n  ".join("{:<7s} {:d} {:.1f}%".format(k.strip(), aligned[k], (aligned[k] / total[k]) * 100.0) for k in ks), "\n  ".join("{:<7s} {}".format(k.strip(), not_aligned.get(k, 0)) for k in ks)))

        ha_sequences = 0
        hi_names = 0
        lineages = collections.defaultdict(int)
        clades = collections.defaultdict(int)
        for e1 in self.data:
            for e2 in e1["s"]:
                if e2.get("g", "HA") == "HA":
                    ha_sequences += 1
                    if e2.get("h"):
                        hi_names += 1
                for clade in e2.get("c", []):
                    clades["{}-{}".format(e1["v"], clade)] += 1
            if e1.get("l"):
                lineages[e1["l"]] += 1

        print("HI names: {} ({:.1f}%)".format(hi_names, hi_names / ha_sequences * 100.0))
        print("Lineages:\n  {}".format("\n  ".join("{:<8s} {:>4d}".format(l, lineages[l]) for l in sorted(lineages))))
        print("Clades:\n  {}".format("\n  ".join("{:<10s} {:>4d}".format(c, clades[c]) for c in sorted(clades))))

        # --------------------------------------------------

    def reset_hi_data(self):
        for e in self._iterate_sequences():
            if "hi_name" in e:
                del e["hi_name"]

    def reset_clade_data(self):
        for e in self._iterate_sequences():
            if "clades" in e:
                del e["clades"]

        # --------------------------------------------------

    def _iterate_sequences(self):
        for e1 in self.data:
            for e2 in e1["s"]:
                yield e2

        # --------------------------------------------------

    def _add_sequence(self, data):
        name = data.get("name")
        if name:
            if name[1] in ["/", "("] and name[0] != data["virus_type"][0]:
                module_logger.warning('Virus type ({}) and name ({}) mismatch'.format(data["virus_type"], name))
            entry = self.find_by_name(name)
            if entry is None:
                if name[:8] not in ["A(H1N1)/", "A(H3N2)/"] and name[:2] != "B/":
                    module_logger.warning('Suspicious name {!r}'.format(name))
                entry = {"N": name, "s": [], "v": data["virus_type"]}
                self._insert_by_name(entry)
            try:
                new = self._update_db_entry(entry, data)
            except Exclude as err:
                module_logger.info('Sequence excluded ({}): {}'.format(err, data["name"]))
                new = False
                if not entry["s"]:        # no sequences in entry, remove the entry
                    index = self.find_by_name(entry["N"], return_index=True)
                    if index is not None:
                        del self.data[index]
        else:
            module_logger.warning('Entry without name: {}'.format(data["lab_id"]))
            new = False
        return new

    def _update_db_entry(self, entry, data):
        if entry["v"] != data["virus_type"]:
            raise RuntimeError("Cannot add {!r} to {!r} db entry\n{}\n{}".format(data["virus_type"], entry["v"], data, entry))
        if data.get("date") and data["date"] not in entry.get("d", []):
            entry.setdefault("d", []).append(data["date"])
            entry["d"].sort()
        sameseq = self._look_for_the_same_sequence(entry, data)
        new = False
        if sameseq is None:
            new_entry = {"l": {}, "g": "HA"}
            self._update_entry_passage(entry_passage=new_entry, data=data, sequence_match="new", db_entry=entry)
            entry["s"].append(new_entry)
            module_logger.debug('new sequence entry added {} {}'.format(data["name"], data.get("passage", "")))
            new = True
        elif sameseq["type"] == "update":
            self._update_entry_passage(entry_passage=entry["s"][sameseq["index"]], data=data, sequence_match=sameseq["sequence_match"], db_entry=entry)
        elif sameseq["type"] == "different-genes":
            module_logger.warning(sameseq["message"])
        else:
            module_logger.error("[INTERNAL] Unrecoginzed sameseq: {}".format(sameseq))
        return new

    def _look_for_the_same_sequence(self, entry, data):
        """If no matching sequence found returns None. Else updates
        passage entry, if passage is the same. Returns string describing
        match type.
        """
        r = None
        data_passage = data.get("passage", "")
        for e_no, e in enumerate(entry["s"]):
            sequence_match = self._sequences_match(master=e["n"], s=data["sequence"])
            if sequence_match:
                if data_passage and data_passage not in e.get("p", []):
                    # if e["p"]:
                    #     module_logger.warning('[SAMESEQ] different passages {!r} {!r} {!r}'.format(data["name"], e["p"], data_passage))
                    r = {"type": "update", "sequence_match": sequence_match, "index": e_no}
                elif e.get("g") and data.get("gene") and e["g"] != data["gene"]:
                    r = {"type": "different-genes", "message": '[SAMESEQ] different genes {!r} {!r} {!r}'.format(data["name"], e["g"], data["gene"])}
                else:
                    r = {"type": "update", "sequence_match": sequence_match, "index": e_no}
                break
        return r

    def _sequences_match(self, master, s):
        if master == s:
            r = "equal"
        elif s in master:
            r = "sub"
        elif master in s:
            r = "super"
        else:
            r = None
        return r

    def _update_entry_passage(self, entry_passage, data, sequence_match, db_entry):
        if data.get("passage") and data["passage"] not in entry_passage.get("p", []):
            entry_passage.setdefault("p", []).append(data["passage"])
        lab_e = entry_passage["l"].setdefault(data["lab"], [])
        if data.get("lab_id") and data["lab_id"] not in lab_e:
            lab_e.append(data["lab_id"])
        if data.get("gene") and not entry_passage.get("g"):
            entry_passage["g"] = data["gene"]
        if sequence_match in ["super", "new"]:   # update sequences with the longer one
            entry_passage["n"] = data["sequence"]
            try:
                align_data = amino_acids.translate_to_aa_and_align(data["sequence"], name=data["name"])
                if align_data:
                    entry_passage["a"] = align_data["aa"]
                    if align_data.get("shift") is not None:
                        entry_passage["s"] = align_data["shift"]
                        entry_passage["t"] = - align_data["offset"] + entry_passage["s"] * 3
                        self._check_aligned(entry_passage, db_entry)
                else:
                    module_logger.warning('Not translated/aligned {}'.format(data["name"]))
            except amino_acids.SequenceIsTooShort as err:
                raise Exclude(str(err))

    sReH1AlignCheck1 = re.compile(r"[HX][AX][NX][NX][SX]")
    sReH1AlignCheck2 = re.compile(r"[SV]YIVE")
    sReH3AlignCheck1 = re.compile(r"T[IX](TN|AN|TS|TD|TH|TY|TK|MN|XN)")
    sReBAlignCheck1 = re.compile(r"T[AITX][AIST](PTK|PIK|LTK|PTR|PAK|PXK|PTQ|PKK|PVK|PTX)")

    def _check_aligned(self, entry_passage, db_entry):
        shift = entry_passage["s"]
        if db_entry["v"] == "A(H1N1)":
            if shift <= 0 and not self.sReH1AlignCheck1.match(entry_passage["a"][-shift + 7: -shift + 12]) and not self.sReH1AlignCheck2.match(entry_passage["a"][-shift + 76: -shift + 81]):
                module_logger.warning('Alignment problem (H1 HANNS or SYIVE) for {}: {} {}'.format(db_entry["N"], entry_passage["a"][-shift + 7: -shift + 12], entry_passage["a"][-shift + 76: -shift + 81]))
        elif db_entry["v"] == "A(H3N2)":
            if shift <= 0 and not self.sReH3AlignCheck1.match(entry_passage["a"][-shift + 27: -shift + 31]):
                module_logger.warning('Alignment problem (H3 TITN) for {}: {}'.format(db_entry["N"], entry_passage["a"][-shift + 27: -shift + 31]))
        elif db_entry["v"] == "B":
            if shift <= 0 and not self.sReBAlignCheck1.match(entry_passage["a"][-shift + 32: -shift + 38]):
                module_logger.warning('Alignment problem (B TTTPTK) for {}: {}'.format(db_entry["N"], entry_passage["a"][-shift + 32: -shift + 38]))

    # def _update_entry_passage(self, entry_passage, data, sequence_match, db_entry):
    #     if data.get("passage") and data["passage"] not in entry_passage["p"]:
    #         entry_passage["p"].append(data["passage"])
    #     lab_e = entry_passage["l"].setdefault(data["lab"], [])
    #     if data.get("lab_id") and data["lab_id"] not in lab_e:
    #         lab_e.append(data["lab_id"])
    #     if data.get("gene") and not entry_passage.get("g"):
    #         entry_passage["g"] = data["gene"]
    #     if sequence_match in ["super", "new"]:   # update sequences with the longer one
    #         aa = amino_acids.translate_sequence_to_amino_acid(data["sequence"], name=data["name"])
    #         entry_passage["n"] = data["sequence"]
    #         entry_passage["a"] = aa.translated
    #         self.align(aa.translated, entry_passage, data, db_entry)   # may raise
    #         if entry_passage.get("s") is not None:
    #             entry_passage["t"] = - aa.offset + entry_passage["s"] * 3
    #         if "BEIJING CHAOYANG/158/2010" in data["name"]:
    #             module_logger.warning('{} ALIGN nuc:{} aa:{} shift:{}'.format(data["name"], data["sequence"][:50], aa.translated[:50], entry_passage["s"]))
    #             #raise NotImplementedError()

    # @classmethod
    # def align(cls, sequence, entry_passage, data, db_entry, verbose=False):
    #     try:
    #         aligment_data = amino_acids.align(sequence, verbose=verbose)
    #     except amino_acids.SequenceIsTooShort as err:
    #         raise Exclude(str(err))
    #     if aligment_data:
    #         module_logger.debug('aligment_data {}'.format(aligment_data))
    #         if aligment_data["virus_type"] != db_entry["v"]:
    #             module_logger.warning('Virus type detection mismatch {} vs. {}'.format(aligment_data, data))
    #         if aligment_data.get("lineage"):
    #             if db_entry.get("l"):
    #                 if db_entry["l"] != aligment_data["lineage"]:
    #                     module_logger.warning('Lineage detection mismatch for {}: {} vs. {}'.format(db_entry["N"], aligment_data, db_entry["l"]))
    #             else:
    #                 db_entry["l"] = aligment_data["lineage"]
    #         if entry_passage.get("g"):
    #             if entry_passage["g"] != aligment_data["gene"]:
    #                 module_logger.warning('Gene detection mismatch for {}: {} vs. {}'.format(db_entry["N"], aligment_data, entry_passage["g"]))
    #         else:
    #             entry_passage["g"] = aligment_data["gene"]
    #         entry_passage["s"] = aligment_data["shift"]
    #     elif verbose: # if db_entry["v"] in ["A(H3N2)", "A(H1N1)"]:
    #         module_logger.warning('Not aligned {:<45s} len:{:3d} {}'.format(db_entry["N"], len(sequence), sequence[:40]))


        # --------------------------------------------------

    # for report
    def _multiple_sequences_per_passage(self, entry):
        passages = {}
        for e in entry["s"]:
            for p in e.get("p", []):
                passages.setdefault(p, []).append(e.get("g"))
        return sorted(p for p, g in passages.items() if len(set(g)) != len(g))

        # --------------------------------------------------

    def normalize(self, data):
        normalized = {k: self._normalize_x(k, data) for k in ["name", "passage", "date"]}
        for entry in data:
            if entry.get("name"):
                norm = normalized["name"][entry["name"]]
                if norm[:2] == "A/" and entry["virus_type"][0] == "A":
                    norm = "{}{}".format(entry["virus_type"], norm[1:])
                elif norm[1] != "/" and norm.count("/") == 2:
                    norm = "{}/{}".format(entry["virus_type"], norm)
                if norm != entry["name"]:
                    entry["raw_name"] = entry["name"]
                    entry["name"] = norm
            for key in ["passage", "date"]:
                if entry.get(key):
                    norm = normalized[key][entry[key]]
                    if norm != entry[key]:
                        entry["raw_" + key] = entry[key]
                        entry[key] = norm

    def _normalize_x(self, key, data):
        source = set(e[key] for e in data if e.get(key))
        module_logger.info('{} {}s to normalize'.format(len(source), key))
        return getattr(acmacs, "normalize_{}s".format(key))(source)

    # ----------------------------------------------------------------------

    sReNamePrefix = re.compile(r"^(B|A|A\(H\d+\)|A\(H\d+N\d+\))", re.I)
    sReYear = re.compile(r"^/(19[0-9][0-9]|20[0-2][0-9])$")

    def _split_name(self, raw_name):
        name = None
        passage = None
        if self.sReNamePrefix.match(raw_name):
            fields = raw_name.split()
            for name_parts in range(len(fields)):
                if self.sReYear.match(fields[name_parts][-5:]):
                    name = " ".join(fields[:name_parts + 1])
                    passage = " ".join(fields[name_parts + 1:])
                    break
        return name, passage

    # ----------------------------------------------------------------------

    @classmethod
    def normalize_virus_type(cls, vt):
        if vt is None:
            pass
        elif vt.upper() in ["B", "BV", "BY", "BVIC", "BYAM", "B/VIC", "B/YAM"]:
            vt = "B"
        elif vt.upper() in ["H1PDM", "H1SEAS", "H1", "A(H1N1)"]:
            vt = "A(H1N1)"
        elif vt.upper() in ["H3", "A(H3N2)"]:
            vt = "A(H3N2)"
        else:
            raise ValueError("Unrecognized virus type {}".format(vt))
        return vt

    sLabs = {"CDC": "CDC", "CNIC": "CNIC", "CRICK": "NIMR", "MELB": "MELB", "NIID": "NIID", "NIMR": "NIMR", "VDRL": "MELB"}

    @classmethod
    def normalize_lab(cls, lab):
        try:
            return cls.sLabs[lab.strip().upper()]
        except:
            raise ValueError("Unrecognized lab: {!r}".format(lab))

    sLineages = {"VICTORIA": "VICTORIA", "YAMAGATA": "YAMAGATA", "VIC": "VICTORIA", "YAM": "YAMAGATA"}

    @classmethod
    def normalize_lineage(cls, lineage):
        try:
            return cls.sLineages[lineage.strip().upper()]
        except:
            raise ValueError("Unrecognized lineage: {!r}".format(lineage))

    # ----------------------------------------------------------------------

    sReDate = re.compile(r"^\s*(?P<year>\d\d\d\d)[-\s]*(?P<month>\d\d)[-\s]*(?P<day>\d\d)\s*$")

    @classmethod
    def _fix_date(cls, date):
        if date:
            m = cls.sReDate.match(date)
            if not m:
                raise ValueError("Cannot parse date: {!r}".format(date))
            date = "-".join((m.group("year"), m.group("month"), m.group("day")))
        return date

    # ----------------------------------------------------------------------
    # Adoption if python 3.5 bisect module functions

    def find_by_name(self, name, return_index=False):
        lo = 0
        hi = len(self.data)
        while lo < hi:
            mid = (lo + hi) // 2
            mid_name = self.data[mid]["N"]
            if mid_name == name:
                return mid if return_index else self.data[mid]
            if mid_name < name:
                lo = mid + 1
            else:
                hi = mid
        return None

    def _insert_by_name(self, entry):
        lo = 0
        hi = len(self.data)
        while lo < hi:
            mid = (lo + hi) // 2
            if self.data[mid]["N"] < entry["N"]:
                lo = mid + 1
            else:
                hi = mid
        self.data.insert(lo, entry)

# ======================================================================
### Local Variables:
### eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
### End:
