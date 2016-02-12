# -*- Python -*-
# license
# license.

"""
Class to access sequence database
"""

import os, sys, json
import logging; module_logger = logging.getLogger(__name__)
from . import open_file, acmacs, amino_acids, utility
from .utility import timeit

# TODO
# export fasta with encoded names
# match HI by cdcid and by name
#      If there are multiple sequences per name/passage, choose the last one (Colin's request of Wed, 4 Sep 2013 20:04:22 +0200)
# convert to aa
# align by aa, mark HA, NA
# detect clades

# ----------------------------------------------------------------------

class Exclude (Exception): pass

# ----------------------------------------------------------------------

# self.names is dict:
#    {<name>:
#        {
#            "data": [
#                {
#                    "passages": [<passage>],
#                    "nuc": <sequence-nucleotides>,
#                    "aa": <sequence-amino-acids>,
#                    "labs": {<lab> :[<lab_id>]},
#                    "gene": <HA|NA>,
#                    "hi_name": "",
#                    "clades": [""],
#                },
#            ],
#            "virus_type": <virus_type>,
#            "lineage": "VICTORIA",  VICTORIA, YAMAGATA, 2009PDM, SEASONAL
#            "dates": [<date>]
#        },
#    }

class SeqDB:

    def __init__(self, path_to_db):
        self.path_to_db = path_to_db
        self.load()

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

    def get(self, entry, amino_acid, aligned):
        if amino_acid:
            if aligned:
                s = self._aligned(entry)
            else:
                s = entry["aa"]
        else:
            s = entry["nuc"]
        return s

    def get_aa(self, source=None, aligned=True):
        """Returns list of aa sequences from source"""
        if source is None:
            source = self.all()
        return [self._aligned(e2) if aligned else e2["aa"] for n, e1 in source.items() for e2 in e1["data"]]

        # --------------------------------------------------

    def load(self):
        if os.path.isfile(self.path_to_db):
            module_logger.info('Reading {}'.format(os.path.realpath(self.path_to_db)))
            data = open_file.read_json(self.path_to_db)
            if data.get("  version") != "sequence-database-v1":
                raise RuntimeError("Unrecognized sequence database version: {}".data.get("  version"))
        else:
            data = {}
        self.names = data.get("1-names", {})
        self.cdcids = data.get("2-cdcids", {})

    def save(self):
        data = {"  version": "sequence-database-v1", "1-names": self.names, "2-cdcids": self.cdcids}
        module_logger.info('Writing {}'.format(self.path_to_db))
        with timeit("Written in"):
            open_file.write_json(self.path_to_db, data, indent=1, sort_keys=True)

        # --------------------------------------------------

    def iterate_sequences(self):
        """Yields {"name":, "virus_type":, "lineage":, "dates":, "seq": <"data" entry>}"""
        for name, db_entry in self.names.items():
            e = {k: v for k,v in db_entry.items() if k != "data"}
            e["name"] = name
            for seq in db_entry["data"]:
                e["seq"] = seq
                yield e

    def select(self, lab, virus_type, lineage, gene):
        data = self.all()
        if lab:
            data = self.select_by("lab", lab.upper(), data)
        virus_type, lineage = utility.fix_virus_type_lineage(virus_type, lineage)
        if virus_type:
            data = self.select_by("virus_type", virus_type, data)
        if lineage:
            data = self.select_by("lineage", lineage, data)
        data = self.select_by("gene", gene.upper(), data)
        return data

    def all(self):
        return self.names

    def select_by(self, field, value, source=None):
        if source is None:
            source = self.all()
        if field in ["virus_type", "lineage"]:
            data = {n: e for n,e in source.items() if e.get(field) == value}
        elif field == "lab":
            def fix_lab(entry):
                e = {f: v for f,v in entry.items() if v != "labs"}
                e["labs"] = {value: entry["labs"][value]}
                return e
            def filter_lab(entry):
                data = [fix_lab(e) for e in entry["data"] if value in e["labs"]]
                if data:
                    r = {f: v for f,v in entry.items() if f != "data"}
                    r["data"] = data
                else:
                    r = None
                return r
            data = {n: e for n, e in ((nn, filter_lab(ee)) for nn, ee in source.items()) if e}
        elif field == "gene":
            def filter_gene(entry):
                data = [e for e in entry["data"] if e.get("gene", "HA") == value]
                if data:
                    r = {f: v for f,v in entry.items() if f != "data"}
                    r["data"] = data
                else:
                    r = None
                return r
            data = {n: e for n, e in ((nn, filter_gene(ee)) for nn, ee in source.items()) if e}
        else:
            raise ValueError("Unsupported field {!r} to select by".format(field))
        return data

    def select_aligned(self, source=None):

        def filter_aligned(entry):
            data = [e for e in entry["data"] if e.get("shift") is not None]
            if data:
                r = {f: v for f,v in entry.items() if f != "data"}
                r["data"] = data
            else:
                r = None
            return r

        if source is None:
            source = self.all()
        data = {n: e for n, e in ((nn, filter_aligned(ee)) for nn, ee in source.items()) if e}
        return data

    def names_sorted_by(self, field, source=None):
        if source is None:
            source = self.all()
        if field == "name":
            r = sorted(source)
        elif field == "date":
            r = sorted(source, key=lambda n: source[n]["dates"][-1] if source[n].get("dates") else "0000") # entries lacking date come first
        else:
            raise ValueError("Unsupported field {!r} to sort by".format(field))
        return r

        # --------------------------------------------------

    def report(self):
        print("SeqDB names:", len(self.names))
        print("SeqDB cdcids:", len(self.cdcids))

        strange_names = [n for n in self.names if n.count("/") < 3]
        if strange_names:
            print("Strange names {}:\n\t{}".format(len(strange_names), "\n\t".join(strange_names)))

        multi_seq_per_passage = [ee for ee in ((n, self._multiple_sequences_per_passage(e)) for n, e in self.names.items()) if ee[1]]
        if multi_seq_per_passage:
            # print("Multiple sequences per passage {}:\n\t{}".format(len(multi_seq_per_passage), "\n\t".join("{!r} {!r}".format(n, p) for n, p in multi_seq_per_passage)))
            print("Multiple sequences per passage {}".format(len(multi_seq_per_passage)))

        not_translated_to_aa = [ee[0] for ee in ((n, all(bool(e2.get("aa")) for e2 in e["data"])) for n, e in self.names.items()) if not ee[1]]
        if not_translated_to_aa:
            print("Not translated to amino-acids {}".format(len(not_translated_to_aa)))

        aligned = {}                      # by virus_type
        not_aligned = {}
        for n, e1 in self.names.items():
            for e2 in e1["data"]:
                if e2.get("shift") is not None:
                    aligned.setdefault(e1["virus_type"], 0)
                    aligned[e1["virus_type"]] += 1
                else:
                    not_aligned.setdefault(e1["virus_type"], 0)
                    not_aligned[e1["virus_type"]] += 1
        aligned[" All"] = sum(aligned.values())
        not_aligned[" All"] = sum(not_aligned.values())
        total = {k: (aligned[k] + not_aligned.get(k, 0)) for k in aligned}
        ks = sorted(aligned)
        print("Aligned:\n  {}\nNot aligned\n  {}".format("\n  ".join("{:<7s} {:d} {:.1f}%".format(k.strip(), aligned[k], (aligned[k] / total[k]) * 100.0) for k in ks), "\n  ".join("{:<7s} {}".format(k.strip(), not_aligned.get(k, 0)) for k in ks)))

        ha_sequences = 0
        hi_names = 0
        for e1 in self.names.values():
            for e2 in e1["data"]:
                if e2.get("gene", "HA") == "HA":
                    ha_sequences += 1
                    if e2.get("hi_name"):
                        hi_names += 1
        print("HI names: {} ({:.1f}%)".format(hi_names, hi_names / ha_sequences * 100.0))

        lineages = {}
        for e1 in self.names.values():
            if e1.get("lineage"):
                lineages.setdefault(e1["lineage"], 0)
                lineages[e1["lineage"]] += 1
        print("Lineages:\n  {}".format("\n  ".join("{:<8s} {:>4d}".format(l, lineages[l]) for l in sorted(lineages))))

        clades = {}
        for e1 in self.names.values():
            for e2 in e1["data"]:
                for clade in e2.get("clades", []):
                    key = "{}-{}".format(e1["virus_type"], clade)
                    if key in clades:
                        clades[key] += 1
                    else:
                        clades[key] = 1
        print("Clades:\n  {}".format("\n  ".join("{:<10s} {:>4d}".format(c, clades[c]) for c in sorted(clades))))

        # --------------------------------------------------

    def reset_hi_data(self):
        for e1 in self.names.values():
            for e2 in e1["data"]:
                if "hi_name" in e2:
                    del e2["hi_name"]

    def reset_clade_data(self):
        for e1 in self.names.values():
            for e2 in e1["data"]:
                if "clades" in e2:
                    del e2["clades"]

        # --------------------------------------------------

    def _add_sequence(self, data):
        name = data.get("name")
        if name:
            if name[1] in ["/", "("] and name[0] != data["virus_type"][0]:
                module_logger.warning('Virus type ({}) and name ({}) mismatch'.format(data["virus_type"], name))
            entry = self.names.get(name)
            if entry is None:
                if name[:8] not in ["A(H1N1)/", "A(H3N2)/"] and name[:2] != "B/":
                    module_logger.warning('Suspicious name {!r}'.format(name))
                entry = {"data": [], "virus_type": data["virus_type"], "dates": []}
                self.names[name] = entry
            new = self._update_db_entry(entry, data)
            self._update_cdcid(data)
        else:
            module_logger.warning('Entry without name: {}'.format(data["lab_id"]))
            new = False
        return new

    def _update_db_entry(self, entry, data):
        if entry["virus_type"] != data["virus_type"]:
            raise RuntimeError("Cannot add {!r} to {!r} db entry\n{}\n{}".format(data["virus_type"], entry["virus_type"], data, entry))
        if data.get("date") and data["date"] not in entry["dates"]:
            entry["dates"].append(data["date"])
            entry["dates"].sort()
        sameseq = self._look_for_the_same_sequence(entry, data)
        new = False
        try:
            if sameseq is None:
                new_entry = {"labs": {}, "passages": []}
                self._update_entry_passage(entry_passage=new_entry, data=data, sequence_match="new", db_entry=entry)
                entry["data"].append(new_entry)
                module_logger.debug('new sequence entry added {} {}'.format(data["name"], data.get("passage", "")))
                new = True
            elif sameseq["type"] == "update":
                self._update_entry_passage(entry_passage=entry["data"][sameseq["index"]], data=data, sequence_match=sameseq["sequence_match"], db_entry=entry)
            elif sameseq["type"] == "different-genes":
                module_logger.warning(sameseq["message"])
            else:
                module_logger.error("[INTERNAL] Unrecoginzed sameseq: {}".format(sameseq))
        except Exclude as err:
            module_logger.info('Sequence excluded ({}): {}'.format(err, data["name"]))
        return new

    def _look_for_the_same_sequence(self, entry, data):
        """If no matching sequence found returns None. Else updates
        passage entry, if passage is the same. Returns string describing
        match type.
        """
        r = None
        data_passage = data.get("passage", "")
        for e_no, e in enumerate(entry["data"]):
            sequence_match = self._sequences_match(master=e["nuc"], s=data["sequence"])
            if sequence_match:
                if data_passage and data_passage not in e["passages"]:
                    # if e["passages"]:
                    #     module_logger.warning('[SAMESEQ] different passages {!r} {!r} {!r}'.format(data["name"], e["passages"], data_passage))
                    r = {"type": "update", "sequence_match": sequence_match, "index": e_no}
                elif e.get("gene") and data.get("gene") and e["gene"] != data["gene"]:
                    r = {"type": "different-genes", "message": '[SAMESEQ] different genes {!r} {!r} {!r}'.format(data["name"], e["gene"], data["gene"])}
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
        if data.get("passage") and data["passage"] not in entry_passage["passages"]:
            entry_passage["passages"].append(data["passage"])
        lab_e = entry_passage["labs"].setdefault(data["lab"], [])
        if data.get("lab_id") and data["lab_id"] not in lab_e:
            lab_e.append(data["lab_id"])
        if data.get("gene") and not entry_passage.get("gene"):
            entry_passage["gene"] = data["gene"]
        if sequence_match in ["super", "new"]:   # update sequences with the longer one
            aa = amino_acids.translate_sequence_to_amino_acid(data["sequence"])
            self._align(aa, entry_passage, data, db_entry)   # may raise
            entry_passage["nuc"] = data["sequence"]
            entry_passage["aa"] = aa

    def _update_cdcid(self, data):
        if data["lab"] == "CDC" and data.get("lab_id"):
            if len(data["lab_id"]) >= 8:
                existing = self.cdcids.get(data["lab_id"])
                if existing is None:
                    self.cdcids[data["lab_id"]] = data["name"]
                elif existing != data["name"]:
                    module_logger.warning('[CDCID] {!r} for another name: {!r} existing: {!r}'.format(data["lab_id"], data["name"], existing))
            else:
                module_logger.warning('[CDCID] too short (ignored) {!r} {!r}'.format(data["lab_id"], data["name"]))

    def _align(self, sequence, entry_passage, data, db_entry):
        try:
            aligment_data = amino_acids.align(sequence)
        except amino_acids.SequenceIsTooShort as err:
            raise Exclude(str(err))
        if aligment_data:
            if aligment_data["virus_type"] != db_entry["virus_type"]:
                module_logger.warning('Virus type detection mismatch {} vs. {}'.format(aligment_data, data))
            if aligment_data.get("lineage"):
                if db_entry.get("lineage"):
                    if db_entry["lineage"] != aligment_data["lineage"]:
                        module_logger.warning('Lineage detection mismatch for {}: {} vs. {}'.format(data["name"], aligment_data, db_entry["lineage"]))
                else:
                    db_entry["lineage"] = aligment_data["lineage"]
            if entry_passage.get("gene"):
                if entry_passage["gene"] != aligment_data["gene"]:
                    module_logger.warning('Gene detection mismatch for {}: {} vs. {}'.format(data["name"], aligment_data, entry_passage["gene"]))
            else:
                entry_passage["gene"] = aligment_data["gene"]
            entry_passage["shift"] = aligment_data["shift"]
        elif db_entry["virus_type"] in ["A(H3N2)", "A(H1N1)"]:
            module_logger.warning('Not aligned {} len:{} {}'.format(data["name"], len(sequence), sequence))

    def _aligned(self, entry):
        s = entry["aa"]
        shift = entry.get("shift")
        if shift is None:
            raise RuntimeError("Not aligned")
        if shift < 0:
            s = s[-shift:]
        elif shift > 0:
            s = ("X" * shift) + s
        return s

        # --------------------------------------------------

    # for report
    def _multiple_sequences_per_passage(self, entry):
        passages = {}
        for e in entry["data"]:
            for p in e["passages"]:
                passages.setdefault(p, []).append(e.get("gene"))
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

# ======================================================================
### Local Variables:
### eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
### End:
