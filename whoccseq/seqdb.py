# -*- Python -*-
# license
# license.

"""
Class to access sequence database
"""

import os, sys, json
import logging; module_logger = logging.getLogger(__name__)
from . import open_file, acmacs

# ----------------------------------------------------------------------

# self.names is dict {<name>: {"data": [{"passages": [<passage>], "sequence": <sequence>, "labs": {<lab> :[<lab_id>]}, "gene": <HA|NA>}, ...], "virus_type": <virus_type>, "dates": [<date>]}, ...}

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
        for entry in data:
            self._add_sequence(entry)

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
        open_file.write_json(self.path_to_db, data, indent=1, sort_keys=True)

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

        # --------------------------------------------------

    def _add_sequence(self, data):
        if data.get("name"):
            entry = self.names.get(data["name"])
            if entry is None:
                entry = {"data": [], "virus_type": data["virus_type"], "dates": []}
                self.names[data["name"]] = entry
            self._update_db_entry(entry, data)
        else:
            module_logger.warning('Entry without name: {}'.format(data["lab_id"]))

    def _update_db_entry(self, entry, data):
        if entry["virus_type"] != data["virus_type"]:
            raise RuntimeError("Cannot add {!r} to {!r} db entry\n{}\n{}".format(data["virus_type"], entry["virus_type"], data, entry))
        if data.get("date") and data["date"] not in entry["dates"]:
            entry["dates"].append(data["date"])
        sameseq = self._look_for_the_same_sequence(entry, data)
        if sameseq is None:
            new_entry = {"labs": {}, "passages": []}
            self._update_entry_passage(entry_passage=new_entry, data=data, sequence_match="new")
            entry["data"].append(new_entry)
        elif sameseq["type"] == "update":
            self._update_entry_passage(entry_passage=entry["data"][sameseq["index"]], data=data, sequence_match=sameseq["sequence_match"])
        elif sameseq["type"] == "different-genes":
            module_logger.warning(sameseq["message"])
        else:
            module_logger.error("[INTERNAL] Unrecoginzed sameseq: {}".format(sameseq))

    def _look_for_the_same_sequence(self, entry, data):
        """If no matching sequence found returns None. Else updates
        passage entry, if passage is the same. Returns string describing
        match type.
        """
        r = None
        data_passage = data.get("passage", "")
        for e_no, e in enumerate(entry["data"]):
            sequence_match = self._sequences_match(master=e["sequence"], s=data["sequence"])
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

    def _update_entry_passage(self, entry_passage, data, sequence_match):
        if data.get("passage") and data["passage"] not in entry_passage["passages"]:
            entry_passage["passages"].append(data["passage"])
        lab_e = entry_passage["labs"].setdefault(data["lab"], [])
        if data.get("lab_id") and data["lab_id"] not in lab_e:
            lab_e.append(data["lab_id"])
        if data.get("gene") and not entry_passage.get("gene"):
            entry_passage["gene"] = data["gene"]
        if sequence_match in ["super", "new"]:   # update sequences with the longer one
            entry_passage["sequence"] = data["sequence"]

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
