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

# self.names is dict {<name>: {"passages": {<passage>: [{"sequence": <sequence>, "labs": [<lab>], "lab_ids": [], "gene": <HA|NA>}, ...], "virus_type": <virus_type>, "dates": []}, ...}
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

        multi_seq_per_passage = [(n, p) for n, e1 in self.names.items() for p, e2 in e1["passages"].items() if len(e2) > 1]
        if multi_seq_per_passage:
            print("Multiple sequences per passage {}:\n\t{}".format(len(multi_seq_per_passage), "\n\t".join("{!r} {!r}".format(n, p) for n, p in multi_seq_per_passage)))

        # --------------------------------------------------

    def _add_sequence(self, data):
        if data.get("name"):
            entry = self.names.get(data["name"])
            if entry is None:
                entry = {"passages": {}, "virus_type": data["virus_type"], "dates": []}
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
        entry_passages = entry["passages"]
        if sameseq is None:
            new_entry_passage = {"sequence": data["sequence"], "labs": [data["lab"]], "lab_ids": []}
            if data.get("gene"):
                new_entry_passage["gene"] = data["gene"]
            if data.get("lab_id"):
                new_entry_passage["lab_ids"].append(data["lab_id"])
            entry_passages.setdefault(data.get("passage", ""), []).append(new_entry_passage)
        elif sameseq["type"] == "update-empty-entry-passage":
            # data has passage and matches seq in entry without passage, update entry
            entry_passage = entry_passages[""][sameseq["index"]]
            del entry_passages[""][sameseq["index"]]
            if not entry_passages[""]:
                del entry_passages[""]
            self._update_entry_passage(entry_passage=entry_passage, data=data, sequence_match=sameseq["sequence_match"])
            entry_passages.setdefault(data["passage"], []).append(entry_passage)

    def _look_for_the_same_sequence(self, entry, data):
        """If no matching sequence found returns None. Else updates
        passage entry, if passage is the same. Returns string describing
        match type.
        """
        r = None
        data_passage = data.get("passage", "")
        data_gene = data.get("gene")
        for passage, e1 in entry["passages"].items():
            for e2_no, e2 in enumerate(e1):
                sequence_match = self._sequences_match(master=e2["sequence"], s=data["sequence"])
                if sequence_match:
                    if data_passage and passage != data_passage:
                        if passage:
                            module_logger.warning('[SAMESEQ] different passages {!r} {!r} {!r}'.format(data["name"], passage, data_passage))
                            r = {"type": "different-passages"}
                        else:
                            r = {"type": "update-empty-entry-passage", "sequence_match": sequence_match, "index": e2_no}
                    elif e2.get("gene") and data_gene and e2["gene"] != data["gene"]:
                        module_logger.warning('[SAMESEQ] different genes {!r} {!r} {!r}'.format(data["name"], e2["gene"], data["gene"]))
                        r = {"type": "different-genes"}
                    else:
                        self._update_entry_passage(entry_passage=e2, data=data, sequence_match=sequence_match)
                        # if data.get("lab_id") and data["lab_id"] not in e2["lab_ids"]:
                        #     e2["lab_ids"].append(data["lab_id"])
                        # if data["lab"] not in e2["labs"]:
                        #     e2["labs"].append(data["lab"])
                        # if data_gene and not e2.get("gene"):
                        #     e2["gene"] = data["gene"]
                        # if m == "super":   # update sequences with the longer one
                        #     e2["sequence"] = data["sequence"]
                        r = {"type": "good"}
                if r:
                    break
            if r:
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
        if data.get("lab_id") and data["lab_id"] not in entry_passage["lab_ids"]:
            entry_passage["lab_ids"].append(data["lab_id"])
        if data["lab"] not in entry_passage["labs"]:
            entry_passage["labs"].append(data["lab"])
        if data.get("gene") and not entry_passage.get("gene"):
            entry_passage["gene"] = data["gene"]
        if sequence_match == "super":   # update sequences with the longer one
            entry_passage["sequence"] = data["sequence"]

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
