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
        print("Strange names:\n\t", "\n\t".join(n for n in self.names if n.count("/") < 3))
        print("Multiple sequences per passage:\n\t", "\n\t".join("{!r} {!r}".format(n, p) for n, e1 in self.names.items() for p, e2 in e1["passages"].items() if len(e2) > 1), sep="")

        # --------------------------------------------------

    def _add_sequence(self, data):
        if data.get("name"):
            entry = self.names.get(data["name"])
            if entry is None:
                entry = {"passages": {}, "virus_type": data["virus_type"], "dates": []} # "lab": data["lab"],
                self.names[data["name"]] = entry
            self._update_db_entry(entry, data)
        else:
            module_logger.warning('Entry without name: {}'.format(data["lab_id"]))

    def _update_db_entry(self, entry, data):
        if entry["virus_type"] != data["virus_type"]:
            raise RuntimeError("Cannot add {!r} to {!r} db entry\n{}\n{}".format(data["virus_type"], entry["virus_type"], data, entry))
        if data.get("date") and data["date"] not in entry["dates"]:
            entry["dates"].append(data["date"])
        passage = data.get("passage", "")
        entry_passage = entry["passages"].setdefault(passage, [])
        ignore = False
        if not ignore:
            new_entry_passage = {"sequence": data["sequence"], "labs": [data["lab"]], "lab_ids": []}
            if data.get("gene"):
                new_entry_passage["gene"] = data["gene"]
            if data.get("lab_id"):
                new_entry_passage["lab_ids"].append(data["lab_id"])
            entry_passage.append(new_entry_passage)

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
