# -*- Python -*-
# license
# license.

"""
Functions for reading and generating fasta files (old formats).
"""

import os, re, csv
import logging; module_logger = logging.getLogger(__name__)
from . import open_file, fasta as fasta_m

# ----------------------------------------------------------------------

def read_fasta_with_csv(fasta_file, csv_file, lab, virus_type, **_):
    if lab == "CDC":
        reader = FastaReaderWithCSV_CDC(lab=lab, virus_type=virus_type)
    else:
        reader = FastaReaderWithCSV(lab=lab, virus_type=virus_type)
    return reader.read(fasta_f=open_file.open_for_reading_text(fasta_file), csv_f=open_file.open_for_reading_text(csv_file))

# ----------------------------------------------------------------------

class FastaReaderWithCSV:
    """Gets sequences from fasta file, extracts names and passages from csv file. First line of the csv file is a list of fields, extracted field names: FASTA_ID, NAME, PASSAGE, DATE, GENE."""

    def __init__(self, lab, virus_type):
        self.lab = lab
        self.virus_type = virus_type

    def read(self, fasta_f, csv_f):
        self.csv_name = os.path.basename(csv_f.name)
        self.csv_fields = None
        self.not_found_in_csv = 0
        name_data = {entry['fasta_id']: entry for entry in (self.read_csv_entry(csv_entry, csv_f.name) for csv_entry in csv.reader(csv_f)) if entry}
        r = [self._update_entry(e) for e in (self.name_extract_from_csv(raw_name=raw_name, sequence=sequence, name_data=name_data) for raw_name, sequence in fasta_m.read_from_string(fasta_f.read(), fasta_f.name)) if e]
        module_logger.debug('{} sequences imported ({} ignored) from {} {}'.format(len(r), self.not_found_in_csv, fasta_f.name, self.csv_name))
        return r

    def name_extract_from_csv(self, raw_name, sequence, name_data):
        entry = fasta_m.name_parser().parse(raw_name.strip())
        name_entry = self.find_name_entry(name_data, entry)
        if name_entry:
            entry.update(name_entry)
            entry["sequence"] = sequence
            entry["source"] = self.csv_name
            # module_logger.info('entry for {}'.format(entry['name']))
        else:
            if raw_name[-2:] not in ["_1", "_2", "_3", "_5", "_6", "_7", "_8"]:
                module_logger.debug('[NOCSVENTRY]: no entry in {!r} for {!r} {!r}'.format(self.csv_name, entry['name'], raw_name))
            self.not_found_in_csv += 1
            entry = None
        return entry

    def find_name_entry(self, name_data, entry):
        return name_data.get(entry['name'])

    sExpectedFields = {"CDC": ['fasta_id', 'name', 'passage'], "CNIC": ['fasta_id', 'name'], None: ['fasta_id', 'name', 'passage', 'date']}

    def read_csv_entry(self, csv_entry, csv_filename):
        result = None
        if not self.csv_fields:
            self.csv_fields = {re.sub(r'^#', '', n.strip().lower()): index for index, n in enumerate(csv_entry)}
            not_found = [expected for expected in (self.sExpectedFields.get(self.lab) or self.sExpectedFields[None]) if expected not in self.csv_fields]
            if not_found:
                raise FastaReaderError('[{}]: CSV file does not have mandatory fields: {}'.format(csv_filename, not_found))
        elif not csv_entry[0].strip() or csv_entry[0].strip()[0] != '#':   # ignore commented out entries
            result = {field: csv_entry[index].strip() if len(csv_entry) > index else None for field, index in self.csv_fields.items()}
            result['fasta_id'] = result['fasta_id'].upper()
        return result

    # overrides
    def name_only(self, raw_name, m, report_prefix):
        if self.lab == 'CDC':
            return {'name': raw_name.upper()}
        else:
            return super().name_only(raw_name, m, report_prefix)

    def _update_entry(self, entry):

        for ft, f in (("lab_id", ["lab num", "lab #", "labnumber", "lab id"]), ("passage", ["passage history", "passage history-of received"])):
            if ft not in entry:
                for ff in f:
                    if ff in entry:
                        entry[ft] = entry[ff]
                        break

        e = {k: v for k, v in entry.items() if k in ["lab_id", "date", "passage", "gene", "sequence", "name", "source"] and v}
        e.update(lab=self.lab, virus_type=self.virus_type)
        return e

# ======================================================================

class FastaReaderWithCSV_CDC (FastaReaderWithCSV):
    """Additonal stuff to match more csv entries for brain-damaged CDC files."""

    def read_csv_entry(self, csv_entry, csv_filename):
        result = super().read_csv_entry(csv_entry, csv_filename)
        if result and not result['fasta_id'] and result['lab_id']:
            result['fasta_id'] = result['lab_id']
        return result

    def find_name_entry(self, name_data, entry):
        def remove_4(s):
            if s[-2:] == '_4':
                return s[:-2]
            else:
                return s

        def remove_eg(s):
            return re.sub(r'EG\d[FR]?_4$', '', s)

        def remove_fr(s):
            return re.sub(r'[FR]_4$', '', s)

        return name_data.get(entry['name']) or name_data.get(remove_4(entry['name'])) or name_data.get(remove_eg(entry['name'])) or name_data.get(remove_fr(entry['name']))

# ----------------------------------------------------------------------

# ======================================================================
### Local Variables:
### eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
### End:
