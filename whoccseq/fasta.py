# -*- Python -*-
# license
# license.

"""
Functions for reading and generating fasta files.
"""

import os, re, csv
import logging; module_logger = logging.getLogger(__name__)
from . import open_file

# ======================================================================

def generate(filename, data, names_order=None, predicate=None, sequence_access=None, split_at=75):
    module_logger.info('generating fasta into {}'.format(filename))
    with open_file.open_for_writing_binary(filename) as f:
        if isinstance(data, dict):
            _generate_from_dict(f=f, data=data, names_order=names_order or sorted(data), predicate=predicate or (lambda n, s: True), sequence_access=sequence_access or (lambda a: a), split_at=split_at)
        elif isinstance(data, list) and len(data) and isinstance(data[0], (list, tuple)) and len(data[0]) == 2:
            _generate_from_list(f, data, split_at=split_at)
        else:
            raise ValueError("Cannot write {} to fasta (dict or list of (name, seq) pairs expected".format(type(data)))

def _write_seq_to_file(f, name, seq, split_at):
    f.write(b">")
    f.write(encode_name(name).encode('utf-8'))
    f.write(b"\n")
    f.write(sequence_split(seq, chunk_len=split_at).encode('utf-8'))
    f.write(b"\n")

def _generate_from_dict(f, data, names_order, predicate, sequence_access, split_at):
    for name in names_order:
        seq = data[name]
        if predicate(name, seq):
            _write_seq_to_file(f, name, sequence_access(seq), split_at)

def _generate_from_list(f, data, split_at):
    for name, seq in data:
        _write_seq_to_file(f, name, seq, split_at)

# ----------------------------------------------------------------------

def generate_phylip(filename, data):
    if not isinstance(data, list) or not isinstance(data[0], (list, tuple)) or len(data[0]) != 2:
        raise ValueError("Invalid data passed to generate_phylip: list of tuples with two elements (name, sequence) expected")
    entries = len(data)
    max_name = max(len(entry[0]) for entry in data)
    max_seq = max(len(entry[1]) for entry in data)
    with open(filename, 'w') as f:
        f.write('{} {}\n'.format(entries, max_seq))
        for entry in data:
            f.write("{:<{}s} {:-<{}s}\n".format(entry[0], max_name, entry[1], max_seq))

# ----------------------------------------------------------------------

# use URL-kind encoding to avoid ambiguity
def encode_name(name):
    for char in "% :()!*';@&=+$,?#[]":
        name = name.replace(char, '%{:02X}'.format(ord(char)))
    return name

sReNameDecoder = re.compile(r'%([0-9A-Fa-f][0-9A-Fa-f])')

def decode_name(name):
    def replace(match):
        return chr(int(match.group(1), 16))
    return sReNameDecoder.sub(replace, name)

# ----------------------------------------------------------------------

# def encode_name_old(name):
#     return name.replace(' ', '^').replace(':', '%').replace('(', '<').replace(')', '>')

# def decode_name_old(name):
#     return name.replace('^', ' ').replace('%', ':').replace('<', '(').replace('>', ')')

# ----------------------------------------------------------------------

def sequence_split(sequence, chunk_len=75, separator="\n"):
    if chunk_len and chunk_len > 0:
        return separator.join(sequence[i:i+chunk_len] for i in range(0, len(sequence), chunk_len))
    else:
        return sequence

# ----------------------------------------------------------------------

sReSubtype = re.compile(r'^(?:B/|A\((H3|H1))')

def detect_subtype(sequences):
    for name in sequences:
        m = sReSubtype.match(name)
        if m:
            subtype = m.group(1) or "B"
            logging.info('Subtype: {} ({})'.format(subtype, name))
            return subtype
    raise RuntimeError("cannot detect subtype by sequence names")

# ----------------------------------------------------------------------

def read(filename, duplicates="merge", report_duplicates=True, decode_names=True, uppercase_sequences=True):
    """Returns dict {name: sequence}"""
    return FastaReaderBasicDict(decode_names=decode_names, uppercase_sequences=uppercase_sequences).read(filename, duplicates=duplicates, report_duplicates=report_duplicates)

def read_as_list(filename, decode_names=True, uppercase_sequences=True):
    """Returns list of tuples (name, sequence) keeping the order from filename."""
    return FastaReaderBasicList(decode_names=decode_names, uppercase_sequences=uppercase_sequences).read(filename)

# ======================================================================

class FastaReaderError (Exception):
    pass

# ======================================================================

class FastaReaderBasic:

    def __init__(self, uppercase_sequences=False):
        self.uppercase_sequences = uppercase_sequences

    def read(self, filename, **kwargs):
        return self.read_entries(self.open(filename), filename=os.path.basename(filename))

    def open(self, filename):
        self.filename = filename
        source = open_file.open_for_reading_binary(filename)
        return source

    def read_entries(self, source, filename):
        sequence = ''
        raw_name = None
        for line_no, line in enumerate(source, start=1):
            try:
                line = line.decode('utf-8').strip()
            except:
                raise FastaReaderError('{filename}:{line_no}: cannot decode line'.format(filename=filename, line_no=line_no))
            if not line or line[0] == ';':               # empty or comment line
                pass
            elif line[0] == '>':
                if raw_name:
                    if not sequence:
                        raise FastaReaderError('{filename}:{line_no}: {raw_name!r} without sequence'.format(raw_name=raw_name, filename=filename, line_no=line_no))
                    sequence = self.normalize_sequence(sequence)
                    self.check_sequence(sequence, line_no)
                    yield {'raw_name': raw_name, 'sequence': sequence, 'source': filename}
                sequence = ''
                raw_name = line[1:]
            else:
                if not raw_name:
                    raise FastaReaderError('{filename}:{line_no}: sequence without name'.format(filename=filename, line_no=line_no))
                sequence += line
        if raw_name:
            if not sequence:
                raise FastaReaderError('{filename}:EOF: {raw_name!r} without sequence'.format(raw_name=raw_name, filename=filename))
            sequence = self.normalize_sequence(sequence)
            self.check_sequence(sequence, line_no)
            yield {'raw_name': raw_name, 'sequence': sequence, 'source': filename}

    sReSequence = re.compile(r"^[A-Za-z\-~:\*\.]+$")

    def check_sequence(self, sequence, line_no):
        if not self.sReSequence.match(sequence):
            raise FastaReaderError('{filename}:{line_no}: invalid sequence read: {sequence}'.format(sequence=sequence, filename=self.filename, line_no=line_no))

    def normalize_sequence(self, sequence):
        sequence = sequence.replace("/", "-")   # / found in H1pdm sequences
        if self.uppercase_sequences:
            sequence = sequence.upper()
        return sequence

# ======================================================================

class FastaReaderBasicDict (FastaReaderBasic):
    """Gets dict of raw_name to sequence entries"""

    def __init__(self, decode_names=True, uppercase_sequences=False):
        super().__init__(uppercase_sequences=uppercase_sequences)
        self.decode_names = decode_names

    def read(self, filename, duplicates="merge", report_duplicates=True, **kwargs):
        result = {}
        for entry in self.read_entries(self.open(filename), filename=os.path.basename(filename)):
            self.add_entry(result, entry, duplicates, report_duplicates)
        return result

    def add_entry(self, result, entry, duplicates, report_duplicates):
        if self.decode_names:
            name = decode_name(entry['raw_name'])
        else:
            name = entry['raw_name']
        if name in result and not self.equal(entry['sequence'], result[name]):
            if duplicates == "error":
                raise FastaReaderError("Duplicating name: {}, different sequences.".format(name))
            elif duplicates == "merge":
                result[name].append(entry['sequence'])
            if report_duplicates:
                module_logger.warning("Duplicating name: {}, different sequences.".format(name))
        else:
            if duplicates == "merge":
                result[name] = [entry['sequence']]
            else:
                result[name] = entry['sequence']

    def equal(self, new, existing):
        if isinstance(existing, str):
            r = new == existing
        elif isinstance(existing, list):
            r = any(new == e for e in existing)
        else:
            raise ValueError("Unsupported {}".format(type(existing)))
        return r

# ======================================================================

class FastaReaderBasicList (FastaReaderBasic):
    """Gets dict of raw_name to sequence entries"""

    def __init__(self, decode_names=True, uppercase_sequences=False):
        super().__init__(uppercase_sequences=uppercase_sequences)
        self.decode_names = decode_names

    def read(self, filename, **kwargs):
        def name(entry):
            if self.decode_names:
                return decode_name(entry['raw_name'])
            else:
                return entry['raw_name']
        return [(name(entry), entry['sequence']) for entry in self.read_entries(self.open(filename), filename=os.path.basename(filename))]

# ======================================================================

class FastaReaderWithNameNormalizing (FastaReaderBasic):

    def __init__(self, lab=None, virus_type=None, uppercase_sequences=False):
        super().__init__(uppercase_sequences=uppercase_sequences)
        self.lab = lab
        # from ..viruses import VirusesAPI
        # self.viruses_api = VirusesAPI(lab=self.lab, virus_type=virus_type)

    # def normalize_name(self, raw_name, name_type='whocc_fasta_match'):
    #     name, location = self.viruses_api.normalize_name(raw_name, name_type=name_type)
    #     return {k: v for k, v in (('name', name), ('country', location and location.country), ('continent', location and location.continent), ('latitude', location and location.latitude), ('longitude', location and location.longitude)) if v}

    def normalize_name(self, raw_name, name_type='whocc_fasta_match'):
        return {"name": raw_name}

    def normalize_date(self, date):
        return date
        # try:
        #     return date and parsers.datetime_parse_to_string(date, include_time=False)
        # except parsers.ParsingError as err:
        #     module_logger.warning(err)
        #     return None

    def normalize_passage(self, passage, report_prefix):
        return passage
        # try:
        #     normalized = passage and parsers.passage_parse_to_string(lab=self.lab, source=passage, report_prefix='[PASSAGE] ' + (report_prefix or ''))
        # except parsers.ParsingError as err:
        #     module_logger.error('Passage {!r} parsing: {}'.format(passage, err))
        #     normalized = passage
        # return normalized

    def normalize_lab_id(self, lab_id):
        if self.lab == 'CDC' and lab_id and len(lab_id) > 10 and lab_id[10] == '_':
            lab_id = lab_id[:10]
        return lab_id

# ======================================================================

class FastaReaderWithNameParsing (FastaReaderWithNameNormalizing):

    sParsers = (
        #(re.compile(r'^(?P<type>B|H3|H1)(?P<location>[A-Z]+)(?P<isolation_number>\d+)(?P<year>\d\d)[A-Z]*\s', re.I), 'nimr_glued'),
        (re.compile(r'^(?P<name>[^\s]+)\s+PileUp\sof', re.I), 'nimr_20090914'),
        (re.compile(r'^(?P<designation>[^|]+)\s+\|\s+(?P<passage>[^\s]+)\s*\|\s+(?P<name>(?:EPI|201)\d+)\s*$', re.I), 'cdc_20100913'), # name | passage | fasta_id
        (re.compile(r'^(?P<name1>EPI\d+)\s+\|\s+(?P<gene>HA|NA)\s+\|\s+(?P<designation>[^\|]+)\s+\|\s+(?P<name>EPI_[A-Z_0-9]+)\s+\|\s*(?P<passage>[^\s]+)?\s*\|\s*(?P<flu_type>.+)?\s*$', re.I), 'gisaid_melb'), # name1 | gene | designation | name | passage | flu_type
        (re.compile(r'^(?P<name>[^|]+)\s+\|\s+(?:(?P<year1>\d+)-(?P<month1>\d+)-(?P<day1>\d+)|(?P<year2>\d+)-(?P<month2>\d+)\s+\(day unknown\)|(?P<year3>\d+)\s+\(month and day unknown\))\s+\|\s+(?P<passage>[^\|]*)\s+\|\s+(?P<lab_id>[^\|]*)?\s+\|\s+(?P<lab>[A-Za-z ]+)\s*$', re.I), 'gisaid'), # name | date | passage | lab_id | lab
        (re.compile(r'^(?P<name>[^|]+)\s+\|\s+(?:(?P<year1>\d+)-(?P<month1>\d+)-(?P<day1>\d+)|(?P<year2>\d+)-(?P<month2>\d+)\s+\(day unknown\)|(?P<year3>\d+)\s+\(month and day unknown\))\s+\|\s+(?P<passage>[^\|]+)\s+\|\s+(?P<lab_id>[^\|]+)?\s+\|.*$', re.I), 'gisaid'), # name | date | passage | lab_id | something-else
        (re.compile(r'^(?P<name>[^|]+)\s+\|\s+(?:(?P<year1>\d+)-(?P<month1>\d+)-(?P<day1>\d+)|(?P<year2>\d+)-(?P<month2>\d+)\s+\(day unknown\)|(?P<year3>\d+)\s+\(month and day unknown\))\s+\|\s+(?P<passage>[^\|]+)?\s*\|\s*(?P<lab_id>.+)?\s*$', re.I), 'gisaid'), # name | date | passage? | lab_id
        (re.compile(r'^(?P<name>[^|]+)\s+\|\s+(?P<passage>[^\s]+)\s+\|\s*(?P<lab_id>.+)?\s*$', re.I), 'gisaid_without_date'), # name | passage | lab_id
        (re.compile(r'^(?P<name>[^\s]+)\s\s+(?P<date>[\d/]+)\s\s+(?P<passage>[^\s]+)\s*$', re.I), 'melb_20100823'), # name  date  passage
        (re.compile(r'^(?P<name>\d+S\d+)\s+\"Contig\s+\d+\"\s+\(\d+,\d+\)$', re.I), 'melb_20110921'),
        (re.compile(r'^(?P<name>[^_]+/\d\d\d\d)[\s_]+[^/]*(?:[\s_]?:(?P<passage>[EC]\d*))[^/]*$', re.I), 'name_passage'), # CNIC
        (re.compile(r'^(?P<name>[^_]+/\d\d\d\d)[\s_]+(?:Jan|Feb|Mar|Apr|may|Jun|Jul|Aug|Sep|Oct|Nov|Dec)?[^/]*$', re.I), 'name_only'), # CNIC
        (re.compile(r'^(?P<name>[^\s]+_4)$', re.I), 'name_only'),   # CDC
        (re.compile(r'^(?P<name>[^_]+)[\s_]+[^/]*$', re.I), 'name_only'), # CNIC
        (re.compile(r'^(?P<name>.+/\d\d\d\d)[\s\-]*(?P<passage>[\.\dA-Z]+)?$', re.I), 'name_passage'), # NIID
        (re.compile(r'.'), 'simple'),
        )

    def read(self, filename, name_type='whocc_fasta_match', **kwargs):
        report_prefix = '[{}]'.format(os.path.basename(filename))
        return (e for e in (self.postprocess(self.name_parse(entry, report_prefix=report_prefix), report_prefix=report_prefix, name_type=name_type) for entry in self.read_entries(self.open(filename), filename=os.path.basename(filename))) if e)

    def postprocess(self, entry, report_prefix, name_type):
        if entry.get('name'):
            entry.update(self.normalize_name(entry['name'], name_type=name_type))
        if entry.get('lab_id'):
            entry['lab_id'] = self.normalize_lab_id(entry['lab_id'])
        # if self.lab == 'CDC' and not entry.get('lab_id'):
        #     module_logger.warning('[NOCDCID] {}: entry without lab_id removed: {}'.format(report_prefix or '', entry['raw_name'].strip()))
        #     entry = None
        return entry

    def name_parse(self, entry, report_prefix=None):
        for rex, func_name in self.sParsers:
            m = rex.match(entry['raw_name'])
            if m:
                # if func_name != 'simple':
                #     module_logger.info('{}'.format(func_name))
                entry.update(getattr(self, func_name)(entry['raw_name'], m, report_prefix=report_prefix))
                break
        # module_logger.debug('name_parse {}'.format(entry))
        return entry

    def simple(self, raw_name, m, report_prefix):
        return {'name': raw_name}

    # def nimr_glued(self, raw_name, m, report_prefix=None):
    #     def convert_type(t):
    #         if re.match(r'^H\d\d?$', t):
    #             t = 'A'
    #         return t
    #     return {'name': '/'.join((convert_type(m.group('type')), m.group('location'), m.group('isolation_number'), m.group('year')))}

    def nimr_20090914(self, raw_name, m, report_prefix):
        return {'name': m.group('name').upper()}

    def gisaid(self, raw_name, m, report_prefix=None, with_date=True):
        # module_logger.info('gisaid with_date:{} {} --> {}'.format(with_date, raw_name, m.groups()))
        year = (with_date and (m.group('year1') or m.group('year2') or m.group('year3'))) or None
        try:
            lab = self._fix_gisaid_lab(m.group('lab'))
        except IndexError:
            lab = None
        return {k: v for k, v in (('name', m.group('name')),
                                  ('date', year and self.normalize_date('-'.join((year, m.group('month1') or m.group('month2') or '01', m.group('day1') or '01')))),
                                  ('passage', self.normalize_passage(m.group('passage'), report_prefix=report_prefix)),
                                  ('lab_id', m.group('lab_id')),
                                  ('lab', lab),
                                  ) if v is not None}

    def _fix_gisaid_lab(self, lab):
        return lab.replace("Centers for Disease Control and Prevention", "CDC").replace("Crick Worldwide Influenza Centre", "NIMR").replace("WHO Collaborating Centre for Reference and Research on Influenza", "MELB")

    def gisaid_without_date(self, raw_name, m, report_prefix):
        return self.gisaid(raw_name=raw_name, m=m, report_prefix=report_prefix, with_date=False)

    def gisaid_melb(self, raw_name, m, report_prefix):
        return {'name': m.group('name'), 'gene': m.group('gene')}

    def melb_20100823(self, raw_name, m, report_prefix):
        return {'name': m.group('name').upper(), 'date': self.normalize_date(m.group('date')), 'passage': self.normalize_passage(m.group('passage'), report_prefix=report_prefix)}

    def melb_20110921(self, raw_name, m, report_prefix):
        return {'name': m.group('name').upper(), '_extract_options': {'shared_fasta': True}}

    def name_only(self, raw_name, m, report_prefix):
        return {'name': m.group('name').upper()}

    def cdc_20100913(self, raw_name, m, report_prefix):
        return {'name': m.group('name').upper(), 'passage': self.normalize_passage(m.group('passage'), report_prefix=report_prefix)}

    def name_passage(self, raw_name, m, report_prefix):
        return {'name': m.group('name').upper(), 'passage': self.normalize_passage(m.group('passage'), report_prefix=report_prefix)}

# ======================================================================

class FastaReaderWithCSV (FastaReaderWithNameParsing):
    """Gets sequences from fasta file, extracts names and passages from csv file. First line of the csv file is a list of fields, extracted field names: FASTA_ID, NAME, PASSAGE, DATE, GENE."""

    def read(self, filename, csv_filename, name_type='whocc_fasta_match', **kwargs):
        self.csv_fields = None
        report_prefix = '[{}]'.format(os.path.basename(filename))
        name_data = {entry['fasta_id']: entry for entry in (self.read_csv_entry(csv_entry, csv_filename) for csv_entry in csv.reader(open_file.open_for_reading_text(csv_filename))) if entry}
        return (e for e in (self.name_extract_from_csv(entry=entry, name_data=name_data, report_prefix=report_prefix, name_type=name_type) for entry in self.read_entries(self.open(filename), filename=os.path.basename(filename)) if entry) if e)

    def name_extract_from_csv(self, entry, name_data, report_prefix=None, name_type='whocc_fasta_match'):
        entry = self.name_parse(entry, report_prefix=report_prefix)
        name_entry = self.find_name_entry(name_data, entry)
        if name_entry:
            entry.update(name_entry)
            entry.update(self.normalize_name(entry['name'], name_type=name_type))
        else:
            entry = self.no_csv_entry(entry, report_prefix=report_prefix)
        return entry

    def find_name_entry(self, name_data, entry):
        return name_data.get(entry['name'])

    def no_csv_entry(self, entry, report_prefix):
        extract_options = entry.pop('_extract_options', {})
        if extract_options.get('shared_fasta'):
            entry = None
        else:
            module_logger.warning('[NOCSVENTRY] {}: no entry in csv file for {!r} {!r}'.format(report_prefix or '', entry['name'], entry['raw_name'].strip()))
        return entry

    def read_csv_entry(self, csv_entry, csv_filename):
        result = None
        if not self.csv_fields:
            if self.lab == 'CDC':
                expected_fields = ('fasta_id', 'name', 'passage')
            elif self.lab == 'CNIC':
                expected_fields = ('fasta_id', 'name')
            else:
                expected_fields = ('fasta_id', 'name', 'passage', 'date')
            self.csv_fields = {re.sub(r'^#', '', n.strip().lower()): index for index, n in enumerate(csv_entry)}
            not_found = [expected for expected in expected_fields if expected not in self.csv_fields]
            if not_found:
                raise FastaReaderError('[{}]: CSV file does not have mandatory fields: {}'.format(csv_filename, not_found))
        elif not csv_entry[0].strip() or csv_entry[0].strip()[0] != '#':   # ignore commented out entries
            result = {field: csv_entry[index].strip() if len(csv_entry) > index else None for field, index in self.csv_fields.items()}
            if result.get('passage'):
                result['passage'] = self.normalize_passage(result['passage'], report_prefix='[{}]'.format(os.path.basename(csv_filename))) or None
            if result.get('date'):
                result['date'] = self.normalize_date(result['date']) or None
            result['fasta_id'] = result['fasta_id'].upper()
        return result

    # overrides
    def name_only(self, raw_name, m, report_prefix):
        if self.lab == 'CDC':
            return {'name': raw_name.upper()}
        else:
            return super().name_only(raw_name, m, report_prefix)

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

    def no_csv_entry(self, entry, report_prefix):
        module_logger.warning('[NOCSVENTRY] {}: no entry in csv file for {!r} {!r}'.format(report_prefix or '', entry['name'], entry['raw_name'].strip()))
        return None

        # if (len(name) in (13, 14) and name[10:12] == 'EG') or (len(name) == 11 and name[10] in ('F', 'R')):
        #     name = name[:10]

# ======================================================================
### Local Variables:
### eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
### End:
