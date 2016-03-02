# -*- Python -*-
# license
# license.
"""
Files opening and reading/writing functions.
"""

import sys, os, bz2, lzma, json
import logging; module_logger = logging.getLogger(__name__)

# ======================================================================

def open_for_reading_binary(filename):
    """Opens binary file for reading. Handles compressed files transparently."""
    return FileBinaryReader(filename)

# ======================================================================

def open_for_reading_text(filename):
    """Opens text file for reading. Handles compressed files transparently."""
    return FileTextReader(filename)

# ======================================================================

def read_binary(filename):
    """Reads and returns entire binary file content as bytes. Handles
    compressed files transparently."""
    with FileBinaryReader(filename) as f:
        return f.read()

# ======================================================================

def read_or_get_binary(data, try_reading_from_file=True):
    """If data is a filename, reads it. Uncompresses data, if necessary."""
    try:
        if try_reading_from_file and isinstance(data, str) and len(data) < 1024:
            data = read_binary(data)
    except Exception as err:
        raise RuntimeError('Unable to open file {!r}'.format(data))
    # check if data is compressed
    try:
        data = lzma.decompress(data)
    except Exception as err:
        try:
            data = bz2.decompress(data)
        except Exception as err:
            pass
    # do not convert to str!
    return data

# ======================================================================

def write_json(filename, data, indent=None, sort_keys=False, backup=True):
    if indent is None:
        s = json.dumps(data, separators=[',', ':'], indent=indent, sort_keys=sort_keys)
    else:
        s = json_dumps(data, indent=indent, indent_increment=indent)
    with open_for_writing_binary(filename, backup=backup) as fd:
        fd.write(s.encode('utf-8'))

# ======================================================================

def read_json(filename):
    return json.loads(open_for_reading_text(filename).read())

# ======================================================================

def json_dumps(data, indent=2, indent_increment=2):
    """More compact dumper with wide lines."""

    def simple(d):
        r = True
        if isinstance(d, dict):
            r = not any(isinstance(v, (list, tuple, set, dict)) for v in d.values())
        elif isinstance(d, (tuple, list)):
            r = not any(isinstance(v, (list, tuple, set, dict)) for v in d)
        return r

    def end(symbol, indent):
        if indent > indent_increment:
            r = "{:{}s}{}".format("", indent - indent_increment, symbol)
        else:
            r = symbol
        return r

    r = []
    if simple(data):
        if isinstance(data, set):
            r.append(json.dumps(sorted(data), sort_keys=True))
        else:
            r.append(json.dumps(data, sort_keys=True))
    else:
        if isinstance(data, dict):
            r.append("{")
            for no, k in enumerate(sorted(data), start=1):
                comma = "," if no < len(data) else ""
                r.append("{:{}s}{}: {}{}".format("", indent, json.dumps(k), json_dumps(data[k], indent + indent_increment, indent_increment), comma))
            r.append(end("}", indent))
        elif isinstance(data, (tuple, list)):
            r.append("[")
            for no, v in enumerate(data, start=1):
                comma = "," if no < len(data) else ""
                r.append("{:{}s}{}{}".format("", indent, json_dumps(v, indent + indent_increment, indent_increment), comma))
            r.append(end("]", indent))
    return "\n".join(r)

# ----------------------------------------------------------------------

def open_for_writing_binary(filename, compressed=None, backup=True, makedirs=True):
    """Opens binary file for writing. If compressed is None, autodetects if data should be compressed by filename suffix."""
    if filename == '-' or filename is None:
        f = sys.stdout.buffer
    else:
        if compressed is None:
            if filename[-4:] == '.bz2':
                compressed = 'bz2'
            elif filename[-3:] == '.xz':
                compressed = 'xz'
        elif compressed is True:
            compressed = 'xz'
        if backup:
            backup_file(filename, backup)
        if makedirs and '/' in filename:
            try:
                os.makedirs(os.path.dirname(filename))
            except:
                pass
        if compressed == 'bz2':
            f = bz2.BZ2File(filename, mode='w')
        elif compressed == 'xz':
            f = lzma.open(filename, mode='wb', preset=9 | lzma.PRESET_EXTREME)
        else:
            f = open(filename, mode='wb')
    return f

# ======================================================================

def write_binary(filename, data, compressed=None, backup=True, makedirs=True):
    """Writes data (bytes) into a binary file. If compressed is None,
    autodetects if data should be compressed by filename suffix."""
    with open_for_writing_binary(filename, compressed=compressed, backup=backup, makedirs=makedirs) as f:
        f.write(data)

# ======================================================================

def backup_file(filename, backup_dir=None):
    """Backup the file, if it exists. Backups versioning is supported."""
    if isinstance(backup_dir, str):
        newname = os.path.join(backup_dir, os.path.basename(filename))
    else:
        newname = filename
    version = 1
    while os.access(newname, os.F_OK):
        newname = '{}.~{:02d}~'.format(filename, version)
        version += 1
    if newname != filename:
        #module_logger.debug("Backing up file: {} --> {}".format(filename, newname))
        try:
            os.rename(filename, newname)
        except Exception as err:
            module_logger.warning('Cannot create backup copy of {}: {}'.format(filename, err), exc_info=True)

# ======================================================================

class FileBinaryReader:

    def __init__(self, filename):
        self.filename = filename
        self.name = filename
        self.open_xz()

    def read(self, size=-1):
        try:
            return self.f.read(size)
        except (IOError, EOFError, lzma.LZMAError):
            self.open()
            return self.read(size)

    def readline(self):
        try:
            return self.f.readline()
        except (IOError, EOFError, lzma.LZMAError):
            self.open()
            return self.readline()

    def open_plain(self):
        self.f = open(self.filename, mode='rb')
        self.open = self.open_fail

    def open_bz2(self):
        self.f = bz2.BZ2File(self.filename)
        self.open = self.open_plain

    def open_xz(self):
        self.f = lzma.LZMAFile(self.filename)
        self.open = self.open_bz2

    def open_fail(self):
        raise IOError('Unable to open/read ' + repr(self.filename))

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        try:
            self.f.close()
        except:
            pass

    def __iter__(self):
        return self

    def __next__(self):
        s = self.readline()
        if not len(s):
            raise StopIteration()
        return s

# ======================================================================

class FileTextReader (FileBinaryReader):

    def read(self, size=-1):
        s = super().read(size=size)
        if isinstance(s, bytes):
            s = s.decode('utf-8')
        return s

    def readline(self):
        s = super().readline()
        if isinstance(s, bytes):
            s = s.decode('utf-8')
        return s

# ======================================================================
### Local Variables:
### eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
### End:
