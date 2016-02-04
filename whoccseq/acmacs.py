# -*- Python -*-
# license
# license.

"""
API to access AcmacsWeb server
"""

import os, urllib.request, json, traceback
import logging; module_logger = logging.getLogger(__name__)

# ----------------------------------------------------------------------

sApi = None

def api(url_prefix=None):
    global sApi
    if sApi is None:
        sApi = API(url_prefix)
    if url_prefix is not None:
        sApi.url_prefix = url_prefix
    return sApi

# ----------------------------------------------------------------------

def normalize_names(names):
    source = list(names)
    response = api().execute({"C": "name_normalize", "names": source})
    mapping = dict(zip(source, response["names"]))
    #module_logger.info('Names:\n{}'.format(json.dumps(mapping, indent=2, sort_keys=True)))
    return mapping

# ----------------------------------------------------------------------

def normalize_passages(passages):
    source = list(passages)
    response = api().execute({"C": "passage_normalize", "passages": source})
    mapping = dict(zip(source, response["passages"]))
    # module_logger.info('Passages:\n{}'.format(json.dumps(mapping, indent=2, sort_keys=True)))
    return mapping

# ----------------------------------------------------------------------

def normalize_dates(dates):
    source = list(dates)
    response = api().execute({"C": "date_normalize", "dates": source})
    mapping = dict(zip(source, response["dates"]))
    # module_logger.info('Dates:\n{}'.format(json.dumps(mapping, indent=2, sort_keys=True)))
    return mapping

# ----------------------------------------------------------------------

class CommandError (Exception):
    """Raised by api._execute if command resposne contains error and raise_error flag is set."""

# ======================================================================

class API:

    def __init__(self, url_prefix):
        self.url_prefix = url_prefix

    def execute(self, command):
        if self.url_prefix:
            response = self._execute_http(command)
        else:
            raise ValueError('No url_prefix')
        if isinstance(response, dict) and response.get('E'):
            raise CommandError(response['E'])
        return response

    def _execute_http(self, command):
        command['F'] = 'json'
        if "localhost" in self.url_prefix:
            import ssl
            context = ssl.create_default_context()
            context.check_hostname = False
            context.verify_mode = ssl.CERT_NONE
        else:
            context = None
        # print(command)
        response = urllib.request.urlopen(url='{}/api'.format(self.url_prefix), data=json.dumps(command).encode('utf-8'), context=context).read()
        return json.loads(response.decode('utf-8'))

# ======================================================================
### Local Variables:
### eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
### End:
