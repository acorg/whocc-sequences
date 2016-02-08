# -*- Python -*-
# license
# license.
# ----------------------------------------------------------------------

import datetime
import logging; module_logger = logging.getLogger(__name__)
from contextlib import contextmanager

@contextmanager
def timeit(name):
    start = datetime.datetime.utcnow()
    yield
    module_logger.info('{} <{}>'.format(name, datetime.datetime.utcnow() - start))

# ======================================================================
### Local Variables:
### eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
### End:
