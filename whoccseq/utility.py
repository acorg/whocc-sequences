# -*- Python -*-
# license
# license.
# ----------------------------------------------------------------------

import datetime
import logging; module_logger = logging.getLogger(__name__)
from contextlib import contextmanager

# ----------------------------------------------------------------------

@contextmanager
def timeit(name):
    start = datetime.datetime.utcnow()
    yield
    module_logger.info('{} <{}>'.format(name, datetime.datetime.utcnow() - start))

# ----------------------------------------------------------------------

def fix_virus_type_lineage(vt, lineage=None):
    if vt:
        vt = vt.upper()
        if vt == "H1":
            vt = "A(H1N1)"
        elif vt == "H3":
            vt = "A(H3N2)"
        if vt not in ["B", "A(H1N1)", "A(H3N2)"]:
            raise ValueError("Unrecognized virus type: {}".format(vt))
        if lineage:
            lineage = lineage.upper()
            if vt == "B":
                if lineage not in ["YAMAGATA", "VICTORIA"]:
                    raise ValueError("Unrecognized {} lineage: {}".format(vt, lineage))
            elif vt == "A(H1N1)":
                if lineage not in ["2009PDM"]:
                    raise ValueError("Unrecognized {} lineage: {}".format(vt, lineage))
            else:
                raise ValueError("{} does not have lineages".format(vt))
    elif lineage:
        raise ValueError("You must provide virus type to filter by lineage")
    return vt, lineage

# ======================================================================
### Local Variables:
### eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
### End:
