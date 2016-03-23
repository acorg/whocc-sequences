# -*- Python -*-
# license
# license.

import sys, os, re
from pathlib import Path
import logging; module_logger = logging.getLogger(__name__)

# ----------------------------------------------------------------------

class Job:

    def __init__(self, clusters :dict):
        """clusters - cluster_id to number of jobs mapping"""
        self.clusters = clusters

# ----------------------------------------------------------------------

sReCondorProc = re.compile(r'\*\*\s+Proc\s+(\d+)\.(\d+):')

def submit(program, program_args :list, description :str, current_dir :str, capture_stdout=False, email="condor@skepner.eu", notification="Error", machines :list = None):
    current_dir = str(Path(current_dir).resolve())
    desc = [
        ["universe", "vanilla"],
        ["executable", str(Path(program).resolve())],
        ["should_transfer_files", "NO"],
        ["notify_user", email],
        ["notification", notification],
        ["Requirements", "({})".format(" || ".join('machine == "{}"'.format(m) for m in machines)) if machines else None],
        ["initialdir", current_dir],
        ["description", "{} {}".format(description, current_dir)],
        [""],
        ]
    for no, args in enumerate(program_args, start=1):
        desc.extend([
            ["arguments", args],
            ["error", "{}/{:03d}.stderr".format(current_dir, no)],
            ["output", "{}/{:03d}.stdout".format(current_dir, no) if capture_stdout else None],
            ["queue"],
            [""],
            ])
    desc_s = "\n".join(" = ".join(e) if len(e) == 2 else e[0] for e in desc if len(e) == 2 and e[1])
    desc_filename = str(Path(current_dir, "garli.desc"))
    with open(desc_filename, "w") as f:
        f.write("\n".join(desc_s))
    output = _run("condor_submit", "-verbose", desc_filename)
    cluster = collections.defaultdict(int)
    for line in output.splitlines():
        m = sReCondorProc.match(line)
        if m:
            cluster[m.group(1)] += 1
    if not cluster:
        logging.error(output)
        raise RuntimeError("cluster id not found in the submission results {}".format(cluster))
    logging.info("Cluster: {}".format(dict(cluster)))
    return Job(cluster)

# ----------------------------------------------------------------------

def _run(*program_and_args):
    return subprocess.check_output(program_and_args), env={"LD_LIBRARY_PATH": ""}).decode("utf-8")

# ======================================================================
### Local Variables:
### eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
### End:
