# -*- Python -*-
# license
# license.

import sys, os, re, collections, subprocess, time, datetime, pprint
from pathlib import Path
import logging; module_logger = logging.getLogger(__name__)

# ----------------------------------------------------------------------

class JobError (Exception): pass

class Job:

    sJobStatus = {'1': 'pending', '2': 'running', '3': 'removed', '4': 'completed', '5': 'held', '6': 'transferring output', '7': 'suspended'}

    def __init__(self, clusters :dict):
        """clusters - cluster_id to number of jobs mapping"""
        self.clusters = clusters

    def __str__(self):
        return str(self.clusters)

    def status(self):
        # in reality needs to look into history file to check for failed jobs but it is slow
        still_running = collections.defaultdict(int)
        errors = []
        output = _run("condor_q", *list(self.clusters), "-autoformat:tl", "ClusterId", "ProcId", "JobStatus", "HoldReason")
        for line in output.splitlines():
            # module_logger.info('Status line {}'.format(line))
            state = dict(f.split(' = ') for f in line.split('\t'))
            state["JobStatus"] = self.sJobStatus[state["JobStatus"]]
            if state["JobStatus"] not in ["pending", "running", "transferring output", "completed", "suspended"]:
                errors.append(state)
            else:
                still_running[state["ClusterId"]] += 1
        if errors:
            raise JobError(errors)
        r = {}
        for cluster, jobs in self.clusters.items():
            cluster_running = still_running.get(cluster, 0)
            if cluster_running:
                r[cluster] = {"state": "running", "done": (jobs - cluster_running) / jobs}
            else:
                r[cluster] = {"state": "completed"}
        all_jobs = sum(self.clusters.values())
        all_running = sum(still_running.values())
        r["all"] = {"state": "running" if all_running else "completed", "done": (all_jobs - all_running) / all_jobs}
        return r

    def wait(self, check_interval_in_seconds=30, verbose=True):
        start = datetime.datetime.now()
        while True:
            time.sleep(check_interval_in_seconds)
            status = self.status()
            if status["all"]["state"] == "completed":
                break
            if verbose:
                module_logger.info('Percent done: {:.1f}%'.format(status["all"]["done"] * 100.0))
        if verbose:
            module_logger.info('All done in {}'.format(datetime.datetime.now() - start))

# ----------------------------------------------------------------------

sReCondorProc = re.compile(r'\*\*\s+Proc\s+(\d+)\.(\d+):')

def submit(program, program_args :list, description :str, current_dir :str, capture_stdout=False, email="condor@skepner.eu", notification="Error", machines :list = None):
    current_dir = str(Path(current_dir).resolve())
    Path(current_dir).chmod(0o777)        # to allow remote processes runinnig under user nobody to write files
    desc = [
        ["universe", "vanilla"],
        ["executable", str(Path(program).resolve())],
        ["should_transfer_files", "NO"],
        ["notify_user", email],
        ["notification", notification],
        ["Requirements", "({})".format(" || ".join('machine == "{}"'.format(fix_machine_name(m)) for m in machines)) if machines else None],
        ["initialdir", current_dir],
        ["description", "{} {}".format(description, current_dir)],
        [""],
        ]
    stderr_files = [Path(current_dir, "{:04d}.stderr".format(no)) for no, args in enumerate(program_args, start=1)]
    for f in stderr_files:
        f.touch()
        f.chmod(0o777)
    if capture_stdout:
        stdout_files = [Path(current_dir, "{:04d}.stdout".format(no)) for no, args in enumerate(program_args, start=1)]
        for f in stdout_files:
            f.touch()
            f.chmod(0o777)
    else:
        stdout_files = [None] * len(program_args)
    for no, args in enumerate(program_args):
        desc.extend([
            ["arguments", " ".join(str(a) for a in args)],
            ["error", stderr_files[no] and str(stderr_files[no])],
            ["output", stdout_files[no] and str(stdout_files[no])],
            ["queue"],
            [""],
            ])
    desc_s = "\n".join((" = ".join(e) if len(e) == 2 else e[0]) for e in desc if len(e) != 2 or e[1])
    desc_filename = str(Path(current_dir, "garli.desc"))
    with open(desc_filename, "w") as f:
        f.write(desc_s)
    output = _run("condor_submit", "-verbose", desc_filename)
    cluster = collections.defaultdict(int)
    for line in output.splitlines():
        m = sReCondorProc.match(line)
        if m:
            cluster[m.group(1)] += 1
    if not cluster:
        logging.error(output)
        raise RuntimeError("cluster id not found in the submission results {}".format(cluster))
    return Job(dict(cluster))

# ----------------------------------------------------------------------

sMachineName = {"albertine": "i19", "odette": "odette.antigenic-cartography.org", "i18": "odette.antigenic-cartography.org"}
def fix_machine_name(machine):
    return sMachineName.get(machine, machine)

# ----------------------------------------------------------------------

def _run(*program_and_args):
    return subprocess.check_output(program_and_args, env={"LD_LIBRARY_PATH": ""}).decode("utf-8")

# ======================================================================
### Local Variables:
### eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
### End:
