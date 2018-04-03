#!/usr/bin/env python
import subprocess
import multiprocessing
import logging
import argparse
import re
import time

logging.basicConfig(level=logging.INFO)

def run(cmd):
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    out = p.communicate()[0]

    result = out.strip()
    return p.returncode, result

def run_slurm(array_id, script, args):
    cmd = ["sbatch", "--exclusive", "-a", str(array_id), script] + list(args)
    rc, result = run(cmd)
    if rc == 0:
        match = re.search(r'Submitted\s+batch\s+job\s+(\d+)', result.strip())
        if match:
            logging.info("Job %s submitted (%s)", match.group(1), cmd)
            return match.group(1) 
        else:
            logging.error("Couldn't parse %s", result.strip())
    else:
        logging.error("Error running %s", cmd)
    return None

def slurm_job_done(array_id, job_id):
    cmd = ["squeue", "-r", "-j", "%s_%s" % (job_id, array_id), "-h", "-o", "%A"]
    rc, result = run(cmd)
    if rc == 0:
        done = not result.strip()
        logging.info("job %s_%s, done = %s (%s)", job_id, array_id, done, result.strip())
        return done
    else:
        logging.error("Error running %s", cmd)
    return False

def run_job(arguments):
    array_id, script, args = arguments
    logging.info("Running %s %s", array_id, script)
    job_id = run_slurm(array_id, script, args)
    if job_id is not None:
        time.sleep(10)
        done = slurm_job_done(array_id, job_id)
        while not done:
            time.sleep(5)
            done = slurm_job_done(array_id, job_id)
        logging.info("%s - done", arguments)
    else:
        logging.error("Didn't get job id: %s %s", array_id, script)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--start', metavar='S', type=int, help='Starting array id', required=True)
    parser.add_argument('-e', '--end', metavar='E', type=int, help='Ending array id', required=True)
    parser.add_argument('-n', '--parallel', metavar='N', type=int, help='Number of jobs to run in parallel', required=True)
    parser.add_argument('--script', type=argparse.FileType('r'), help='The script to execute', required=True)
    parser.add_argument('extra_args', nargs='*', default=[], help='Any additional arguments')

    args = parser.parse_args()

    arguments = [ (array_id, args.script.name, args.extra_args) for array_id in range(args.start, args.end+1) ]
    pool = multiprocessing.Pool(processes=min(args.parallel, len(arguments)))

    results = pool.map(run_job, arguments)

    pool.close()
    pool.join()
