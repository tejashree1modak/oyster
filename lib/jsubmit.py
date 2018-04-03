#!/usr/bin/env python
import multiprocessing
import logging
import argparse
import os
import time
import random
logging.basicConfig(level=logging.INFO)

def slack_message(msg):
    logging.info("%s", msg)
    #os.system(r"source %s/slack.sh;post_slack_message cluster-jobs "'%s'" jsubmit" % 
            #(os.path.dirname(__file__), msg))

def run_job(array_id, script):
    cmd = "sbatch -W -x n042 --exclusive -a %s %s" % (array_id, script)
    slack_message("Running '%s'" % cmd)
    rc = os.system(cmd)
    if rc != 0: # if cmd runs successfully, rc is 0
        slack_message("Error running %s" % cmd)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', metavar='1,2,3-5,7,19-30', help='The array IDs', required=True)
    parser.add_argument('-n', '--parallel', metavar='N', type=int, help='Number of jobs to run in parallel', required=True)
    parser.add_argument('--script', type=argparse.FileType('r'), help='The script to execute', required=True)
    args = parser.parse_args()

    # ids is a list of arrayids that we will generate based on -a specified
    # on the commandline. The value is comma sperated ids like 1,2,3-5,10 etc
    # we split on , and then if there is a '-' then split on it to get start,end
    # else just the number that is specified and append into ids
    # If we split to get a start and end, then we use range(start,end+1) to create
    # an array - e.g. 3-5 = range(3,6) -> [3,4,5] which is then added to ids list
    # finally we sort the ids list in ascending order.
    ids = []
    for k in args.a.split(','):
        if '-' in k:
            start, end = map(int, k.split('-'))
            ids += range(start, end+1)
        else:
            ids.append(int(k))
    ids.sort()

    # generate a list of tuples e.g. -a 1,3,4-6 will result in =
    # [(1,<script>),(3,<Script>),(4,<script>),(5,<script>),(6,<Script>)]
    arguments = [ (array_id, args.script.name) for array_id in ids ]

    slots = {}  # a dictionary of <arrayid> -> <running process that is running run_job() function>
    nprocesses = args.parallel

    while arguments or slots:
        # if we have the len of slots dict is < nprocesses, then we can stand to add more
        # array ids and start processes until we have nprocesses running
        while len(slots) < nprocesses and arguments:
            array_id, script = arguments.pop(0) # pop removes it from arguments, so eventually it will be empty and 
            slots[array_id] = multiprocessing.Process(target=run_job, args=(array_id,script))
            slots[array_id].start()

        time.sleep(10)                 

        # look a each process in slots dictionary. If it is finished, then
        # delete it from dict. Then when we go around in the loop, we will start
        # another one in its lieu .. and so on until we exhaust arguments list
        for array_id, process in slots.items():
            if not process.is_alive():
                logging.info("%s is done with %s", array_id, process.exitcode)
                del slots[array_id]
