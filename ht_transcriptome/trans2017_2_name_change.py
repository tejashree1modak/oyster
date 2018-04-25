import argparse
import os
import sys
import re
import glob

def change_name(f):
    '''
    f is a name like C-K-0_S60_L007_R1_001.fastq
    we will return a new name like C-K-0_L007_1.fastq
    '''
    name = os.path.basename(f)
    dir = os.path.dirname(f)
    name = re.sub(r'(.*)_S\d+_(L\d+)_R(\d)_\d+.fastq', r'\1_\2_\3.fastq', name)
    print "mv %s %s" % (f, os.path.join(dir, name))
    os.rename(f, os.path.join(dir, name))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=["trans2017_2"], required=True, help="What kind of name change do we want")
    parser.add_argument("dir", nargs=1, help="The directory with all the files")

    args = parser.parse_args()

    if not os.path.isdir(args.dir[0]):
        print "ERROR: '%s' not a directory" % args.dir[0]
        sys.exit(1)

    for f in glob.glob(args.dir[0] + '/*.fastq'):
        change_name(f)
