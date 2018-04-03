import re
import sys
import collections
import csv
import argparse

def get_transcript_id(col):
    for opts in col.split(";"):
        opts = opts.strip()
        if opts.startswith("transcript_id"):
            return opts.split()[-1].strip('" ')
    return None

def do_search(needle, haystack, write_stats):
    table = { id.strip() for id in needle }
    stats = collections.Counter()
    for line in csv.reader(haystack, dialect=csv.excel_tab):
        if len(line) > 2 and line[2] == "transcript":
            tid = get_transcript_id(line[-1])
            if tid is not None and tid in table:
                print "\t".join(line)
                stats[tid] += 1

    if write_stats:
        write_stats = csv.writer(write_stats, dialect=csv.excel)
        write_stats.writerows(stats.most_common())

def do_uniq(search_results):
    stats = collections.Counter()
    for line in csv.reader(search_results, dialect=csv.excel_tab):
        tid = get_transcript_id(line[-1])
        if tid is not None:
            stats[tid] += 1

    write_stats = csv.writer(sys.stdout, dialect=csv.excel)
    write_stats.writerows(stats.most_common())

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    sub = parser.add_subparsers(dest="command")

    search = sub.add_parser("search", description="Search for transcript ids in gtf")
    search.add_argument('--needle', type=argparse.FileType('r'), help='The file with ids to search for, one per line', required=True)
    search.add_argument('--haystack', type=argparse.FileType('r'), default=sys.stdin, help='The gtf file to search "needles" in')
    search.add_argument('--write-stats', type=argparse.FileType('w'), help='Write statistics of the matches')

    uniq = sub.add_parser("unique-matches", description="Get all unique matches from the output of search")
    uniq.add_argument("--search-results", type=argparse.FileType('r'), default=sys.stdin, help='result of search operation')

    args = parser.parse_args()

    if args.command == 'search':
        do_search(args.needle, args.haystack, args.write_stats)
    elif args.command == 'unique-matches':
        do_uniq(args.search_results)
