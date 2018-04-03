import xml.etree.ElementTree as ET
import argparse
import sys
import csv

HSP_KEYS = ('num', 'bit-score', 'score', 'evalue', 'query-from',
            'query-to', 'hit-from', 'hit-to', 'query-frame', 'hit-frame',
            'identity', 'positive', 'gaps', 'align-len')

HIT_KEYS = ('num', 'id', 'def', 'accession', 'len')

ITER_KEYS = ('iter-num', 'query-def', 'query-len')

def construct_empty_result(prefix, keys):
    result = { (prefix + '_' + key) : "" for key in keys }
    return result

def get_child(parent, child_name):
    child = parent.find("./" + child_name)
    if child is not None:
        return child.text.strip()
    return ""

def get_hsp_info(hsp):
    result = construct_empty_result("Hsp", HSP_KEYS)
    if hsp is not None:
        for key in result.keys():
            result[key] = get_child(hsp, key)
    return result

def get_hit_info(it_hit):
    result = construct_empty_result("Hit", HIT_KEYS)
    hsp = it_hit.find("./Hit_hsps/Hsp") if it_hit is not None else None
    if it_hit is not None:
        for key in result.keys():
            result[key] = get_child(it_hit, key)

    result.update(get_hsp_info(hsp))
    return result

def get_top_n(iteration, n):
    result_rows = []
    base = construct_empty_result("Iteration", ITER_KEYS)
    for key in base.keys():
        base[key] = get_child(iteration, key)

    hits = iteration.findall("./Iteration_hits/Hit")[:n] if iteration is not None else []
    if hits:
        for hit in hits:
            result = base.copy()
            result.update(get_hit_info(hit))
            result_rows.append(result)
    else:
        base.update(get_hit_info(None))
        result_rows.append(base)
    
    return result_rows

def get_columns():
    return [ ('Iteration_' + k) for k in ITER_KEYS ] + [
             ('Hit_' + k) for k in HIT_KEYS] + [ 
             ('Hsp_' + k) for k in HSP_KEYS]

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--xml', type=argparse.FileType('r'), default=sys.stdin, help='The combined xml file')
    parser.add_argument('-o', '--out', type=argparse.FileType('w'), default=sys.stdout, help='The out csv file')
    parser.add_argument('-n', '--top', type=int, default=1, help='The top N hits to include (default:1)')

    args = parser.parse_args()
    if args.top <= 0:
        print "ERROR: --top needs to be >= 1"
    else:
        document = ET.parse(args.xml)
        rows = []
        for iteration in document.getroot().findall("./BlastOutput_iterations/Iteration"):
            rows.extend(get_top_n(iteration, args.top))

        writer = csv.DictWriter(args.out, get_columns(), dialect=csv.excel)
        writer.writeheader()
        if rows:
            writer.writerows(rows)
