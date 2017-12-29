import xml.etree.ElementTree as ET
import argparse
import re
import sys

# assuming output files end in '-N' (e.g. blast.xml.1, blast.xml.10) get the suffix to
# sort the files in right order
suffix_re = re.compile('.*\.xml\.(\d+)$')
def get_sort_key(file_obj):
    m = suffix_re.match(file_obj.name)
    if m:
        return int(m.group(1))
    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Concatenate multiple blast out xml files')
    parser.add_argument('-o', '--out', type=argparse.FileType('w'), default=sys.stdout, help='Output file')
    parser.add_argument('files', nargs='+', type=argparse.FileType('r'), help='Input xml files to concatenate')

    args = parser.parse_args()

    first = bi = None
    for fl in sorted(args.files, key=get_sort_key):
        tree = ET.parse(fl)
        if first is None:
            first = tree
            bi = first.getroot().find("./BlastOutput_iterations")
        else:
            # add <BlastOutput_iterations><Iteration> section to first"
            n = len(bi)
            for it in tree.getroot().findall("./BlastOutput_iterations/Iteration"):
                # update the number in <Iteration_iter-num> tag
                iter_num = it.find("./Iteration_iter-num")
                if iter_num is not None:
                    iter_num.text = str(int(iter_num.text.strip()) + n)

                # update the number in <Iteration_iter-num> tag
                qid = it.find("./Iteration_query-ID")
                if qid is not None:
                    qnum = qid.text.strip()
                    qnum = qnum[len("Query_"):] # take the suffix number, like 1 from Query_1, 2 from Query_2
                    qid.text = "Query_{0}".format(int(qnum) + n)

                bi.append(it)   # add it to the first 

    first.write(args.out)
