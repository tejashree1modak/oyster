import xml.etree.ElementTree as ET
import argparse
import collections
import os

def read_fa(fa_file):
    fa = collections.OrderedDict()
    last = None
    keys = []

    for line in fa_file:
        line = line.strip()
        if line:
            if line.startswith('>'):
                last = line[1:]
                fa[last] = 0
                keys.append(last)
            else:
                fa[last] += len(line)
    return keys, fa

def get_info_from_xml(xml_file):
    tree = ET.parse(xml_file)
    df = tree.getroot().find("BlastOutput_query-def")
    if df is not None:
        anno = df.text.strip()
        qlen = tree.getroot().find("BlastOutput_query-len")
        return anno, int(qlen.text.strip()) if qlen is not None else None

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--xml-dir', help='The directory with all xmls', required=True)
    parser.add_argument('--xml-prefix', help='The prefix of each xml file', required=True)
    parser.add_argument('--fa-file', type=argparse.FileType('r'), help='The input fa file', required=True)

    args = parser.parse_args()

    k, fa = read_fa(args.fa_file)
    #print '\n'.join("%03d %s" % (i, v) for i, v in enumerate(k))

    ok = True
    
    for ix, anno in enumerate(k):
        ix += 1
        xml_file = os.path.join(args.xml_dir, '%s.xml.%s' % (args.xml_prefix, ix))
        xanno, xlen = None, None
        if os.path.isfile(xml_file):
            xanno, xlen = get_info_from_xml(xml_file)

        if xanno is None:
            print "Step %s not run" % ix
            ok = False
        else:
            if xanno != anno:# or xlen != fa[anno]:
                print "Step %s, expected(anno=%s,len=%s), actual(anno=%s,len=%s)" % (
                    ix, anno, fa[anno], xanno, xlen)
                ok = False

    if ok:
        print "OK"
