#!/usr/bin/awk -f 
BEGIN {
    FS = OFS = ","
    if (ARGC != 3 || top < 0) {
        printf "ERROR: topreg.awk -v top=<topN> -v header=<0|1> <reg txt file> <top n csv file>\n"
        exit 1
    }
}

FILENAME == ARGV[1] && NR > header {
    if (top > 0 && NR > header + top) {
        nextfile
    }
    reg[$1] = $2
}

FILENAME == ARGV[2]  {
    split($2, id, /_/)
    if (id[2] in reg)  {
        print id[2], reg[id[2]], $6
    }
}

