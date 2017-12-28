#!/usr/bin/awk -f

# and it expects input to be of the form <offset>:<line>,
# basically output of grep -b ^> fasta

BEGIN {
    FS = ":"
    start = -1
    end = -1
}

NR == n { 
    start = $1
}

start >= 0 && NR == n+1 {
    end = $1
    exit
} 

END {
    if (start >= 0 && end > 0)  {
        print start, end-start
    } else if (start >= 0) {
        print start
    }
}
