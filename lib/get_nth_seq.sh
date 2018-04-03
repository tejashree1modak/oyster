#!/bin/bash
PROG=${0##*/}
LIB_DIR=${0%/*}

usage() {
    echo -e "$PROG: [options] [fa file]\n"
    echo "-n <N>            Extract Nth sequence"
    echo "-o <out-file>     Write the sequence to this out file"
    echo "-s <seq/id name>  Search this sequence/id '>..' in fa file and then extract that seq"
    exit
}

while getopts "n:o:s:h" arg; do
    case "$arg" in
        n)  N="$OPTARG" ;;
        o)  out_file="$OPTARG" ;;
        s)  seq="$OPTARG" ;;
        h)  usage ;;
        *)  usage ;;
    esac
done

shift $((OPTIND-1))

if [ -n "$N" -o -n "$seq" ] || [ -n "$out_file" -a -f "$1" ]; then
    query_file="$1"

    if [ -z "$N" ]; then
        if ! N=$(grep '^>' $query_file | awk -v s="$seq" '$1 == ">" s { print NR; exit }' ) ; then
            echo "ERROR: cannot find '$seq'"
            exit 1
        else
            echo "'$seq' is ${N}th sequence"
        fi
    fi

    # if this is an array job, then we will take Nth sequence from the query_file and
    # write it to a temporary file as the input file
    # the output file will ahve array id as its suffix and will have to be merged later on by hand
    #out_file="/tmp/${query_file##*/}.${N}.$$"

    offset=( $(grep -b '^>' $query_file | awk -v FS=':' -v n=$N -f $LIB_DIR/nfasta.awk) )

    if [ ${#offset[*]} -eq 1 ] ; then
        # last offset
        echo dd if=$query_file of=$out_file bs=1 skip=${offset[0]}
        dd if=$query_file of=$out_file bs=1 skip=${offset[0]}
    elif [ ${#offset[*]} -eq 2 ]; then
        echo dd if=$query_file of=$out_file bs=1 skip=${offset[0]} count=${offset[1]}
        dd if=$query_file of=$out_file bs=1 skip=${offset[0]} count=${offset[1]}
    else
        echo "ERROR: Couldnt get $n offset from $infile (${offset[*]})"
        exit 1
    fi
else
    echo "ERROR: (-n or -s), -o and input file required"
    usage
fi
