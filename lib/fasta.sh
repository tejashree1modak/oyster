# --------------------------------------------------------------------------------
# this file has some shell routines that help common tasks with fasta files
# --------------------------------------------------------------------------------
#

project_home=/data3/marine_diseases_lab/tejashree/Bio_project_SRA/
git_lib=$project_home/github/oyster/lib

# given a fasta file, N and name of a tmp file, extract the Nth sequence
# and its annotation from the fasta file and write it to tmp file
function get_nth_sequence()  {
    local infile=$1 n=$2 outfile=$3
    echo $*
    local offset=( $(grep -b '^>' $infile | awk -v n=$n -f $git_lib/nfasta.awk) )

    if [ ${#offset[*]} -eq 1 ] ; then
        # last offset
        echo dd if=$infile of=$outfile bs=1 skip=${offset[0]} status=noxfer
        dd if=$infile of=$outfile bs=1 skip=${offset[0]} status=noxfer
    elif [ ${#offset[*]} -eq 2 ]; then
        echo dd if=$infile of=$outfile bs=1 skip=${offset[0]} count=${offset[1]} status=noxfer
        dd if=$infile of=$outfile bs=1 skip=${offset[0]} count=${offset[1]} status=noxfer
    else
        echo "ERROR: Couldnt get $n offset from $infile (${offset[*]})"
        return 1
    fi
}
