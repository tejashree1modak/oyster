#!/bin/bash
#SBATCH --job-name="blast"
#SBATCH --time=9999:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --output="des4_5d_trans_blast_out.%A-%a"
#SBATCH --error="des4_5d_trans_blast_out.%A-%a"

set -e 
echo "START $(date)"

ncbi_nr_db=/data3/shared/ncbi-nr/nr 
project_home=/data3/marine_diseases_lab/tejashree/Bio_project_SRA/
blast_dir=$project_home/blast 
stringtie=$project_home/stringtie

source $project_home/github/oyster/lib/slack.sh

cd $blast_dir
#change query:
#query_file=$stringtie/test/rightway_unchar_na.fa
#query_file=$project_home/github/oyster/ht_transcriptome/in.fa
query_file=$stringtie/test/design4/des4_5d_trans/des4_5d_trans_unique.fa
#query_file=$stringtie/test/design1/trans/des1_trans_unique.fa
#out_file=test.xml
#CHange:
out_file=des4_5d_trans.xml

num_threads=20

echo $HOSTNAME
module load BLAST+/2.5.0-foss-2016b-Python-2.7.12

if [ -z "$radix" ]; then
    echo "No radix"
    exit 1
fi

SLURM_ARRAY_TASK_ID=$(( ($radix * 1000) + $SLURM_ARRAY_TASK_ID ))
echo "SLURM_ARRAY_TASK_ID = $SLURM_ARRAY_TASK_ID"

# tell system to call post_slack_message with last 50 lines of output
trap "tail -n 50 blast_out-${SLURM_ARRAY_TASK_ID} | post_slack_message cluster-jobs - '$me::log'" ERR

if [ -n "$SLURM_ARRAY_TASK_ID" ]; then
    # if this is an array job, then we will take Nth sequence from the query_file and
    # write it to a temporary file as the input file
    # the output file will ahve array id as its suffix and will have to be merged later on by hand
    tmp_file="/tmp/${query_file##*/}.${SLURM_ARRAY_TASK_ID}.$$"

    offset=( $(grep -b '^>' $query_file | awk -v FS=':' -v n=$SLURM_ARRAY_TASK_ID -f $project_home/github/oyster/lib/nfasta.awk) )

    if [ ${#offset[*]} -eq 1 ] ; then
        # last offset
        echo dd if=$query_file of=$tmp_file bs=1 skip=${offset[0]} status=noxfer
        dd if=$query_file of=$tmp_file bs=1 skip=${offset[0]} status=noxfer
    elif [ ${#offset[*]} -eq 2 ]; then
        echo dd if=$query_file of=$tmp_file bs=1 skip=${offset[0]} count=${offset[1]} status=noxfer
        dd if=$query_file of=$tmp_file bs=1 skip=${offset[0]} count=${offset[1]} status=noxfer
    else
        echo "ERROR: Couldnt get $n offset from $infile (${offset[*]})"
        exit 1
    fi

    query_file=$tmp_file 
    out_file=${out_file}.${SLURM_ARRAY_TASK_ID}
fi

me="blast::${PBS_JOBNAME}::$SLURM_ARRAY_TASK_ID on $HOSTNAME"

echo blastx -db $ncbi_nr_db -outfmt 5 -evalue 1e-3 -word_size 3 -show_gis -num_alignments 20 -max_hsps 20 -num_threads $num_threads -out $out_file -query $query_file
post_slack_message cluster-jobs "$me : blastx -db $ncbi_nr_db -outfmt 5 -evalue 1e-3 -word_size 3 -show_gis -num_alignments 20 -max_hsps 20 -num_threads $num_threads -out $out_file -query $query_file"
blastx -db $ncbi_nr_db -outfmt 5 -evalue 1e-3 -word_size 3 -show_gis -num_alignments 20 -max_hsps 20 -num_threads $num_threads -out $out_file -query $query_file
post_slack_message cluster-jobs "$me : done"

echo "test done" $(date)
