#!/bin/bash
#PBS-l nodes=2
#PBS-l walltime=1000:00:00
#PBS -j oe
#PBS -q default
#PBS -o mstrg_out

module load bio/BEDTools/2.26.0-foss-2016b

#This is a script to grab the lines that match the MSTRG ids obtained from the deseq script as significant or nonsignificant genes

F=/data3/marine_diseases_lab/tejashree/Bio_project_SRA/stringtie/test 
echo $F
#set -e
#echo "START" $(date)
#
#args=($(cat $F/mstrgs_sig.txt))
#while read -r line
#do
#    for i in ${args[@]}
#    do
#        case "$line" in
#            *"$i"*) echo "$line";;
#        esac
#    done
#done <$F/"stringtie_merged.gtf"

# I did this step manually on the head node since I got an error permission denied when I ran it through a script. 
#awk '$3=="transcript"' $F/mstrg_sig.gtf > $F/mstrg_sig_noexons.gtf

bedtools getfasta -fi /data3/marine_diseases_lab/tejashree/Bio_project_SRA/cvir_genome/cvir_edited.fa \
        -bed $F/mstrg_sig_noexon.bed -s -name -fo $F/mstrg_sig_noexon.fa.out 

