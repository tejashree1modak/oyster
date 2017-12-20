#!/bin/bash
#PBS-l nodes=2
#PBS-l walltime=1000:00:00
#PBS -j oe
#PBS -q default
#PBS -o mstrg_out

module load bio/BEDTools/2.26.0-foss-2016b
module load genometools/1.5.9-foss-2016b 

#This is a script to grab the lines that match the MSTRG ids obtained from the deseq script as significant or nonsignificant genes
#Set RUN_STEP prior to running the script.

F=/data3/marine_diseases_lab/tejashree/Bio_project_SRA/stringtie/test 
echo $F

#RUN_STEP=mstrg_getinfo
RUN_STEP=gene_getinfo
#RUN_STEP=transcripts 
#RUN_STEP=formatgtf
#RUN_STEP=bedtools 
#RUN_STEP=editfasta
#RUN_STEP=unique

set -e
echo "START" $(date)

#This part of code extracts those lines from the original stringtie_merged.gtf file that 
#match the MSTRG#. It will extract all fields for that id and print into a new file in the same way as in gtf file. 
#The output of this was written to the output file from which I deleted everything else but keep the results and saved it
# as mstrg_sig.gtf that is fed into the second step.
if [ "${RUN_STEP}" == "mstrg_getinfo" ] ; then
    args=($(cat $F/gene_sig.txt))
    while read -r line
    do
        for i in ${args[@]}
        do
            case "$line" in
                *"$i"*) echo "$line";;
            esac
        done
    done <$F/"stringtie_merged.gtf"
fi

#for getting lines that match the geneids and not mstrg
if [ "${RUN_STEP}" == "gene_getinfo" ] ; then
    while read name; do 
    grep "name" $F/"stringtie_merged.gtf"
    done <$F/"gene_sig.txt"
fi 


# I did this step manually on the head node since I got an error permission denied when I ran it through a script. 
if [ "${RUN_STEP}" == "transcripts" ] ; then
    awk '$3=="transcript"' $F/mstrg_sig.gtf > $F/mstrg_sig_noexons.gtf
fi

#To format gtf file for bedtools. 
if [ "${RUN_STEP}" == "formatgtf" ] ; then
    awk -v OFS='\t' -F'\t' '{print $1,$4,$5,$9,$6,$7}' mstrg_sig_noexons.gtf > mstrg_sig_noexons.bed 
fi

#bedtools getfasta extracts sequences from a FASTA file for each of the intervals defined in a BED/GFF/VCF file.
#(http://bedtools.readthedocs.io/en/latest/content/tools/getfasta.html)
#TIPS: 1. The headers in the input FASTA file must exactly match the chromosome column in the BED file.
#      2. BED files containing a single region require a newline character at the end of the line, otherwise a blank output file is produced.
#      3. Using the -name option, one can set the FASTA header for each extracted sequence to be the “name” columns from the BED feature.
if [ "${RUN_STEP}" == "bedtools" ] ; then
    bedtools getfasta -fi /data3/marine_diseases_lab/tejashree/Bio_project_SRA/cvir_genome/cvir_edited.fa \
                        -bed $F/mstrg_sig_noexons.bed -s -name -fo $F/mstrg_sig_noexons.fa.out 
fi

#Edit the output fasta file to only keep the MSTRG#. (can also use the names.awk script to edit .bed file before running bedtools
if [ "${RUN_STEP}" == "editfasta" ] ; then
    awk -F '"' '{if ($1 ~ /^>/) {print ">" $2} else {print}}' mstrg_sig_noexons.fa.out > mstrg_sig_noexons.fa
fi

#Remove duplicate sequences with genometools
#http://genometools.org/tools/gt_sequniq.html
## for mstrg: For the gene_count significant genes: 164 out of 783 sequences have been removed (20.945%0)
## for geneids: 3 out of 24 sequences have been removed (12.500%) 
if [ "${RUN_STEP}" == "unique" ] ; then
    gt sequniq -o $F/mstrg_sig_noexons_unique.fa $F/mstrg_sig_noexons.fa
fi

