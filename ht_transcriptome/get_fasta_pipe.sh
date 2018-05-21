#!/bin/bash
#SBATCH --job-name="get_fasta"
#SBATCH --time=9999:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --output="get_fasta_out.%A-%a"
#SBATCH --error="get_fasta_out.%A-%a"

set -e
echo "START" $(date)

project_home=/data3/marine_diseases_lab/tejashree/Bio_project_SRA/
fasta=${project_home}/stringtie/test/trans2017_final/fasta/trans_1

cd $fasta

array1=($(ls -1 ${fasta}/*.txt))
echo "Starting: Get fasta files "
for f in ${array1[@]}; do
    echo "working on $f "
    file_basename=$(basename $f)
    file=${file_basename%.txt}
    /data3/marine_diseases_lab/tejashree/Bio_project_SRA/github/oyster/ht_transcriptome/mstrg_getfasta.sh -i $f -s /data3/marine_diseases_lab/tejashree/Bio_project_SRA/stringtie/test/trans2017_final/stringtie_merged.gtf -o /data3/marine_diseases_lab/tejashree/Bio_project_SRA/stringtie/test/trans2017_final/fasta/trans_1/${file}
    echo "$f done"
done




