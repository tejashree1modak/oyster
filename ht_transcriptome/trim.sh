#!/bin/bash
#PBS-l nodes=2
#PBS-l walltime=1000:00:00
#PBS -j oe
#PBS -o trim.sh_out

#Trimmomatic script to quality trim illumina reads. The script takes input as paired reads and 
#1.removes adaptors
#2.Clips first ten bases for EVERY read
#3.Removes reads that are less than 25 bp in length

set -e 
echo "START $(date)"
module load Trimmomatic/0.32-Java-1.7.0_80 

F=/data3/marine_diseases_lab/tejashree/Bio_project_SRA/trinity
G=/data3/marine_diseases_lab/tejashree/Bio_project_SRA/trinity/trimmomatic

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.32.jar PE -threads 6 -phred33 -trimlog /data3/marine_diseases_lab/tejashree/Bio_project_SRA/trinity/trimmomatic/trimlog \
/data3/marine_diseases_lab/tejashree/Bio_project_SRA/trinity/SRR5357619_1.fq.gz /data3/marine_diseases_lab/tejashree/Bio_project_SRA/trinity/SRR5357619_2.fq.gz \
/data3/marine_diseases_lab/tejashree/Bio_project_SRA/trinity/trimmomatic/forward_paired_out.fq.gz /data3/marine_diseases_lab/tejashree/Bio_project_SRA/trinity/trimmomatic/forward_unpaired_out.fq.gz /data3/marine_diseases_lab/tejashree/Bio_project_SRA/trinity/trimmomatic/reverse_paired_out.fq.gz /data3/marine_diseases_lab/tejashree/Bio_project_SRA/trinity/trimmomatic/reverse_unpaired_out.fq.gz \
ILLUMINACLIP:/opt/software/Trimmomatic/0.32-Java-1.7.0_80/adapters/TruSeq3-PE.fa:2:30:10 HEADCROP:10 MINLEN:25 
done 
