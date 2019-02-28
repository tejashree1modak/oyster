#!/bin/bash
#PBS-l nodes=2
#PBS-l walltime=1000:00:00
#PBS -j oe
#PBS -o trinity.sh_o

# Trinity de novo assembly post trimmomatic quality filtering

set -e 
echo "START $(date)"
module load Trinity/2.4.0-foss-2016b 

Trinity --seqType fq --max_memory 100G \
--left /data3/marine_diseases_lab/tejashree/Bio_project_SRA/trinity/trimmomatic/forward_paired_out.fq.gz /data3/marine_diseases_lab/tejashree/Bio_project_SRA/trinity/trimmomatic/forward_unpaired_out.fq.gz \
--right /data3/marine_diseases_lab/tejashree/Bio_project_SRA/trinity/trimmomatic/reverse_paired_out.fq.gz /data3/marine_diseases_lab/tejashree/Bio_project_SRA/trinity/trimmomatic/reverse_unpaired_out.fq.gz \
--output "/data3/marine_diseases_lab/tejashree/Bio_project_SRA/trinity/trinityout" --CPU 6 
done 
echo "STOP" $(date)

