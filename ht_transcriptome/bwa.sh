#!/bin/bash
#SBATCH --job-name="bwa"
#SBATCH --time=9999:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --output="bwa_trans2017_final_out.%A-%a"
#SBATCH --error="bwa_trans2017_final_out.%A-%a"

set -e
echo "START" $(date)

#module load BWA/0.7.15-foss-2016b 
module load SAMtools/1.5-foss-2017a
module load bowtie2/2.2.4
module load metabat/0.32.4-1-foss-2016b
#bwa code: 
#

project_home=/data3/marine_diseases_lab/tejashree/Bio_project_SRA/
bac_genomes=${project_home}/bac_genomes
bwa_index=${bac_genomes}/bwa_index/allgenomes_bac_bwa.fna
bac_reference_fasta=${bac_genomes}/ALLgenomes_bac.fna
cvir_index=${project_home}/cvir_genome/cvir_edited.fa
bwa=${project_home}/bwa
trans2017=${bwa}/trans2017_final/bac
#trans2017=${bwa}/trans2017_final/cvir
qc=$project_home/qc/trans2017_final

cd $trans2017
source $project_home/github/oyster/lib/slack.sh

#FILES=(
#    _
#    C_K_0_TACAGC
#    C_M_0_CGGAAT
#    C_V_0_CACGAT
#    RE_K_6_TCCCGA
#    RE_M_6_CTCAGA
#    RE_V_6_CATGGC
#    RI_K_24_TCGGCA
#    RI_K_6_TCATTC
#    RI_M_24_TAATCG
#    RI_M_6_CTATAC
#    RI_V_24_CCAACA
#    RI_V_6_CAGGCG
#    S4_K_24_TCGAAG
#    S4_K_6_TATAAT
#    S4_M_24_GACGAC
#    S4_M_6_CTAGCT
#    S4_V_24_CATTTT
#    S4_V_6_CACTCA
#    )
    
FILES=( 
    _ 
    C_K_0
    C_M_0
    C_V_0
    C_R1
    C_R2
    C_R3
    RE_K_6
    RE_M_6
    RE_V_6
    RE_R1
    RE_R2
    RE_R3
    RI_K_24
    RI_K_6
    RI_M_24
    RI_M_6
    RI_V_24
    RI_V_6
    RIplusRE_R1
    RIplusRE_R2
    RIplusRE_R3
    S4_K_24
    S4_K_6
    S4_M_24
    S4_M_6
    S4_V_24
    S4_V_6
    S4plusRE_R1
    S4plusRE_R2
    S4plusRE_R3
    )

me="bwa::${SBATCH_JOB_NAME}-${SLURM_ARRAY_TASK_ID} on $HOSTNAME"

#bwa mem -t 4 -a /net/fs03/data3/marine_diseases_lab/tejashree/Bio_project_SRA/bac_genomes/bwa_index/allgenomes_bac_bwa.fna SRR5357619_1.fastq SRR5357619_2.fastq > SRR5357619.fastq.sam

# ------------------------------------------------------------------------------------------
#bwa align paired end reads to reference genomes of all bacteria 
#array1=($(ls -1 ${qc}/${FILES[$SLURM_ARRAY_TASK_ID]}_1.fq))
#if [ "${#array1[@]}" == "0" ]; then
#    echo "($SLURM_ARRAY_TASK_ID) ERROR: No input files matching '${qc}/${FILES[$SLURM_ARRAY_TASK_ID]}_1.fq'"
#    post_slack_message cluster-jobs "ERROR: No input files matching '${qc}/${FILES[$SLURM_ARRAY_TASK_ID]}_1.fq'" "$me"
#    exit 1
#fi

#echo "Starting bwa alignment on ${array1[@]}"
#post_slack_message cluster-jobs "Starting bwa alignment on ${array1[@]}" "$me"
#for left_file in ${array1[@]}; do
#    right_file=$(echo ${left_file}|sed s/_1/_2/)
#    file=$(basename ${left_file%_1.fq})
#
#    echo "starting bwa mem -t4 -a ${cvir_index} ${left_file} ${right_file} > ${trans2017}/${file}.sam"
#    post_slack_message cluster-jobs "bwa mem -t4 -a ${cvir_index} ${left_file} ${right_file} > ${trans2017}/${file}.sam" "${me}.${file}"
#    if [ "$debug" ] ; then
#        touch ${trans2017}/${file}.sam
#    else
#        #bwa mem -t4 -a ${bwa_index} ${left_file} ${right_file} > ${trans2017}/${file}.sam
#        bwa mem -t4 -a ${cvir_index} ${left_file} ${right_file} > ${trans2017}/${file}.sam
#    fi
#    echo done
#    post_slack_message cluster-jobs "done" "${me}.${file}"
#done
#echo "Done bwa alignment"
#post_slack_message cluster-jobs "Done bwa alignment" "$me"
#echo "STOP" $(date)

#Samtools sort and convert to bam
#samtools flagtsat
#jgi_summarize_bam_contig_depths

array2=($(ls -1 ${trans2017}/${FILES[$SLURM_ARRAY_TASK_ID]}.sam))
echo "Starting coverage estmation "
for f in ${array2[@]}; do
    echo "working on $f "
    #samtools sort -m 50G -o $f.bam $f
    #samtools flagstat ${f}.bam > ${f}.stat
    #samtools depth ${f}.bam > ${f}_bamtoolsdepth.txt
    echo "running: jgi_summarize_bam_contig_depths  --outputDepth ${f}.jgiDepth.txt --showDepth --referenceFasta ${bac_reference_fasta} $f.bam"
    jgi_summarize_bam_contig_depths  --outputDepth ${f}.jgiDepth.txt --showDepth --referenceFasta ${bac_reference_fasta} $f.bam 
    echo "$f done"
done

