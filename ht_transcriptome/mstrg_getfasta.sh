#!/bin/bash
PROG=${0#*/}
module load bio/BEDTools/2.26.0-foss-2016b
module load genometools/1.5.9-foss-2016b 
module load python/2.7.10

function usage {
    echo "${PROG} [options]"
    echo 
    echo "-i <id-file>          The id file, with one id at in a line"
    echo "-s <stringtie-file>   The stringtie-merged.gtf file"
    echo "-o <output-prefix>    The prefix for the generated outfiles"
}

# https://sookocheff.com/post/bash/parsing-bash-script-arguments-with-shopts/
while getopts "i:s:o:h" arg ; do 
    case "${arg}" in
    i) ID_FILE=$OPTARG ;;
    s) STRINGTIE_MERGED=$OPTARG ;;
    o) OUT_PFX=$OPTARG ;;
    *) usage
       exit 0
       ;;
  esac
done
shift $((OPTIND -1))

if [ -z "${ID_FILE}" -o -z "${STRINGTIE_MERGED}" -o -z "${OUT_PFX}" ]; then
    echo "ERROR: Required options missing"
    usage
    exit 1
fi

GIT=/data3/marine_diseases_lab/tejashree/Bio_project_SRA/github/oyster/ht_transcriptome/

# Validation
# 
# 1) make sure that input id file itself doesn't have any duplicates
# sort -u <input> | uniq -c | awk '$1 > 1' | wc -l => 0
# 

echo -n "Validating id file $ID_FILE ... "
if [ $(sort -u $ID_FILE | uniq -c | awk '$1 > 1' | wc -l | awk '{print $1}') -gt 0 ]; then
    echo "ERROR: Duplicates in $ID_FILE : $(sort -u $ID_FILE | uniq -c | awk '$1 > 1' | wc -l)"
    exit 1
else
    echo done
fi


# 3) search only the transcript sections of the stringtie merged file
OUT=${OUT_PFX}.gtf
awk -F '\t' '$3 == "transcript"' $STRINGTIE_MERGED  | python ${GIT}/get_ids.py search --needle $ID_FILE --haystack - > $OUT
echo "Found $(wc -l $OUT) entries matching $ID_FILE"

# 4) convert the format of $OUT to a suitable format for bedtools
awk -v OFS='\t' -F'\t' '{print $1,$4,$5,$9,$6,$7}' $OUT > ${OUT_PFX}.bed 

# 5) load the bedtools module and run its getfasta
bedtools getfasta -fi /data3/marine_diseases_lab/tejashree/Bio_project_SRA/cvir_genome/cvir_edited.fa -bed ${OUT_PFX}.bed -s -name -fo ${OUT_PFX}.fa.out
echo -n "Validating ${OUT_PFX}.fa.out ... "

if [ "$(grep -c '^>' ${OUT_PFX}.fa.out)" -ne "$(wc -l ${OUT_PFX}.bed | awk '{print $1}')" ]; then
    echo "ERROR: Mismatch in number of sequences in ${OUT_PFX}.fa.out $(grep -c '^>' ${OUT_PFX}.fa.out) vs $(wc -l ${OUT_PFX}.bed | awk '{print $1}')"
    exit 1
else
    echo done
fi

# 6) update the annotations
# and then validate the unique annotations
awk -F '"' ' { if ($1 ~ /^>/) { a=gensub(/.*:([0-9]+)-([0-9]+).*/, "\\1_\\2", $NF); printf ">%s_%s_%s\n", $2, $4, a } else { print } } ' ${OUT_PFX}.fa.out > ${OUT_PFX}.fa
echo -n "Validating ${OUT_PFX}.fa ... "
if [ "$( grep '^>' ${OUT_PFX}.fa | sort | uniq -c | awk '$1 > 1' | wc -l | awk '{print $1}' )" -gt 0 ]; then
    echo "WARNING: Duplicate entries in ${OUT_PFX}.fa"
else
    echo done
fi
    
# 7) use genome tools
gt sequniq --force -o ${OUT_PFX}_unique.fa ${OUT_PFX}.fa
