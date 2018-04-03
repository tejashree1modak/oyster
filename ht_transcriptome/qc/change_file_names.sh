#!/bin/bash

# the purpose of this script it to change trans2017 and any future files names
# from (...R1.fastq, ...R2.fastq) to the corresponding (..._1.fastq, ..._2.fastq)
#                                                                 ^
# Usage:                                                          | 
# change_file_names.sh <directory>                                |
#       where directory has R1, R2.fastq which will be changed to /


if [ $# -ne 1 -o ! -d "$1" ] ; then
    echo "Usage: $0 <dir>"
fi

for r1_file in "$1"/*.R1.fastq ; do
    r2_file="${r1_file%1.fastq}2.fastq"

    if [ -f "$r2_file" ] ; then 
        # create symlinks so that we don't modify the original file name
        r1="${r1_file%.R1.fastq}_1.fastq"
        r2="${r2_file%.R2.fastq}_2.fastq"

        ln -f -s "$r1_file" "$r1" && echo "Created link $r1 -> $r1_file"
        ln -f -s "$r2_file" "$r2" && echo "Created link $r2 -> $r2_file"
    fi
done
