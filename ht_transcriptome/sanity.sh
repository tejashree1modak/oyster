#!/bin/bash

:> good
:> incomplete
:> norun

for ((i=1; i <= 400 ;i++)); do
    if ls -1 de5_rif12vs5_blast_out.*-$i  > /dev/null 2>&1 ; then
        # an out file exists for $i, check if it has blastx
        if grep -qw blastx de5_rif12vs5_blast_out.*-$i; then
            if egrep -q "^oktest done" de5_rif12vs5_blast_out.*-$i; then
                echo de5_rif12vs5_blast_out.*-$i | tr ' ' '\n' >> good
            else
                echo de5_rif12vs5_blast_out.*-$i started blastx but did not complete
                echo de5_rif12vs5_blast_out.*-$i | tr ' ' '\n' >> incomplete
            fi
        else
            echo de5_rif12vs5_blast_out.*-$i did not run blastx
            echo de5_rif12vs5_blast_out.*-$i | tr ' ' '\n' >> norun
        fi
    else
        echo $i was not submitted
    fi
done
