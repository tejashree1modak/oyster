BEGIN {
    OFS = "\t"
    FS = "\t"
    print "gene_id", "transcript_id", "ref_gene_id"
}

$0 !~ /^#/ && $3 == "transcript" {
    split($9, annotations, /;[[:space:]]*/)
    row["gene_id"] = row["transcript_id"] = row["ref_gene_id"] = ""
    for (i = 1; i <= length(annotations) ;i++) {
        split(annotations[i], values, /[[:space:]]+/)
        if (values[1] ~ /(gene_id|transcript_id|ref_gene_id)/)  {
            row[values[1]] = values[2]
        }
    }

    print row["gene_id"], row["transcript_id"], row["ref_gene_id"]
}