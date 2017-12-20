BEGIN {
    FS = "\t"
    OFS = "\t"
    print "type", "ID", "Parent", "product"
}

NF == 9  {
    typ = $3
    split($9, annotations, /;/)
    row["ID"] = row["Parent"] = row["product"] = ""
    for (i=1; i <= length(annotations) ;i++) {
        split(annotations[i], values, /=/)
        if (length(values) == 2 && values[1] ~ /(ID|Parent|product)/) {
            row[values[1]] = values[2]
        }
    }
    if ("product" in row && length(row["product"]) > 0) {
        print typ, row["ID"], row["Parent"], row["product"]
    }
}