import csv

# ----
# parsing ; sperated annotations
# gtf and gff files have annotations of kind
#   (gtf) - gene_id "MSTRG.1"; transcript_id "gene39493"; exon_number "1"; gene_name "COX1"; ref_gene_id "gene39493";
#   (gff) - ID=rna1654;Parent=gene957;Dbxref=GeneID:111131527,Genbank:XM_022479118.1;Name=XM_022479118.1;gbkey=mRNA;gene=LOC111131527;model_evidence=Supporting evidence includes similarity to: 4 Proteins;product=tRNA wybutosine-synthesizing protein 2 homolog;transcript_id=XM_022479118.1
# 
def parse_annotations(annotation_line, kind):
    ret = {}
    if kind == 'gtf':
        for value in annotation_line.split(';'):
            value = value.strip().split()
            if len(value) == 2:
                ret[value[0]] = value[1].strip('"')
    elif kind == 'gff':
        for value in annotation_line.split(';'):
            value = value.strip().split('=')
            if len(value) == 2:
                ret[value[0].strip()] = ret[value[1].strip()]
    return ret

# --
# parsing gtf files
# NC_007175.2     StringTie       transcript      1       1623    1000    +       .       gene_id "MSTRG.1"; transcript_id "gene39493"; gene_name "COX1"; ref_gene_id "gene39493";
# NC_007175.2     StringTie       exon    1       1623    1000    +       .       gene_id "MSTRG.1"; transcript_id "gene39493"; exon_number "1"; gene_name "COX1"; ref_gene_id "gene39493";
# NC_007175.2     StringTie       transcript      80      358     1000    -       .       gene_id "MSTRG.2"; transcript_id "MSTRG.2.1";

def parse_gtf(gtf_file):
    gtf_dict = {}
    for line in csv.DictReader(gtf_file, ["chromosome", "col1", "type", "start", "stop", "col5", "strand", "col7", "annotation"]):
        if line["type"] == "transcript":
            annotation = parse_annotations(line['annotation'], "gtf")
            gtf_dict.setdefault(annotation['gene_id'], []).append({'transcript_id' : annotation['transcript_id'],
                                                              'gene_id' : annotation['gene_id'], 'start' : int(line['start']),
                                                              'stop' : int(line['stop'])})
    return gtf_dict

# --
# parsing gff files
# NC_035780.1     RefSeq  region  1       65668440        .       +       .       ID=id0;Dbxref=taxon:6565;Name=1;chromosome=1;collection-date=22-Mar-2015;country=USA;gbkey=Src;genome=chromosome;isolate=RU13XGHG1-28;isolation-source=Rutgers Haskin Shellfish Research Laboratory inbred lines (NJ);mol_type=genomic DNA;tissue-type=whole sample
# NC_035780.1     Gnomon  gene    13578   14594   .       +       .       ID=gene0;Dbxref=GeneID:111116054;Name=LOC111116054;gbkey=Gene;gene=LOC111116054;gene_biotype=lncRNA

def parse_gff(gff_file):
    gff_dict = {}
    for line in csv.DictReader(gff_file, ["chromosome", "col1", "type", "start", "stop", "col5", "strand", "col7", "annotation"]):
        annotation = parse_annotations(line, "gff")
        gff_dict.setdefault(annotation['Parent'], []).append({'product' : annotation['product'],
                                                            'Parent' : annotation['Parent'], 'gff_start' : int(line['start']),
                                                            'gff_stop' : int(line['stop'])})
    return gff_dict

# --
# parsing blast2go files
# Tags    SeqName Description     Length  #Hits   e-Value sim mean        #GO     GO IDs  GO Names        Enzyme Codes    Enzyme Names    InterPro IDs    InterPro GO IDs InterPro GO Names
# [BLASTED, MAPPED, ANNOTATED]    MSTRG.15820     spire homolog 1-like isoform X7 24121   20      0.0     95.4    3       P:GO:0016192; F:GO:0003779; P:GO:0045010        P:vesicle-mediated transport; F:actin binding; P:actin nucleation 

def parse_bast2go(b2g_file):
    b2g_dict = {}
    for line in csv.DictReader(b2g_file):
        b2g_dict.setdefault(line['SeqName'], []).append(line['Description'])
    return b2g_dict

def verify(b2g_file, gtf_file, gff_file):
    b2g_dict = parse_bast2go(b2g_file)
    gtf_dict = parse_gtf(gtf_file)
    gff_dict = parse_gff(gff_file)

    # join b2g_dict -> gtf_dict
    join1 = {} 
    for seq_name, descriptions in b2g_dict.items():
        if seq_name in gtf_dict:
            for description in descriptions:
                for gtf_value in gtf_dict[seq_name]:
                    d = gtf_value.copy()
                    d['description'] = description
                    join1.setdefault(seq_name, []).append(d)


    # join gtf_dict -> gff_dict
    join2 = {}
    for gene_id, values in gtf_dict.items():
        if gene_id in gff_dict:
            for value in values:
                for gff_value in gff_dict[gene_id]:
                    d = value.copy()
                    d.update(gff_value)
                    join2.setdefault(gene_id, []).append(result)