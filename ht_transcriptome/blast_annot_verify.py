import csv
import subprocess
import tempfile
import os
import json
import sys

# getting transcripts from gff using awk

#Read in the blast2go_go_table_20171020_1653.txt file as dict. Has headers and is tab separated.
#blast_dict = Seqname , Description
# parsing blast2go files
# Tags    SeqName Description     Length  #Hits   e-Value sim mean        #GO     GO IDs  GO Names        Enzyme Codes    Enzyme Names    InterPro IDs    InterPro GO IDs InterPro GO Names
# [BLASTED, MAPPED, ANNOTATED]    MSTRG.15820     spire homolog 1-like isoform X7 24121   20      0.0     95.4    3       P:GO:0016192; F:GO:0003779; P:GO:0045010        P:vesicle-mediated transport; F:actin binding; P:actin nucleation 
blast_dict = {}
with open('/Users/tejashree/b2gWorkspace/blast2go_go_table_20171020_1653.txt') as blastfile:
    reader = csv.DictReader(blastfile, delimiter = '\t')
    for line in reader:
        seqname = line["SeqName"]
        desc = line["Description"]
        blast_dict[seqname] = desc
print "Parsed", len(blast_dict), "entries from", blastfile.name

#Read in the stringtie_merged.gtf file as a dict
# First few lines of the file will be removed by using cut before giving the file as input
#No header to make keys? Yes 
# parsing gtf files
# NC_007175.2     StringTie       transcript      1       1623    1000    +       .       gene_id "MSTRG.1"; transcript_id "gene39493"; gene_name "COX1"; ref_gene_id "gene39493";
# NC_007175.2     StringTie       exon    1       1623    1000    +       .       gene_id "MSTRG.1"; transcript_id "gene39493"; exon_number "1"; gene_name "COX1"; ref_gene_id "gene39493";
# NC_007175.2     StringTie       transcript      80      358     1000    -       .       gene_id "MSTRG.2"; transcript_id "MSTRG.2.1";
gtf_dict = {}  # { 'gene_id' : [{ 'ref_gene_id':"", 'transcript_id':""}]}

with open("/tmp/sanitized_gtf", "w") as sanitized_gtf:
    gtf_file = '/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2016/genome/transcriptome/transcriptome/design1/stringtie_merged.gtf'

    get_transcripts_script = os.path.join(os.path.dirname(__file__), "get_transcripts.awk")
    subprocess.check_call(["awk", "-f", get_transcripts_script, gtf_file], stdout=sanitized_gtf)


reader = csv.DictReader(open("/tmp/sanitized_gtf"), dialect="excel-tab")
# sanitized_gtf is a file of 3 columns
# gene_id transcript_id ref_gene_id .. ref_gene_id can be empty
for line in reader:
    gene_id = line["gene_id"].strip('"') 
    transcript_id = line["transcript_id"].strip('"') 
    ref_gene_id = line["ref_gene_id"].strip('"') 

    if gene_id in gtf_dict:
        gtf_dict[gene_id].append({'ref_gene_id' : ref_gene_id, 'transcript_id': transcript_id})
    else:
        gtf_dict[gene_id] = [{'ref_gene_id' : ref_gene_id, 'transcript_id': transcript_id}]
print "Created sanitized gtf data genes =", len(gtf_dict), "entries =", sum(len(entries) for entries in gtf_dict.values())

for gene_id in gtf_dict:
    gene_id_rows = [ row for row in gtf_dict[gene_id] if row['ref_gene_id'] ]
    if gene_id_rows:
        ref_gene_ids = { row["ref_gene_id"] for row in gene_id_rows }
        gtf_dict[gene_id] = []
        for rgi in ref_gene_ids:
            gtf_dict[gene_id].append({'ref_gene_id': rgi, "transcript_id" : ''})
    else:
        transcript_id_rows = [ row for row in gtf_dict[gene_id]
                                if row['transcript_id'] and not row['transcript_id'].startswith("MSTRG") ]
        if transcript_id_rows:
            gtf_dict[gene_id] = transcript_id_rows
        else:
            gtf_dict[gene_id] = [ {'transcript_id' : '', 'ref_gene_id' : ''}]
print "Cross referenced sanitized gtf data genes =", len(gtf_dict), "entries =", sum(len(entries) for entries in gtf_dict.values())
        
# gtf dictionary
# { "gene_id" : [  {"transcript_id" : T, ref_gene_id : R}
#                       where if R is "na" there could be one or more T entries that are valid
#                             if R is not "na" there could be one or more R entries for a gene_id (transcript_id will be "na" for such)
#                             else there will be 1 entry in this list, with both "na"
#                ]}

with open("/tmp/sanitized_gff", "w") as sanitized_gff:
    gff_file = '/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2016/genome/transcriptome/transcriptome/files/ref_C_virginica-3.0_top_level.gff3'

    get_transcripts_script = os.path.join(os.path.dirname(__file__), "get_gff_annotations.awk")
    subprocess.check_call(["awk", "-f", get_transcripts_script, gff_file], stdout=sanitized_gff)


gff_ID_dict = {}
gff_Parent_dict = {}

reader = csv.DictReader(open("/tmp/sanitized_gff"), dialect="excel-tab")
for line in reader:
    ID, Parent, product, type = line["ID"], line["Parent"], line["product"], line["type"]
    if Parent:
        gff_Parent_dict.setdefault(Parent, []).append(line)
    else:
        gff_ID_dict[ID] = [line]
print "Created cross referenced gff genes ID entries =", len(gff_ID_dict), ", Parent entries =", len(gff_Parent_dict)
     
# gff_ID_dict
# { "ID" : [ {"Parent": "", "ID": I , type": T, "product": pd} ]}
# gff_Parent_dict
# { "Parent" : [ {"Parent": P, "ID": I, "type": T, "product": pd} ]}

# combine / join gff and gtf

# at this point, gtf_dict will look like
# {
#    'gene1' : [
#                   {'ref_gene_id" : "r1", "transcript_id": ""}      # if ref_gene_id is set, transcript_id will always be empty
#                   {'ref_gene_id" : "r2", "transcript_id": ""}      # if ref_gene_id is set, transcript_id will always be empty
#                   ... more of the same ...
#               ],
#    'gene2' : [
#                   {'ref_gene_id" : "", "transcript_id": "t1"}      # if transcript_id is set, ref_gene_id will always be empty
#                   {'ref_gene_id" : "", "transcript_id": "t2"}      
#                   ... more of the same ...
#               ],
#    'gene3' : [
#                   {'ref_gene_id" : "", "transcript_id": ""}      # Exactly 1 item in this list with both ref_gene_id and transcript_id empty
#               ]
# }
gtf_join_gff = {}
for gene_id, row in gtf_dict.items():
    for gene_anno in row:
        transcript_id, ref_gene_id = gene_anno["transcript_id"], gene_anno["ref_gene_id"]
        infos = []
        if ref_gene_id:
            if ref_gene_id in gff_Parent_dict:
                infos = gff_Parent_dict[ref_gene_id]
        elif transcript_id:
            if transcript_id in gff_ID_dict:
                infos = gff_ID_dict[transcript_id]
        else:
            infos = [{"Parent" : None, "ID":None, "type":None, "product":None}]

        for info in infos:
            gtf_join_gff.setdefault(gene_id, []).append({"ref_gene_id":ref_gene_id,
                            "transcript_id":transcript_id, "type":info["type"], "product":info["product"]})

writer = csv.DictWriter(file("/tmp/gtf_join_gff.csv", "w"), ["gene_id", "ref_gene_id", "transcript_id", "type", "product"])
writer.writeheader()
ix = 0
for gene_id, rows in gtf_join_gff.items():
    ix += len(rows)
    for row in rows:
        line = row.copy()
        line['gene_id'] = gene_id
        writer.writerow(line)
print "Generated gtf-join-gff data base of", len(gtf_join_gff), "genes, with", ix, "total entries"

blast_join = {}
ix = 0
for bgene_id in blast_dict:
    if bgene_id in gtf_join_gff:
        ix += len(gtf_join_gff[bgene_id])
        for row in gtf_join_gff[bgene_id]:
            row = row.copy()
            row["description"] = blast_dict[bgene_id]
            blast_join.setdefault(bgene_id, []).append(row)
    else:
        ix += 1
        blast_join[bgene_id] = [{"ref_gene_id":None, "transcript_id":None, "type":None, "product":None, "description":blast_dict[bgene_id]}]
print "Generated gtf-join-gff-join-blast data base of", len(blast_join), "genes, with", ix, "total entries"

writer = csv.DictWriter(file("/tmp/blast_join.csv", "w"), ["gene_id", "ref_gene_id", "transcript_id", "type", "product", "description"])
writer.writeheader()
for gene_id, rows in blast_join.items():
    for row in rows:
        line = row.copy()
        line['gene_id'] = gene_id
        writer.writerow(line)