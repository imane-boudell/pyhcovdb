from Bio import Entrez
from Bio import SeqIO
import xlsxwriter

Entrez.email = "parse@example.com"

gene_file_fasta = open('mers_genes_sequences.fasta', 'w')
protein_file_fasta = open('mers_proteins_sequences.fasta', 'w')
workbook = xlsxwriter.Workbook('mers.xlsx')
worksheet = workbook.add_worksheet()

print("Started processing of MERS sequences...")

list_count = 0
list_retmax = 0

try:
    handle = Entrez.esearch(db="nucleotide", term="MERS-CoV[Organism] AND genomic_dna[filter]")
    rec_list = Entrez.read(handle)
    list_retmax = int(rec_list['RetMax'])
    list_count = int(rec_list['Count'])
except Exception as e:
    pass

try:
    if list_retmax < list_count:
        handle = Entrez.esearch(db="nucleotide", term="MERS-CoV[Organism] AND genomic_dna[filter]", retmax=list_count)
        rec_list = Entrez.read(handle)
except Exception as e:
    pass

id_list = rec_list['IdList']
handle.close()

print("Fetching %s records..." % len(id_list))

'''
id_file = open('id_mers.txt', 'w')
id_file.write("\n".join(id_list))
id_file.close()
'''

import math
import numpy as np

myrecs = list()
n = len(id_list)

if n > 9999: # number of records largers than fetch max
    #split up the list
    c = n / 2000.0
    n_chunks = math.ceil(c)
    for sub_list in np.array_split(id_list, n_chunks):
        subl = sub_list.tolist()
        handle_gb = Entrez.efetch(db="nucleotide", rettype="gb", id=subl)
        recs = list(SeqIO.parse(handle_gb, "gb"))
        myrecs.extend(recs)

else:
    # get all records at once
    handle_gb = Entrez.efetch(db="nucleotide", rettype="gb", id=id_list)
    myrecs = list(SeqIO.parse(handle_gb, "gb"))

handle_gb.close()

print("Collecting %s records..." % len(myrecs))

name = "Unknown"
id = "Unknown"
gene_sequence = "Unknown"
prot_sequence = "Unknown"
product = "Unknown"
protein_id = "Unknown"
collection_date = "Unknown"
release_date = "Unknown"
host = "Unknown"
isolate = "Unknown"
country = "Unknown"
organism = "Unknown"
strain = "Unknown"
isolation_source = "Unknown"

print("Iterating over records...")

i = 0
worksheet.write(0, i ,"Gene ID")
i += 1
worksheet.write(0, i ,"Gene Symbol")
i += 1
worksheet.write(0, i ,"Gene Product")
i += 1
worksheet.write(0, i ,"Protein ID")
i += 1
worksheet.write(0, i ,"Isolate")
i += 1
worksheet.write(0, i ,"Strain")
i += 1
worksheet.write(0, i ,"Isolation Source")
i += 1
worksheet.write(0, i ,"Organism")
i += 1
worksheet.write(0, i ,"Host")
i += 1
worksheet.write(0, i ,"Collection Date")
i += 1
worksheet.write(0, i ,"Release Date")
i += 1
worksheet.write(0, i ,"Country")

row = 1
gene_names  = ["ORF1a", "ORF1ab", "ORF3", "ORF4a", "ORF4b", "ORF5", "ORF8b", "Spike", "NS3", "NS4a", "NS4b", "NS5", "NSP3"]


def word_in (word, phrase):
    return word in phrase.split()


for rec in myrecs:
    #print(rec.name)  #for genebank id
    id = rec.name
    gene_sequence = str(rec.seq)
    gene_flag = False
    source_flag = False
    i = 0

    for feature in rec.features:
        if feature.type == 'CDS' and gene_flag == False:
            if 'gene' in feature.qualifiers:
                name = feature.qualifiers['gene'][0]
                name = name.capitalize()
                if 'orf' in name.lower():
                    tmp = name.lower()
                    name = tmp.replace('orf','ORF')
                elif 'nsp' in name.lower():
                    tmp = name.lower()
                    name = tmp.replace("nsp","NSP")
                elif 'ns' in name.lower():
                    tmp = name.lower()
                    name = tmp.replace("ns","NS")
            if 'product' in feature.qualifiers:
                product = feature.qualifiers['product'][0]
                product = product.capitalize()
                if 'orf' in product.lower():
                    tmp = product.lower()
                    product = tmp.replace("orf","ORF")
                elif product.lower().startswith("nsp"):
                    tmp = product.lower()
                    product = tmp.replace("nsp","NSP")
                elif product.lower().startswith("ns"):
                    tmp = product.lower()
                    product = tmp.replace("ns","NS")
            if 'protein_id' in feature.qualifiers:
                protein_id = feature.qualifiers['protein_id'][0]
            if 'translation' in feature.qualifiers:
                prot_sequence = feature.qualifiers['translation'][0]
            gene_flag = True
            continue
        if feature.type == 'source' and source_flag == False:
            if 'collection_date' in feature.qualifiers:
                collection_date = feature.qualifiers['collection_date'][0]
            if 'host' in feature.qualifiers:
                host = feature.qualifiers['host'][0]
            if 'country' in feature.qualifiers:
                country = feature.qualifiers['country'][0].split(':')[0]
            if 'isolate' in feature.qualifiers:
                isolate = feature.qualifiers['isolate'][0]
            if 'strain' in feature.qualifiers:
                strain = feature.qualifiers['strain'][0]
            if 'isolation_source' in feature.qualifiers:
                isolation_source = feature.qualifiers['isolation_source'][0]
            source_flag = True
            continue

    for label, value in rec.annotations.items():
        if label == 'date': # release date
            release_date = value
        if label == 'organism':
            organism = value
    #print "============================================"
    if "-" in collection_date:
        tmp = collection_date.rsplit('-', 1)[1]
        if len(tmp) != 4: #strip from left instead
            collection_date = collection_date.split('-', 1)[0]
        else:
            collection_date = tmp
    #map missing gene symbols
    if name == 'NA':
        tmp = product.split()[0]
        if len(tmp) == 1:
            name = tmp
        else:
            if 'nucleo' in tmp.lower():
                name = "N"
            else: # map name according to predefined list
                for mystr in gene_names:
                    if word_in(mystr.lower(), tmp.lower()):
                        name = mystr
                        break

    # write to file
    if prot_sequence != "Unknown" and gene_sequence != "Unknown":
        gene_file_fasta.write(">%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s\n" % (id, name, product, protein_id, isolate, strain, isolation_source, organism, host, collection_date, release_date, country,))
        gene_file_fasta.write(gene_sequence)
        gene_file_fasta.write("\n")
        protein_file_fasta.write(">%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s\n" % (id, name, product, protein_id, isolate, strain, isolation_source, organism, host, collection_date, release_date, country,))
        protein_file_fasta.write(prot_sequence)
        protein_file_fasta.write("\n")
        # write to Excel sheet
        worksheet.write(row, i, id)
        i += 1
        worksheet.write(row, i, name)
        i += 1
        worksheet.write(row, i, product)
        i += 1
        worksheet.write(row, i, protein_id)
        i += 1
        worksheet.write(row, i, isolate)
        i += 1
        worksheet.write(row, i, strain)
        i += 1
        worksheet.write(row, i, isolation_source)
        i += 1
        worksheet.write(row, i, organism)
        i += 1
        worksheet.write(row, i, host)
        i += 1
        worksheet.write(row, i, collection_date)
        i += 1
        worksheet.write(row, i, release_date)
        i += 1
        worksheet.write(row, i, country)
        row += 1

    name = "Unknown"
    id = "Unknown"
    gene_sequence = "Unknown"
    prot_sequence = "Unknown"
    product = "Unknown"
    protein_id = "Unknown"
    collection_date = "Unknown"
    release_date = "Unknown"
    host = "Unknown"
    isolate = "Unknown"
    country = "Unknown"
    organism = "Unknown"
    strain = "Unknown"
    isolation_source = "Unknown"

print("Finished processing of MERS sequences...")

gene_file_fasta.close()
protein_file_fasta.close()
workbook.close()
