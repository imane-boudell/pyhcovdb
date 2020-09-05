from Bio import SeqIO
import pandas as pd


file_names = ['mers', 'sars', 'sarscov2']
all_data = [] # creating a list of lists

header_list = ['genbank_genome_accession','gene_symbol','gene_product_name','genbank_protein_accession','strain_name','isolate','isolation_source','virus_specimen','host','collection_date','country','fasta','sequence_type']

for my_gene_file in file_names:
    # get gene sequences
    print("Processing genes sequences...")
    fasta_sequences = SeqIO.parse(open(my_gene_file+"_genes_sequences.fasta"), 'fasta')
    for fasta in fasta_sequences:
        mylist = []
        sequence = str(fasta.seq)
        header = fasta.description
        split_list = header.split("|")
        id = split_list[0]
        mylist.append(id)
        name = split_list[1]
        mylist.append(name)
        product = split_list[2]
        mylist.append(product)
        protein_id = split_list[3]
        mylist.append(protein_id)
        strain = split_list[5]
        mylist.append(strain)
        isolate = split_list[4]
        mylist.append(isolate)
        isolation_source = split_list[6]
        mylist.append(isolation_source)
        #organism = split_list[7]
        organism = my_gene_file+'_virus'
        mylist.append(organism)
        host = split_list[8]
        mylist.append(host)
        collection_date = split_list[9]
        mylist.append(collection_date)
        country = split_list[11]
        mylist.append(country)
        #manipulate sequence then add it
        seq_format = ">" + id + "|" + name + "|" + product + "|" + protein_id + "|" + strain + "|" + isolate + "|" + isolation_source + "|" + organism + "|" + host + "|" + collection_date + "|" + country + "\n" + sequence + "\n"
        mylist.append(seq_format)
        mylist.append("nucl")
        #print mylist
        all_data.append(mylist)

    # get prot  sequences
    print("Processing protein sequences...")
    fasta_sequences = SeqIO.parse(open(my_gene_file + "_proteins_sequences.fasta"), 'fasta')
    for fasta in fasta_sequences:
        mylist = []
        sequence = str(fasta.seq)
        header = fasta.description
        split_list = header.split("|")
        id = split_list[0]
        mylist.append(id)
        name = split_list[1]
        mylist.append(name)
        product = split_list[2]
        mylist.append(product)
        protein_id = split_list[3]
        mylist.append(protein_id)
        strain = split_list[5]
        mylist.append(strain)
        isolate = split_list[4]
        mylist.append(isolate)
        isolation_source = split_list[6]
        mylist.append(isolation_source)
        #organism = split_list[7]
        organism = my_gene_file+'_virus'
        mylist.append(organism)
        host = split_list[8]
        mylist.append(host)
        collection_date = split_list[9]
        mylist.append(collection_date)
        country = split_list[11]
        mylist.append(country)
        # manipulate sequence then add it
        seq_format = ">" + id + "|" + name + "|" + product + "|" + protein_id + "|" + strain + "|" + isolate + "|" + isolation_source + "|" + organism + "|" + host + "|" + collection_date + "|" + country + "\n" + sequence + "\n"
        mylist.append(seq_format)
        mylist.append("prot")
        # print mylist
        all_data.append(mylist)
        
df = pd.DataFrame(all_data, columns=header_list)
#print(df.head(5))
print("Writing to file...")
df.to_csv("all_viruses.csv", encoding='utf-8', index=False, line_terminator='\n')
