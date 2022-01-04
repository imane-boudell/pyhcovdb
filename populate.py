import sqlalchemy
from sqlalchemy import String, Integer, Text
import pandas as pd
import numpy as np
import sys
from sqlalchemy_utils import create_database, database_exists

database_username = 'root'
database_password = 'rootpassword'
database_ip       = 'localhost'
database_name     = 'hcov'
url = 'mysql+pymysql://{0}:{1}@{2}/{3}'.format(database_username, database_password, database_ip, database_name)
if not database_exists(url):
    create_database(url)
database_connection = sqlalchemy.create_engine('mysql+pymysql://{0}:{1}@{2}/{3}'.
                                               format(database_username, database_password, 
                                                      database_ip, database_name))

# Epitopse
epitopes = ['mers', 'sars', 'sarscov2']


for ep in epitopes:
    df = pd.read_csv("epitopes/%s.csv" % ep, keep_default_na=False)
    #df.drop(df.columns[[0,3,4,7,8,17,18,19,20,21,22,23,24,25,26,27]], axis=1, inplace=True)
    df.drop(df.columns[[0]], axis=1, inplace=True)
    df.columns = map(str.lower, df.columns)
    df.columns = map(lambda s: s.replace(" ", "_"), df.columns)

    
    # map types to.. strings (who cares)
    types = {}
    for c in df.columns:
        #print(c)
        types[c] = String(512)
        #print('=========================')

    df.to_sql(con=database_connection, name="%s_epitopes" % ep, index_label='id', chunksize=10, if_exists='replace', dtype=types)

    print(df.head(5))


#df = pd.read_csv("mytest.csv")  # , keep_default_na=False)
df = pd.read_csv("all_viruses.csv")  # , keep_default_na=False)

#df.drop(df.columns[[0]], axis=1, inplace=True)
#df = df[df["gene_product_name"].str.startswith("CHECK_") == False]
print(df.head(5))

types = {"id": Integer(), 'genbank_genome_accession': String(255), 'gene_symbol': String(255),
         'gene_product_name': String(255),
         'genbank_protein_accession': String(255), 'strain_name': String(255), 'isolate': String(255),
         'isolation_source': String(255), 'virus_specimen': String(255), 'host': String(255),
         'collection_date': String(8),
         'country': String(255), 'sequence_type': String(4),'fasta': Text()}



df.to_sql(con=database_connection, name='virus', index_label='id', chunksize=10, if_exists='replace', dtype=types)

# Add unique keys:

with database_connection.connect() as con:
    con.execute('ALTER TABLE `virus` ADD PRIMARY KEY (`id`);')
    con.execute('ALTER TABLE `virus` ADD INDEX (`genbank_genome_accession`);')
    con.execute('ALTER TABLE `virus` ADD INDEX (`genbank_protein_accession`);')
    con.execute('ALTER TABLE `virus` ADD INDEX (`virus_specimen`);')
    con.execute('ALTER TABLE `virus` ADD INDEX (`sequence_type`);')
    con.execute('ALTER TABLE `virus` ADD INDEX (`gene_symbol`);')
    con.execute('ALTER TABLE `virus` ADD INDEX (`strain_name`);')
    con.execute('ALTER TABLE `virus` ADD INDEX (`collection_date`);')
    con.execute('ALTER TABLE `virus` ADD INDEX (`host`);')
    con.execute('ALTER TABLE `virus` ADD INDEX (`country`);')
    con.execute('ALTER TABLE `mers_epitopes` ADD PRIMARY KEY (`id`);')
    con.execute('ALTER TABLE `sars_epitopes` ADD PRIMARY KEY (`id`);')
    con.execute('ALTER TABLE `sarscov2_epitopes` ADD PRIMARY KEY (`id`);')

