from fastapi import FastAPI
from db import Virus, MersEp, SarsEp, Sarscov2Ep, row2dict, result2dict
from sqlalchemy import and_, or_, func
from typing import List
from starlette.middleware.cors import CORSMiddleware
from collections import defaultdict
from Bio.Align.Applications import MuscleCommandline, ClustalOmegaCommandline
import uuid
import os
import subprocess

app = FastAPI()

#cors
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_credentials=True, allow_methods=["*"], allow_headers=["*"])

def runCommand(command):
    return subprocess.run(command, stdout=subprocess.PIPE).stdout.decode('utf-8')

def getUniqeValues(field_name, filter_):
    col = getattr(Virus, field_name)
    result = Virus.query.with_entities(col).filter(filter_).distinct().all()
    buff = [i[0] for i in result]
    buff = list(sorted(buff, key=lambda x: (x is None or x == "N/A", x)))
    
    return buff

def getUniqeValues2(field_name1, field_name2, filter_):
    col1 = getattr(Virus, field_name1)
    col2 = getattr(Virus, field_name2)
    result = Virus.query.with_entities(col1,col2).filter(filter_).distinct().all()

    return result

@app.get("/")
def read_root():
    return row2dict(Sarscov2Ep.query.first())

@app.post("/algo/msa/{msa_type}")
def algo_msa(msa_type: str, seq_id: List[int], consensus: bool = None):
    if len(seq_id) > 10:
        return "Cannot process more than 10 sequences for MSA. Operation aborted."

    result = Virus.query.with_entities("id", "fasta").filter(Virus.id.in_(seq_id))
    
    result_dict = {}
    for r in result:
        result_dict[r[0]] = r[1]

    fasta_file = "/root/tmp/%s" % str(uuid.uuid4())
    with open(fasta_file, "w") as fasta:
        # Ensure ordering of sequences based on input
        for i in seq_id:
            fasta.write(result_dict[i]+ "\n\n")
    msa_command = None
    if msa_type == "muscle":
        msa_command = MuscleCommandline("muscle", input=fasta_file, html=True, quiet=True)
        ret = msa_command()
    elif msa_type == "clustalo":
        msa_command = ClustalOmegaCommandline(infile=fasta_file)
        ret = msa_command()
    else: # if msa_type == "mview":
        clustal_file = "/root/tmp/%s" % str(uuid.uuid4())
        msa_command = ClustalOmegaCommandline(infile=fasta_file, outfile=clustal_file)
        msa_command()
        con = "on" if consensus else "off"
        ret = runCommand(["mview", "--css", "on", "--pcid", "aligned", "--ruler", "on", "--width", "80", 
            "-coloring", "mismatch", "-colormap", "pink", "-consensus", con, "-con_threshold", "100", 
            "-html", "head", "-in", "fasta", clustal_file])
        os.remove(clustal_file)

    os.remove(fasta_file)

    return ret 

@app.post("/map/by_criteria/{virus_specimen}")
def read_map_criteria(virus_specimen: str, gene_symbol: str = None, host: str = None, 
        country: str = None, collection_date: int = None):

    ftr = Virus.virus_specimen == virus_specimen

    if gene_symbol != None:
        ftr = and_(ftr, Virus.gene_symbol == gene_symbol)
    
    if host != None:
        ftr = and_(ftr, Virus.host == host)

    if country != None:
        ftr = and_(ftr, Virus.country == country)
    
    if collection_date != None:
        ftr = and_(ftr, Virus.collection_date == collection_date)

    result = Virus.query.with_entities("country").filter(ftr).all()
    ret = defaultdict(int)
    
    for i in result:
        ret[i[0]] += 1

    return [[k,v] for k,v in ret.items()]

@app.get("/viruses/search/by_accession/{accession_number}")
def read_virus(accession_number: str):
    result = Virus.query.filter(or_(Virus.genbank_genome_accession == accession_number, 
        Virus.genbank_protein_accession == accession_number)).all()
    return result2dict(result)


@app.get("/epitopes/{specimen}")
def read_epitopes(specimen: str):
    EpTable = None
    if specimen == "sars":
        EpTable = SarsEp
    elif specimen == "mers":
        EpTable = MersEp
    else:
        EpTable = Sarscov2Ep

    result = EpTable.query.all()
    return result2dict(result) 

@app.post("/viruses/search/by_criteria/{sequence_type}/{virus_specimen}")
def read_virus_by_criteria(virus_specimen: str, sequence_type: str, gene_symbol: List[str] = None, host: List[str] = None, 
        country: List[str] = None, collection_date: List[str] = None):
    
    in_gene_symbol = True if gene_symbol is None else Virus.gene_symbol.in_(gene_symbol)
    in_host = True if host is None else Virus.host.in_(host)
    in_country = True if country is None else Virus.country.in_(country)
    in_collection_date = True if collection_date is None else Virus.collection_date.in_(collection_date)

    and_filter = and_(Virus.sequence_type == sequence_type, Virus.virus_specimen == virus_specimen, 
        in_gene_symbol, in_host, in_country, in_collection_date)
    
    result = Virus.query.filter(and_filter)
    ret = result2dict(result.all())
    return ret

@app.post("/viruses/search_criteria/result_count/{sequence_type}/{virus_specimen}")
def read_search_criteria_count(virus_specimen: str, sequence_type: str, gene_symbol: List[str] = None, host: List[str] = None, 
        country: List[str] = None, collection_date: List[str] = None):
    

    in_gene_symbol = True if gene_symbol is None else Virus.gene_symbol.in_(gene_symbol)
    in_host = True if host is None else Virus.host.in_(host)
    in_country = True if country is None else Virus.country.in_(country)
    in_collection_date = True if collection_date is None else Virus.collection_date.in_(collection_date)

    and_filter = and_(Virus.sequence_type == sequence_type, Virus.virus_specimen == virus_specimen, in_gene_symbol, in_host, in_country, in_collection_date)

    result = Virus.query.filter(and_filter)
    ret = result.count()

    return ret

@app.post("/viruses/search_criteria/{sequence_type}/{virus_specimen}")
def read_search_criteria_ex(virus_specimen: str, sequence_type: str, gene_symbol: List[str] = None, host: List[str] = None, 
        country: List[str] = None, collection_date: List[str] = None):

    in_gene_symbol = True if gene_symbol is None else Virus.gene_symbol.in_(gene_symbol)
    in_host = True if host is None else Virus.host.in_(host)
    in_country = True if country is None else Virus.country.in_(country)
    in_collection_date = True if collection_date is None else Virus.collection_date.in_(collection_date)

    and_filter = and_(Virus.sequence_type == sequence_type, Virus.virus_specimen == virus_specimen, 
        in_gene_symbol, in_host, in_country, in_collection_date)

    ret = {}
    #ret["genes"] = getUniqeValues2("gene_symbol", "gene_product_name", and_filter)
    ret["gene_symbol"] = getUniqeValues("gene_symbol", and_filter)
    ret["gene_product_name"] = getUniqeValues("gene_product_name", and_filter)
    ret["host"] = getUniqeValues("host", and_filter)
    ret["country"] = getUniqeValues("country", and_filter)
    ret["collection_date"] = getUniqeValues("collection_date", and_filter)
    

    return ret
