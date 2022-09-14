from Bio import Entrez
import os

Entrez.email = os.getenv("email")

accession_no = ["LN871587.1", "FM207520.1", "HF564650.1", "FM207545.1", "AM157442.1", "U01332.1", "X07714.1", "AF177667.1", "AF177666.1","Z47544.1", "X77334.1", "NR_112172.1", "LT599799.2", "KT378441.1"]

for accession in accession_no:
    val = Entrez.efetch(db="nucleotide", id= accession, rettype="fasta", retmode="text")
    name = str(accession) + str(".fasta")
    with open(name, mode="a") as file:
        file.write(str(val.readlines()))
