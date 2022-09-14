from Bio import Entrez
import os

Entrez.email = os.getenv("email")

for accession in accession_no:
    val = Entrez.efetch(db="nucleotide", id= accession, rettype="fasta", retmode="text")
    all_vals = str(val.readlines())
    with open("16S_rRNA.fasta", mode="w") as file:
        file.write(all_vals)

