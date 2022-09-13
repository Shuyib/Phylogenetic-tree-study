from Bio import Entrez
import os

Entrez.email = os.getenv("email")

handle = Entrez.efetch(db="nucleotide", id="X07714.1", rettype="fasta", retmode="text")

handle.readlines()

