'''
This script calculates a couple of common bioinformatic metrics. It takes a multifasta file as input and returns a pandas dataframe with the metrics.
'''

# Importing the libraries and modules that we are going to use  to retrieve the data from NCBI
# Required python library: Biopython
# Required modules:
# 1.Entrez for retrieving sequences
# 2.SeqIO for exploring the sequences
# 3.os provides functions for interacting with the operating system

import os
from Bio import Entrez
from Bio import SeqIO


# To use NCBI'S E-Utilities you must specify your email address with each request
# This function grabs your email from your .bashrc file
# it should look like this to work export email="johndoe@provider.com"
Entrez.email = os.getenv("email")

# change working directory:
os.chdir("updated_data")


def multifastaloader(*accession_numbers):
    """This function takes a list of accession numbers as input,
    fetches the fasta sequences from NCBI and then outputs them in a multi fasta file ;sequences.fasta

    Parameters
    ----------
    *accession_numbers : a list of strings with accession numbers used in NCBI to retrieve sequences with those identifiers.


    Returns
    -------
    type
        Nonetype
    A single multifasta file is returned to the working directory.

    """

    try:
        for i in accession_numbers:
            query = Entrez.efetch(
                db="nucleotide", id=i, rettype="fasta", retmode="text"
            )
            record = SeqIO.parse(
                query, "fasta"
            )  # converting the query to a sequence record in biopython inorder to explore the sequences

            # writing the query to a file
            with open("sequences.fasta", "a", encoding="UTF-8") as multi_fasta:
                SeqIO.write(record, multi_fasta, "fasta")
    except Exception as error_message:
        print(f"An error occurred: {error_message}")



accession_no = [
    "LN871587.1",
    "FM207520.1",
    "HF564650.1",
    "FM207545.1",
    "AM157442.1",
    "U01332.1",
    "LN871587.1",
    "FM207520.1",
    "HF564650.1",
    "FM207545.1",
    "AM157442.1",
    "U01332.1",
    "X07714.1",
    "AF177667.1",
    "AF177666.1",
    "Z47544.1",
    "X77334.1",
    "NR_112172.1",
    "LT599799.2",
    "KT378441.1",
]


multifastaloader(accession_no)
