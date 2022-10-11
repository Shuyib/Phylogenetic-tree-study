# IPython log file

get_ipython().run_line_magic('logstart', 'load_multifasta.py')
sequences = []
identifiers = []
from Bio import SeqIO

with open("sequences.fasta") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        identifiers.append(record.id)
        sequences.append(record.seq)
        
exit()
