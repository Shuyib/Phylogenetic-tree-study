# IPython log file

get_ipython().run_line_magic('logstart', 'load_multifasta.py')

from Bio import SeqIO
from Bio.SeqUtils import GC

def GA(seq):
    ga = sum(seq.count(x) for x in ["G", "A", "g", "a", "S", "s"])
    try:
        return ga * 100.0 / len(seq)
    except ZeroDivisionError:
        return 0.0


sequences = []
identifiers = []
gc_content = []
ga_content = [] 
length_seq= []


with open("sequences.fasta") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        gc_content.append(GC(str(record)))
        ga_content.append(GA(str(record)))
        identifiers.append(record.id)
        sequences.append(record.seq)
        length_seq.append(len(record))
        dataframe = pd.DataFrame(zip(identifiers, length_seq, gc_content, ga_content), columns=["accession_id", "seq_length", "gc_content", "ga_content"])
        
exit()
