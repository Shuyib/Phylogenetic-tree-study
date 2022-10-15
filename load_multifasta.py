"""
Python utility module to calculate a couple of sequence identification metrics.
"""

from Bio import SeqIO # manipulation of sequences
from Bio.SeqUtils import GC # sequence metric module Bio.SeqUtils
import pandas as pd  # data manipulation library


def GA(seq):
    """Calculate G+C content, return percentage (as float between 0 and 100)

    Parameters
    ----------
    seq : DNA sequence (ATCG)


    Returns
    -------
    Float64
        % of Guanine and Adenine in the given sequence
    
    Example
    -------
    >>> from Bio.SeqUtils import GC
    >>> GA("ACTGN")
    40.0


    """
    g_a = sum(seq.count(x) for x in ["G", "A", "g", "a", "S", "s"])
    try:
        return g_a * 100.0 / len(seq)
    except ZeroDivisionError:
        return 0.0


def calculate_seq_metrics(multifasta):
    """this function calculates a couple of common bioinformatic metrics.

    Parameters
    ----------
    multifasta : UTF-8 formatted multifasta file


    Returns
    -------
    NoneType
        Pandas.DataFrame

    """
    # define empty lists to be populated with data
    sequences = []
    identifiers = []
    gc_content = []
    ga_content = []
    length_seq = []
    
    
    # use a filehandler to open the file
    with open(multifasta, encoding="UTF-8") as handle:
        # foreach record find individual gc content,
        # ga_content, store identifiers, sequences
        # finding the length of each sequence
        # pair the prior stored data using zip
        # append a column name to the pandas dataframe
        for record in SeqIO.parse(handle, "fasta"):
            gc_content.append(GC(str(record)))
            ga_content.append(GA(str(record)))
            identifiers.append(record.id)
            sequences.append(str(record.seq))
            length_seq.append(len(record))
            dataframe = pd.DataFrame(
                zip(identifiers, length_seq, gc_content, ga_content, sequence),
                columns=["accession_id", "seq_length", "gc_content", "ga_content","sequence"],
            )
        return dataframe


if __name__ == "__main__":
    store_df = calculate_seq_metrics("sequences.fasta")
    store_df.to_csv("sequence_metrics.csv")
