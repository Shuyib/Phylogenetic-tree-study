"""
Python utility module to calculate a couple of sequence identification metrics.
"""
import os  # operating system module
from Bio import SeqIO  # manipulation of sequences
from Bio.SeqUtils import GC  # sequence metric module Bio.SeqUtils
import pandas as pd  # data manipulation library

# specify working directory
os.chdir("updated_data")


def GA(seq):
    """Calculate G+A content, return percentage (as float between 0 and 100)

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


def kmer_frequency(sequence, k=4):
    """This function finds the substring of length k in a fasta sequence.

    Parameters
    ----------
    sequence : takes a multifasta sequence in UTF-8 format

    k :
         (Default value = 4) number of substrings to generate

    Returns
    -------
    List with different permutations of substrings based on k set in the function.

    """
    return [sequence[x : x + k] for x in range(len(sequence) - k + 1)]


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
                zip(identifiers, length_seq, gc_content, ga_content, sequences),
                columns=[
                    "accession_id",
                    "seq_length",
                    "gc_content",
                    "ga_content",
                    "sequence",
                ],
            )
            dataframe["4mers"] = dataframe["sequence"].map(kmer_frequency)
            dataframe["5mers"] = dataframe["sequence"].map(
                lambda x: kmer_frequency(sequence=x, k=5)
            )
            dataframe["6mers"] = dataframe["sequence"].map(
                lambda x: kmer_frequency(sequence=x, k=6)
            )
        return dataframe


if __name__ == "__main__":
    store_df = calculate_seq_metrics("sequences.fasta")
    store_df.to_csv("sequence_metrics.csv")
