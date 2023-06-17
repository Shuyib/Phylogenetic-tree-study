'''
This script tokenizes the sequences in the sequence_metrics.csv file and compares the 
sequences to each other using cosine similarity. 

'''

import os
import pandas as pd
import numpy as np
from sklearn.preprocessing import LabelEncoder
from tensorflow.keras.preprocessing.sequence import pad_sequences
from sklearn.metrics.pairwise import cosine_similarity

# if not installed, install tensorflow, keras and scikit-learn
# with os module, you can run terminal commands from python
# os.system("pip install tensorflow")
# os.system("pip install keras")
# os.system("pip install scikit-learn")


# change working directory:
os.chdir("updated_data")

# load the sequence_metrics.csv file
df = pd.read_csv("sequence_metrics.csv")

# label encode the sequences in the sequence column & identifier column
# label encoding is used to convert the characters in the sequence column to integers
# copy the original accession id column to a new column called accession_id "]
df2 = df["accession_id"].copy()
df["accession_id"] = LabelEncoder().fit_transform(df["accession_id"])


# label encoding the sequences in the sequence column
df["sequence_split"] = df["sequence"].apply(lambda x: LabelEncoder().fit_transform(list(x)))


# padding the onehot encoded characters to the same length as the longest sequence in the dataset
padded = pad_sequences(df["sequence_split"], padding="post")

# adding the padded sequences to the dataframe
df["padded"] = padded.tolist()

# cosine similarity: a measure of similarity between two non-zero vectors of an inner product space
cosine_similarity = cosine_similarity(df["padded"].tolist())

# get the average cosine similarity score for each sequence
# why? because cosine similarity returns a matrix where each row represents the similarity of one sequence to all other sequences
cosine_similarity = cosine_similarity.mean(axis=1)

# converting the cosine similarity matrix to a dataframe and adding the accession id as the index
cosine_similarity_df = pd.DataFrame(cosine_similarity, index=df["accession_id"])

# resetting the index
cosine_similarity_df = cosine_similarity_df.reset_index()

# renaming column labeled 0 to cosine_similarity average
cosine_similarity_df = cosine_similarity_df.rename(columns={0: "cosine_similarity_average"})

# adding the accession id column back to the dataframe
cosine_similarity_df["accession_id"] = df2

# add encoded column back to the dataframe: 
# call the column label encoded column
cosine_similarity_df["label_encoded"] = df["accession_id"]
cosine_similarity_df = cosine_similarity_df[["accession_id", "label_encoded", "cosine_similarity_average"]]

# remove duplicate rows
cosine_similarity_df = cosine_similarity_df.drop_duplicates("accession_id")

# saving the cosine similarity matrix as a csv file
cosine_similarity_df.to_csv("cosine_similarity.csv")

# how to interpret the cosine similarity scores
# 0.0 - 0.2 : low similarity
# 0.2 - 0.4 : moderate similarity
# 0.4 - 0.6 : high similarity
# 0.6 - 0.8 : very high similarity
# 0.8 - 1.0 : identical sequences

# next steps:
# you could use this for the clustering analysis to see how the sequences are related to each other
# you could also use this for the classification analysis to see if the sequences are similar to each other