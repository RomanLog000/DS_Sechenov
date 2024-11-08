# Step 1: Compute alignments for T cell receptor sequences using a weighted BLOSUM62 matrix

# Install Biopython library for sequence alignment
# !pip install biopython

# Import necessary libraries
from Bio import pairwise2
import pandas as pd


# Function to compute alignment score using BLOSUM62 matrix
def compute_alignment_score(seq1, seq2):
    blosum62 = {
        ('A', 'A'): 4, ('R', 'R'): 5, ('N', 'N'): 6, ('D', 'D'): 6, ('C', 'C'): 9,
        ('Q', 'Q'): 5, ('E', 'E'): 5, ('G', 'G'): 6, ('H', 'H'): 8, ('I', 'I'): 4,
        ('L', 'L'): 4, ('K', 'K'): 5, ('M', 'M'): 5, ('F', 'F'): 6, ('P', 'P'): 7,
        ('S', 'S'): 5, ('T', 'T'): 5, ('W', 'W'): 11, ('Y', 'Y'): 7, ('V', 'V'): 4,
        ('B', 'B'): 6, ('Z', 'Z'): 5, ('X', 'X'): 5, ('*', '*'): 1,
        ('A', 'R'): -1, ('A', 'N'): -2, ('A', 'D'): -2, ('A', 'C'): 0, ('A', 'Q'): -1,
        ('A', 'E'): -1, ('A', 'G'): 0, ('A', 'H'): -2, ('A', 'I'): -1, ('A', 'L'): -1,
        ('A', 'K'): -1, ('A', 'M'): -1, ('A', 'F'): -2, ('A', 'P'): -1, ('A', 'S'): 1,
        ('A', 'T'): 0, ('A', 'W'): -3, ('A', 'Y'): -2, ('A', 'V'): 0,
        # Continue with remaining entries...
    }

    # Initialize alignment score
    alignment_score = 0

    # Compute alignment score for each aligned pair of residues
    for residue1, residue2 in zip(seq1, seq2):
        if (residue1, residue2) in blosum62:
            alignment_score += blosum62[(residue1, residue2)]
        else:
            alignment_score += blosum62[(residue2, residue1)]  # Use reverse order if pair not found

    return alignment_score


# Read T cell receptor sequences from table
# Replace 'your_table.csv' with the path to your table file
tcr_data = pd.read_csv('CMV_2exp_cleaned.csv', sep='\t')

# Iterate over pairs of sequences and compute alignment scores
alignment_scores = []
for i in range(len(tcr_data)):
    for j in range(i + 1, len(tcr_data)):
        for k in range(1, 6):  # Loop over technical replicas
            seq1_pos = tcr_data.loc[i, f'cdr3aa_pos_rep{k}']
            seq2_pos = tcr_data.loc[j, f'cdr3aa_pos_rep{k}']
            seq1_neg = tcr_data.loc[i, f'cdr3aa_neg_rep{k}']
            seq2_neg = tcr_data.loc[j, f'cdr3aa_neg_rep{k}']

            alignment_score_pos = compute_alignment_score(seq1_pos, seq2_pos)
            alignment_score_neg = compute_alignment_score(seq1_neg, seq2_neg)

            alignment_scores.append((i, j, k, alignment_score_pos, alignment_score_neg))

# Convert alignment scores to DataFrame
alignment_df = pd.DataFrame(alignment_scores,
                            columns=['index1', 'index2', 'replica', 'alignment_score_pos', 'alignment_score_neg'])

# Save alignment scores to a CSV file
alignment_df.to_csv('alignment_scores.csv', index=False)
