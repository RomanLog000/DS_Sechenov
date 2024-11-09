import pandas as pd

# Read the CSV file
# Replace 'input_file.csv' with the path to your input CSV file
df = pd.read_csv('Book1.csv', sep='\t')

# Filter out sequences with gaps
filtered_df = df.dropna(subset=['cdr3aa_pos_rep1', 'cdr3aa_pos_rep2', 'cdr3aa_pos_rep3', 'cdr3aa_pos_rep4', 'cdr3aa_pos_rep5',
                                 'cdr3aa_neg_rep1', 'cdr3aa_neg_rep2', 'cdr3aa_neg_rep3', 'cdr3aa_neg_rep4', 'cdr3aa_neg_rep5'])

# Write the filtered data to a new CSV file
# Replace 'output_file.csv' with the desired path for the output CSV file
filtered_df.to_csv('output_file.csv', index=False)
