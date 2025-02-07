
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import numpy as np

# Read sequences from FASTA files
def read_fasta(file_path):
    with open(file_path, "r") as file:
        record = SeqIO.read(file, "fasta")
        return record.seq

seq1 = read_fasta("seq1.fasta")
seq2 = read_fasta("seq2.fasta")

print("Sequence 1:", seq1)
print("Sequence 2:", seq2)

# Define the substitution matrix
substitution_matrix = {
    'A': {'A': 10, 'C': 2, 'G': 5, 'T': 2},
    'C': {'A': 2, 'C': 10, 'G': 2, 'T': 5},
    'G': {'A': 5, 'C': 2, 'G': 10, 'T': 2},
    'T': {'A': 2, 'C': 5, 'G': 2, 'T': 10}
}

gap_penalty = -5

def sequence_alignment(seq1, seq2, substitution_matrix, gap_penalty):
    m, n = len(seq1), len(seq2)
    dp = np.zeros((m + 1, n + 1))
    
    # Initialize the DP table
    for i in range(m + 1):
        dp[i][0] = i * gap_penalty
    for j in range(n + 1):
        dp[0][j] = j * gap_penalty
    
    # Fill the DP table
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = dp[i - 1][j - 1] + substitution_matrix[seq1[i - 1]][seq2[j - 1]]
            delete = dp[i - 1][j] + gap_penalty
            insert = dp[i][j - 1] + gap_penalty
            dp[i][j] = max(match, delete, insert)
    print (dp)
    return dp[m][n]



seq11 = "AATAAT"
seq22= "AAGG"

# Compute the optimal alignment score
optimal_cost = sequence_alignment(seq1, seq2, substitution_matrix, gap_penalty)
print("Optimal Alignment Cost:", optimal_cost)


