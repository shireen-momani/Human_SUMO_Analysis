import os
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
pd.set_option("expand_frame_repr", False)

"""
Performs in silico digestion of proteins using Lys-C protease specificity rules.

Lys-C cleavage specificity:
- Cleaves peptide bonds at the C-terminal side of lysine (K) residues (i.e., cleaves after K).
- Exception: does NOT cleave when lysine is followed by proline ("proline rule").

Input:
- Example FASTA file included: protein_database.fasta
"""

def generate_peptides(proseq, cut_sites, max_missed_cleavages):
    """
    Parameters:
    -----------
    proseq (str): Protein sequence to be digested.  
    cut_sites (list of int): Positions where Lys-C cleaves, including start and end positions.  
    max_missed_cleavages (int): Maximum number of missed cleavage sites allowed.
    """
    peptides = []
    
    # i: starting cut site index
    # j: number of cut sites to skip (1 = no missed cleavage, 2 = 1 missed, etc.)
    for i in range(len(cut_sites)):
        for j in range(1, min(max_missed_cleavages + 2, len(cut_sites) - i + 1)):
            if i + j < len(cut_sites):
                start_pos = cut_sites[i]
                end_pos = cut_sites[i + j]
                peptide_seq = proseq[start_pos:end_pos]
                peptides.append(f"{peptide_seq}, {start_pos + 1}, {end_pos}")
    
    return peptides


def lys_c(proseq):
    """
    Parameters:
    -----------
    proseq (str): The protein sequence to digest
    """
    peptides = []
    cut_sites = [0]  # Start with position 0 (N-terminus)
    
    # Search sequence for Lys-C cleavage sites
    for i in range(0, len(proseq) - 1):
        if proseq[i] == 'K' and proseq[i + 1] != 'P':
            cut_sites.append(i + 1)

    # Add C-terminus as final cut site
    if cut_sites[-1] != len(proseq):
        cut_sites.append(len(proseq))
    
    # Dynamically determine max missed cleavages to prevent index errors
    # Proteins with fewer lysines (cut sites) get proportionally fewer missed cleavages
    num_cuts = len(cut_sites)
    
    if num_cuts >= 10:
        max_missed = 8  # Large proteins: up to 8 missed cleavages
    elif num_cuts == 9:
        max_missed = 7
    elif num_cuts == 8:
        max_missed = 6
    elif num_cuts == 7:
        max_missed = 5
    elif num_cuts == 6:
        max_missed = 4
    elif num_cuts == 5:
        max_missed = 3
    elif num_cuts == 4:
        max_missed = 2
    elif num_cuts == 3:
        max_missed = 1
    else:  
        # No Lys-C sites found - return entire protein as single peptide
        peptides.append(f"{proseq}, {1}, {len(proseq)}")
        return peptides
    
    # Generate all possible peptides
    peptides = generate_peptides(proseq, cut_sites, max_missed)
    
    return peptides


output = open('lys_c_digested_peptides_temp.txt', 'w')
# Process each protein in the FASTA database
for record in SeqIO.parse('example_input/protein_database.fasta', 'fasta'):
    proseq = str(record.seq)
    peptide_list = lys_c(proseq)
    
    for peptide in peptide_list:
        pep_split = peptide.split(', ')
        output.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (
            record.id,           # Protein accession
            pep_split[0],        # Peptide sequence
            pep_split[1],        # Start position on protein
            pep_split[2],        # End position on protein
            len(record.seq),     # Protein length
            record.description   # Protein description
        ))

output.close()

# Handle N-terminal methionine removal
# Many proteins undergo co-translational removal of the initiator methionine
# This creates a mass shift of 131 Da that must be accounted for in MS analysis

df_digest = pd.read_csv('lys_c_digested_peptides_temp.txt', sep='\t', header=None)

df_digest.columns = ['id', 'peptide', 'start', 'end', 'prot len', 'description']

dict_digest = df_digest.to_dict('records')

id = []
peptide = []
start = []
end = []
prot_len = []
description = []

# Generate peptides with N-terminal Met removed
for index in range(len(dict_digest)):
    pep_str = str(dict_digest[index]['peptide'])
    if pep_str[0] == 'M':  # Check if peptide starts with methionine
        id.append(dict_digest[index]['id'])
        peptide.append(dict_digest[index]['peptide'][1:])  # Remove first Met
        start.append(dict_digest[index]['start'] + 1)      # Adjust start position
        end.append(dict_digest[index]['end'])
        prot_len.append(dict_digest[index]['prot len'])
        description.append(dict_digest[index]['description'])

# Create dataframe with Met-cleaved peptides
df = pd.DataFrame()
df['id'] = id
df['peptide'] = peptide
df['start'] = start
df['end'] = end
df['prot len'] = prot_len
df['description'] = description

# Combine original and Met-cleaved peptides
df_ = pd.concat([df_digest, df], axis=0)

# Filter out very short peptides (Adjust based on MS instrument sensitivity)
df_ = df_[(df_['end'] - df_['start'] + 1) >= 7] 


Record_pep = []
dict_df_ = df_.to_dict('records')

for index in range(len(dict_df_)):
    pep_seq = str(dict_df_[index]['peptide'])
    
    # Create unique identifier
    acc = f"{dict_df_[index]['id']}_{dict_df_[index]['prot len']}_{dict_df_[index]['start']}_{dict_df_[index]['end']}"
    
    descr = dict_df_[index]['description']
    record = SeqRecord(Seq(pep_seq), id=acc, description=descr)
    Record_pep.append(record)

if not os.path.exists('output'):
        os.makedirs('output')  
        
# Write peptide database in FASTA format
SeqIO.write(Record_pep, "output/lys_c_digested_peptides.fasta", "fasta")
print("Writing to FASTA file ... Completed")

# Remove temporary files
os.remove("lys_c_digested_peptides_temp.txt")

