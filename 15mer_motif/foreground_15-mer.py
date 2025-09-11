import os
import pandas as pd
from Bio import SeqIO
pd.set_option("expand_frame_repr", False)


"""The process for creating the foreground 15-mer sequences was as follows: Pre-aligned 15-mer peptides 
were centred on each SUMOylated lysine and extended seven residues toward the N-terminus and seven 
residues toward the C-terminus. If the site was located less than seven residues from the N/C-terminus, 
the 15-mer was completed with the letter "X" to reach a length of fifteen residues. 
The 15-mer peptides were then filtered for redundancy to only include unique sequences."""

def get_15_mers(input_csv, fasta_file, output_txt):
    pd.set_option("expand_frame_repr", False)

    df = pd.read_csv(input_csv, sep=',')
    df = df.dropna(axis=0, subset=['Protein position'])

    df_txt = df.drop_duplicates(subset=['Protein'], keep='last')
    id_file = df_txt['Protein']
    wanted = id_file.tolist()
    output_file = "UP000005640_9606_less_seqs.fasta"

    records = (r for r in SeqIO.parse(fasta_file, format="fasta") if r.id in wanted)
    count = SeqIO.write(records, output_file, "fasta")

    if count < len(wanted):
        print("Warning %i IDs not found in %s" % (len(wanted) - count, fasta_file))

    # Parse sequences once and store in dictionary for efficiency
    seq_dict = {}
    for record in SeqIO.parse(output_file, format="fasta"):
        seq_dict[record.id] = record.seq
    
    Mer = []
    dict_df = df.to_dict('records')
    for i in range(len(dict_df)):
        acc = dict_df[i]['Protein']
        pos = int(dict_df[i]['Protein position'])
        
        if acc in seq_dict:
            sequence = seq_dict[acc]
            sR = sequence[pos:pos + 7]
            if pos - 8 >= 0:
                sL = sequence[pos - 8:pos - 1]
            else:
                sL = sequence[0:pos - 1]
            if len(sL) < 7:
                sL = 'X' * (7 - len(sL)) + sL
            if len(sR) < 7:
                sR = sR + 'X' * (7 - len(sR))
            Mer.append(sL + 'K' + sR)

    df_15_mers = pd.DataFrame({'15_mer': Mer})
    print('SUMO site count ' + str(len(df_15_mers)))
    df_15_mers = df_15_mers.drop_duplicates(subset=['15_mer'], keep='first')
    print('15-mer count ' + str(len(df_15_mers)))
    if not os.path.exists('output'):
        os.makedirs('output')
    df_15_mers.to_csv(output_txt, sep=' ', index=False, header=False)
    os.remove('UP000005640_9606_less_seqs.fasta')


get_15_mers('example_input/sumo_site_pos.csv', 'example_input/UP000005640_9606.fasta', 
            'output/15mers.txt')

