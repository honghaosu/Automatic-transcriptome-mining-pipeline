#!/usr/bin/env python
# coding: utf-8

# In[30]:


import os
import subprocess
import Bio
import sys
import argparse
import pathlib
from pathlib import Path
import tempfile
from Bio import SeqIO, SearchIO
import pandas as pd


# In[6]:


#Define what a directory path is

#Parse arguments       
parser = argparse.ArgumentParser(description='Database named after .protein.fa')
parser.add_argument('-i', '--query', type=str, help='Query sequences')
parser.add_argument('-c', '--db', type=pathlib.Path, help='Directory of database')
parser.add_argument('-e', '--eval', type=float, help='Blast e-value')
parser.add_argument('-q', '--qcov', type=str, help='Query coverage')
parser.add_argument('-bo', '--blasto', type=str, help='Blast output')
parser.add_argument('-o', '--fasta', type=str, help='Name and sequence from Blast')

args = parser.parse_args()


# In[10]:


#Concatenate all peps into the same variable
print("Just for the first round of mining")
fasta_dir = args.db
fasta_files = [f for f in os.listdir(fasta_dir) if f.endswith(".protein.fa")]

database = ''

for fasta_file in fasta_files:
    with open(os.path.join(fasta_dir, fasta_file), 'r') as f:
        database += f.read()


# In[ ]:


#makeblastdb
with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp:
    tmp.write(database)
    tmp_path = tmp.name

def make_blastdb(db_input, db_output_dir):
    db_output_dir_str = str(db_output_dir)
    cmd = ["makeblastdb", "-in", db_input, "-dbtype", "prot", "-out", f"{db_output_dir_str}/blastpdb"]
    subprocess.run(cmd)
    print("makeblastdb completed")

make_blastdb(tmp_path, args.db)

os.remove(tmp_path)


# In[ ]:


#Blastp
query_file = args.query
db_file = f"{args.db}/blastpdb"
output_file = args.blasto
e_value = args.eval
qcov_hsp_perc = args.qcov

cmd = ["blastp", "-query", query_file, "-db", db_file, "-out", output_file, "-evalue", str(e_value), "-qcov_hsp_perc", str(qcov_hsp_perc), "-outfmt", "7"]
subprocess.run(cmd, check=True)


# In[52]:


#Filter

  
process_cmd=['sed', '-e', '/^#/d', output_file] 
filtered=subprocess.check_output(process_cmd)

# Write the output to a new file
tmp="temp.tsv" 
with open(tmp, 'wb') as f:
    f.write(filtered)

blast_df=pd.read_table("temp.tsv", header=None)
print(blast_df.head(10))
blast_df_filtered=blast_df[blast_df[3]>600]#May make this an argument
print("Sequences with amino acid length lower than 600 are removed")
#print(blast_df_filtered.head(10))

blast_sequence=blast_df_filtered[1]
#print(blast_sequence.head(10))
blast_sequence_unique=list(set(blast_sequence))
#print(blast_sequence_unique)

os.remove(tmp)


# In[50]:


#Make a dictionary
with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp:
    tmp.write(database)
    tmp_path = tmp.name

    
def parse_fasta(fasta_file):
    """
    Parses a FASTA file and returns a dictionary of sequences.
    """
    sequences_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    return sequences_dict

sequences_dict = parse_fasta(tmp_path)

sequences = {}
for seq_id, seq_record in sequences_dict.items():
    sequences[seq_id] = str(seq_record.seq)
#The name of the dictionary is called sequences



# In[69]:


#Search against dictionary
Pep_sequence=[]
for key in blast_sequence_unique:
    Pep_sequence.append(sequences[key])

#Combine sequence name and sequence

data_dict = {'Name': blast_sequence_unique, 'Sequence': Pep_sequence}

# Create the DataFrame
Final_result = pd.DataFrame(data_dict)
Final_result.to_csv('temp.tsv', sep='\t', index=False, header=False)
tsv=SeqIO.parse("temp.tsv", "tab")
count = SeqIO.write(tsv, args.fasta, "fasta")
print("Converted %i records" % count)

os.remove("temp.tsv")


# In[58]:


query=args.query
new=args.fasta
cat_cmd=['cat', query, new]
cat_process=subprocess.run(cat_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
with open("merged.fasta", "w") as f:
    f.write(cat_process.stdout.decode())

mafft_cmd=['mafft', '--maxiterate', '1000','--globalpair', 'merged.fasta']
mafft_process = subprocess.run(mafft_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
mafft_output=args.fasta+".mafft"
if mafft_process.returncode != 0:
    print(f"MAFFT exited with an error: {mafft_process.stderr.decode()}")
else:
    # write the output to a file
    with open(mafft_output, "w") as f:
        f.write(mafft_process.stdout.decode())
    print(f"MAFFT finished successfully. Output written to {mafft_output}")
    
os.remove("merged.fasta")


# In[59]:


fasttree_cmd=['FastTree', mafft_output]
fasttree_process=subprocess.run(fasttree_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
fasttree_output=mafft_output+".fst"
if fasttree_process.returncode != 0:
    print(f"FastTree exited with an error: {fasttree_process.stderr.decode()}")
else:
    with open(fasttree_output, "w") as f:
        f.write(fasttree_process.stdout.decode())
    print(f"FastTree finished successfully. Output written to {fasttree_output}")

