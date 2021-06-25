import sys
from Bio import SeqIO
import pandas as pd
import numpy as np
import os



df = pd.read_csv(sys.argv[1].strip(), sep="\t")


####select high quality genomes

df2=df[(df["Is complete?"] == True) & (df["Is high coverage?"] == True)]


df3=df2[(df2["N-Content"] <= 0.001) | (df2["N-Content"].isna()) ]


LL=df3["Virus name"].to_list()

df3.set_index("Virus name",inplace=True)
print(df3)
outF=open("selected_for_phylogeny.fasta","w")
for i in SeqIO.parse(sys.argv[2].strip(),"fasta"):
    if str(i.id).strip().split("|")[0] in LL:
        i.id=str(i.id).strip().split("|")[0]+"|"+str(df3.loc[str(i.id).strip().split("|")[0]]["Pango lineage"])
        #print(i.id)
        SeqIO.write(i,outF,"fasta")
outF.close()

os.system("cat Reference_SARS-CoV-2.fasta selected_for_phylogeny.fasta NEW.fasta > selected+reference.fasta")

