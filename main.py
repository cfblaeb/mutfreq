from sam2df import bam2df
import plots

import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO


fasta_fn = "/home/laeb/data/storage/NGS/others/Florian/ACC1.fa"
fasta = str(SeqIO.read(open(fasta_fn), 'fasta').seq)

try:
    df = pd.read_pickle("lowmut")
except FileNotFoundError:
    align_fn = "/home/laeb/data/storage/NGS/others/Florian/bow/lowmut.bam"
    df = bam2df(align_fn, fasta)
    df.to_pickle("lowmut")


plt.style.use(['ggplot'])
plots.coverage_plot(df)
plt.show()
#plots.seq_plot(df, fasta, 100, 200)
#plots.mut_plot(df, fasta)
#plt.show()

