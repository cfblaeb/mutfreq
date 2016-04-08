from sam2df import sam2df
import plots

import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO


#fasta_fn = "/home/laeb/data/storage/NGS/others/Florian/ACC1.fa"
fasta_fn = "/home/laeb/Downloads/ACC1.fa"
fasta = str(SeqIO.read(open(fasta_fn), 'fasta').seq)

try:
    df = pd.read_pickle("lowmut")
except FileNotFoundError:
    #align_fn = "/home/laeb/data/storage/NGS/others/Florian/bow/lowmut.bam"
    align_fn = "/home/laeb/Downloads/lowmut.sam"
    df = sam2df(align_fn, len(fasta))
    df.to_pickle("lowmut")

plt.style.use(['ggplot'])
plots.coverage_plot(df)
plt.show()
#plots.seq_plot(df, fasta, 100, 200)
plots.mut_plot(df, fasta)
plt.show()

