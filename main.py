from sam2df import sam2df
import plots

import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO

target = "highmut"

fasta_fn = "/home/laeb/data/storage/NGS/others/Florian/ACC1.fa"
fasta = str(SeqIO.read(open(fasta_fn), 'fasta').seq)

try:
    df = pd.read_pickle(target)
except FileNotFoundError:
    align_fn = "/home/laeb/data/storage/NGS/others/Florian/bow/{}.sam".format(target)
    df = sam2df(align_fn, len(fasta))
    df.to_pickle(target)

plt.style.use(['ggplot'])
plots.coverage_plot(df)
plt.show()
#plots.seq_plot(df, fasta, 100, 200)
plots.mut_plot(df, fasta)
plt.show()

