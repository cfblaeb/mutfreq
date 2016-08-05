from sam2df import sam2df
from functions import identify_mutations
import plots

import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO

fasta_fn = "/home/laeb/data/storage/NGS/others/Florian/genome/ACC1.fa"
fasta = str(SeqIO.read(open(fasta_fn), 'fasta').seq)

plt.style.use(['ggplot'])
working_dir = "/home/laeb/data/storage/NGS/others/Florian/second"
sid_to_name = {1:'M_9',2:'M_10',3:'M_9_8',4:'M_9_9',5:'H_9',6:'H_10',7:'H_9_8',8:'H_9_9',9:'M_8',10:'M_8_8',11:'M_8_9',12:'M_8_8_8',13:'M_8_8_9',14:'H_8',15:'H_8_8',16:'H_8_9',17:'H_8_8_8',18:'H_8_8_9',19:'H',20:'M'}

for x in range(20):
    target = x+1
    print("{}:{}".format(target, sid_to_name[target]))
    pickle_fn = "{}/pickles/{}.pickle".format(working_dir, target)
    sam_fn = "{}/sams/{}.sam".format(working_dir, target)

    try:
        df = pd.read_pickle(pickle_fn)
    except FileNotFoundError:
        df = sam2df(sam_fn, len(fasta))
        df.to_pickle(pickle_fn)

    f, axs = plt.subplots(nrows=2, sharex=True, squeeze=True)
    plots.coverage_plot(df, axs[0])
    plots.mut_plot(df, axs[1], fasta)
    f.tight_layout()
    f.savefig("{}/figs/svg/{}.svg".format(working_dir, sid_to_name[target]))

    with open('{}/csvs/{}.csv'.format(working_dir, sid_to_name[target]), 'w') as f:
        f.write("Position, Percentage, WT nucleotide, Mutated nucleotide\n")
        for mut in identify_mutations(df, fasta):
            f.write("{}, {}, {}, {}\n".format(mut['pos'], mut['perc'], mut['wt'], mut['mut']))