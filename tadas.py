from glob import glob
import pandas as pd
from matplotlib import pyplot as plt
from align import run_bowtie2, parse_bowtie2_output

data_folder = "/home/laeb/data/bioinf/Tadas"
ref_seq = "{}/wt.fa".format(data_folder)
files = sorted(glob("{}/AmpliconSeq_220716/*.gz".format(data_folder)))
mutcount = parse_bowtie2_output(run_bowtie2(ref_seq, [(files[0], files[1])]), ref_seq)

mutlens = []
for muts, count in mutcount.items():
    mutlens += [len(muts)] * count

pd.Series(mutlens).plot(kind='hist', bins=20)
plt.show()
