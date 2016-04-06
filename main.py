from bam2df import bam2df

align_fn = "/home/laeb/data/storage/NGS/others/Florian/aligned1.bam"

df = bam2df(align_fn)
print(df.head())