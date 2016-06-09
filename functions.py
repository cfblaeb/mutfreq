def set_fasta_index(df, fasta):
    df['fasta'] = [nuc for nuc in fasta]
    return df.set_index('fasta')


def identify_mutations(df, fasta, threshold=0.05, print_it=False):
    df = set_fasta_index(df, fasta)

    df2 = (df.T / df.T.sum()).T
    df2.loc['A', 'A'] = 0
    df2.loc['T', 'T'] = 0
    df2.loc['G', 'G'] = 0
    df2.loc['C', 'C'] = 0

    df2['i'] = range(len(fasta))
    df3 = df2.set_index('i')

    mutspots = df3.loc[(df3.T > threshold).sum()[(df3.T > threshold).sum() > 0].index]
    muts = []
    for i, row in mutspots.iterrows():
        for nuc, perc in row.iteritems():
            if perc>threshold:
                muts.append({'pos': i, 'wt': fasta[i], 'mut': nuc, 'perc': perc})
                if print_it:
                    print("Pos:{}, Perc:{:.1%}, {}-->{}".format(i, perc, fasta[i], nuc))

    return muts
