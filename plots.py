def set_fasta_index(df, fasta):
    df['fasta'] = [nuc for nuc in fasta]
    return df.set_index('fasta')


def coverage_plot(df):
    ax = df.T.sum().plot(title='Read depth per position')
    ax.set_ylabel("Read depth")
    ax.set_xlabel("Position")
    return ax


def seq_plot(df, fasta, start, end):
    df = set_fasta_index(df, fasta)

    ax = df[start:end].plot(kind='bar', stacked=True, title='seq plot')
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax.set_ylabel("Read depth")
    return ax


def mut_plot(df, fasta):
    df = set_fasta_index(df, fasta)

    df2 = (df.T / df.T.sum()).T
    df2.loc['A', 'A'] = 0
    df2.loc['T', 'T'] = 0
    df2.loc['G', 'G'] = 0
    df2.loc['C', 'C'] = 0

    df2['i'] = range(len(fasta))
    df3 = df2.set_index('i')
    ax = df3.T.sum().plot(logy=True, title='Mutation frequency per position')
    ax.set_xlabel('Position')
    ax.set_ylabel('Frequency of mutations (0 to 1, log scale)')
    return ax
