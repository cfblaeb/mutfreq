from functions import set_fasta_index

#OLD STYLE PLOTS THAT REQUIRE OLD STYLE SUMMARY DF. USE (not yet made) FUNCTION TO CONVERT MuTATION COUNTER TO SUMMARY_DF
def coverage_plot(df, ax):
    df.T.sum().plot(ax=ax, title='Read depth per position')
    ax.set_ylabel("Read depth")
    ax.set_xlabel("Position")
    ax.set_ylim(0)

#OLD STYLE PLOTS THAT REQUIRE OLD STYLE SUMMARY DF. USE (not yet made) FUNCTION TO CONVERT MuTATION COUNTER TO SUMMARY_DF
def seq_plot(df, ax, fasta, start, end):
    df = set_fasta_index(df, fasta)

    df[start:end].plot(ax=ax, kind='bar', stacked=True, title='seq plot')
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax.set_ylabel("Read depth")

#OLD STYLE PLOTS THAT REQUIRE OLD STYLE SUMMARY DF. USE (not yet made) FUNCTION TO CONVERT MuTATION COUNTER TO SUMMARY_DF
def mut_plot(df, ax, fasta):
    df = set_fasta_index(df, fasta)

    df2 = (df.T / df.T.sum()).T
    df2.loc['A', 'A'] = 0
    df2.loc['T', 'T'] = 0
    df2.loc['G', 'G'] = 0
    df2.loc['C', 'C'] = 0

    df2['i'] = range(len(fasta))
    df3 = df2.set_index('i')
    df3.T.sum().plot(ax=ax, logy=True, title='Mutation frequency per position')
    ax.set_xlabel('Position')
    ax.set_ylabel('Frequency of mutations (0 to 1, log scale)')
