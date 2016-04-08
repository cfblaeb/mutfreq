from pysam import AlignmentFile
from pandas import DataFrame


def bam2df(bam_fn, fasta):
    df = DataFrame(index=range(len(fasta)), columns=['A', 'T', 'C', 'G', '^', '+', 'N']).fillna(0)

    bamfile = AlignmentFile(bam_fn, "rb")

    reads = 0
    for read in bamfile.fetch():
        if not read.is_secondary and not read.is_unmapped:
            reads += 1
            if reads % 10 == 0:
                print(reads)
            #Go over cigar
            df_pos = read.reference_start
            read_pos = 0
            for code, l in read.cigartuples:
                # M 0 alignment match (can be a sequence match or mismatch)"
                # I 1 insertion to the reference
                # D 2 deletion from the reference
                # N 3 skipped region from the reference
                # S 4 soft clipping (clipped sequences present in SEQ)
                # H 5 hard clipping (clipped sequences NOT present in SEQ)
                # P 6 padding (silent deletion from padded reference)
                # = 7 sequence match
                # X 8 sequence mismatch

                if code not in [0, 1, 2]:
                    print("WARNING, unhandled cigar code: {}, {}".format(code, read.cigarstring))
                else:
                    if code == 0:  # match, mark df_pos with base at read_pos and move both ahead
                        for tracker in range(l):
                            df.loc[df_pos, read.seq[read_pos]] += 1
                            df_pos += 1
                            read_pos += 1
                    elif code == 1:  # insertion, mark df_pos as insert and move read_pos ahead
                        df.loc[df_pos, '+'] += 1
                        for tracker in range(l):
                            read_pos += 1
                    elif code == 2:  # deletion, mark df_pos with deletion at each pos and move df_pos ahead
                        for tracker in range(l):
                            df.loc[df_pos, '^'] += 1
                            df_pos += 1

    bamfile.close()
    return DataFrame(df)
