from pysam import AlignmentFile
from pandas import DataFrame


def bam2df(bam_fn):
    pos_data = []
    bamfile = AlignmentFile(bam_fn, "rb")
    for pileupcolumn in bamfile.pileup():
        nucs = {'A': 0, 'C': 0, 'T': 0, 'G': 0, 'N': 0}
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                nuc = pileupread.alignment.query_sequence[pileupread.query_position]
                if nuc in nucs:
                    nucs[nuc] += 1
                else:
                    print("warning, unknown nuc: {}\nIn read {}\nAt pos {}\n".format(nuc,
                                                                                     pileupread.alignment.query_name,
                                                                                     pileupread.query_position))
        pos_data.append(nucs)
    bamfile.close()
    return DataFrame(pos_data)
