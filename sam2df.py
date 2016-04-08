from pandas import DataFrame


def sam2df(sam_fn, fasta_len):
    reads = 0
    DI_reads = 0
    dna_data = []
    for x in range(fasta_len):
        dna_data.append({'A': 0, 'G': 0, 'C': 0, 'T': 0, 'N': 0})

    with open(sam_fn) as f:
        for read in f:
            if not read[0] == '@':
                parts = read.split("\t")
                flag = int(parts[1])
                reads += 1
                if reads % 100000 == 0:
                    print(reads)

                if flag in [83, 99, 147, 163]:
                    cigar = parts[5]
                    if "I" in cigar or "D" in cigar:
                        DI_reads += 1
                    else:
                        pos = int(parts[3]) - 1
                        seq = parts[9]

                        for x in range(len(seq)):  # go over the positions in the dna
                            dna_data[x+pos][seq[x]] += 1  # write WT

    print("{} reads, {} DI, {:.2%}".format(reads, DI_reads, DI_reads/reads))
    return DataFrame(dna_data)
