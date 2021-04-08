from subprocess import run
from shlex import split as shsplit
from collections import defaultdict, Counter
from Bio import SeqIO
from re import findall
import pandas as pd
from gzip import open as gzopen


def format_mut(kind: str, **kwargs):
    """
    format the mutations in a systematic way
    :param kind: SNP, DEL or INS
    :param kwargs: depending on kind, but all require pos, SNP and INS require seq, and DEL requires len and SNP requires ref
    :return: a systematically formatted mutation string
    """
    if kind == 'SNP':
        return "SNP_{}_{}_{}".format(kwargs['pos'], kwargs['ref'], kwargs['seq'])
    elif kind == 'DEL':
        return "DEL_{}_{}".format(kwargs['pos'], kwargs['len'])
    elif kind == 'INS':
        return "INS_{}_{}".format(kwargs['pos'], kwargs['seq'])


def find_snps(ref_seq, seq):
    """
    Will compare 2 sequences and return a tuple for every mismatch (pos, seq1, seq2)
    :param ref_seq: a str
    :param seq: a str
    :return: list of tuples with info on mismatches
    """
    for ipos, iseq in enumerate(zip(ref_seq, seq)):
        if iseq[0] != iseq[1]:
            yield (ipos, iseq[0], iseq[1])


def run_bowtie2(ref_index: str, fastqs: list, output_name: str, bowtie2_flags="-p 8 -X 700 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50"):
    """
    Runs bowtie2. You can also manually run a mapper, but the the parser might not like the output.
    :param ref_index: name of the bowtie2 index to use
    :param fastqs: list of read pairs [(read1_1_fwd.fastq, read1_2_fwd.fastq),(read2_1_rev.fastq, read2_2_rev.fastq)]
    :param bowtie2_flags: Bowtie2 will be run with "--end-to-end --no-unal --no-hd --no-sq -x {ref_index} -1 {fwd reads} -2 {rev reads} -S {output_name}
    :return: return code of run, 0 = successful
    """

    command = f"bowtie2 {bowtie2_flags} --end-to-end --no-unal --no-hd --no-sq -x {ref_index} -1 {','.join([x[0] for x in fastqs])} -2 {','.join([x[1] for x in fastqs])} -S {output_name}"
    return run(shsplit(command), stderr=open(f"bowtie2_{output_name}.log", 'w')).returncode


def parse_bowtie2_output(bt2_stdout, ref_genome):
    """
    this parses the bowtie2 output into 2 different data structure. A mutation Counter, and a ref_seq mutation dataframe
    :param bt2_stdout: sam file generated by run_bowtie2
    :param ref_genome: path to ref genome in fasta format. Only first entry in fasta file will be used
    :return: (mutation Counter, ref_seq_mutation dataframe)
    """

    read_pairs = defaultdict(set)  # this is for the mutation Counter

    ref_seq = str(next(SeqIO.parse(open(ref_genome), 'fasta')).seq)
    # this is for the ref_seq mutation dataframe
    dna_data = []
    for x in range(len(ref_seq)):
        dna_data.append({'A': 0, 'G': 0, 'C': 0, 'T': 0, 'N': 0, '+': 0, '^': 0})

    with open(bt2_stdout) as f:
        for i, out in enumerate(f.readlines()):
            if i % 100000 == 0:
                print(i)
            # if i > 100000:  # for quicker testing
            #    break
            if out[0] == '@':  # skip headers
                continue
            parts = out.split("\t")
            alignment_flag = int(parts[1].strip())
            if alignment_flag in [83, 99, 147, 163]:  # good alignment
                refpos = int(parts[3].strip()) - 1  # position of alignment (sam is 1 based!!!)
                cigar = parts[5].strip()  # cigar which will tell about indels
                read = parts[9].strip()  # the aligned sequence
                MD = parts[17].strip()[5:]  # MD which will tell of snps and deletions

                # create a list of SNPs, inserts and deletions.
                if str(cigar[:-1]) == str(MD) == str(len(read)):  # if WT
                    read_pairs[parts[0].strip()].update([])  # just update mut counter with a WT (empty set)
                    for ipos, inuc in enumerate(read):
                        dna_data[refpos + ipos][inuc] += 1  # and dna_data
                else:
                    muts = []  # list of mutations
                    # figure out gross things from cigar and SNPs from MD

                    readpos = 0
                    # should probably test the cigar string for characters other than MDI....
                    for segment in findall("\d+[MDI]", cigar):  # split on m's, d's and i's. Better hope theres no soft clipping....
                        # go over each segment and split the number from the letter
                        letter = segment[-1]
                        number = int(segment[:-1])
                        # if letter = M, then look for SNPS and increment refpos and readpos by number
                        if letter is 'M':
                            for spos, sfrom, sto in find_snps(ref_seq[refpos:], read[readpos:number]):
                                muts.append(format_mut('SNP', pos=refpos + spos, ref=sfrom, seq=sto))

                            for ipos, inuc in enumerate(read[readpos:number]):  # go over the positions in the dna
                                dna_data[refpos + ipos][inuc] += 1

                            refpos += number  # ref and read pos moves forward
                            readpos += number
                        # if letter = D, then just report deletion and increment refpos
                        elif letter is 'D':
                            muts.append(format_mut('DEL', pos=refpos, len=number))

                            for x in range(number):
                                # if refpos + x < len(dna_data):
                                dna_data[refpos + x]['^'] += 1

                            refpos += number  # only reference moves forward
                        # if letter = I, then report insertion and increment readpos
                        elif letter is 'I':
                            muts.append(format_mut('INS', pos=refpos, seq=read[readpos:readpos + number]))

                            dna_data[refpos]['+'] += 1

                            readpos += number  # only read post moves forward

                    # finally save the mutations
                    read_pairs[parts[0].strip()].update(muts)

        # Now you have a mutation set for each read, lets count them
    return Counter([frozenset(x) for x in read_pairs.values()]), pd.DataFrame(dna_data)


def read_flash_merged(fastq: str, ref_seq: str, qual_threshold: float, allow_N: bool=False) -> Counter:
    """

    :param fastq: the path in str format to a merged fastq
    :param ref_seq: the reference sequence in str format
    :param qual_threshold: the minimum mean quality threshold
    :param allow_N: keep (true) or remove (false) all sequences with N
    :return: counter
    """
    variants = Counter()
    qual_threshold += 33  # phred33
    ref_len = len(ref_seq)
    blocksize = 4  # read each block in fastq
    with gzopen(fastq) as f:
        for i, line in enumerate(f.readlines()):
            if i % 100000 == 0:
                print(i/4)
            ln = i % blocksize
            if ln == 1:
                seq = line.strip()
                seq_len = len(seq)
            elif ln == 3:
                if seq_len == ref_len:  # skip if length != len(ref_fasta)
                    # per base qual threshold (removed because reverse reads always have at least 1 bad base)
                    # qual_good = True
                    #for q in line.strip():
                    #    if q < qual_threshold:
                    #        qual_good = False
                    #        print(line)
                    #        break

                    # average threshold
                    if sum(seq)/seq_len >= qual_threshold:
                        seq = seq.decode()
                        if seq == ref_seq:  # check if WT
                            variants.update(("WT", ))
                        elif allow_N is False and 'N' in seq:  # check for N in seq
                            variants.update(("NinSeq", ))
                        else:  # else discover SNPs
                            variants.update((frozenset(f"{i}{x}>{y}" for i, (x, y) in enumerate(zip(ref_seq, seq)) if x != y),))
                    else:
                        variants.update(("low_qual", ))
                else:
                    variants.update(("indel", ))
    return variants
