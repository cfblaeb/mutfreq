from subprocess import run, PIPE
from shlex import split as shsplit
from os.path import basename, splitext, dirname
from collections import defaultdict, Counter
from Bio import SeqIO
from re import findall


def format_mut(kind, **kwargs):
    """
    format the mutations in a systematic way
    :param kind: SNP, DEL or INS
    :param kwargs: depending on kind, but all require pos, SNP and INS require seq, and DEL requires len and SNP requires ref
    :return: a systematically formated mutation string
    """
    if kind == 'SNP':
        return "SNP_{}_{}_{}".format(kwargs['pos'], kwargs['ref'], kwargs['seq'])
    elif kind == 'DEL':
        return "DEL_{}_{}".format(kwargs['pos'], kwargs['len'])
    elif kind == 'INS':
        return "INS_{}_{}".format(kwargs['pos'], kwargs['seq'])


def find_snps(ref_seq, seq):
    """
    Will compare 2 sequences and return a tubple for every mismatch (pos, seq1, seq2)
    :param ref_seq: a str
    :param seq: a str
    :return: list of tubles with info on mismatches
    """
    for ipos, iseq in enumerate(zip(ref_seq, seq)):
        if iseq[0] != iseq[1]:
            yield(ipos, iseq[0], iseq[1])


def run_bowtie2(ref_genome, fastqs):
    """
    Runs bowtie2-build and align
    :param ref_genome: path to reference fasta, should only contain 1 sequence
    :param fastqs: list of read pairs [(read1_1.fastq, read1_2.fastq),(read2_1.fastq, read2_2.fastq)]
    :return: returns output of bowtie stdout
    """
    # create genome index in genome folder
    run(shsplit("bowtie2-build {} {}".format(ref_genome, splitext(basename(ref_genome))[0])), cwd=dirname(ref_genome), stderr=open("bowtie2-build.log", 'w'))
    # align reads
    # You may want:
    #  -N 1 (mismatches in seed alignment, default 0)
    command = "bowtie2 -p 8 --very-sensitive --no-unal --no-head -x {} -1 {} -2 {}".format(
                splitext(ref_genome)[0],
                ",".join([x[0] for x in fastqs]),
                ",".join([x[1] for x in fastqs]))

    return run(shsplit(command), stdout=PIPE, stderr=open("bowtie2.log", 'w'), cwd=dirname(ref_genome), universal_newlines=True)


def parse_bowtie2_output(bt2_stdout, ref_genome):
    read_pairs = defaultdict(set)
    ref_seq = str(next(SeqIO.parse(open(ref_genome), 'fasta')).seq)
    for i, out in enumerate(bt2_stdout.stdout.splitlines()):
        if i % 100000 == 0:
            print(i)
        #if i > 100000:
        #    break
        parts = out.split("\t")
        alignment_flag = int(parts[1].strip())
        if alignment_flag in [83, 99, 147, 163]:  # good alignment
            cigar = parts[5].strip()         # cigar which will tell about indels
            read = parts[9].strip()           # the aligned sequence
            MD = parts[17].strip()[5:]           # MD which will tell of snps and deletions

            # create a list of SNPs, inserts and deletions.
            if str(cigar[:-1]) == str(MD) == str(len(read)):  # if WT
                read_pairs[parts[0].strip()].update([])
            else:
                muts = []  # list of mutations
                # figure out gross things from cigar and SNPs from MD
                refpos = int(parts[3].strip())-1  # position of alignment (sam is 1 based!!!)
                readpos = 0
                # should probably test the cigar string for characters other than MDI....
                for segment in findall("\d+[MDI]", cigar): #split on m's, d's and i's. Better hope theres no soft clipping....
                    #go over each segment and split the number from the letter
                    letter = segment[-1]
                    number = int(segment[:-1])
                    #if letter = M, then look for SNPS and increment refpos and readpos by number
                    if letter is 'M':
                        for spos, sfrom, sto in find_snps(ref_seq[refpos:], read[readpos:number]):
                            muts.append(format_mut('SNP', pos=refpos + spos, ref=sfrom, seq=sto))
                        refpos += number
                        readpos += number
                    #if letter = D, then just report deletion and increment refpos
                    elif letter is 'D':
                        muts.append(format_mut('DEL', pos=refpos, len=number))
                        refpos += number
                    #if letter = I, then report insertion and increment readpos
                    elif letter is 'I':
                        muts.append(format_mut('INS', pos=refpos, seq=read[readpos:readpos+number]))
                        readpos += number

                # finally save the mutations
                read_pairs[parts[0].strip()].update(muts)

    # Now you have a mutation set for each read, lets count them
    return Counter([frozenset(x) for x in read_pairs.values()])

