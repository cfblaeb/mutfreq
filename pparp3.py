__author__ = 'laeb'

#I am trying trimmomatic in order to reduce error rates
#java -jar ~/data/programs/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 8 -trimlog tomlog.txt -basein PARP1-Lib2-control_S1_L001_R1_001.fastq -baseout trimcontrol.fq.gz LEADING:30 TRAILING:30 SLIDINGWINDOW:4:30 MINLEN:80

#reads aligned to PARP1template.fa using either Bowtie2 or other aligner. If the aligner doesnt provide proper MD field, samtools can generate it.
#bowtie2 --very-sensitive -p 8 -x PARP1template -1 ../PARP1-Lib2-1-P6_S1_L001_R1_001.fastq -2 ../PARP1-Lib2-1-P6_S1_L001_R2_001.fastq | samtools view -b - | samtools sort -o selected.bam -T prefix
#Bowtie2 alignments extracted using the following: (other aligners will need change to the field numbers)
#samtools view control.bam | awk '{if ($2==83 || $2==99 || $2==147 || $2==163) {print $2 "\t" $4 "\t" $6 "\t" substr($18,6) "\t" $10} }' > control.txt
#2 is alignment flag (we select only properly paired, both sequences aligned)
#4 is position of alignment
#6 is cigar which will tell about indels
#18 is MD which will tell of snps and deletions
#10 is the aligned sequence

#a given sequence position can mutate to a different residue or a deletion/insertion can appear
#a different res can simply be recorded as such on dna and protein level
#a deletion can simply be reported as that on the dna level but on the protein level this can cause changes to all downstream AA's
#an insert can in principal also be reported the same way...but the inserted material would then need to be on a seperate graph and only the fact that an insert occurs here will be recorded
#and on the protein level they suffer from the same issue as deletions.
#indels on the protein level could be translated..but lets see

from Bio.Seq import translate
import sys

error_file = "/home/laeb/data/storage/NGS/MiSeq/pparp-folding/exp2_2/error.txt"
error_fasta = "gctatcgcactttaacgtttcgtgctgccccctcagtctatgcaatagaccataaactgcaaaaaaaagtccgctgataaggcttgaaaagttcatttccagacccatttttacatcgtagccgatgaggacgcgcctgatgggtgttctggctacctgacctgtccattgtggaaggtcttacattctcgctgatttca"

control_file = "/home/laeb/data/storage/NGS/MiSeq/pparp-folding/exp2_2/control.txt"
sorted_file = "/home/laeb/data/storage/NGS/MiSeq/pparp-folding/exp2_2/sorted.txt"
template_fasta = "AACTTTAAGAAGGAGATATACATATGGTGAACTCCTCTGCTTCAGCAGATAAGCCATTATCCAACATGAAGATCCTGACTCTCGGGAAGCTGTCCCGGAACAAGGATGAAGTGAAGGCCATGATTGAGAAACTCGGGGGGAAGTTGACGGGGACGGCCAACAAGGCTTCCCTGTGCATCAGCACCAAAAAGGAGGTGGAAAAGATGAATAAGAAGATGGAGGAAGTAAAGGAAGCCAACATCCGAGTTGTGTCTGAGGACTTCCTCCAGGACGTCTCCGCCTCCACCAAGAGCCTTCAGGAGTTGTTCTTAGCGCACATCTTGTCCCCTTGGGGGGCAGAGGTGAAGGCAGAGCCTGTTGAAGTTGTGGCCCCAACTAGTCATCATCACCACCATCAT"
template_protein = "MVNSSASADKPLSNMKILTLGKLSRNKDEVKAMIEKLGGKLTGTANKASLCISTKKEVEKMNKKMEEVKEANIRVVSEDFLQDVSASTKSLQELFLAHILSPWGAEVKAEPVEVVAPTSHHHHH"


def analyser_fil(fil, nuc_seq, aa_seq="", orfs=0):
    '''
    fil is a txt file made as described in header
    nuc_seq is nucleotide fasta sequence
    aa_seq is amino acid fasta sequence
    orfs is the nuc at which the open reading frame starts
    '''

    nuc_seq = nuc_seq.upper()
    aa_seq = aa_seq.upper()
    nuc_wt = 0
    AA_wt = 0
    reads = 0

    dna_data = []
    for x in nuc_seq:
        dna_data.append({'A': 0, 'G': 0, 'C': 0, 'T': 0, '^': 0, '+': 0, 'N':0})

    prot_data = []
    if aa_seq:
        for x in aa_seq:
            prot_data.append({'G':0, 'P':0, 'A':0, 'V':0, 'L':0, 'I':0, 'M':0, 'C':0, 'F':0, 'Y':0, 'W':0, 'H':0, 'K':0, 'R':0, 'Q':0, 'N':0, 'E':0, 'D':0, 'S':0, 'T':0, '*': 0, '^': 0, '+': 0, 'X': 0})

    silent_dna = []
    for x in nuc_seq:
        silent_dna.append({'A': 0, 'G': 0, 'C': 0, 'T': 0, '^': 0, '+': 0})

    reads = 0
    with open(fil) as f:
        for alignment in f:
            reads += 1
            if reads%100000==0:
                sys.stdout.write("{}\r".format(reads))
                sys.stdout.flush()
            #if reads > 100000:
            #	break
            parts = alignment.split("\t")
            alignment_flag = int(parts[0].strip())
            pos = int(parts[1].strip())-1 #4 is position of alignment (sam is 1 based so we subtract 1 to make it zero based)
            cigar = parts[2].strip() #6 is cigar which will tell about indels
            MD = parts[3].strip() #18 is MD which will tell of snps and deletions
            seq = parts[4].strip() #10 is the aligned sequence

            if aa_seq:
                #make a translation of the sequence
                inframe_seq = seq[(orfs-pos)%3:-(len(seq[(orfs-pos)%3:])%3) or None] #the in-frame part of the sequence (trims the start and end to fit into protein translation frame
                transed = translate(inframe_seq)

                prot_pos = ((pos+((orfs-pos)%3))-orfs)/3 #how many AA's into the protein sequence does this start
                if prot_pos<0: #our in_frame dna seq starts before the actual protein starts
                    transed = transed[-prot_pos:]
                    prot_pos = 0

                if prot_pos + len(transed) > len(aa_seq): #the translated protein extends beyond the known universe..let truncate it..
                    transed = transed[:len(aa_seq)-(prot_pos + len(transed))]

            #ok now look at the actual sequence alignment and figure out if its just wt or something mutated
            if str(cigar[:-1]) == str(MD) == str(len(seq)): #we got a WT
                #wt counting for fun
                nuc_wt += 1
                AA_wt += 1

                for x in range(pos,pos+len(seq)): #go over the positions in the dna
                    dna_data[x][nuc_seq[x]] += 1 #write WT
                if aa_seq:
                    for x in range(prot_pos, prot_pos+len(transed)):
                        prot_data[x][aa_seq[x]] += 1 #write WT
            else: #something is not WT
                if cigar==str(len(seq))+'M': #its just snp(s)
                    #nucleotides
                    for x, nuc in enumerate(seq):
                        try:
                            dna_data[pos+x][nuc] += 1
                        except:
                            print(sys.exc_info()[0])
                            print("unknown nuc: {}, at pos+x: {}+{}, MD: {}".format(nuc, pos, x, MD))
                            print(seq)

                    if aa_seq:
                        for x, AA in enumerate(transed):
                            try:
                                prot_data[prot_pos+x][AA] += 1
                            except:
                                print(sys.exc_info()[0])
                                print("unknown AA: {}, at pos+x: {}+{}, MD: {}".format(AA, prot_pos, x, MD))
                                print(transed)

                        #if it is WT protein then I should save the nucs in a silent_nucs thing
                        if transed in aa_seq:
                            AA_wt += 1
                            for x, AA in enumerate(transed):
                                cp = orfs + ((prot_pos+x) * 3)
                                #check if nucs are silently mutated
                                for qqq, nuc in enumerate(seq[cp-pos:cp-pos+3]): #go over the codon in the read
                                    if nuc != 'N' and nuc != nuc_seq[cp+qqq]:
                                        try:
                                            #print("silent mutation {},{},{},{},{}".format(MD,cigar,AA, nuc_seq[cp:cp+3], seq[cp-pos:cp-pos+3]))
                                            silent_dna[cp+qqq][nuc] += 1
                                        except:
                                            print(sys.exc_info()[0])
                                            print("silent: unknown nuc: {}, at pos+x: {}+{}, MD: {}".format(nuc, cp, qqq, MD))
                else: #we got an indel
                    #we handle indels in a simple way, We dont count nucs/AA's from indel sequences, but simply mark the position of the indel starting position
                    #the cigar always starts with a number of matched nucs (e.g. 132M)
                    #this is then followed by a number and a D or an I
                    #if we split on M, the zeroth element will be the number of matches and the zeroth plus 1 will be the indel
                    deletion = False
                    insertion = False
                    position = int(cigar.split("M")[0])
                    if cigar.find("D")!=-1:
                        deletion = True
                    if cigar.find("I")!=-1:
                        insertion = True
                    if deletion and insertion:
                        if cigar.find("D")>cigar.find("I"): #the insert was first
                            deletion = False
                        else: #the deletion was first
                            insertion = False
                    if deletion:
                        dna_data[pos+position]['^'] += 1
                    else:
                        dna_data[pos+position]['+'] += 1

                    if aa_seq:
                        if position>=orfs and position<orfs+prot_pos*3+len(transed)*3: #the indel is present in the open reading frame
                            if deletion:
                                indel_prot_pos = ((pos+((orfs-pos)%3))-orfs)/3
                                prot_data[int((pos+position-orfs)/3)]['^'] += 1
                            else:
                                prot_data[int((pos+position-orfs)/3)]['+'] += 1

    print('\t' + "\t".join([str(x) for x in range(len(nuc_seq))]))
    print('\t' + "\t".join(x for x in nuc_seq))
    for nuc in dna_data[0].keys():
        print(str(nuc) + '\t' + "\t".join([str(x[nuc]) for x in dna_data]))

    print("\n")
    if aa_seq:
        print('\t' + "\t".join([str(x) for x in range(len(aa_seq))]))
        print('\t' + "\t".join(x for x in aa_seq))
        for slAA in prot_data[0].keys():
            print(str(slAA) + '\t' + "\t".join([str(x[slAA]) for x in prot_data]))

        print("\n")

        print('\t' + "\t".join([str(x) for x in range(len(nuc_seq))]))
        print('\t' + "\t".join(x for x in nuc_seq))
        for nuc in silent_dna[0].keys():
            print(str(nuc) + '\t' + "\t".join([str(x[nuc]) for x in silent_dna]))

    print("\n")
    afile = open(fil+"_processed.csv", 'w')
    afile.write('\t' + "\t".join([str(x) for x in range(len(nuc_seq))]))
    afile.write("\n")
    afile.write('\t' + "\t".join(x for x in nuc_seq))
    afile.write("\n")
    for nuc in dna_data[0].keys():
        afile.write(str(nuc) + '\t' + "\t".join([str(x[nuc]) for x in dna_data]))
        afile.write("\n")

    afile.write("\n")

    if aa_seq:
        afile.write('\t' + "\t".join([str(x) for x in range(len(aa_seq))]))
        afile.write("\n")
        afile.write('\t' + "\t".join(x for x in aa_seq))
        afile.write("\n")
        for slAA in prot_data[0].keys():
            afile.write(str(slAA) + '\t' + "\t".join([str(x[slAA]) for x in prot_data]))
            afile.write("\n")

        afile.write("\n")

        afile.write('\t' + "\t".join([str(x) for x in range(len(nuc_seq))]))
        afile.write("\n")
        afile.write('\t' + "\t".join(x for x in nuc_seq))
        afile.write("\n")
        for nuc in silent_dna[0].keys():
            afile.write(str(nuc) + '\t' + "\t".join([str(x[nuc]) for x in silent_dna]))
            afile.write("\n")

        afile.write("\n")

    afile.close()

    print("the below is not meaningful. There are 2 pcr products to a sequence. Each are read by a fwd and rev read. Here I am looking individually at the 2 reads for the 2 pcrs individually. I.e. I cannot say if a sequence is truly wt. just that, that one part is...")
    print("there were {} reads".format(reads))
    print(nuc_wt)
    print(AA_wt)


analyser_fil(error_file, error_fasta)
analyser_fil(control_file, template_fasta, template_protein, 23)
analyser_fil(sorted_file, template_fasta, template_protein, 23)
