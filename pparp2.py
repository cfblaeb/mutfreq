__author__ = 'laeb'

from Bio.SeqIO import parse
import sys

template_fasta = "AACTTTAAGAAGGAGATATACATATGGTGAACTCCTCTGCTTCAGCAGATAAGCCATTATCCAACATGAAGATCCTGACTCTCGGGAAGCTGTCCCGGAACAAGGATGAAGTGAAGGCCATGATTGAGAAACTCGGGGGGAAGTTGACGGGGACGGCCAACAAGGCTTCCCTGTGCATCAGCACCAAAAAGGAGGTGGAAAAGATGAATAAGAAGATGGAGGAAGTAAAGGAAGCCAACATCCGAGTTGTGTCTGAGGACTTCCTCCAGGACGTCTCCGCCTCCACCAAGAGCCTTCAGGAGTTGTTCTTAGCGCACATCTTGTCCCCTTGGGGGGCAGAGGTGAAGGCAGAGCCTGTTGAAGTTGTGGCCCCAACTAGTCATCATCACCACCATCAT"
template_protein = "MVNSSASADKPLSNMKILTLGKLSRNKDEVKAMIEKLGGKLTGTANKASLCISTKKEVEKMNKKMEEVKEANIRVVSEDFLQDVSASTKSLQELFLAHILSPWGAEVKAEPVEVVAPTSHHHHH"

orfs = 23 #open reading frame start


prot_data = []
for x in template_protein:
	prot_data.append({'G':0, 'P':0, 'A':0, 'V':0, 'L':0, 'I':0, 'M':0, 'C':0, 'F':0, 'Y':0, 'W':0, 'H':0, 'K':0, 'R':0, 'Q':0, 'N':0, 'E':0, 'D':0, 'S':0, 'T':0, '*': 0})

dna_data = []
for x in template_fasta:
	dna_data.append({'A': 0, 'G': 0, 'C': 0, 'T': 0})

files = ["/home/laeb/temp/control/out.extendedFrags.fastq", "/home/laeb/temp/selected/out.extendedFrags.fastq"]

def ham_dist(frag, temp, max_distance, pos):  # Calculate Hamming distance between 2 sequences
	anyi = 0
	num_mismatches = 0

	while num_mismatches <= max_distance and anyi < len(frag):
		if frag[anyi] != temp[anyi+pos]:
			num_mismatches += 1
		anyi += 1

	return num_mismatches


def visual_align(s1, s2):
	vis = ''
	for x in range(len(s1)):
		if s1[x] == s2[x]:
			vis += '-'
		else:
			vis += 'X'
	return vis


def align(frag, temp, max_mm):
	for nypos in range(len(temp)-len(frag)+1):
		save_ham = ham_dist(frag, temp, max_mm, nypos)
		if save_ham <= max_mm:
			return (nypos, save_ham)
	return (-1, save_ham)

aps = {}
failed = [0,0]
unique_control_seqs = set()
ham_scores = {}

for fi in files:
	i = 0
	print("Counting file: {}".format(fi))
	for sp in parse(fi, "fastq"):
		i+=1
		if i%100000==0:
			sys.stdout.write("{}\r".format(i))
			sys.stdout.flush()

		if len(sp.seq)==298:
			score = ham_dist(sp.seq,template_fasta,298,0)
			if score>10:
				print(sp.seq)
				print(visual_align(sp.seq,template_fasta[:298]))
				print(template_fasta[:298])

fi = open("control", 'w')
fi.write("alignemnt position\n")
for key,value in aps.iteritems():
	fi.write(str(key) + "\t" + str(value) + "\n")
fi.write("ham scores\n")
for key,value in ham_scores.iteritems():
	fi.write(str(key) + "\t" + str(value) + "\n")
fi.close()