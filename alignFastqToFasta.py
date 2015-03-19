from Bio.SeqIO import parse


__author__ = 'laeb'

# read fastq reads one by one
# (option) maybe trim low quality ends
# align to fasta by finding lowest hamming distance
#
# (special case) if running with paired end reads that dont actually overlap, try mapping paired read as well
#
# output a list of the number of a's,g's,c's and t's per pos in template
# output a list aa changes

debug = False
template_fasta = "AACTTTAAGAAGGAGATATACATATGGTGAACTCCTCTGCTTCAGCAGATAAGCCATTATCCAACATGAAGATCCTGACTCTCGGGAAGCTGTCCCGGAACAAGGATGAAGTGAAGGCCATGATTGAGAAACTCGGGGGGAAGTTGACGGGGACGGCCAACAAGGCTTCCCTGTGCATCAGCACCAAAAAGGAGGTGGAAAAGATGAATAAGAAGATGGAGGAAGTAAAGGAAGCCAACATCCGAGTTGTGTCTGAGGACTTCCTCCAGGACGTCTCCGCCTCCACCAAGAGCCTTCAGGAGTTGTTCTTAGCGCACATCTTGTCCCCTTGGGGGGCAGAGGTGAAGGCAGAGCCTGTTGAAGTTGTGGCCCCAACTAGTCATCATCACCACCATCAT"
#template_protein = "VNSSASADKPLSNMKILTLGKLSRNKDEVKAMIEKLGGKLTGTANKASLCISTKKEVEKMNKKMEEVKEANIRVVSEDFLQDVSASTKSLQELFLAHILSPWG"
template_protein = "MVNSSASADKPLSNMKILTLGKLSRNKDEVKAMIEKLGGKLTGTANKASLCISTKKEVEKMNKKMEEVKEANIRVVSEDFLQDVSASTKSLQELFLAHILSPWGAEVKAEPVEVVAPTSHHHHH"

# 26-334 is the in-frame protein
orfs = 23 #open reading frame start
#orfe = 334 # open reading frame end

prot_data = []
for x in template_protein:
	prot_data.append({'G':0, 'P':0, 'A':0, 'V':0, 'L':0, 'I':0, 'M':0, 'C':0, 'F':0, 'Y':0, 'W':0, 'H':0, 'K':0, 'R':0, 'Q':0, 'N':0, 'E':0, 'D':0, 'S':0, 'T':0, '*': 0})

dna_data = []
for x in template_fasta:
	dna_data.append({'A': 0, 'G': 0, 'C': 0, 'T': 0})

files = [["/home/laeb/temp/PARP1-Lib2-control_S1_L001_R1_001.fastq",1],["/home/laeb/temp/PARP1-Lib2-control_S1_L001_R2_001.fastq",0]]

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


def align(frag, temp):
	for nypos in range(len(temp)-len(frag)+1):
		save_ham = ham_dist(frag, temp, 10, nypos)
		if save_ham <= 10:
			return (nypos, save_ham)
	return (-1, save_ham)

aps = {}
failed = [0,0]
unique_control_seqs = set()
ham_scores = {}

for fi in files:
	i = 0
	print("Counting file: {}".format(fi[0]))
	for sp in parse(fi[0], "fastq"):
		i+=1
		if i%10000==0:
			print(i)

		#if i>100000:
		#	break

		if fi[1]==1:
			seq = sp.seq
		else:
			seq = sp.seq.reverse_complement()

		ap_ham = align(str(seq), template_fasta)
		ap = ap_ham[0]
		if ap >= 0: #if we have a good match and (as is expected) it is forward then add it to the data
		#if ah < 5 and ap >= 240: #if we have a good match and (as is expected) it is forward then add it to the data
			#save info regarding alignment position for qual control
			if ap in aps.keys():
				aps[ap]+=1
			else:
				aps[ap]=1
			#save hams scores
			if ap_ham[1] in ham_scores.keys():
				ham_scores[ap_ham[1]]+=1
			else:
				ham_scores[ap_ham[1]]=1
			#save unique seqs for qual control
			unique_control_seqs.add(seq)
			#count AA's for prot_data
			inframe_seq = seq[(orfs-ap)%3:-(len(seq[(orfs-ap)%3:])%3) or None] #the in-frame part of the sequence (trims the start and end to fit into protein translation frame

			transed = inframe_seq.translate()

			prot_pos = ((ap+((orfs-ap)%3))-orfs)/3 #how many AA's into the protein sequence does this start
			if prot_pos<0: #our in_frame dna seq starts before the actual protein starts
				transed = transed[-prot_pos:]
				prot_pos = 0

			if prot_pos + len(transed) > len(template_protein): #the translated protein extends beyond the known universe..let truncate it..
				transed = transed[:len(template_protein)-(prot_pos + len(transed))]

			#print(template_protein)
			#print(' ' * prot_pos + visual_align(str(transed), template_protein[prot_pos:prot_pos + len(str(transed))]))
			#print(' ' * prot_pos + str(transed))

			for pos, let in enumerate(transed):
				try:
					prot_data[prot_pos+pos][let] += 1
				except:
					print("odd...AA is not in my standard list of AA's. Maybe its an X caused by an N in the read. anyway im not counting this one")
					print(i)
					print(ap)
					print(pos)
					print(let)
					print(inframe_seq)

			#count nuc bases in dna_data
			for pos, let in enumerate(seq):
				try:
					dna_data[ap+pos][let] += 1
				except:
					print("odd...nuc is not a,t,c or g. Is it perhaps an N? anyway im not counting this one")
					print(i)
					print(ap)
					print(pos)
					print(let)
					print(seq)
		else:
			failed[fi[1]]+=1

fi = open("control_data", 'w')
fi.write('\t' + "\t".join([str(x) for x in range(len(template_fasta))]))
fi.write("\n")
fi.write('\t' + "\t".join(x for x in template_fasta))
fi.write("\n")
for nuc in dna_data[0].keys():
	fi.write(str(nuc) + '\t' + "\t".join([str(x[nuc]) for x in dna_data]) + "\n")

fi.write('\n\n')

fi.write('\t' + "\t".join([str(x) for x in range(len(template_protein))]))
fi.write("\n")
fi.write('\t' + "\t".join(x for x in template_protein))
fi.write("\n")
for slAA in prot_data[0].keys():
	fi.write(str(slAA) + '\t' + "\t".join([str(x[slAA]) for x in prot_data]) + "\n")
fi.close()


#qual control
fi = open("control", 'w')
fi.write("failed\n")
fi.write(str(failed[0]) + "\t" + str(failed[1]) + "\n")
fi.write("alignemnt position\n")
for key,value in aps.iteritems():
	fi.write(str(key) + "\t" + str(value) + "\n")
fi.write("ham scores\n")
for key,value in ham_scores.iteritems():
	fi.write(str(key) + "\t" + str(value) + "\n")
fi.write("unique aligned seqs\n")
fi.write(str(len(unique_control_seqs)) +"\n")
fi.close()
print('\npart2\n')
#part2
prot_data = []
for x in template_protein:
	prot_data.append({'G':0, 'P':0, 'A':0, 'V':0, 'L':0, 'I':0, 'M':0, 'C':0, 'F':0, 'Y':0, 'W':0, 'H':0, 'K':0, 'R':0, 'Q':0, 'N':0, 'E':0, 'D':0, 'S':0, 'T':0, '*': 0})

dna_data = []
for x in template_fasta:
	dna_data.append({'A': 0, 'G': 0, 'C': 0, 'T': 0})

files = [["/home/laeb/temp/PARP1-Lib2-1-P6_S1_L001_R1_001.fastq",1],["/home/laeb/temp/PARP1-Lib2-1-P6_S1_L001_R2_001.fastq",0]]

aps = {}
failed = [0,0]
unique_seqs = set()
ham_scores = {}

for fi in files:
	i = 0
	print("Counting file: {}".format(fi[0]))
	for sp in parse(fi[0], "fastq"):
		i+=1
		if i%10000==0:
			print(i)

		#if i>100000:
		#	break

		if fi[1]==1:
			seq = sp.seq
		else:
			seq = sp.seq.reverse_complement()

		ap_ham = align(str(seq), template_fasta)
		ap = ap_ham[0]
		if ap >= 0: #if we have a good match and (as is expected) it is forward then add it to the data
		#if ah < 5 and ap >= 240: #if we have a good match and (as is expected) it is forward then add it to the data
			#save info regarding alignment position for qual control
			if ap in aps.keys():
				aps[ap]+=1
			else:
				aps[ap]=1
			#save hams scores
			if ap_ham[1] in ham_scores.keys():
				ham_scores[ap_ham[1]]+=1
			else:
				ham_scores[ap_ham[1]]=1
			#save unique seqs for qual control
			unique_seqs.add(seq)
			#count AA's for prot_data
			inframe_seq = seq[(orfs-ap)%3:-(len(seq[(orfs-ap)%3:])%3) or None] #the in-frame part of the sequence (trims the start and end to fit into protein translation frame

			transed = inframe_seq.translate()

			prot_pos = ((ap+((orfs-ap)%3))-orfs)/3 #how many AA's into the protein sequence does this start
			if prot_pos<0: #our in_frame dna seq starts before the actual protein starts
				transed = transed[-prot_pos:]
				prot_pos = 0

			if prot_pos + len(transed) > len(template_protein): #the translated protein extends beyond the known universe..let truncate it..
				transed = transed[:len(template_protein)-(prot_pos + len(transed))]

			#print(template_protein)
			#print(' ' * prot_pos + visual_align(str(transed), template_protein[prot_pos:prot_pos + len(str(transed))]))
			#print(' ' * prot_pos + str(transed))

			for pos, let in enumerate(transed):
				try:
					prot_data[prot_pos+pos][let] += 1
				except:
					print("odd...AA is not in my standard list of AA's. Maybe its an X caused by an N in the read. anyway im not counting this one")
					print(i)
					print(ap)
					print(pos)
					print(let)
					print(inframe_seq)

			#count nuc bases in dna_data
			for pos, let in enumerate(seq):
				try:
					dna_data[ap+pos][let] += 1
				except:
					print("odd...nuc is not a,t,c or g. Is it perhaps an N? anyway im not counting this one")
					print(i)
					print(ap)
					print(pos)
					print(let)
					print(seq)
		else:
			failed[fi[1]]+=1

fi = open("data", 'w')
fi.write('\t' + "\t".join([str(x) for x in range(len(template_fasta))]))
fi.write("\n")
fi.write('\t' + "\t".join(x for x in template_fasta))
fi.write("\n")
for nuc in dna_data[0].keys():
	fi.write(str(nuc) + '\t' + "\t".join([str(x[nuc]) for x in dna_data]) + "\n")

fi.write('\n\n')

fi.write('\t' + "\t".join([str(x) for x in range(len(template_protein))]))
fi.write("\n")
fi.write('\t' + "\t".join(x for x in template_protein))
fi.write("\n")
for slAA in prot_data[0].keys():
	fi.write(str(slAA) + '\t' + "\t".join([str(x[slAA]) for x in prot_data]) + "\n")
fi.close()

fi = open("selected", 'w')
fi.write("failed\n")
fi.write(str(failed[0]) + "\t" + str(failed[1]) + "\n")
fi.write("alignemnt position\n")
for key,value in aps.iteritems():
	fi.write(str(key) + "\t" + str(value) + "\n")
fi.write("ham scores\n")
for key,value in ham_scores.iteritems():
	fi.write(str(key) + "\t" + str(value) + "\n")
fi.write("unique aligned seqs\n")
fi.write(str(len(unique_seqs)) +"\n")
fi.close()

stuff_in_control_but_not_in_selected = unique_control_seqs.difference(unique_seqs)
stuff_in_selected_but_not_in_control = unique_seqs.difference(unique_control_seqs)

print(len(stuff_in_control_but_not_in_selected))
print(len(stuff_in_selected_but_not_in_control))
print(stuff_in_selected_but_not_in_control)