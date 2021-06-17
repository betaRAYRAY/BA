from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
from random import *

def read_sequence(sequence):
	# read input data
	sequence = sequence.upper()
	return sequence

def read_primers(prim_raw):	

	prim_raw = prim_raw.upper()
	primers = prim_raw.split(",")
	# remove empty strings
	primers = [i for i in primers if i]
	# remove duplicates
	primers = list(dict.fromkeys(primers))
	return primers

def hamming_distance(s1, s2): # string 1 and 2
	if len(s1) != len(s2):
		print("error: sequences for hamming distance have to be of same length.")
	else:
		distance = 0
		for i in range(0,len(s1)):
			if s1[i] != s2[i]:
				distance += 1
		return distance

def list_replace(string, list_before, list_after):
	k = 0
	for i in list_before:
		string = string.replace(i, list_after[k])
		k = k+1
	return(string)

def align_hamming(s, p, hamming_max): # sequence, primers, max dist
	primer_places = []	
	for k in p:
		k_len = len(k)
		for i in range(0,len(s)-k_len+1):
			if (hamming_distance(k, s[i:i+k_len])) <= hamming_max:
				primer_places.append([i,i+k_len])
	return(primer_places)

def parse_fasta(path):
	data = []
	with open(path) as handle:
		for values in SimpleFastaParser(handle):
			data.append(values)
	return data 

def build_rerverse_primers(p):
	k = 0
	for i in p:
		prim_seq = Seq(i)
		p[k] =  str(prim_seq.reverse_complement()) # for rna: reverse_complement_rna()
		k = k+1
	return(p)
		
def chimera_model_1(primer_products):
	proportion = random()	
	product1 = choice(primer_products)
	product2 = choice(primer_products)
	chimerism_1 = [product1[0][0], (product1[1][1]-product1[0][0])*proportion]
	chimerism_2 = [(product2[1][1]-product2[0][0])*proportion, product2[1][1]]	# perhaps 1 char less? (1 char overlap)
	print(str([chimerism_1,chimerism_2]))
	return [chimerism_1,chimerism_2]
