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
	
	for i in range(len(p)):
		primer_length = len(p[i])
		for j in range(len(s) - primer_length + 1):
			if (hamming_distance(p[i], s[j:j+primer_length])) <= hamming_max[i]:
				primer_places.append([j,j+primer_length])
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

def randomize_products(primer_products, termination_probability, product_number):
	print("pn")
	print(product_number)
	print("tp")
	print(termination_probability)
	if (termination_probability == "" and product_number == ""):
		return primer_products
	if (termination_probability == ""):
		prod_array = []
		for i in range(int(product_number)):
			element = random.choice(primer_products)
			prod_array.append(element)
		return prod_array
	if (product_number == ""):
		for i in primer_products:
			terminated = False
			# length = 
			while(terminated == False):
				print("todo")

	else:
		for i in range(product_number):
			# choose a random primer product
			print("TODO")
			# shorten it using termination_probability
			print("TODO")
	return 1

# primer products
def primer_products(forward, reverse, t_p, p_n):
    products = []
    
    for i in forward:
        for j in reverse:
            if i[1] < j[0]:
                products.append([i,j])
    TODO = randomize_products(products, t_p, p_n)
    return (products)   # form: [[f_start,fw_end],[rv_s,rv_e]], ...

def chimera_model_1(primer_products):
	proportion = random()	
	product1 = choice(primer_products)
	product2 = choice(primer_products)
	chimerism_1 = [product1[0][0], (product1[1][1]-product1[0][0])*proportion]
	chimerism_2 = [(product2[1][1]-product2[0][0])*proportion, product2[1][1]]	# perhaps 1 char less? (1 char overlap)
	print(str([chimerism_1,chimerism_2]))
	return [chimerism_1,chimerism_2]
