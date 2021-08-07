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
	if (len(primer_products) == 0):
		print("error: no primer products possible for the used parameters")
		return(primer_products)
	if (termination_probability == "" and product_number == ""):
		return primer_products
	if (termination_probability == ""):
		prod_array = []
		for i in range(int(product_number)):
			element = choice(primer_products)
			prod_array.append(element)
		return prod_array
	if (product_number == ""):
		# fw primer (and rv primer if reached) always preserved as a whole
		resultset = []
		for i in primer_products:
			terminated = False
			length = i[1][0] - i[0][1] # start of rv primer - end of fw primer
			j = 0
			while(terminated == False and j <= length):
				randnum = random()
				if (randnum <= float(termination_probability)):
					shortened_product = [[i[0][0], i[0][1]], [i[0][1]+j, i[0][1]+j]] #no rv primer, start and end rv primer same pos
					resultset.append(shortened_product)	
					terminated = True
				j = j+1
			if(terminated == False):
				resultset.append(i)	# i is left unchanged
		return(resultset)	
	else:

		resultset = []
		for i in range(int(product_number)):
			# choose a random primer product
			element = choice(primer_products)
			# shorten it using termination_probability
			terminated = False
			length = element[1][0] - element[0][1] # start of rv primer - end of fw primer
			j = 0
			while(terminated == False and j <= length):
				randnum = random()
				if (randnum <= float(termination_probability)):
					shortened_product = [[element[0][0], element[0][1]], [element[0][1]+j, element[0][1]+j]] #no rv primer, start and end rv primer same pos
					resultset.append(shortened_product)	
					terminated = True
				j = j+1
			if(terminated == False):
				resultset.append(element)	# element is left unchanged
		return(resultset)	

# primer products
def primer_products(forward, reverse, t_p, p_n):
    products = []
    
    for i in forward:
        for j in reverse:
            if i[1] < j[0]:
                products.append([i,j])
    products = randomize_products(products, t_p, p_n)
    return (products)   # form: [[fw_start,fw_end],[rv_s,rv_e]], ...

def chimera_model_1(primer_products):
	proportion = random()	
	product1 = choice(primer_products)
	product2 = choice(primer_products)
	chimerism_1 = [product1[0][0], (product1[1][1]-product1[0][0])*proportion]
	chimerism_2 = [(product2[1][1]-product2[0][0])*proportion, product2[1][1]]	# perhaps 1 char less? (1 char overlap)
	#print(str([chimerism_1,chimerism_2]))
	return [chimerism_1,chimerism_2]
