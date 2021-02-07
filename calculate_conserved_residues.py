from Bio import SeqIO
from Bio import AlignIO
import glob, os 
import numpy as np

from collections import Counter 

files = glob.glob("domains_sequences/*.aln")

aa_elements = {}
aa_elements_file = open("aa_props.txt","r").read().splitlines()

for line in aa_elements_file[1:]:
	splitted = line.split("\t")
	aa = splitted[0]
	n,o,s = int(splitted[1]), int(splitted[2]), int(splitted[3])
	aa_elements[aa] = {"n":n, "o":o, "s":s}


def calc_N(position):
	"""
	Calculates nitrogen content per alignment column
	"""
	return aa_elements[position]["n"]

def calc_O(position):
        """
        Calculates oxygen content per alignment column
        """
        return aa_elements[position]["o"]

def calc_S(position):
        """
        Calculates sulfur content per alignment column
        """
        return aa_elements[position]["s"]

def find_majority(column, majority=0.5):
	"""
	Find the most common amino acid in a given alignment position (column)
	returns the element or 'None'
	"""
	countered = Counter(column)
	total_counts = sum(countered.values())

	top_element = countered.most_common(1)[0]
	fraction = float(top_element[1])/total_counts

	if fraction >= majority and top_element[0] != "-":
		return top_element[0]

	else:
		return None


def calc_NOS_enrichment(align_array, majority=0.5):
	conserved_positions = []
	num_conserved_positions = 0

        Ns = []
        Os = []
        Ss = []
        for column in align_array.T:
                majority_aa = find_majority(column, majority)

                if majority_aa != None:
                        conserved_positions.append(column)
                        num_conserved_positions += 1
                        ns = calc_N(majority_aa) # number of N atoms in sidechain of the dominating aa in this column
                        Ns.append(ns)
                        oxs = calc_O(majority_aa)
                        Os.append(oxs)
                        ss = calc_S(majority_aa)
                        Ss.append(ss)
	if num_conserved_positions == 0:
		Ns_enrichment = 0
		Os_enrichment = 0
		Ss_enrichment = 0
	else:
        	Ns_enrichment = round(sum(Ns)/float(num_conserved_positions),3)
       		Os_enrichment = round(sum(Os)/float(num_conserved_positions),3)
        	Ss_enrichment = round(sum(Ss)/float(num_conserved_positions),3)

	return Ns_enrichment, Os_enrichment, Ss_enrichment, num_conserved_positions

enrichment_dict = {}
for file in files:
	pfam_name = os.path.basename(file).split(".")[0]
	try:
		align = AlignIO.read(file, "fasta")
	except: continue

	align.get_alignment_length()
	num_conserved_positions = 0
	align_array = np.array([list(rec) for rec in align], np.character)
	Ns_enrichment50, Os_enrichment50, Ss_enrichment50, num_cons50 = calc_NOS_enrichment(align_array, 0.5)
	Ns_enrichment75, Os_enrichment75, Ss_enrichment75, num_cons75 = calc_NOS_enrichment(align_array, 0.75)
	enrichment_dict[pfam_name] = [Ns_enrichment50, Ns_enrichment75, Os_enrichment50, Os_enrichment75, Ss_enrichment50, Ss_enrichment75, num_cons50, num_cons75]  
	#print(pfam_name, Ns_enrichment50, Ns_enrichment75)

## merge with taxonomy_dict

taxonomy_file = open("taxonomy_distribution.tab", "r").read().splitlines()
tax_dict = {}
outfile = open("elements_abundance.tab","w")
outfile.write("\t".join(["pfam","num_bact_phyla","bact_phyla", "num_arch_phyla", "arch_phyla", "Ns_enrichment50", "Ns_enrichment75", "Os_enrichment50", "Os_enrichment75", "Ss_enrichment50", "Ss_enrichment75", "cons_pos50", "cons_pos75"]))
for line in taxonomy_file:
	splitted = line.split("\t")
	tax_dict[splitted[0]] = "\t".join(splitted[1:])
	
	if splitted[0] not in enrichment_dict:
		continue

	enrichments = [str(x) for x in enrichment_dict[splitted[0]]]

	outfile.write(splitted[0] + "\t" + tax_dict[splitted[0]] + "\t" + "\t".join(enrichments) + "\n" )
