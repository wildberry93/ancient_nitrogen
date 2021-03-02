from Bio.SeqUtils.ProtParam import ProteinAnalysis as PA
from Bio import SeqIO
from Bio import AlignIO
import glob, os 
import numpy as np

import statistics 
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

def calculate_pi_full_sequence(sequence):
	"""
	Calculate isoelectric point for the individual sequence (the contemporary sequences proxy)
	"""
	pi = PA(sequence)

	return pi.isoelectric_point()


def calculate_avg_pi_per_sequence(sequence_alignment):
	ips = []
	for rec in sequence_alignment:
		sequence = str(rec.seq).replace("-","")
		ips.append(calculate_pi_full_sequence(sequence))

	return(statistics.mean(ips))

def calculate_pi_per_conserved_pos(sequences, positions):
	"""
	Calculate PI for the conserved positions in all the proteins separately and avg the value
	"""
	pis = []
	for rec in sequences:
		sequence = str(rec.seq).split()
		sequence_positions = "".join([sequence[i] for i in positions])
		sequence_positions_pi = calculate_pi_full_sequence(sequence_positions)
		pis.append(sequence_positions_pi)

	return statistics.mean(pis)


def avg_NOS_per_pos(column):
	"""
	Calculate the average NOS content per column of alignment 
	(simulates the content of NOS in contemporary proteins, as opposed to the conserved acnestors
	"""
	countered = Counter(column)
	total_counts = sum(countered.values())

	N = 0
	O = 0
	S = 0

	for element in countered:
		if element not in aa_elements:
			N += 0
			O += 0
			S += 0
			continue
		elements = aa_elements[element]
		numb_occur = countered[element]

		N += numb_occur*elements["n"]
		O += numb_occur*elements["o"]
		S += numb_occur*elements["s"]

	return N,O,S

def avg_NOS(align_array):
	total_N = 0
	total_O = 0
	total_S = 0
	col_num = 0
        for column in align_array.T:

		col_num += 1
		num_proteins = len(column)
                N,O,S = avg_NOS_per_pos(column)
		N_avg = float(N)/num_proteins
		O_avg = float(O)/num_proteins
		S_avg = float(S)/num_proteins
		total_N += N_avg
		total_O += O_avg
		total_S += S_avg

	final_N = total_N/col_num
	final_O = total_O/col_num
	final_S = total_S/col_num

	return final_N, final_O, final_S

def calc_NOS_enrichment(align_array, majority=0.5):
	conserved_positions = []
	num_conserved_positions = 0
	consensus_sequence = []
        Ns = []
        Os = []
        Ss = []
        for column_id,column in enumerate(align_array.T):
                majority_aa = find_majority(column, majority)
		
                if majority_aa != None:
			consensus_sequence.append(majority_aa)
                        conserved_positions.append(column_id)
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

	consensus_sequence_str = "".join(consensus_sequence)
	return Ns_enrichment, Os_enrichment, Ss_enrichment, num_conserved_positions, consensus_sequence_str, conserved_positions

consensus_sequences_50_file = open("consensus_sequences50.fasta","w")
consensus_sequences_75_file = open("consensus_sequences75.fasta","w")

enrichment_dict = {}
for file in files:
	pfam_name = os.path.basename(file).split(".")[0]
	try:
		align = AlignIO.read(file, "fasta")
	except: continue

	align.get_alignment_length()

	num_conserved_positions = 0
	align_array = np.array([list(rec) for rec in align], np.character)

	avg_pi_per_sequence = calculate_avg_pi_per_sequence(align)

	Ns_enrichment50, Os_enrichment50, Ss_enrichment50, num_cons50, cons_seqs50, cons_pos50 = calc_NOS_enrichment(align_array, 0.5)
	Ns_enrichment75, Os_enrichment75, Ss_enrichment75, num_cons75, cons_seqs75, cons_pos75 = calc_NOS_enrichment(align_array, 0.75)
	N_cont, O_cont, S_cont = avg_NOS(align_array)


	header = ">" + pfam_name
	consensus_sequences_50_file.write(header + "\n" + cons_seqs50 + "\n")
	consensus_sequences_75_file.write(header + "\n" + cons_seqs75 + "\n")
	if cons_seqs50:
		pi_50 = PA(cons_seqs50).isoelectric_point()
		pi_conserved_pos50 = calculate_pi_per_conserved_pos(align, cons_pos50)
	else: 
		pi_50 = "none"
		pi_conserved_pos50 = "none"
	if cons_seqs75:
		pi_75 = PA(cons_seqs75).isoelectric_point()
		pi_conserved_pos75 = calculate_pi_per_conserved_pos(align, cons_pos75)
	else:
		pi_75 = "none"
		pi_conserved_pos75 = "none"

	enrichment_dict[pfam_name] = [Ns_enrichment50, Ns_enrichment75, Os_enrichment50, Os_enrichment75, 
					Ss_enrichment50, Ss_enrichment75, num_cons50, num_cons75,N_cont, O_cont, S_cont, 
					pi_50, pi_75, avg_pi_per_sequence, pi_conserved_pos50, pi_conserved_pos75]  

## merge with taxonomy_dict

taxonomy_file = open("taxonomy_distribution_clan.tab", "r").read().splitlines()
tax_dict = {}
outfile = open("elements_abundance_with_contemporary.tab","w")
outfile.write("\t".join(["pfam","num_bact_phyla","bact_phyla", "num_arch_phyla", "arch_phyla", "Ns_enrichment50", "Ns_enrichment75", "Os_enrichment50", "Os_enrichment75", "Ss_enrichment50", "Ss_enrichment75", "cons_pos50", "cons_pos75", "N_cont", "O_cont", "S_cont", "pi_50", "pi_75","\n"]))
for line in taxonomy_file:
	splitted = line.split("\t")
	tax_dict[splitted[0]] = "\t".join(splitted[1:])
	## dating

	#print(tax_dict[splitted[0]])
	if splitted[0] not in enrichment_dict:
		continue

	enrichments = [str(x) for x in enrichment_dict[splitted[0]]]

	outfile.write(splitted[0] + "\t" + tax_dict[splitted[0]] + "\t" + "\t".join(enrichments) + "\n" )
