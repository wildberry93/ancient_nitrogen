from Bio import SeqIO
from Bio import AlignIO
import glob 
import numpy as np

files = glob.glob("domains_sequences/*.aln")

for file in files[1:2]:
	print(file)
	align = AlignIO.read(file, "fasta")
	align.get_alignment_length()

	align_array = np.array([list(rec) for rec in align], np.character)

	print(align_array[:,0])
