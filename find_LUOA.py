intax = open("taxonomy_distribution.tab","r").read().splitlines()
tax_dict = {}

for line in intax:
	splitted = line.split("\t")
	pfam = splitted[0]
	phyla = splitted[2]

	tax_dict[pfam] = phyla.split(",")
