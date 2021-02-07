infile = open("pfam2accession", "r").read().splitlines()

taxonomy = open("/mnt/sda1/proteomes/obsolete/selected_ref_proteomes_bact_arch_fulltaxonomy_oneperfam", "r").read().splitlines()
p2p = open("/mnt/sda1/proteomes/obsolete/p2p2taxon_prokaryotes", "r").read().splitlines()

tax_dict = {}
for line in taxonomy[1:]:
        splitted = line.split("\t")
        species = splitted[2]
        phylum = splitted[1]
	kingdom = splitted[-2]
        tax_dict[species] = (phylum, kingdom)


p2p_dict = {}

for line in p2p:
        splitted = line.split("\t")

        protein = splitted[1]
        species = splitted[0]

        p2p_dict[protein] = species

pfam_dict = {}
for pfam in infile:
	splitted = pfam.split("\t")
	sequences = splitted[1].split(",")
	pfam = splitted[0]
	pfam_dict[pfam] = {"phyla":[], "kingdoms":[]}
	
	for sequence in sequences:
		sequence = sequence.split(" ")[0]
		species = p2p_dict[sequence]
		taxonomy = tax_dict[species]

		phylum = taxonomy[0]
		kingdom = taxonomy[1]

		if phylum not in pfam_dict[pfam]["phyla"]:

			pfam_dict[pfam]["phyla"].append(phylum)


		if kingdom not in pfam_dict[pfam]["kingdoms"]:
			pfam_dict[pfam]["kingdoms"].append(kingdom)

#print(pfam_dict)
outfile = open("taxonomy_distribution.tab", "w")
for pfam in pfam_dict:
	num_phyla = str(len(pfam_dict[pfam]["phyla"]))
	num_kingdoms = str(len(pfam_dict[pfam]["kingdoms"]))

	outfile.write(pfam+ "\t" + num_phyla + "\t" + ",".join(pfam_dict[pfam]["phyla"]) + "\t" + num_kingdoms + "\t" + ",".join(pfam_dict[pfam]["kingdoms"]) +"\n")
