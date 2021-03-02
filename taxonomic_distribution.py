infile = open("pfam2accession", "r").read().splitlines()

taxonomy = open("/mnt/sda1/proteomes/obsolete/selected_ref_proteomes_bact_arch_fulltaxonomy_oneperfam", "r").read().splitlines()
p2p = open("/mnt/sda1/proteomes/obsolete/p2p2taxon_prokaryotes", "r").read().splitlines()
clan_dict = {}
clan_file = open("clan_membership.txt","r").read().splitlines()
pfam_names_file = open("Pfam-A.clans.tsv","r").read().splitlines()
pfam_names_dict = {}

outgroup = ["Aquificae","Thermotogae","Fusobacteria"]

def is_LUBA_perm(b_phyla, a_phyla):
	if a_phyla: 
		return False

	if len(b_phyla) >= 13 and set(outgroup).intersection(b_phyla):
		return True 

def is_LUOA(b_phyla, a_phyla):
	"""
	Calculate LUOA
	"""
	if a_phyla:
		return False

	marine = ["Spirochaetes","Chlamydiae","Planctomycetes","Bacteroidetes",
                "Chlorobi","Ignavibacteriae","Rhodothermaeota","Proteobacteria","Acidobacteria"]

	terrestrial = ["Tenericutes","Firmicutes","DeinococcusThermus","Cyanobacteria",
                "Actinobacteria","Chloroflexi"]

	num_marine = []
	num_terrestrial = []

	for b_p in b_phyla:
		if b_p in marine:
			num_marine.append(b_p)
		elif b_p in terrestrial:
			num_terrestrial.append(b_p)

	if len(num_marine) >= 2 and len(num_terrestrial) >= 2:
		return True

	else:
		return False

def date(b_phyla, a_phyla):
	"""
	Assign as LUCA, LUOA, LUBA, LUAA and less stringent of those
	"""
	total_bact = 18
	total_arch = 8

	if len(b_phyla) == total_bact and len(a_phyla) == total_arch:
		return "LUCA"

	if len(b_phyla) == total_bact and len(a_phyla) == 0:
		return "LUBA"

	if len(b_phyla) == 0 and len(a_phyla) == 8:
		return "LUAA"

	### lenient criteria

	if len(b_phyla) >= 13 and len(a_phyla) >= 6:
		return "LUCA_perm"

	if is_LUBA_perm(b_phyla, a_phyla): #len(b_phyla) >= 13 and len(a_phyla) == 0:
		return "LUBA_perm"

	if len(b_phyla) == 0 and len(a_phyla) >= 6:
		return "LUAA_perm"

	if is_LUOA(b_phyla, a_phyla): 
                return "LUOA"

	else:
		return "Unassigned"

for line in pfam_names_file:
	splitted = line.split("\t")
	name = splitted[3]
	pfid = splitted[0]

	pfam_names_dict[name] = pfid

tax_dict = {}

for fam in clan_file:
	splitted = fam.split("\t")
	clan = splitted[0]
	family = splitted[1]

	clan_dict[family] = clan

phylum2kingdom = {}
for line in taxonomy[1:]:
        splitted = line.split("\t")
        species = splitted[2]
        phylum = splitted[1]
	kingdom = splitted[-2]
        tax_dict[species] = (phylum, kingdom)
	phylum2kingdom[phylum] = kingdom

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
outfile = open("taxonomy_distribution_clan.tab", "w")
for pfam in pfam_dict:
	num_phyla = str(len(pfam_dict[pfam]["phyla"]))
	num_kingdoms = str(len(pfam_dict[pfam]["kingdoms"]))

	pfid = pfam_names_dict[pfam]

	if pfid not in clan_dict:
		clan = "no clan"
	else:
		clan  = clan_dict[pfid]

	phyla = pfam_dict[pfam]["phyla"]

	num_bact = 0
	num_arch = 0
	bacteria_phyla = []
	archaea_phyla = []

	for phylum in phyla:
		if phylum2kingdom[phylum] == "Archaea":
			num_arch += 1
			archaea_phyla.append(phylum)
		if phylum2kingdom[phylum] == "Bacteria":
			num_bact += 1
			bacteria_phyla.append(phylum)
	pfam_date = date(bacteria_phyla, archaea_phyla) 
	outfile.write(pfam+ "\t" + pfid + "\t" + clan + "\t" + pfam_date + "\t" + num_phyla + "\t" + str(num_bact) + "\t" + str(num_arch)  + "\t" + ",".join(pfam_dict[pfam]["phyla"]) + "\t" + num_kingdoms + "\t" + ",".join(pfam_dict[pfam]["kingdoms"]) +"\n")
