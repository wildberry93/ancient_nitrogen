infile = open("tRNA-synt_1c_renamed_trimmed.fasta.phy","r").read().splitlines()
outfile = open("tRNA-synt_1c_renamed_trimmed_2.fasta.phy","w")
outfile.write(infile[0] + "\n")

for line in infile[1:]:
	splitted = line.split()
	new_name = splitted[0].replace("_303_bp","")
	outfile.write(new_name + "                    " + splitted[1] + "\n")
