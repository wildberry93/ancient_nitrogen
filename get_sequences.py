import os

infile = open("pfam2accession").read().splitlines()

for pfam in infile:
	pfam_id = pfam.split("\t")[0]
	sequences =  pfam.split("\t")[1]

	temp_file = open("temp", "w")

	for seq in sequences.split(","):
		#print(seq)
		seq = seq.split()[0] + "\t" + seq.split()[1] + "-" + seq.split()[2]
		temp_file.write(seq + "\n")
	
	temp_file.close()


	os.system("blastdbcmd -db /mnt/sda1/proteomes/dating/proteomes_perorder/all_proteomes/all_proteomes.fasta -dbtype prot -entry_batch %s -out %s" % ("temp", pfam_id+".fasta"))


	os.remove("temp")
