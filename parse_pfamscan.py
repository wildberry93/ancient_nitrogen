def get_pfam_results(result_file, evalue_threshold=1e-2):
        """
        Parse hmmscan tabular results. Return a dictionary mapping ids to Pfam accessions.
        """

        results = {}
        for line in open(result_file):
                if line[:1] == '#':
                        continue
                col = line.split()
                if len(col) < 14:
                        continue
 
                gi = col[0]
                accession = col[6]
		start = col[1]
		end = col[2]

                if '|' in gi:
                        gi = gi.split('|')[1]
                if '_' in gi:
                        gi = gi.split('_')[0]

		#print(gi)
		#print(accession)
                if accession not in results:
                        results[accession] = []
                results[accession].append(gi + " " +start + " " + end)
                #if evalue <= evalue_threshold:
        return results		
if __name__ == '__main__':
	infile = "all_proteomes_small.pfamscan"
	results = get_pfam_results(infile)
	file = open("pfam2accession", "w")
	for pfam in results:
		file.write(pfam + "\t" + ",".join(results[pfam])+"\n")
