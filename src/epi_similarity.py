#Calculate amino acid similarity of sequences at each year to 2016-2017 vaccine sequence at site A and site B
import os

beginY = 1968
endY = 2010
sites = ['a', 'b', 'c', 'd']

def read_epitope_AB(epifName):

	epitopes = [[] for i in range(len(sites)) ]
	
	epif = open(epifName, "r")
	for line in epif:
		each = line.split("\n")[0].split("\t")
		idx = sites.index(each[1])
		epitopes[idx].append(int(each[0]))
			
	epif.close()
	
	return epitopes
	
def read_vaccineSeq(vacfName):
	vacSeqs = []
	vacf = open(vacfName, "r")
	for line in vacf:
		if line.find(">") >= 0:
			seq = line.split("\n")[0] + " "
		else:
			seq += line.split("\n")[0]
			vacSeqs.append(seq)
			
	vacf.close()
	return vacSeqs
	
def vacSeq_dif(vacSeqs, epitope):
	vs_dif_all = 0
	num_comp = 0
	for v1 in range(len(vacSeqs)-1):
		for v2 in range(v1+1, len(vacSeqs)):
			vac1 = vacSeqs[v1].split(" ")[1]
			vac2 = vacSeqs[v2].split(" ")[1]
			
			num_comp += 1
			vs_dif = 0
			for s in epitope:
				if vac1[s-1] != vac2[s-1]:
					print (s, v1, v2, vac1[s-1], vac2[s-2])
					vs_dif += 1

			vs_dif_all += vs_dif
	vs_dif_avg = 1.0*vs_dif_all/(len(epitope)*num_comp)
	
	return vs_dif_avg
	
			
	
#When there are 2 or more identical samples (same virus name and same sequence), leave only 1.
#Remove samples with deletion or ambiguity at epitope site.
def remove_dupl_ambig(sequences, epitope_mixed):
	years = []

	for year in sequences:

		oneYear = []
		for seq in year:
			if seq in oneYear:
				continue
			else:
				sequence = seq.split(" ")[1]
				
				remove = 0
				for s in range(len(sequence)):
					if s+1 in epitope_mixed:
						if sequence[s+1] == "-" or sequence[s+1] == "?" or sequence[s+1] == "*":
							remove = 1
							break
							
				if remove == 1:
					continue

				oneYear.append(seq)
				
		years.append(oneYear)
	return years

def seqs_byYear(infName):
	inf = open(infName, "r")

	byYear_r = [[] for i in range(beginY, endY+1)]
	USbyYear_r = [[] for i in range(beginY, endY+1)]
	for line in inf:
		if line.find(">") >= 0:
			year = int(line.split("|")[2].split("/")[0])
			country = line.split("|")[3]
			seq = line.split("|")[1] + " "
		else:
			seq += line.split("\n")[0]
			byYear_r[year-beginY].append(seq)
			if country == "USA":
				USbyYear_r[year-beginY].append(seq)
	inf.close()
	
	return byYear_r, USbyYear_r

def calc_epi_similarity(byYear, vacSeqs, epitope):
	epi_similarities = []
	
	for yearSeq in byYear:
		mean_simil = 0
		
		for seq in yearSeq:
			seq = seq.split(" ")[1]
			similarity = 0
			for vacSeq in vacSeqs:
				vacSeq = vacSeq.split(" ")[1]
				for s in epitope:
					if seq[s-1] == vacSeq[s-1]:
						similarity += 1
			
			similarity = 1.0*similarity/(len(epitope)*len(vacSeqs)) #per site
			mean_simil += similarity
			
		mean_simil = mean_simil/len(yearSeq) #per seq
		epi_similarities.append(mean_simil)
	
	return epi_similarities

infName = os.path.normpath("../data/H3_6810AA.fasta")
shifName = os.path.normpath("../data/shih_epitope.txt")
vacfName = os.path.normpath("../data/HongKong4801_AA.fas")
similfName = os.path.normpath("../result/similarities_")

byYear_r, USbyYear_r = seqs_byYear(infName)

epitopes = read_epitope_AB(shifName)
epitope_mixed = []
for epitope in epitopes[:2]:
	epitope_mixed += epitope
	
byYear = remove_dupl_ambig(byYear_r, epitope_mixed)
USbyYear = remove_dupl_ambig(USbyYear_r, epitope_mixed)

vacSeqs = read_vaccineSeq(vacfName)
#for epitope in epitopes:
#	vac_dif = vacSeq_dif(vacSeqs, epitope)
#	print (vac_dif)
vacSeqs = vacSeqs[0:1]
	
for e in range(len(epitopes[:2])):
	epif = open(similfName+sites[e]+".csv", "w")
	epif.write("year,similarity\n")
	similarities = calc_epi_similarity(byYear, vacSeqs, epitopes[e])
	for y in range(len(similarities)):
		epif.write(str(y+beginY)+","+str(similarities[y])+"\n")
	epif.close()

	
	
	
	
	
	
	
	
	
	