#Calculate amino acid similarity of sequences at each year to 2016-2017 vaccine sequence at site A and site B
import os

beginY = 1968
endY = 2010


def read_epitope_AB(epifName):
	epiA = []
	epiB = []
	epif = open(epifName, "r")
	for line in epif:
		each = line.split("\n")[0].split("\t")
		if each[1] == 'a':
			epiA.append(int(each[0]))
		elif each[1] == 'b':
			epiB.append(int(each[0]))
	epif.close()
	
	return epiA, epiB
	
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
			
	
#When there are 2 or more identical samples (same virus name and same sequence), leave only 1.
#Remove samples with deletion or ambiguity at epitope site.
def remove_dupl_ambig(sequences, epiA, epiB):
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
					if s+1 in epiA or s+1 in epiB:
						if sequence[s] == "-" or sequence[s] == "?" or sequence[s] == "*":
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

def calc_epi_similarity(byYear, vacSeqs, epiB):
	B_similarities = []
	
	for yearSeq in byYear:
		mean_simil = 0
		
		for seq in yearSeq:
			similarity = 0
			for s in epiB:
				if seq[s-1] == vacSeqs[0][s-1]:
					similarity += 1
			
			similarity = 1.0*similarity/len(epiB) #per site
			mean_simil += similarity
			
		mean_simil = mean_simil/len(yearSeq) #per seq
		B_similarities.append(mean_simil)
	
	return B_similarities

infName = os.path.normpath("../data/H3_6810AA.fasta")
shifName = os.path.normpath("../data/shih_epitope.txt")
vacfName = os.path.normpath("../data/HongKong4801_AA.fas")

byYear_r, USbyYear_r = seqs_byYear(infName)

epiA, epiB = read_epitope_AB(shifName)
print (epiB)

byYear = remove_dupl_ambig(byYear_r, epiA, epiB)
USbyYear = remove_dupl_ambig(USbyYear_r, epiA, epiB)

vacSeqs = read_vaccineSeq(vacfName)

B_similarities = calc_epi_similarity(byYear, vacSeqs, epiB)
for y in range(len(B_similarities)):
	print (y+beginY, B_similarities[y])
	
	
	
	
	
	
	