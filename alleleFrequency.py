#Calculate allele frequencies at 160 from 1968 to 2010

inf = open("6810AA.fasta", "r")

beginY = 1968
endY = 2010

#When there are 2 or more identical samples (same virus name), leave only 1.
#Remove samples with deletion or ambiguity at site 160.
def remove_dupl_ambig(sequences):
	years = []
	for year in sequences:
		oneYear = []
		for seq in year:
			if seq in oneYear:
				continue
			else:
				site = seq.split(" ")[1][159]
				if site == "-" or site == "?":
					continue
				oneYear.append(seq)
		years.append(oneYear)
	return years

#Calculate allele frequency separately for global samples and US samples
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

byYear = remove_dupl_ambig(byYear_r)
USbyYear = remove_dupl_ambig(USbyYear_r)

frequency = []
useUS = []
for i in range(len(byYear)):
	alleles = [0,0,0,0,0]
	if len(USbyYear[i]) >= 70:
		year = USbyYear[i]
		useUS.append(1)
	else:
		year = byYear[i]
		useUS.append(0)
		
	for virus in year:
		seq = virus.split(" ")[1]
		if seq[159] == "T": 
			alleles[0] += 1
		elif seq[159] == "K":
			alleles[1] += 1
		elif seq[159] == "A":
			alleles[2] += 1
		elif seq[159] == "R":
			alleles[3] += 1
		else:
			alleles[4] += 1
			
	s = sum(alleles)
	freq = []
	for a in alleles:
		f = round( a/s, 3)
		freq.append(f)
	frequency.append(freq)
	

y = 1968
for idx in range(len(frequency)):
	i = frequency[idx]
	print (y, i[0], i[1], i[2], i[3], i[4], useUS[idx])
	y += 1
	
	
	
	
	
	
	
	
	
	
	
	
	
	