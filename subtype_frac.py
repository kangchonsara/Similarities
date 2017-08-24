#Calculate relative influenza subtype frequency H3 vs. H1+B from 1968 to 2010
beginY = 1968
endY = 2010

#Put samples in a list by year
def byYear(infName):
	inf = open(infName, "r")
	years = [[] for i in range(beginY, endY+1)]
	USyears = [[] for i in range(beginY, endY+1)]
	for line in inf:
		if line.find(">") >= 0:
			each = line.split("|")
			year = int(each[2].split("/")[0])
			#If there are 2 or more identical samples (same virus name), leave only 1
			if each[1] in years[year-beginY]:
				continue
			else:
				years[year-beginY].append(each[1])
			
			country = each[3]
			if country == "USA":
				USyears[year-beginY].append(each[1])
	inf.close()
	return years, USyears
	
	
h3n2f = "H3_6810AA.fasta"
h1n1f = "H1.fasta"
bf = "B.fasta"
	
h3n2Y, USh3n2Y = byYear(h3n2f)
h1n1Y, USh1n1Y = byYear(h1n1f)
bY, USbY = byYear(bf)

#Get frequencies of influenza subtype separately for glboal samples and US samples
for i in range(endY-beginY+1):
	UStotal = len(USh3n2Y[i])+len(USh1n1Y[i])+len(USbY[i])
	if UStotal >= 70:	
		h3n2 = len(USh3n2Y[i])
		h1n1_b = len(USh1n1Y[i]) + len(USbY[i])
		total = UStotal
		USonly = 1
	else:
		h3n2 = len(h3n2Y[i])
		h1n1_b = len(h1n1Y[i]) + len(bY[i])
		total = h3n2 + h1n1_b
		USonly = 0
		
	frac_h3n2 = round( h3n2/total, 3)
	frac_h1n1_b = round( h1n1_b/total, 3)
	print (beginY+i, frac_h3n2, frac_h1n1_b, USonly)

	
	
	
	
	
	
	
	
	