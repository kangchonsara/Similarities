import os

infName = os.path.normpath("../raw/HongKong2014_AA.fas")
inf = open(infName, "r")
outfName = os.path.normpath("../data/HongKong4801_AA.fas")
outf = open(outfName, "w")

variants = []

for line in inf:
	if line.find(">") >= 0:
		each = line.split("_|_")
		name = each[1]
		if name == "A/Hong_Kong/4801/2014":
			seq = each[0] + "|" + each[1] + "|" + each[2] + "\n"
		else:
			seq = ''
	else:
		if line.find("?") >= 0:
			continue
			
		if seq != '':
			if line in variants:
				continue
			else:
				variants.append(line)
				seq += line
			outf.write(seq)
inf.close()
outf.close()
			