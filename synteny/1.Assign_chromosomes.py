import sys, os, glob


flist =glob.glob("input/chr.assign/*.csv")
for fNAME in flist:
	sID = os.path.basename(fNAME).split(".")[0]

	# Reading chr.assignment
	nScaffoldName_nChrName_dic = {}
	fpin = open(fNAME,'r')
	for line in fpin:
		part = line.strip('\n').split(',')
		nScaffoldName = part[0]
		try:
			nChrName = "chr"+part[1]
		except:
			print(fNAME,part)
			sys.exit()
		nAssignQuality = part[2][0]
		if nAssignQuality.lower() == "y":
			nScaffoldName_nChrName_dic.setdefault(nScaffoldName,nChrName)

	# Updating chr.assignment on BUSCO results
	fNAME = "input/BUSCO/BUSCO_"+sID+".txt"
	fpin = open(fNAME,'r')
	tmpline_list = []
	for line in fpin:
		tmpline = line
		part = line.strip('\n').split("\t")
		try:
			tmpline = part[0]+'\t'+part[1]+'\t'+nScaffoldName_nChrName_dic[part[2]]+'\t'+'\t'.join(part[3:])+'\n'
		except:
			pass
		tmpline_list.append(tmpline)
	fpin.close()
	fpout = open(fNAME,'w')
	for tmpline in tmpline_list:
		if tmpline.startswith('#'):
			fpout.write(tmpline)
		else:
			line = tmpline.strip().split()
			len_l = len(line)
			if len_l > 2:
				seq = line[2].replace('chr', '')
				try:
					if seq == 'Y' or seq == 'X' or isinstance(int(seq), int):
						fpout.write(tmpline)
				except ValueError:
					pass
			else:
				fpout.write(tmpline)
	fpout.close()

	# Updating chr_assignment on chrsize
	fNAME = "input/chrsize/chrsize_"+sID+".txt"
	fpin = open(fNAME,'r')
	tmpline_list = []
	for line in fpin:
		tmpline = line
		part = line.strip('\n').split("\t")
		try:
			tmpline = nScaffoldName_nChrName_dic[part[0]]+'\t'+part[1]+'\n'
		except:
			pass
		tmpline_list.append(tmpline)
	fpin.close()
	fpout = open(fNAME,'w')
	for tmpline in tmpline_list:
		fpout.write(tmpline)
	fpout.close()

