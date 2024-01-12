##############################################################################################
##
## this script prepares downloadable data for PON-mtRNA predictions
## And also creates a dictionary of position and amino aicds for the mmr proteins
## Abhishek Niroula 2015-09-13
##
###############################################################################################
#!/usr/bin/python

import os
import pickle


wd="/home/structure/pontrna/"

myl=[]

#writefile=open(wd+"ponmmr_site/PON-MMR2_predictions.txt","w")
f=open(wd+"pon-mtrna_site/mysql_data.txt")
for line in f:
	if line.startswith("tRNA"):
		continue
	splt=line[:-1].split("\t")
	tRNA=splt[0]
        pos=splt[1]
	ori=splt[2]
	if pos not in myl:
		myl.append(pos)

#writefile.close()
f.close()

f=open(wd+"pon-mtrna_site/trna_position_list.pkl","wb")
pickle.dump(myl,f)
f.close()

