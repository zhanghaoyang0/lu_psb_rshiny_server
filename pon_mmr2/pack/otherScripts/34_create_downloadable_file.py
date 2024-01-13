##############################################################################################
##
## this script prepares downloadable data for PON-MMR2 predictions
## And also creates a dictionary of position and amino aicds for the mmr proteins
## Abhishek Niroula 2015-04-13
##
###############################################################################################
#!/usr/bin/python

import os
import pickle


wd="/home/structure/ponmmr/"

myd={}

writefile=open(wd+"ponmmr_site/PON-MMR2_predictions.txt","w")
f=open(wd+"ponmmr_site/mysql_data.txt")
for line in f:
	if line.startswith("gene"):
		continue
	splt=line[:-1].split("\t")
	gene=splt[0]
        if gene not in ["MLH1","MSH2","MSH6","PMS2"]:
            continue
	writefile.write(line)
	pos=splt[3]
	ori=splt[4]
	if gene+"_"+pos not in myd.keys():
		myd[gene+"_"+pos]=ori

writefile.close()
f.close()

f=open(wd+"ponmmr_site/mmr_position_dict.pkl","wb")
pickle.dump(myd,f)
f.close()

