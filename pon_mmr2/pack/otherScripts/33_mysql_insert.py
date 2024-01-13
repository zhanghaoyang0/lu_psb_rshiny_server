##############################################################################################
##
## this script inserts the mmr predictions in the mysql table
## Abhishek Niroula 2015-04-13
##
###############################################################################################
#!/usr/bin/python

import MySQLdb as mdb
import os



wd="/home/structure/ponmmr/"


db=mdb.connect("localhost","abhishek","mysqlpassword","ponmmr")
cur=db.cursor()


f=open(wd+"ponmmr_site/mysql_data.txt")
for line in f:
	if line.startswith("gene"):
		continue
	splt=line[:-1].split("\t")
	gene=splt[0]
        if gene not in ["MLH1","MSH2","MSH6","PMS2"]:
            continue
	uniprot=splt[1]
	aas=splt[2]
	pos=splt[3]
	ori=splt[4]
	new=splt[5]
	prob=splt[6]
	clas=splt[7]
	if len(splt)>8 and splt[8]!="":
		insight="Class "+splt[8]
	else:
		insight=""
	cur.execute("""INSERT INTO predictionTable (geneName, UniProtKB,AAS,Position,OriginalAA,NewAA,probabilityPathogenicity,Classification,InSiGHT_class) VALUES(%s, %s, %s, %s, %s, %s, %s, %s, %s)""",(gene,uniprot,aas,int(pos),ori,new,float(prob),clas,insight))
	db.commit()

db.close()
f.close()

