##############################################################################################
##
## this script inserts the PON-mtRNA predictions in the mysql table
## Abhishek Niroula 2015-09-13
##
###############################################################################################
#!/usr/bin/python

import MySQLdb as mdb
import os



wd="/home/structure/pontrna/"


db=mdb.connect("localhost","abhishek","mysqlpassword","ponmtRNA")
cur=db.cursor()


f=open(wd+"pon-mtrna_site/mysql_data.txt")
for line in f:
	if line.startswith("tRNA"):
		continue
	splt=line[:-1].split("\t")
	trna=splt[0]
        pos=splt[1]
	ori=splt[2]
	new=splt[3]
	prob=splt[5]
	clas=splt[4]
	cur.execute("""INSERT INTO predictionTable (mtRNA,position,reference,altered,class,probability) VALUES(%s, %s, %s, %s, %s, %s)""",(trna,int(pos),ori,new,clas,float(prob)))
	db.commit()

db.close()
f.close()

