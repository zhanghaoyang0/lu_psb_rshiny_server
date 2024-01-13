#!/usr/bin/python
# the old script use sql, haoyang rewrite this script without sql at 202401
import sys
import pandas as pd

jobid=sys.argv[1]
# jobid='pmmr24932722'

## Remove hash
def removeHash(text):
    if " " in text:
        text=text.replace(" ","")
    if "#" not in text:
        return text
    splittext=text.split("#")
    return " ".join(splittext[0].split())


# pre-calculated prediction, from 
# db=mdb.connect("localhost","abhishek","mysqlpassword","ponmtRNA") in structure
db = pd.read_csv("data/pon_mmr2.txt", delimiter="\t") 

#################################################################################
# query
#################################################################################
## Collecting the variations
protGene={"P40692":"MLH1","P43246":"MSH2","P52701":"MSH6","P54278":"PMS2"}
fil=open('temp/'+jobid+".query","r")
writefile=open('temp/'+jobid+"_results.txt","w")
writefile.write("Gene\tUniProtKB Accession\tAmino acid substitution\tPosition\tOriginal amino acid\tNew amino acid\tProbability of pathogenicity\tClassification\tInSiGHT class\n")

readf=fil.read()
fil.close()
if readf[-1]=="\n":
    readf=readf[:-1]

entries=readf[1:].replace("\r","").split("\n>")


for entry in entries:
    t=""
    query={}
    text=""
    splitentry=entry.split("\n")
    name=removeHash(splitentry[0])
    for var in splitentry[1:]:
        if var=="" or var=="\n" or var=="\r\n":
            continue
        variation=removeHash(var)
        ref=variation[0]
        alt=variation[-1]
        pos=variation[1:-1]
        if name not in query.keys():
            query[name]=[variation]
        else:
            query[name].append(variation)
    keys=query.keys()
    # q="""select * from predictionTable where """
    dat = pd.DataFrame()
    for key in keys:
        variations=query[key]
        if key in protGene.values():
            # q+="""(geneName='%s' AND (""" % key
            sub = db[(db["geneName"] == key) & db["AAS"].isin(variations)]
        else:
            # q+="""(UniProtKB='%s' AND (""" % key
            sub = db[(db["UniProtKB"] == key) & db["AAS"].isin(variations)]
        dat = pd.concat([dat, sub], ignore_index=True)
        # for v in variations:
        #     q+="""AAS='%s' OR """ % v
        # q=q[:-4]+""")) OR """
    # db=mdb.connect("localhost","abhishek","mysqlpassword","ponmmr")
    # cur=db.cursor()
    # cur.execute(q)
    # dat=cur.fetchall()
    dat = dat.values.tolist()
    if not dat or dat=="None":
        t="Send email!!"
    else:
        for i in dat:
            for j in i:
                text+=str(j)+"\t"
            text=text[:-1]+"\n"
        writefile.write(text)

writefile.close()



