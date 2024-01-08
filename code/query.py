#!/usr/bin/python
# /home/ha0214zh/website/pack/ponp/PONP_site/queryFiles/trna1052227
# the old script use sql, haoyang rewrite this script without sql at 202401

import sys
import os.path
import pickle
import pandas as pd
import re

jobid=sys.argv[1]
queryFilesPath='../temp/'
dbPath='../data/'

## test
# jobid='3ed8495fe587d64df9f530e8e08a77de'

## Remove hash
def removeHash(text):
    if " " in text:
        text=text.replace(" ","")
    if "#" not in text:
        return text
    splittext=text.split("#")
    return " ".join(splittext[0].split())

## this function checks if the text is in correct format
nucs={"A":"T","C":"G","G":"C","T":"A"}
nkeys=nucs.keys()

## position list of the human mt-tRNAs
f=open(dbPath+"trna_position_list.pkl","rb")
poslist=pickle.load(f)
f.close()

## the following values comes from the training dataset.
## these are present in the local file 
## /media/abhishek/DATA/tRNA/data/analysis/unknown_cases_classified.xlsx
##

segPositive=2.25
segNegative=0.1668
biochemPositive=8.07495
biochemNegative=0.1336
histochemPositive=5.9667
histochemNegative=0.0282

# pre-calculated prediction, from 
# db=mdb.connect("localhost","abhishek","mysqlpassword","ponmtRNA") in structure
db = pd.read_csv(dbPath+"pon_mt_trna.txt", delimiter="\t") 
#################################################################################
# make proc_evid
# /home/ha0214zh/website/pack/ponp/PONP_site/trna_scripts.py
#################################################################################

def trnaformatCheck(jobid):
    position=[]
    refer=[]
    alter=[]
    seg=[]
    biochem=[]
    histochem=[]
    fil=open(queryFilesPath+jobid+".query","r")
    readf=removeHash(fil.read()).replace("\n\n","\n")
    fil.close()
    if len(readf)<2:
        return "There are no variations in the submission."
    if readf[-1]=="\n":
        readf=readf[:-1]
    entries=readf.replace("\r","").split("\n")
    for entry in entries:
        splitentry=entry.split(",")
        if len(splitentry)<3:
            return "The submitted variation (%s) is not in the correct format." %entry
        pos=splitentry[0]
        ref=splitentry[1]
        alt=splitentry[2]
        if not pos.isdigit():
            continue
        if ref not in nkeys or alt not in nkeys:
            continue
        position.append(pos)
        refer.append(ref)
        alter.append(alt)
        # if not pos.isdigit():
        #     return "The position in the submitted variation (%s) is not acceptable." %entry
        # else:
        #     position.append(pos)
        # if ref not in nkeys or alt not in nkeys:
        #     return "The nucleotide in the submitted variation (%s) could not be recognized." %entry
        # else:
        #     refer.append(ref)
        #     alter.append(alt)
        if len(splitentry)>3:
            seg.append(splitentry[3])
            if len(splitentry)>4:
                biochem.append(splitentry[4])
                if len(splitentry)>5:
                    histochem.append(splitentry[5])
                else:
                    histochem.append("a")
            else:
                biochem.append("a")
                histochem.append("a")
        else:
            seg.append("a")
            biochem.append("a")
            histochem.append("a")
    evidText=""
    for k in range(0,len(position)):
        t=position[k]+"\t"+refer[k]+"\t"+alter[k]
        evidText+=t+"\t"+seg[k]+"\t"+biochem[k]+"\t"+histochem[k]+"\n"
    writefile=open(queryFilesPath+jobid+"_proc_evid.txt","w")
    writefile.write(evidText)
    writefile.close()
    return "Right evid"


trnaformatCheck(jobid)

#################################################################################
# query
#################################################################################
## Now look for the presence of evidence and classification
evidence_based=""
evid_f=queryFilesPath+jobid+"_proc_evid.txt"
fil=open(evid_f,"r") # input
evid_varlist=[]

if os.path.exists(evid_f):
    fil=open(evid_f,"r")
    for entry in fil:
        splitentry=entry[:-1].split("\t")
        ref=splitentry[1]
        alt=splitentry[2]
        pos=splitentry[0]
        seg_evid=splitentry[3]
        biochem_evid=splitentry[4]
        histochem_evid=splitentry[5]
        if pos not in poslist:
            print("We could not find the position for variation %s in our database. Please check that the variation position is within mitochondrial tRNA. It the position is correct, please contact us at proteinbioinformatik@med.lu.se.\n" % entry.replace("\n",""))
            continue
        # q="""select * from predictionTable where """
        # q+="""(position='%s' AND (""" % pos
        # q+="""(reference='%s' AND (""" % ref
        # q+="""altered='%s'))))""" % alt
        # db=mdb.connect("localhost","abhishek","mysqlpassword","ponmtRNA")
        # cur=db.cursor()
        # cur.execute(q)
        # dat=cur.fetchall()
        dat = db[(db["position"] == int(pos)) & (db["reference"] == ref) & (db["altered"] == alt)] # rewrite by haoyang
        dat = dat.values.tolist() # rewrite by haoyang
        if not dat or dat=="None":
            # q="""select * from predictionTable where """
            # q+="""(position='%s' AND (""" % pos
            # q+="""(reference='%s' AND (""" % nucs[ref]
            # q+="""altered='%s'))))""" % nucs[alt]
            # db=mdb.connect("localhost","abhishek","mysqlpassword","ponmtRNA")
            # cur=db.cursor()
            # cur.execute(q)
            # dat=cur.fetchall()
            dat = db[(db["position"] == int(pos)) & (db["reference"] == nucs[ref]) & (db["altered"] == nucs[alt])]
            dat = dat.values.tolist()
            if not dat or dat=="None":
                error.write("The variation %s could not be mapped to our reference sequence" % entry)
                continue
        if len(dat)<1:
            error.write("The vaiation %s could not be mapped to our reference sequence." % entry)
        for i in dat:
            text=[]
            for j in i:
                text.append(str(j))
            prior_prob=float(text[5])
            if seg_evid=="1":
                if biochem_evid=="1":
                    if histochem_evid=="1":
                        logg_ratio=segPositive*biochemPositive*histochemPositive
                    elif histochem_evid=="0":
                        logg_ratio=segPositive*biochemPositive*histochemNegative
                    else:
                        logg_ratio=segPositive*biochemPositive
                elif biochem_evid=="0":
                    if histochem_evid=="1":
                        logg_ratio=segPositive*biochemNegative*histochemPositive
                    elif histochem_evid=="0":
                        logg_ratio=segPositive*biochemNegative*histochemNegative
                    else:
                        logg_ratio=segPositive*biochemNegative
                else:
                    if histochem_evid=="1":
                        logg_ratio=segPositive*histochemPositive
                    elif histochem_evid=="0":
                        logg_ratio=segPositive*histochemNegative
                    else:
                        logg_ratio=segPositive
            elif seg_evid=="0":
                if biochem_evid=="1":
                    if histochem_evid=="1":
                        logg_ratio=segNegative*biochemPositive*histochemPositive
                    elif histochem_evid=="0":
                        logg_ratio=segNegative*biochemPositive*histochemNegative
                    else:
                        logg_ratio=segNegative*biochemPositive
                elif biochem_evid=="0":
                    if histochem_evid=="1":
                        logg_ratio=segNegative*biochemNegative*histochemPositive
                    elif histochem_evid=="0":
                        logg_ratio=segNegative*biochemNegative*histochemNegative
                    else:
                        logg_ratio=segNegative*biochemNegative
                else:
                    if histochem_evid=="1":
                        logg_ratio=segNegative*histochemPositive
                    elif histochem_evid=="0":
                        logg_ratio=segNegative*histochemNegative
                    else:
                        logg_ratio=segNegative
            else:
                if biochem_evid=="1":
                    if histochem_evid=="1":
                        logg_ratio=biochemPositive*histochemPositive
                    elif histochem_evid=="0":
                        logg_ratio=biochemPositive*histochemNegative
                    else:
                        logg_ratio=biochemPositive
                elif biochem_evid=="0":
                    if histochem_evid=="1":
                        logg_ratio=biochemNegative*histochemPositive
                    elif histochem_evid=="0":
                        logg_ratio=biochemNegative*histochemNegative
                    else:
                        logg_ratio=biochemNegative
                else:
                    if histochem_evid=="1":
                        logg_ratio=histochemPositive
                    elif histochem_evid=="0":
                        logg_ratio=histochemNegative
                    else:
                        logg_ratio="lol"
            if logg_ratio!="lol":
                p_odds=logg_ratio*(prior_prob/(1-prior_prob))
                post_prob=round(p_odds/(1+p_odds),4)
                if post_prob>0.99:
                    clas="Pathogenic"
                elif post_prob>=0.95:
                    clas="Likely pathogenic"
                elif post_prob<0.001:
                    clas="Neutral"
                elif post_prob<0.05:
                    clas="Likely neutral"
                else:
                    clas="VUS"
            else:
                post_prob="NA"
                clas="NA"
            text.append((seg_evid+", "+biochem_evid+", "+histochem_evid).replace("a","NA"))
            text.append(post_prob)
            text.append(clas)
            evid_varlist.append(text)
    fil.close()
else:
    evidence_based="Not applicable"


nVar=0

## Now the predictions to the files
if len(evid_varlist)>0:
    writefile=open(queryFilesPath+jobid+"_results.txt","w")
    writefile.write("mt-tRNA\tVariation\tML_probability_of_pathogenicity\tEvidence\tPosterior_probability_of_pathogenicity\tClassification\n")
    if evidence_based!="Not applicable":
        nVar=len(evid_varlist)
        for entry in evid_varlist:
            mtRNA=entry[0]
            variation=entry[1]+","+entry[2]+","+entry[3]
            prior_prob=round(float(entry[5]),2)
            evidence_av=entry[6]
            if entry[8]=="NA":
                post_prob=""
                classification=entry[4]
            else:
                post_prob=str(round(entry[7],2))
                classification=entry[8]
            writefile.write(mtRNA+"\t"+variation+"\t"+str(prior_prob)+"\t"+evidence_av+"\t"+post_prob+"\t"+classification+"\n")
