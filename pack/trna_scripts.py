#!/usr/bin/env python
from Bio import SeqIO
import re
import pickle

## this function whether the submitted email is in a valid format
from django.core.validators import email_re
def is_valid_email(email):
    if email_re.match(email):
        return True
    return False



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

queryFilesPath="/home/structure/ponp/PONP_site/queryFiles/"
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
            return "The position in the submitted variation (%s) is not acceptable." %entry
        else:
            position.append(pos)
        if ref not in nkeys or alt not in nkeys:
            return "The nucleotide in the submitted variation (%s) could not be recognized." %entry
        else:
            refer.append(ref)
            alter.append(alt)
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

