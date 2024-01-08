#!/usr/bin/python

import MySQLdb as mdb
import sys
import os.path
import subprocess
import smtplib
from email.MIMEMultipart import MIMEMultipart
from email.MIMEBase import MIMEBase
from email.MIMEText import MIMEText
from email.Utils import COMMASPACE, formatdate
from email import Encoders
import pickle
import datetime
import random

job=sys.argv[2]
email=sys.argv[1]


#################################################################################
##
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

##################################################################################

## Citation
cite="""Abhishek Niroula and Mauno Vihinen\n
PON-mt-tRNA: a multifactorial probability-based method for classification of mitochondrial tRNA variations\n
Nucleic Acids Res. 2016 44(5):2020-2027.\n"""

## Disclaimer notice
dis="""Disclaimer\nThis non-profit server, its associated data and services are for research purposes only. The responsibility of Protein Structure and Bioinformatics Group, Lund University is limited to applying the best efforts in providing and publishing good programs and data. The developers have no responsibility for the usage of results, data or information which this server has provided.\n"""
liab="""Liability\nIn preparation of this site and service, every effort has been made to offer the most current and correct information possible. We render no warranty, express or implied, as to its accuracy or that the information is fit for a particular purpose, and will not be held responsible for any direct, indirect, putative, special, or consequential damages arising out of any inaccuracies or omissions. In utilizing this service, individuals, organizations, and companies absolve Lund University or any of their employees or agents from responsibility for the effect of any process, method or product or that may be produced or adopted by any part, notwithstanding that the formulation of such process, method or product may be based upon information provided here.\n"""


## resultFile
column="""Column description:\nmt-tRNA: Mitochondrial tRNA
Variation: Variation in the mt-tRNA
ML_probability_of_pathogenicity: Average probability of pathogenicity predicted by 20 machine leanrning (random forests) predictors
Evidence: Evidence submitted by the user
Posterior_probability_of_pathogenicity: Posterior probability after integrating ML-probability of pathogenicity and evidence. This value is available only when evidence from at least one source is submitted.
Classification: Pathogenicity class of variation. If posterior probability is not empty, the classification is based on posterior probability otherwise, the classification is based on prior probability."""


## Remove hash
def removeHash(text):
    if " " in text:
        text=text.replace(" ","")
    if "#" not in text:
        return text
    splittext=text.split("#")
    return " ".join(splittext[0].split())

## Function to send an email with attachment
def attachEmail(receiver,message,sender,subject,resultFile,path):
    msg=MIMEMultipart()
    msg['Subject']=subject
    msg['From']="proteinbioinformatik@med.lu.se"
    msg['To']=receiver
    msg.attach(MIMEText(message))
    for f in resultFile:
        part=MIMEBase('application',"octet-stream")
        part.set_payload(open(path+f,"rb").read())
        Encoders.encode_base64(part)
        part.add_header('Content-Disposition','attachment',filename='%s' %f)
        msg.attach(part)
    server=smtplib.SMTP('mail.lu.se',25)
    server.sendmail(sender,receiver,msg.as_string())


## An email to the developer
def sendEmail(receiver,message,sender,subject):
    server=smtplib.SMTP('mail.lu.se',25)
    msg=MIMEText(message)
    msg['Subject']=subject
    msg['From']="proteinbioinformatik@med.lu.se"
    msg['To']=receiver
    server.sendmail(sender,receiver,msg.as_string())
    server.quit()


## position list of the human mt-tRNAs
f=open("/home/structure/pontrna/pon-mtrna_site/trna_position_list.pkl","rb")
poslist=pickle.load(f)
f.close()

## Collecting the variations
queryFilesPath="/home/structure/ponp/PONP_site/queryFiles/"

nucs={"A":"T","C":"G","G":"C","T":"A"}

#fil=open(queryFilesPath+job+"_processed.txt","r")

#writefile=open(queryFilesPath+jobid+"_results.txt","w")
#writefile.write("mt-tRNA\tPosition\tOriginal_nucleotide\tNew_nucleotide\tClassification\tPrior_probability_of_pathogenicity\n")

error=open(queryFilesPath+job+"_error.txt","w")
#varlist=[]
#for entry in fil:
#    splitentry=entry[:-1].split("\t")
#    ref=splitentry[1]
#    alt=splitentry[2]
#    pos=splitentry[0]
#    if pos not in poslist:
#        error.write("We could not find the position for variation %s in our database. Please check that the variation position is within mitochondrial tRNA. It the position is correct, please contact us at proteinbioinformatik@med.lu.se" % entry)
#        continue
#    q="""select * from predictionTable where """
#    q+="""(position='%s' AND (""" % pos
#    q+="""(reference='%s' AND (""" % ref
#    q+="""altered='%s'))))""" % alt
#    db=mdb.connect("localhost","abhishek","mysqlpassword","ponmtRNA")
#    cur=db.cursor()
#    cur.execute(q)
#    dat=cur.fetchall()
#    if not dat or dat=="None":
#        q="""select * from predictionTable where """
#        q+="""(position='%s' AND (""" % pos
#        q+="""(reference='%s' AND (""" % nucs[ref]
#        q+="""altered='%s'))))""" % nucs[alt]
#        db=mdb.connect("localhost","abhishek","mysqlpassword","ponmtRNA")
#        cur=db.cursor()
#        cur.execute(q)
#        dat=cur.fetchall()
#        if not data or dat=="None":
#            error.write("The variation %s could not be mapped to our reference sequence" % entry)
#    for i in dat:
#        text=[]
#        for j in i:
#            text.append(str(j))
#        varlist.append(text)

#fil.close()
#db.close()

## Now look for the presence of evidence and classification
evidence_based=""
evid_f=queryFilesPath+job+"_proc_evid.txt"
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
            error.write("We could not find the position for variation %s in our database. Please check that the variation position is within mitochondrial tRNA. It the position is correct, please contact us at proteinbioinformatik@med.lu.se.\n" % entry.replace("\n",""))
            continue
        q="""select * from predictionTable where """
        q+="""(position='%s' AND (""" % pos
        q+="""(reference='%s' AND (""" % ref
        q+="""altered='%s'))))""" % alt
        db=mdb.connect("localhost","abhishek","mysqlpassword","ponmtRNA")
        cur=db.cursor()
        cur.execute(q)
        dat=cur.fetchall()
        if not dat or dat=="None":
            q="""select * from predictionTable where """
            q+="""(position='%s' AND (""" % pos
            q+="""(reference='%s' AND (""" % nucs[ref]
            q+="""altered='%s'))))""" % nucs[alt]
            db=mdb.connect("localhost","abhishek","mysqlpassword","ponmtRNA")
            cur=db.cursor()
            cur.execute(q)
            dat=cur.fetchall()
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

#db.close()

nVar=0

## Now the predictions to the files
if len(evid_varlist)>0:
    writefile=open(queryFilesPath+job+"_results.txt","w")
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

#for entry in varlist:
#    writefile.write(entry[0]+"\t"+entry[1]+","+entry[2]+","+entry[3]+"\t"+str(round(float(entry[5]),2))+"\tNA\t\t"+entry[4]+"\n")


## Write the data in the mysql database
#nVar+=len(varlist)
db=mdb.connect("localhost","abhishek","mysqlpassword","ponmtRNA")
cur=db.cursor()
try:
    now=datetime.datetime.now()
    message="Ok"
    cur.execute("""INSERT INTO userTable (Email, Date, Cases) VALUES(%s,%s, %s)""",(email,now.date(),str(nVar)))
    serial=str(cur.lastrowid)
    jobid="trna"+str(random.randint(0,999))+str(serial)
    cur.execute("""UPDATE userTable SET jobID=%s WHERE Serial=%s""",(jobid,serial))
    db.commit()
    os.makedirs(queryFilesPath+jobid)
    os.chmod(queryFilesPath+jobid,0777)
    os.rename(queryFilesPath+job+".query",queryFilesPath+jobid+"/"+jobid+".query")
#    os.rename(queryFilesPath+job+"_processed.txt",queryFilesPath+jobid+"/"+jobid+"_processed.txt")
    if os.path.exists(queryFilesPath+job+"_proc_evid.txt"):
        os.rename(queryFilesPath+job+"_proc_evid.txt",queryFilesPath+jobid+"/"+jobid+"_proc_evid.txt")
    if os.path.exists(queryFilesPath+job+"_results.txt"):
        os.rename(queryFilesPath+job+"_results.txt",queryFilesPath+jobid+"/"+jobid+"_results.txt")
    if os.path.exists(queryFilesPath+job+"_error.txt"):
        os.rename(queryFilesPath+job+"_error.txt",queryFilesPath+jobid+"/"+jobid+"_error.txt")
    subprocess.call(['chmod', "0777","%s*" %queryFilesPath+jobid+"/"])
except:
    db.rollback()
    receiver="abhishek.niroula@med.lu.se"
    subject="PON-mt-tRNA script breaks!!!"
    message="PON-mt-tRNA script breaks while trying to write in the userTable, please check it for jobid: %s" % email

db.close()

########################


## Email messages and settings
sender="proteinbioinformatik@med.lu.se"
if message=="Ok":
    emailMessage="""Dear PON-mt-tRNA user,
The results for your submission are ready. PON-mt-tRNA server (http://structure.bmc.lu.se/PON-mt-tRNA/)
received a request to send results for job : %s\n
to your e-mail address : %s\n\n
Please find the results files in the attachment.\n\n
If you have any questions, or suggestions, please feel free to contact us!\n
Email : proteinbioinformatik@med.lu.se

The PON-mt-tRNA team""" %(jobid,email)
    ## Now make a list of files to be attached
    atFile=[]
    errF=queryFilesPath+jobid+"/"+jobid+"_error.txt"
    if os.path.exists(errF):
        error=open(errF,"r")
        readError=error.read()
        error.close()
        if len(readError)>5:
            atFile.append(jobid+"_error.txt")
    if os.path.exists(queryFilesPath+jobid+"/"+jobid+"_results.txt"):
        writefile.write("\n\n"+column+"\n\n"+cite+"\n\n"+dis+"\n\n"+liab)
        writefile.close()
        atFile.append(jobid+"_results.txt")
    if len(atFile)>0:
        resultFile=atFile
        subject="PON-mt-tRNA results for jobid: %s" % jobid
        attachEmail(email,emailMessage,sender,subject,resultFile,queryFilesPath+jobid+"/")
    else:
        receiver="abhishek.niroula@med.lu.se"
        subject="PON-mt-tRNA script breaks!!!"
        message="PON-mt-tRNA script breaks while trying to send an email, please check it for jobid: %s" % jobid
        sendEmail(receiver,message,sender,subject)
        writefile.close()
else:
    receiver="abhishek.niroula@med.lu.se"
    subject="PON-mt-tRNA script breaks!!!"
    message="PON-mt-tRNA script breaks because there is no jobid created"
    sendEmail(receiver,message,sender,subject)

