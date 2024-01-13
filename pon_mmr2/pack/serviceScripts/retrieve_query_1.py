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


jobid=sys.argv[2]
email=sys.argv[1]

## Citation
cite="""How to cite PON-MMR2?\nNiroula A, Vihinen M. 2015. Classification of amino acid substitutions in mismatch repair proteins using PON-MMR2. Hum Mutat 36(12):1128-34.\n"""

## Disclaimer notice
dis="""Disclaimer\nThis non-profit server, its associated data and services are for research purposes only. The responsibility of Protein Structure and Bioinformatics Group, Lund University is limited to applying the best efforts in providing and publishing good programs and data. The developers have no responsibility for the usage of results, data or information which this server has provided.\n"""
liab="""Liability\nIn preparation of this site and service, every effort has been made to offer the most current and correct information possible. We render no warranty, express or implied, as to its accuracy or that the information is fit for a particular purpose, and will not be held responsible for any direct, indirect, putative, special, or consequential damages arising out of any inaccuracies or omissions. In utilizing this service, individuals, organizations, and companies absolve Lund University or any of their employees or agents from responsibility for the effect of any process, method or product or that may be produced or adopted by any part, notwithstanding that the formulation of such process, method or product may be based upon information provided here.\n"""


## Email messages and settings
sender="proteinbioinformatik@med.lu.se"
emailMessage="""Dear PON-MMR2 user,
The results for your submission are ready. PON-MMR2 server (http://structure.bmc.lu.se/PON-MMR2/)
received a request to send results for job : %s\n
to your e-mail address : %s\n\n
Please find the results files in the attachment.\n\n
If you have any questions, or suggestions, please feel free to contact us!\n
Email : proteinbioinformatik@med.lu.se

The PON-MMR2 team""" %(jobid,email)




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
    part=MIMEBase('application',"octet-stream")
    part.set_payload(open(path+resultFile,"rb").read())
    Encoders.encode_base64(part)
    part.add_header('Content-Disposition','attachment',filename='%s' %resultFile)
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


## Collecting the variations
queryFilesPath="/home/structure/ponp/PONP_site/queryFiles/"+jobid+"/"
protGene={"P40692":"MLH1","P43246":"MSH2","P52701":"MSH6","P54278":"PMS2"}
fil=open(queryFilesPath+jobid+".query","r")
writefile=open(queryFilesPath+jobid+"_results.txt","w")
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
    q="""select * from predictionTable where """
    for key in keys:
        if key in protGene.values():
            q+="""(geneName='%s' AND (""" % key
        else:
            q+="""(UniProtKB='%s' AND (""" % key
        variations=query[key]
        for v in variations:
            q+="""AAS='%s' OR """ % v
        q=q[:-4]+""")) OR """
    q=q[:-4]
    db=mdb.connect("localhost","abhishek","mysqlpassword","ponmmr")
    cur=db.cursor()
    cur.execute(q)
    dat=cur.fetchall()
    if not dat or dat=="None":
        t="Send email!!"
    else:
        for i in dat:
            for j in i:
                text+=str(j)+"\t"
            text=text[:-1]+"\n"
        writefile.write(text)

if t!="Send email!!":
    writefile.write("\n\n"+cite+"\n\n"+dis+"\n\n"+liab)
    writefile.close()
    resultFile="%s_results.txt" % jobid
    subject="PON-MMR2 results for jobid: %s" % jobid
    attachEmail(email,emailMessage,sender,subject,resultFile,queryFilesPath)
else:
    receiver="abhishek.niroula@med.lu.se"
    subject="PON-MMR2 script breaks!!!"
    message="PON-MMR2 script breaks somewhere, please check it for jobid: %s" % jobid
    sendEmail(receiver,message,sender,subject)
    writefile.close()

