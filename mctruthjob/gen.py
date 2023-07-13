#!/bin/python

import sys
import os

df="data"
rootfolder="root"
IsoLepPDG=13

stdhepfolder="/cefs/higgs/liangh/data/stdhep/162.5GeV/WhizardAis/data/"
data=[]
data.append( (stdhepfolder+"/4fermions/E162.5.Pww_sl.e0.p0.whizard195", ("ww_sl0muq", "ww_sl0tauq") ) )
data.append( (stdhepfolder+"/4fermions/E162.5.Pzz_sl.e0.p0.whizard195", ("zz_sl0mu", "zz_sl0tau") ) )
data.append( (stdhepfolder+"/2fermions/E162.5.Pqq.e0.p0.whizard195", ("qq",) ) )

def list_files_in_folder(folder_path):
    files = []
    for file_name in os.listdir(folder_path):
        file_path = os.path.join(folder_path, file_name)
        if os.path.isfile(file_path):
            files.append(file_path)
    files = sorted(files)
    print(files)
    return files

def readfile(filename):
    with open(filename, 'r') as file:
        contents = file.read()
    return contents

def writefile(filename, c):
    with open(filename, 'w') as file:
        file.write(c)


pwd = os.getcwd()

def genfor(proc_folder, ch):
    print(proc_folder, ch)
    os.system("mkdir -p %s/%s"%(df,ch) )    

    filelist = list_files_in_folder(proc_folder)
    filelist = [file for file in filelist if ch in file]
    listlist = [filelist[i:i + 5] for i in range(0, len(filelist), 5)]

    lastjobidx = 0
    for i, group in enumerate(listlist):
        lastjobidx = i

        filenames = "\n".join(group)


        rootfilename = os.getcwd() + "/" + "%s/%s/%s_%d.root"%(df,ch,ch, i)
        xmlfilename = os.getcwd() + "/" + "%s/%s/%s_%d.xml"%(df,ch,ch, i)
        jobfilename = os.getcwd() + "/" + "%s/%s/%s_%d.sh"%(df,ch,ch, i)


        xml = readfile("SteerTemplate.xml")
        #xml = xml.replace("PLACEHOLDER_SCLIO", filenames)
        xml = xml.replace("PLACEHODER_STDHEP", filenames)
        xml = xml.replace("PLACEHOLDER_ROOT", rootfilename)
        xml = xml.replace("PLACEHOLDER_CenterOfMassEnergy", "%f"%240.0)
        xml = xml.replace("PLACEHOLDER_IsoLepPDG", "%d"%IsoLepPDG)
        writefile(xmlfilename, xml)

        job = readfile("JobTemplate.sh")
        job = job.replace("PLACEHOLDER_XML",   xmlfilename  )
        writefile(jobfilename, job)
        os.system("chmod +x %s"%jobfilename)
       
        if lastjobidx > 10 and False:
            break
    
    subjobsh = "%s/%s/SubJobs.sh"%(df,ch)
    mergesh = "%s/%s/Merge.sh"%(df,ch)
    rootfile = "%s/%s/%s/%s.root"%(pwd, df, ch,ch)

    writefile(subjobsh, "hep_sub %s -n %d"%(   os.getcwd() + "/%s/%s/%s_%%{ProcId}.sh" % (df,ch,ch)         ,  lastjobidx  ))
    os.system("chmod +x %s"%subjobsh)
    writefile(mergesh, "rm -f %s; hadd -f6 %s %s/%s/%s/%s_*.root"%(rootfile, rootfile,pwd,df,ch,ch)  )
    os.system("chmod +x %s"%mergesh)
    return subjobsh, mergesh, rootfile

os.system("mkdir -p %s"%df)
os.system("mkdir -p %s/%s"% (df, rootfolder) )
mergeall = ""
subjoball = ""
copyroot=""


for (proc_folder, chs) in data:
    for ch in chs:
        subjobsh, mergesh, rootfile = genfor(proc_folder, ch)
        mergeall +=  pwd + "/%s\n"%mergesh 
        subjoball +=  pwd + "/%s\n"%subjobsh
        copyroot += "cp -rf %s %s/%s/%s\n"%(rootfile, pwd, df, rootfolder)

mergeallsh = "%s/mergeall.sh"%df
writefile(mergeallsh, mergeall)
os.system("chmod +x %s"%mergeallsh)

subjoballsh = "%s/suball.sh"%df
writefile(subjoballsh, subjoball)
os.system("chmod +x %s"%subjoballsh)

copyrootsh = "%s/copyroot.sh"%df
writefile(copyrootsh, copyroot)
os.system("chmod +x %s"%copyrootsh)








