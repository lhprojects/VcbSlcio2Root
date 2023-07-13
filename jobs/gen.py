#!/bin/python

prefix = "/cefs/data/DstData/CEPC240/CEPC_v4/"
#IsoLepPDG = 11
#df="data_e1_20230618"
IsoLepPDG = 13
df="data_m_20230701"
rootfolder ="root"

data=[]
enable=True
def append(d):
    if enable:
        data.append(d)

enable=True
append(('/cefs/data/DstData/CEPC240/CEPC_v4/4fermions/E240.Pzz_l.e0.p0.whizard195', [ 'zz_l0mumu',  'zz_l04mu',  'zz_l0taumu',  'zz_l04tau',  'zz_l0tautau',  ] ,492))
append(('/cefs/data/DstData/CEPC240/CEPC_v4/4fermions/E240.Psznu_l.e0.p0.whizard195', [ 'sznu_l0mumu',  'sznu_l0tautau',  ] ,310))
append(('/cefs/data/DstData/CEPC240/CEPC_v4/4fermions/E240.Pzzorww_h.e0.p0.whizard195', [ 'zzorww_h0cscs',  'zzorww_h0udud',  ] ,14118))
append(('/cefs/data/DstData/CEPC240/CEPC_v4/4fermions/E240.Pww_l.e0.p0.whizard195', [ 'ww_l0ll',  ] ,1994))

enable=True
append(('/cefs/data/DstData/CEPC240/CEPC_v4/4fermions/E240.Pww_sl.e0.p0.whizard195', [ 'ww_sl0muq',  'ww_sl0tauq',  ] ,23886))

enable=True
append(('/cefs/data/DstData/CEPC240/CEPC_v4/4fermions/E240.Psw_sl.e0.p0.whizard195', [ 'sw_sl0qq',  ] ,12861))
append(('/cefs/data/DstData/CEPC240/CEPC_v4/4fermions/E240.Pww_h.e0.p0.whizard195', [ 'ww_h0cuxx',  'ww_h0ccds',  'ww_h0uusd',  'ww_h0uubd',  'ww_h0ccbs',  ] ,18689))
append(('/cefs/data/DstData/CEPC240/CEPC_v4/4fermions/E240.Psw_l.e0.p0.whizard195', [ 'sw_l0mu',  'sw_l0tau',  ] ,4316))
append(('/cefs/data/DstData/CEPC240/CEPC_v4/4fermions/E240.Psznu_sl.e0.p0.whizard195', [ 'sznu_sl0nu_down',  'sznu_sl0nu_up',  ] ,720))

enable=True
append(('/cefs/data/DstData/CEPC240/CEPC_v4/4fermions/E240.Pzz_sl.e0.p0.whizard195', [ 'zz_sl0nu_up',  'zz_sl0tau_up',  'zz_sl0nu_down',  'zz_sl0mu_up',  'zz_sl0mu_down',  'zz_sl0tau_down',  ] ,2736))

enable=True
append(('/cefs/data/DstData/CEPC240/CEPC_v4/4fermions/E240.Pzz_h.e0.p0.whizard195', [ 'zz_h0cc_nots',  'zz_h0uu_notd',  'zz_h0dtdt',  'zz_h0utut',  ] ,2546))
append(('/cefs/data/DstData/CEPC240/CEPC_v4/4fermions/E240.Psze_l.e0.p0.whizard195', [ 'sze_l0mu',  'sze_l0nunu',  'sze_l0tau',  'sze_l0e',  ] ,5428))
append(('/cefs/data/DstData/CEPC240/CEPC_v4/4fermions/E240.Pszeorsw_l.e0.p0.whizard195', [ 'szeorsw_l0l',  ] ,1230))
append(('/cefs/data/DstData/CEPC240/CEPC_v4/4fermions/E240.Pzzorww_l.e0.p0.whizard195', [ 'zzorww_l0mumu',  'zzorww_l0tautau',  ] ,2131))
append(('/cefs/data/DstData/CEPC240/CEPC_v4/4fermions/E240.Psze_sl.e0.p0.whizard195', [ 'sze_sl0dd',  'sze_sl0uu',  ] ,1575))

#data.append('/cefs/data/DstData/CEPC240/CEPC_v4/2fermions/E240.Pnn.e0.p0.whizard195', [ 'nn',  ] ,3049))
append(('/cefs/data/DstData/CEPC240/CEPC_v4/2fermions/E240.Pe3e3.e0.p0.whizard195', [ 'e3e3',  ] ,2967))
#data.append('/cefs/data/DstData/CEPC240/CEPC_v4/2fermions/E240.Pn3n3.e0.p0.whizard195', [ 'n3n3',  ] ,1647))
#data.append('/cefs/data/DstData/CEPC240/CEPC_v4/2fermions/E240.Pn1n1.e0.p0.whizard195', [ 'n1n1',  ] ,3171))
append(('/cefs/data/DstData/CEPC240/CEPC_v4/2fermions/E240.Pbhabha.e0.p0.whizard195', [ 'bhabha',  ] ,3108))
append(('/cefs/data/DstData/CEPC240/CEPC_v4/2fermions/E240.Pqq.e0.p0.whizard195', [ 'qq',  ] ,8326))
append(('/cefs/data/DstData/CEPC240/CEPC_v4/2fermions/E240.Pe2e2.e0.p0.whizard195', [ 'e2e2',  ] ,3014))
#data.append('/cefs/data/DstData/CEPC240/CEPC_v4/2fermions/E240.Pn2n2.e0.p0.whizard195', [ 'n2n2',  ] ,1588))

append(('/cefs/data/DstData/CEPC240/CEPC_v4/higgs/E240.Pqqh_X.e0.p0.whizard195', [ 'qqh_X',  ] ,569))
append(('/cefs/data/DstData/CEPC240/CEPC_v4/higgs/E240.Pe2e2h_X.e0.p0.whizard195', [ 'e2e2h_X',  ] ,81))
append(('/cefs/data/DstData/CEPC240/CEPC_v4/higgs/E240.Pe3e3h_X.e0.p0.whizard195', [ 'e3e3h_X',  ] ,88))
append(('/cefs/data/DstData/CEPC240/CEPC_v4/higgs/E240.Pe1e1h_X.e0.p0.whizard195', [ 'e1e1h_X',  ] ,77))
append(('/cefs/data/DstData/CEPC240/CEPC_v4/higgs/E240.Pnnh_X.e0.p0.whizard195', [ 'nnh_X',  ] ,439))

allch=[]
for (proc, chs, files) in data:
    for ch in chs:
        allch.append(ch)

allch = sorted(allch)
print(allch)
print(len(allch))
import sys
import os

def list_files_in_folder(folder_path):
    files = []
    for file_name in os.listdir(folder_path):
        file_path = os.path.join(folder_path, file_name)
        if os.path.isfile(file_path):
            files.append(file_path)
    return sorted(files)

def readfile(filename):
    with open(filename, 'r') as file:
        contents = file.read()
    return contents

def writefile(filename, c):
    with open(filename, 'w') as file:
        file.write(c)


pwd = os.getcwd()

def genfor(proc, ch):
    print(proc, ch)
    os.system("mkdir -p %s/%s"%(df,ch) )    

    filelist = list_files_in_folder(proc)
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
        xml = xml.replace("PLACEHOLDER_SCLIO", filenames)
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
for (proc, chs, files) in data:
    for ch in chs:
        subjobsh, mergesh, rootfile = genfor(proc, ch)
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









