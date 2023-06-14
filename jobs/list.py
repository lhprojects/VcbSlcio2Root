#!/bin/python

import os



import os

def get_headnames(directory):
    if True:
        for root, dirs, files in os.walk(directory):
            print "data.append(('%s', [" % str(root),
            headnames_set = set()
            for file in files:
                headname = file.split('.', 1)[0]  # get the part of the file name before the first '.'
                if headname not in headnames_set:
                    print "'%s', "%headname,
                    headnames_set.add(headname)
            print "]",
            print ",%d)"%len(files)

    if False:
        os.system("echo >> eventnumber")
        for root, dirs, files in os.walk(directory):
            headnames_set = set()
            for file in files:
                headname = file.split('.', 1)[0]  # get the part of the file name before the first '.'
                if headname not in headnames_set:
                    headnames_set.add(headname)
                    os.system("lcio_event_counter %s >> eventnumber"%file)


get_headnames('/cefs/data/DstData/CEPC240/CEPC_v4/4fermions')
get_headnames('/cefs/data/DstData/CEPC240/CEPC_v4/2fermions')
get_headnames('/cefs/data/DstData/CEPC240/CEPC_v4/higgs')
