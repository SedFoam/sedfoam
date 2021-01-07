#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 14:01:08 2021

@author: chauchat
""" 
import os,sys,time,subprocess
import numpy as np

#
#caseList = ['1DSedim','1DBedLoad',...
#            '1DAvalancheMuI','1DSimpleShear','1DBoundaryLayer',...
#            '2DChannel','1DSheetFlow']
caseList = ['1DSedim','1DBedLoad']

action = "run"
#action = "clean"

basepath = '../'


#
# Reading SedFoam results
#
for case in caseList:
    try:
        print("Running case:",case)
        print('cd to directory: '+basepath+case)
        os.chdir(basepath+case)
        if (action == "run"):
            proc = subprocess.Popen(
               ['./Allrun'], stdout=subprocess.PIPE)
        elif (action == "clean"):
            proc = subprocess.Popen(
               ['./Allclean'], stdout=subprocess.PIPE)
        time.sleep(3)
    except:
        print("There is a problem with this case")
        sys.exit(0)
