#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 16:40:59 2019

@author: vla18041
"""

#######################################
# importing libraries
#######################################

import sys
import time
import numpy
sys.path.append("/home/karlik/sergiole/Sergio/Articles/2020/TMD_fits/DataProcessor")
PathToHarpy="/home/karlik/sergiole/Sergio/Articles/2020/TMD_fits/artemide-development-master/harpy"
sys.path.append(PathToHarpy)
import DataProcessor.harpyInterface
import DataProcessor.DataMultiSet




MAINPATH="/home/karlik/sergiole/Sergio/Articles/2020/TMD_fits/DataProcessor/"
#%%
#######################################
#Initialize artemide with a replica
#######################################
import harpy
path_to_constants=MAINPATH+"FittingPrograms/SV19/Constants-files/"
harpy.initialize(path_to_constants+"DY_mT_n3lo/NNPDF31/const-NNPDF31_model_4.0")
#harpy.initialize(path_to_constants+"DY_n3lo/const-NNPDF31_n3lo")
#harpy.initialize(path_to_constants+"DY_nnlo/const-HERA20_NNLO")
#harpy.initialize(path_to_constants+"DY_n3lo/const-HERA20_n3lo")
#harpy.initialize(path_to_constants+"DY_nnlo/const-MMHT14_NNLO")
#harpy.initialize(path_to_constants+"DY_nnlo/const-CT14_NNLO")
#harpy.initialize(path_to_constants+"DY_nnlo/const-PDF4LHC_NNLO")
#harpy.setNPparameters([1.93000, 0.04270, 0.22400, 9.24000, 375.00000, 2.15000, -4.97000, 0.00000]) #NNPDF31 NNLO+N3LL DY+SIDIS
#harpy.setNPparameters_TMDR(-2)
#harpy.setNPparameters_uTMDPDF(-2)
#harpy.setNPparameters_uTMDPDF_f([0.2240,0.2240,0.2240,0.2240,0.2240,0.2240,0.2240,0.2240,0.2240,0.2240,
                  #9.24,9.24,9.24,9.24,9.24,9.24,9.24,9.24,9.24,9.24,9.24,
                  #375.0,375.0,375.0,375.0,375.0,375.0,375.0,375.0,375.0,375.0,375.0,
                  #2.15,2.15,2.15,2.15,2.15,2.15,2.15,2.15,2.15,2.15,2.15,
                  #-4.97,-4.97,-4.97,-4.97,-4.97,-4.97,-4.97,-4.97,-4.97,-4.97,-4.97,
                  #0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])  # SV19 NNPDF31
#harpy.setNPparameters_uTMDPDF_f([0.37,0.37,0.37,0.145,0.33,0.37,0.33,0.145,0.37,0.37,0.37,
                  #10.,10.,10.,10.,10.4,0.,10.,0.,10.4,10.,10.,10.,
                  #680.0,680.0,680.0,680.0,680.0,680.0,680.0,680.0,680.0,680.0,680.0,
                  #2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,
                  #0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                  #0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                  #0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])  # MODEL 2.2 HERA20

#%%
#######################################
#Set PDF replica
#######################################
#harpy.setPDFreplica(int(sys.argv[1]))


#%%
### read the list of files and return the list of DataSets
def loadThisData(listOfNames):
    import DataProcessor.DataSet

    path_to_data="/home/karlik/sergiole/Sergio/Articles/2020/TMD_fits/DataProcessor/DataLib/unpolDY/"


    dataCollection=[]
    for name in listOfNames:
        loadedData=DataProcessor.DataSet.LoadCSV(path_to_data+name+".csv")
        dataCollection.append(loadedData)

    return dataCollection

#%%
#################### LOG save function
LOGPATH=MAINPATH+"FittingPrograms/SV19/W_boson/CT18NNLO/"+"W_MODEL_2.2_["+time.ctime()+"].log"
def SaveToLog(logTitle,text):
    with open(LOGPATH, 'a+', encoding='utf-8') as file:
        file.write(time.ctime())
        file.write(' --> '+logTitle+'\n')
        file.write(text)
        file.write('\n \n \n')

#%%
##################Cut function
def cutFunc(p):
    par=1.0
    if p["type"]=="DY":
        if(p["xSec"]>0):
            err=numpy.sqrt(sum([i**2 for i in p["uncorrErr"]]))/p["xSec"]
        else:
            err=100.
        delta=p["<qT>"]/p["<Q>"]

        if(p["id"][0] == "E"):
            delta=p["<qT>"]/p["Q"][1]


        if(p["id"][0:4] == "E605"):
            if(p["Q"][0]==10.5):#UPSILON resonance-bin
                return False , p
        elif(p["id"][0:4] == "E772"):
            if(p["Q"][0]<10):#these bins seems broken
                return False , p
        elif(p["id"][0:4] == "E615"):
            if(9<p["<Q>"]<11.2):#UPSILON resonance-bin
                return False , p
        elif(p["id"][0:4] == "E228"):
            if(9<p["<Q>"]<11):#UPSILON resonance-bin
                return False , p
        else:
            if(9<p["<Q>"]<11):#UPSILON resonance-bin
                return False , p

    if p["type"]=="SIDIS":
        if p["<z>"]>0.8:
            return False , p
        ## bins with low z drop
        if p["<z>"]<0.2:
            return False , p

        par=1.0
        if p["xSec"]<0.00000001:
            err=1
            delta=1
        else:
            ##############3 I MULTIPLY THE ERROR BY 100 (so it does not affect the cuts)
            err=10000#*numpy.sqrt(p.uncorrErrorsSquare)/p.xSec
            gamma2=(2.0*p["M_target"]*p["<x>"]/p["<Q>"])**2
            rho2=(p["M_product"]/p["<z>"]/(p["<Q>"]))**2
            qT=p["<pT>"]/p["<z>"]*numpy.sqrt((1+gamma2)/(1-gamma2*rho2))
            delta=qT/(p["<Q>"])

            ### compute the largest possible qT (approximate)
            gamma2WORST=(2.0*p["M_target"]*p["x"][1]/p["<Q>"])**2
            # it is definitely not a TMD point
            if gamma2WORST*rho2>1:
                return False , p
            qTWORST=p["pT"][1]/p["z"][0]*numpy.sqrt((1+gamma2WORST)/(1-gamma2WORST*rho2))

            ## drop if qT>Q/2
            if qTWORST>p["<Q>"]/2:
                return False , p

        ### drop Q<2
        if p["<Q>"]<2 :
            return False , p
    par=1.0
    if p["type"]=="DY_mT":
        if(p["xSec"]>0):
            err=numpy.sqrt(sum([i**2 for i in p["uncorrErr"]]))/p["xSec"]
        else:
            err=100.
        delta=p["<qT>"]/p["<Q>"]

#    return delta<0.5 and p.qT_avarage<80
    return (delta<0.1 or (delta<0.25 and par/err*delta**2<1)) , p

#%%
### Loading the data set
# theData=DataProcessor.DataMultiSet.DataMultiSet("DYmTset",loadThisData(['A7_mT', 'A7_mT_e', 'A7_mT_m',
#                                                                         'CMS8_mT_e', 'CMS8_mT_m',
#                                                                         'D01.8_mT', 'CDF1.8_mT',
#                                                                         '13_high_mT_minus','13_high_mT_plus',
#                                                                         '13_middle_mT_minus','13_middle_mT_plus',
#                                                                         '13_low_mT_minus','13_low_mT_plus']))

theData=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisData(['CMS8_mT_e', 'CMS8_mT_m']))
# theData=DataProcessor.DataMultiSet.DataMultiSet("DYset",loadThisData(['Z_13_TeV_low', 'Z_13_TeV_middle']))
setDY=theData.CutData(cutFunc)

print('Loaded ', setDY.numberOfSets, 'data sets with', sum([i.numberOfPoints for i in setDY.sets]), 'points.')
print('Loaded experiments are', [i.name for i in setDY.sets])
#print('Transverse mass ranges are', [i.points[0]["mT"] for i in setDY.sets])

#for i in setDY.sets:
#    for j in range(i.numberOfPoints):
#        i.points[j]["mT"][1]=450.0

#print('Transverse mass ranges are', [i.points[0]["mT"] for i in setDY.sets])

#%%
#######################################
# Test the main replica
#######################################

#print("We do not set the mass of the lepton to zero")
#print("PDF: MSHT20 MODEL 4.0")
#print("Theoretical uncertainties: Scale variation")
# DataProcessor.harpyInterface.PrintChi2Table(setDY,printDecomposedChi2=True)
#Chi2=DataProcessor.harpyInterface.ComputeChi2(setDY)
#print(Chi2)

#%%
##########################
# Spectrum qT
##########################

#for i in range(25):
    #MW=60.+5.*float(i)
    #harpy.setEWparameters([MW,2.0850000000000000])
# text="""
# Z boson
# We do set the mass of leptons to zero
# We use the MODEL 4.0
# PDF: NNPDF31
# Data: Z_13_TeV_low.csv, Z_13_TeV_middle.csv
# """
# print(text)
# startT=time.time()
XSec=DataProcessor.harpyInterface.ComputeXSec(setDY)
print(XSec)
#
# ########################### c1=1.0, c2=0.5, c3=1.0, c4=1.0 ###############################
#
# harpy.varyScales(1.0,0.5,1.0,1.0)
# XSecc205=DataProcessor.harpyInterface.ComputeXSec(setDY)
#
#
# ########################### c1=1.0, c2=2.0, c3=1.0, c4=1.0 ###############################
#
# harpy.varyScales(1.0,2.0,1.0,1.0)
# XSecc22=DataProcessor.harpyInterface.ComputeXSec(setDY)
#
#
# ########################### c1=1.0, c2=1.0, c3=1.0, c4=0.5 ###############################
#
# harpy.varyScales(1.0,1.0,1.0,0.5)
# XSecc405=DataProcessor.harpyInterface.ComputeXSec(setDY)
#
#
# ########################### c1=1.0, c2=1.0, c3=1.0, c4=2.0 ###############################
#
# harpy.varyScales(1.0,1.0,1.0,2.0)
# XSecc42=DataProcessor.harpyInterface.ComputeXSec(setDY)
#
# k=0
# for i in setDY.sets:
#     output="{"
#     print(i.name)
#     result={}
#     result={"id":i.name}
#     result["qT-XSec"]=[]
#     print(i.numberOfPoints)
#     for j in range(i.numberOfPoints):
#         result["qT-XSec"].append([i.points[j]["qT"][0],XSec[k],
#                                     XSecc205[k],XSecc22[k],XSecc405[k],XSecc42[k]])
#         scale=str(XSecc205[k])+","+str(XSecc22[k])+","+str(XSecc405[k])+","+str(XSecc42[k])
#         output=output+"{"+str(i.points[j]["qT"][0])+","+str(XSec[k])+","+scale+"},"
#         k=k+1
#     #print(result)
#     output=output+"}"
#     print(output)
# endT=time.time()
#
# print(endT-startT)

#%%
#######################################
# Minimisation
#######################################
#totalN=setDY.numberOfPoints

#def chi_2(x):
   #startT=time.time()
   #harpy.setNPparameters(x)
   ##harpy.setEWparameters([x[0],x[1]])
   ##NPparametersf=[0.324,0.324,0.324,x[1],x[0],0.324,x[0],x[1],0.324,0.324,0.324,
                  ##13.200,13.200,13.200,x[3],x[2],13.200,x[2],x[3],13.200,13.200,13.200,
                  ##356.,356.,356.,x[5],x[4],356.,x[4],x[5],356.,356.,356.,
                  ##2.050,2.050,2.050,2.050,2.050,2.050,2.050,2.050,2.050,2.050,2.050,
                  ##-10.400,-10.400,-10.400,-10.400,-10.400,-10.400,-10.400,-10.400,-10.400,-10.400,-10.400,
                  ##0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
   ##harpy.setNPparameters_uTMDPDF_f(NPparametersf)
   #print('np set =',["{:8.4f}".format(i) for i in x], end =" ")

   #ccDY2,cc3=DataProcessor.harpyInterface.ComputeChi2(setDY)

   #cc=ccDY2/totalN
   #endT=time.time()
   #print(':->',cc,'       t=',endT-startT)
   #return ccDY2

#%%

#from iminuit import Minuit


#######################################
# Set Model
#######################################

################################ Model 1.0 #####################################

################################ HERA20 ########################################
#initialValues=(2.,0.0331,0.23,5.6,0.252,8.,570.,0.,0.)
#initialErrors=(0.1,0.0016,0.0022,0.5,0.015,0.6,13.,0.1, 0.1)
#searchLimits=((1.4,4.5),(0.0001,5.),(0.,2.),(0.,25.),(0.,2.),(0.,25.),(0.,1000.),None,None)
#parametersToMinimize=(True,False,False,False,False,False,False,True,True)

################################ NNPDF31 #######################################
#initialValues=(2.,0.029,0.35,2.6,0.152,7.6,230.,0.,0.)
#initialErrors=(0.1,0.007,0.07,1.1,0.031,2.8,180.,0.1, 0.1)
#searchLimits=((1.4,4.5),(0.0001,5.),(0.,2.),(0.,25.),(0.,2.),(0.,25.),(0.,1000.),None,None)
#parametersToMinimize=(True,False,False,False,False,False,False,True,True)

################################ CT18NNLO ######################################
#initialValues=(2.,0.0391,0.16,7.9,0.212,5.3,700.0,0.,0.)
#initialErrors=(0.1,0.0024,0.04,1.,0.030,1.,130.,0.1, 0.1)
#searchLimits=((1.4,4.5),(0.0001,5.),(0.,2.),(0.,25.),(0.,2.),(0.,25.),(0.,1000.),None,None)
#parametersToMinimize=(True,False,False,False,False,False,False,True,True)

################################ Model 2.0 #####################################

################################ HERA20 ########################################
#initialValues=(2.,0.0377,0.,13.3,0.45,0.,0.239,8.1,1000.)
#initialErrors=(0.1,0.0021,0.2,0.9,0.09,0.4,0.027,1.1,90)
#searchLimits=((1.4,4.5),(0.0001,5.),(0.,2.),(0.,25.),(0.,2.),(0.,25.),(0.,2.),(0.,25.),(0.,1000.))
#parametersToMinimize=(True,False,False,False,False,False,False,False,False)

################################ NNPDF31 #######################################
#initialValues=(2.,0.02993,0.3173,2.189,0.37,3.8,0.1371,8.067,248.6)
#initialErrors=(0.1,0.014,0.0021,0.017,0.10,1.7,0.0018,0.14, 0.07)
#searchLimits=((1.4,4.5),(0.0001,5.),(0.,2.),(0.,25.),(0.,2.),(0.,25.),(0.,2.),(0.,25.),(0.,1000.))
#parametersToMinimize=(True,False,False,False,False,False,False,False,False)

################################# CT18NNLO ######################################
#initialValues=(2.,0.0422,0.04,13.5,0.36,0.,0.208,3.1,730.)
#initialErrors=(0.1,0.0017,0.16,2.3,0.28,1.,0.034,0.6,170)
#searchLimits=((1.4,4.5),(0.0001,5.),(0.,2.),(0.,25.),(0.,2.),(0.,25.),(0.,2.),(0.,25.),(0.,1000.))
#parametersToMinimize=(True,False,False,False,False,False,False,False,False)

################################# Model 2.1 #####################################

################################# HERA20 ########################################
#initialValues=(2.,0.0369,0.15,9.27,0.16,8.6,0.64,0.0,1030.)
#initialErrors=(0.1,0.0030,0.05,0.18,0.05,1.3,0.14,2.1,100)
#searchLimits=((1.4,4.5),(0.0001,5.),(0.,10.),(0.,50.),(0.,10.),(0.,50.),(0.,10.),(0.,50.),(0.,2500.))
#parametersToMinimize=(True,False,False,False,False,False,False,False,False)

################################# NNPDF31 #######################################
#initialValues=(2.,0.0289,0.347,3.1,0.109,7.8,0.254,6.1,252.)
#initialErrors=(0.1,0.0015,0.011,1.9,0.011,2.0,0.026,1.8,15)
#searchLimits=((1.4,4.5),(0.0001,5.),(0.,10.),(0.,50.),(0.,10.),(0.,50.),(0.,10.),(0.,50.),(0.,2500.))
#parametersToMinimize=(True,False,False,False,False,False,False,False,False)

################################# CT18NNLO ######################################
#initialValues=(2.,0.0389,0.159,7.9,0.210,5.3,0.210,5.,700.)
#initialErrors=(0.1,0.0034,0.027,0.5,0.018,1.6,0.04,11,600)
#searchLimits=((1.4,4.5),(0.0001,5.),(0.,10.),(0.,50.),(0.,10.),(0.,50.),(0.,10.),(0.,50.),(0.,2500.))
#parametersToMinimize=(True,False,False,False,False,False,False,False,False)

################################# Model 2.2 #####################################

################################# HERA20 ########################################
#initialValues=(2.,0.0341,0.145,10.4,0.33,0.00,0.37,10.,680.)
#initialErrors=(0.1,0.0011,0.013,0.5,0.05,0.35,0.12,6.,80.)
#searchLimits=((1.4,4.5),(0.0001,5.),(0.,10.),(0.,50.),(0.,10.),(0.,50.),(0.,10.),(0.,50.),(0.,2500.))
#parametersToMinimize=(True,False,False,False,False,False,False,False,False)

################################# NNPDF31 #######################################
#initialValues=(2.,0.0301,0.188,5.54,0.2004,4.375,0.4861,0.00,232.8)
#initialErrors=(0.1,0.0006,0.012,0.26,0.0009,0.012,0.0026,0.06,0.6)
#searchLimits=((1.4,4.5),(0.0001,5.),(0.,10.),(0.,50.),(0.,10.),(0.,50.),(0.,10.),(0.,50.),(0.,2500.))
#parametersToMinimize=(True,True,False,False,False,False,True,True,True)

################################# CT18NNLO ######################################
#initialValues=(2.,0.0418,0.09,12.5,0.29,0.00,0.00,17.,820.)
#initialErrors=(0.1,0.0026,0.05,1.7,0.05,0.15,0.07,6.,200.)
#searchLimits=((1.4,4.5),(0.0001,5.),(0.,10.),(0.,50.),(0.,10.),(0.,50.),(0.,10.),(0.,50.),(0.,2500.))
#parametersToMinimize=(True,False,False,False,False,False,False,False,False)

################################# Model 3.0 #####################################

################################ HERA20 ########################################
#initialValues=(2.,0.0340,0.7244,9.700,0.0035,11.800,1000.00,0.0430,-0.270)
#initialErrors=(0.1,0.0028,0.0400,1.200,0.0600,2.4000,90.00,0.0180,0.050)
#searchLimits=((1.40,4.50),(0.0001,5.0),(0.000,2.000),(0.,25.00),(0.000,2.000),(0.,25.),(0.,1000.),None,(-0.5,10.))
#parametersToMinimize=(True,True,False,True,False,True,True,True,True)

################################# NNPDF31 #######################################
#initialValues=(2.,0.023,0.489,4.7,0.390,12.89,320.,0.136,-0.141)
#initialErrors=(0.1,0.005,0.018,1.2,0.032,0.26,40.,0.025,0.022)
#searchLimits=((1.4,4.5),(0.0001,5.),(0.,2.),(0.,25.),(0.,2.),(0.,25.),(0.,1000.),None,(-0.5,10.))
#parametersToMinimize=(True,False,False,False,False,False,False,False,False)

################################# CT18NNLO ######################################
#initialValues=(2.,0.0329,0.31,11.0,0.482,12.8,740.,0.168,-0.130)
#initialErrors=(0.1,0.0027,0.12,0.7,0.024,1.1,140.,0.011,0.013)
#searchLimits=((1.4,4.5),(0.0001,5.),(0.,2.),(0.,25.),(0.,2.),(0.,25.),(0.,1000.),None,(-0.5,10.))
#parametersToMinimize=(True,False,False,False,False,False,False,False,False)



################################# W mass fit ####################################
#initialValues=(80.379,2.085)
#initialErrors=(0.001,0.001)
#searchLimits=((78.,82.),(2.,3.))
#parametersToMinimize=(False,True)



#%%
#text="""
#W boson
#We do not set the mass of leptons to zero
#We use the MODEL 2.2
#PDF: CT18NNLO
#Non-perturbative parameters fitted: c_0, lambda_1, lambda_2, lambda_3, lambda_4,
#lambda_5, lambda_6, lambda_7
#Data: A7_mT.csv, A7_mT_e.csv, A7_mT_m.csv, CMS8_mT_e.csv, CMS8_mT_m.csv,
#D01.8_mT.csv, CDF1.8_mT.csv
#"""

#SaveToLog("MODEL 2.2; PDF: CT18NNLO",text=text)

#Chi2Total, Chi2=DataProcessor.harpyInterface.ComputeChi2(setDY)
#Chi2Exps={"Chi^2 Total:":[Chi2Total/totalN,totalN]}
#k=0
#for i in setDY.sets:
    #Chi2Exps[i.name]=[Chi2[k]/i.numberOfPoints,i.numberOfPoints]
    #k=k+1
#SaveToLog("Chi^2/d.o.f=",str(Chi2Exps))


#m = Minuit.from_array_func(chi_2, initialValues,
      #error=initialErrors, limit=searchLimits, fix=parametersToMinimize, errordef=1)


#m.tol=0.0001*totalN*10000 ### the last 0.0001 is to compensate MINUIT def
#m.strategy=1

#SaveToLog("MINIMIZATION STARTED",str(m.params))
#%%


#m.tol=0.0001*totalN*10000 ### the last 0.0001 is to compensate MINUIT def
#m.strategy=1


#m.migrad()




#SaveToLog("MINIMIZATION FINISHED",str(m.params))
#SaveToLog("CORRELATION MATRIX",str(m.matrix(correlation=True)))
#%%

#m.hesse()


#SaveToLog("HESSE FINISHED",str(m.params))
#SaveToLog("CORRELATION MATRIX",str(m.matrix(correlation=True)))
#%%

#m.minos()


#SaveToLog("MINOS FINISHED",str(m.params))
#SaveToLog("CORRELATION MATRIX",str(m.matrix(correlation=True)))


########################
#New Chi Squares
########################


#harpy.setNPparameters([m.values["x0"],m.values["x1"],m.values["x2"],m.values["x3"],
                       #m.values["x4"],m.values["x5"],m.values["x6"],m.values["x7"],
                       #m.values["x8"]])


#Chi2Total, Chi2=DataProcessor.harpyInterface.ComputeChi2(setDY)
#Chi2Exps={"Chi^2 Total:":[Chi2Total/totalN,totalN]}
#k=0
#for i in setDY.sets:
    #Chi2Exps[i.name]=[Chi2[k]/i.numberOfPoints,i.numberOfPoints]
    #k=k+1
#SaveToLog("Chi^2/d.o.f=",str(Chi2Exps))
