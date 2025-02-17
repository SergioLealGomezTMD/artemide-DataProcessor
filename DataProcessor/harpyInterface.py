#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  8 14:52:27 2020

@author: vla18041
"""

import harpy
import numpy
from . import DataSet
from . import DataMultiSet


def ComputeXSec(data,method="default"):
    """ Computes the cross-section values for the given data (matched to data)
    data can be: DataSet or DataMultiSet or a point-like dictionary
    methods are:
        default = usual one
        binless = evaluated with avarage values of bins only
    """
    if isinstance(data,DataSet.DataSet) or isinstance(data,DataMultiSet.DataMultiSet):
        YY=_ComputeXSec_Data(data,method)
    elif  isinstance(data,dict):
        YY=_ComputeXSec_Point(data,method)
    else:
        raise Exception('ComputeXSec works only with DataSet or DataMultiSet or a point-like dictionary')
        
    return YY

def _ComputeXSec_Data(data,method="default"):
    """ Computes the cross-section values for the given data (matched to data)
    data can be: DataSet or DataMultiSet
    methods are:
        default = usual one
        binless = evaluated with avarage values of bins only. Multiplied by area of the bin
        central = evaluated with avarage values of bins only. No further modifications.
    """
    
    if method=="default":        
        if data.processType == "DY":
            XX=harpy.DY.xSecList([d["process"] for d in data.points],
                                    [d["s"] for d in data.points],
                                    [d["qT"] for d in data.points],
                                    [d["Q"] for d in data.points],
                                    [d["y"] for d in data.points],
                                    [d["includeCuts"] for d in data.points],
                                    [d["cutParams"] for d in data.points])
        
        if data.processType == "DY_mT":
            XX=harpy.DY.xSecmTList([d["process"] for d in data.points],
                                    [d["s"] for d in data.points],
                                    [d["qT"] for d in data.points],
                                    [d["Q"] for d in data.points],
                                    [d["y"] for d in data.points],
                                    [d["mT"] for d in data.points],
                                    [d["includeCuts"] for d in data.points],
                                    [d["cutParams"] for d in data.points])

            
        elif data.processType=="SIDIS":
            XX=harpy.SIDIS.xSecList([d["process"] for d in data.points],
                                    [d["s"] for d in data.points],
                                    [d["pT"] for d in data.points],
                                    [d["z"] for d in data.points],
                                    [d["x"] for d in data.points],
                                    [d["Q"] for d in data.points],
                                    [d["includeCuts"] for d in data.points],
                                    [d["cutParams"] for d in data.points],
                                    [[d["M_target"],d["M_product"]] for d in data.points])
    elif method=="binless":
        if data.processType == "DY":
            XX1=harpy.DY.xSecListBINLESS([d["process"] for d in data.points],
                                    [d["s"] for d in data.points],
                                    [d["<qT>"] for d in data.points],
                                    [d["<Q>"] for d in data.points],
                                    [d["<y>"] for d in data.points],
                                    [d["includeCuts"] for d in data.points],
                                    [d["cutParams"] for d in data.points])
            XX=[]
            ## I weight the binless computation by the width of the bin
            for i in range(len(XX1)):
                XX.append(XX1[i]*
                          (data.points[i]["Q"][1]**2-data.points[i]["Q"][0]**2)*
                          (data.points[i]["qT"][1]**2-data.points[i]["qT"][0]**2)*
                          (data.points[i]["y"][1]-data.points[i]["y"][0]))
            
        elif data.processType=="SIDIS":
            XX1=harpy.SIDIS.xSecListBINLESS([d["process"] for d in data.points],
                                    [d["s"] for d in data.points],
                                    [d["<pT>"] for d in data.points],
                                    [d["<z>"] for d in data.points],
                                    [d["<x>"] for d in data.points],
                                    [d["<Q>"] for d in data.points],
                                    [[d["M_target"],d["M_product"]] for d in data.points])
            XX=[]
            for i in range(len(XX1)):
                XX.append(XX1[i]*
                          (data.points[i]["Q"][1]**2-data.points[i]["Q"][0]**2)*
                          (data.points[i]["pT"][1]**2-data.points[i]["pT"][0]**2)*
                          (data.points[i]["x"][1]-data.points[i]["x"][0])*
                          (data.points[i]["z"][1]-data.points[i]["z"][0]))
    elif method=="central":
        if data.processType == "DY":
            XX=harpy.DY.xSecListBINLESS([d["process"] for d in data.points],
                                    [d["s"] for d in data.points],
                                    [d["<qT>"] for d in data.points],
                                    [d["<Q>"] for d in data.points],
                                    [d["<y>"] for d in data.points],
                                    [d["includeCuts"] for d in data.points],
                                    [d["cutParams"] for d in data.points])
            
        elif data.processType=="SIDIS":
            XX=harpy.SIDIS.xSecListBINLESS([d["process"] for d in data.points],
                                    [d["s"] for d in data.points],
                                    [d["<pT>"] for d in data.points],
                                    [d["<z>"] for d in data.points],
                                    [d["<x>"] for d in data.points],
                                    [d["<Q>"] for d in data.points],
                                    [[d["M_target"],d["M_product"]] for d in data.points])
     
    YY=data.MatchWithData(XX)
    for i in range(len(YY)):
        if YY[i]<0.:
            YY[i]=100000.

    return YY

def _ComputeXSec_Point(p,method="default"):
    """ Computes the cross-section values for the given point
    methods are:
        default = usual one
        binless = evaluated with avarage values of bins only. Multiplied by area of the bin
        central = evaluated with avarage values of bins only. No further modifications.
    """
    
    if method=="default":        
        if p["type"] == "DY":
            XX1=harpy.DY.xSecList([p["process"]],[p["s"]],[p["qT"]],[p["Q"]],
                                    [p["y"]],[p["includeCuts"]],[p["cutParams"]])
        elif p["type"] == "DY_mT":
            XX1=harpy.DY.xSecmTList([p["process"]],[p["s"]],[p["qT"]],[p["Q"]],
                                    [p["y"]],[p["mT"]],[p["includeCuts"]],[p["cutParams"]])
            
        elif p["type"]=="SIDIS":
            XX1=harpy.SIDIS.xSecList([p["process"]],[p["s"]],[p["pT"]],[p["z"]],
                                    [p["x"]],[p["Q"]],[p["includeCuts"]],
                                    [p["cutParams"]],[[p["M_target"],p["M_product"]]])
        else:
            raise Exception('The dictionary is not a point.')
        
        XX=XX1[0]*p["thFactor"]
        
    elif method=="binless":
        if p["type"] == "DY":
            XX1=harpy.DY.xSecListBINLESS([p["process"]],[p["s"]],[p["<qT>"]],[p["<Q>"]],
                                    [p["<y>"]],[p["includeCuts"]],[p["cutParams"]])
            ### Matching with data (note that the normalization is ignored)
            XX=XX1[0]*(p["Q"][1]**2-p["Q"][0]**2)*(p["qT"][1]**2-p["qT"][0]**2)* \
                          (p["y"][1]-p["y"][0])*p["thFactor"]
            
        elif p["type"]=="SIDIS":
            XX1=harpy.SIDIS.xSecListBINLESS([p["process"]],[p["s"]],[p["<pT>"]],[p["<z>"]],
                                    [p["<x>"]],[p["<Q>"]],[[p["M_target"],p["M_product"]]])            

            ### Matching with data (note that the normalization is ignored)
            XX=XX1[0]*(p["Q"][1]**2-p["Q"][0]**2)*(p["pT"][1]**2-p["pT"][0]**2)* \
                (p["x"][1]-p["x"][0])*(p["z"][1]-p["z"][0])*p["thFactor"]
        else:
            raise Exception('The dictionary is not a point.')
    elif method=="central":
        if p["type"] == "DY":
            XX1=harpy.DY.xSecListBINLESS([p["process"]],[p["s"]],[p["<qT>"]],[p["<Q>"]],
                                    [p["<y>"]],[p["includeCuts"]],[p["cutParams"]])
            ### Matching with data (note that the normalization is ignored)
            XX=XX1[0]*p["thFactor"]
            
        elif p["type"]=="SIDIS":
            XX1=harpy.SIDIS.xSecListBINLESS([p["process"]],[p["s"]],[p["<pT>"]],[p["<z>"]],
                                    [p["<x>"]],[p["<Q>"]],[[p["M_target"],p["M_product"]]])            

            ### Matching with data (note that the normalization is ignored)
            XX=XX1[0]*p["thFactor"]
        else:
            raise Exception('The dictionary is not a point.')
    
       
    return XX

def ComputeChi2(data,method="default"):
    """ Computes the chi^2 values for the given data    
    data can be: DataSet or DataMultiSet
    """
    YY=ComputeXSec(data,method)
    
    ZZ=data.chi2(YY)
    
    if isinstance(data,DataSet.DataSet):
        return ZZ, [ZZ]
    else:
        return ZZ

def PrintChi2Table(data,method="default",printSysShift=True,printDecomposedChi2=False):
    """
    Compute and print the values of chi^2/Npt for experiments

    Parameters
    ----------
    data : DataSet or DataMultiSet
        The data for which the computation is made
    method : string, optional
        The method of computation of xSection.
        The default is "default".
    printSysShift: bool, optional
        If True print the list of systematic shifts (determined by nuisance parameter)
        The default is True
    printDecomposedChi2: bool, optional
        If True the chi62 presented in the decomposed form
        The default is False

    Returns
    -------
    None.

    """
    import time
    
    startT=time.time()
    
    YY=ComputeXSec(data,method)
    
    ZZ=data.chi2(YY)
    
    if isinstance(data,DataSet.DataSet):
        chi2T, chi2Part = ZZ, [ZZ]
    else:
        chi2T, chi2Part = ZZ
    
    #chi2T, chi2Part = ComputeChi2(data,method)
    
    if printSysShift:
        shift=data.DetermineAvarageSystematicShift(YY)
    
    if printDecomposedChi2:
        decChi2=data.DecomposeChi2(YY)
    
    endT=time.time()
    
    
    
    if isinstance(data,DataSet.DataSet):
        maxLength=len(data.name)
    elif isinstance(data,DataMultiSet.DataMultiSet):
        maxLength=max([len(d.name) for d in data.sets])
    else:
        maxLength=10
    
    line="{:{width}}".format("name",width=maxLength)+' | '+'  N  '+' | '
    line2="{:-<{width}}".format("",width=maxLength)+'-|-'+'-----'+'-|-'
    if printDecomposedChi2:
        line+=' chiD^2/N '+' | '+' chiL^2/N '+' | '+' chi^2/N  '+' | '
        line2+='----------'+'-|-'+'----------'+'-|-'+'----------'+'-|-'
    else:
        line+=' chi^2/N  '+' | '
        line2+='----------'+'-|-'
    
    if printSysShift:
        line+='sys.shift%'+' | '        
        line2+='----------'+'-|-'
    print(line)
    print(line2)
   
    #Only one set in MultiSet of just Set
    if len(chi2Part)==1:
        line="{:{width}} | {:5d} |".format(data.name,data.numberOfPoints,width=maxLength)
        if printDecomposedChi2:
            line+=" {:10.3f} | {:10.3f} | {:10.3f} |".format(
                decChi2[0]/data.numberOfPoints,decChi2[1]/data.numberOfPoints,decChi2[2]/data.numberOfPoints)
        else:
            line+=" {:10.3f} |".format(chi2T/data.numberOfPoints)
        if printSysShift:
            line+=" {:10.3f} |".format(shift*100)
        print(line)
    else:
        for i in range(len(chi2Part)):
            line="{:{width}} | {:5d} |".format(data.sets[i].name,data.sets[i].numberOfPoints,width=maxLength)
            if printDecomposedChi2:
                line+=" {:10.3f} | {:10.3f} | {:10.3f} |".format(
                    decChi2[i][0]/data.sets[i].numberOfPoints,
                    decChi2[i][1]/data.sets[i].numberOfPoints,
                    decChi2[i][2]/data.sets[i].numberOfPoints)
            else:
                line+=" {:10.3f} |".format(chi2Part[i]/data.sets[i].numberOfPoints)
            if printSysShift:
                line+=" {:10.3f} |".format(shift[i]*100)
            print(line)
            
        print(line2)
        line="{:{width}} | {:5d} |".format('Total',data.numberOfPoints, width=maxLength)
        if printDecomposedChi2:
            sumdChi=numpy.sum(decChi2, axis=0)
            line+=" {:10.3f} | {:10.3f} | {:10.3f} |".format(
                    sumdChi[0]/data.numberOfPoints,
                    sumdChi[1]/data.numberOfPoints,
                    sumdChi[2]/data.numberOfPoints)
        else:
            line+=" {:10.3f} |".format(chi2T/data.numberOfPoints)
        print(line)
        
    print("Computation time = ",endT-startT,' sec.')
