###########################################################################
##                                                                        ##
##  Blaze - a volume rendering and analytics program                      ##
##  Copyright (C) 2016-2018 Graphics Research Group, IIIT Delhi           ##
##                                                                        ##
##  This program is free software: you can redistribute it and/or modify  ##
##  it under the terms of the GNU General Public License as published by  ##
##  the Free Software Foundation, either version 3 of the License, or     ##
##  (at your option) any later version.                                   ##
##                                                                        ##
##  This program is distributed in the hope that it will be useful,       ##
##  but WITHOUT ANY WARRANTY; without even the implied warranty of        ##
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         ##
##  GNU General Public License for more details.                          ##
##                                                                        ##
##  You should have received a copy of the GNU General Public License     ##
##  along with this program.  If not, see http://www.gnu.org/licenses/.   ##
##                                                                        ##
############################################################################
##           Author: Tushar Arora                                         ##
##           E-mail: tushar15107@iiitd.ac.in                              ##
##           Date  : 14.12.2016                                           ##
############################################################################
import numpy as np
import matplotlib.pyplot as plot
import os
from collections import namedtuple
import re
import array
from sklearn.cluster import DBSCAN
from sklearn import datasets
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import MeanShift, estimate_bandwidth
import math
import random
import operator
import sys
import imp
from multiprocessing import Pool

sampleLHIdxs=[]
data_array=[]
actualData=[]
labelwisemeanl=[]
labelwisemeanh=[]
sigmahlfinaldet=[]
sigmahlfinalinv=[]
verticeSetAlpha=[]
verticeSetAlphaDeviation=[]
# dataLabelsFinals=[]

def gaussianSamplingEdge(i):
    if(i%100000==0):
        print(i)
    # global dataLabelsFinals
    if(actualData[sampleLHIdxs[i]]>0 and data_array[0][sampleLHIdxs[i]]<data_array[1][sampleLHIdxs[i]]):
        # print("here")
        probabilityTemp=[]
        for j in range(len(sigmahlfinaldet)):
            xtmp=np.matrix([[data_array[0][sampleLHIdxs[i]]],[data_array[1][sampleLHIdxs[i]]]])
            mutemp=np.matrix([[labelwisemeanh[j]],[labelwisemeanl[j]]])
            prob=(1/(2*math.pi*(sigmahlfinaldet[j]**(1/2))))*np.exp((-1/2)*(xtmp-mutemp).transpose()*sigmahlfinalinv[j]*(xtmp-mutemp))
            probabilityTemp.append(prob)
        
        max=0   
        maxidx=-1
        for x in range(len(probabilityTemp)):
            if(probabilityTemp[x]>=max):
                max=probabilityTemp[x]
                maxidx=x     
        # dataLabelsFinals[sampleLHIdxs[i]]=maxidx+1  
        # putGlobalVal(sampleLHIdxs[i],maxidx+1)
        t=[]
        t.append(sampleLHIdxs[i])
        t.append(maxidx+1)
        return t        
    else:
        t=[]
        t.append(sampleLHIdxs[i])
        t.append(0)
        return t

def gaussianSamplingMaterial(i):        
    if(i%100000==0):
        print(i)
    # global dataLabelsFinals
    if(actualData[i]>0):
        # print("here")
        probabilityTemp=[]
        for j in range(len(verticeSetAlpha)):
            xMu=verticeSetAlpha[j]
            xStd=verticeSetAlphaDeviation[j]
            prob=(1/((2*math.pi)**(1/2)*(xStd**(1/2))))*np.exp((-1/2)*(((actualData[i]-xMu)**2)/(xStd**2)))
            probabilityTemp.append(prob)
        
        max=0   
        maxidx=-1
        for x in range(len(probabilityTemp)):
            if(probabilityTemp[x]>=max):
                max=probabilityTemp[x]
                maxidx=x     
        t=[]
        t.append(i)
        t.append(maxidx+1)
        return t        
    else:
        t=[]
        t.append(i)
        t.append(0)
        return t

def putGlobalVal(i,x):
	global dataLabelsFinals
	dataLabelsFinals[i]=x

def getGlobalVal():
	global dataLabelsFinals
	print(sum(dataLabelsFinals[1:1000000]))

def multiCoreEdge(sampleLHIdxsip,data_arrayip,actualDataip,labelwisemeanlip,labelwisemeanhip,sigmahlfinalinvip,sigmahlfinaldetip):
	global sampleLHIdxs 
	sampleLHIdxs=sampleLHIdxsip
	global data_array
	data_array=data_arrayip
	global actualData
	actualData=actualDataip
	global labelwisemeanh
	labelwisemeanh=labelwisemeanhip
	global labelwisemeanl
	labelwisemeanl=labelwisemeanlip
	global sigmahlfinaldet
	sigmahlfinaldet=sigmahlfinaldetip
	global sigmahlfinalinv
	sigmahlfinalinv=sigmahlfinalinvip
	# global dataLabelsFinals
	dataLabelsFinals=[0 for i in range(len(data_array[0]))]
	ptemparray=[]
	pool=Pool(processes=3)
	inputs=range(int(len(sampleLHIdxs)))
	t=pool.map(gaussianSamplingEdge,inputs)
	for i in range(len(t)):
		a=t[i][0]
		b=t[i][1]
		dataLabelsFinals[a]=b
	# print(sum(dataLabelsFinals))
	return dataLabelsFinals

def multiCoreMaterial(data_arrayip,actualDataip,verticeSetAlphaip,verticeSetAlphaDeviationip):
	global data_array
	data_array=data_arrayip
	global actualData
	actualData=actualDataip
	global verticeSetAlpha
	verticeSetAlpha=verticeSetAlphaip
	global verticeSetAlphaDeviation
	verticeSetAlphaDeviation=verticeSetAlphaDeviationip
	# global dataLabelsFinals
	dataLabelsFinals=[0 for i in range(len(data_array[0]))]
	pool=Pool(processes=3)
	inputs=range(len(data_array[0]))
	t=pool.map(gaussianSamplingMaterial,inputs)
	for i in range(len(t)):
		a=t[i][0]
		b=t[i][1]
		dataLabelsFinals[a]=b
	print(sum(dataLabelsFinals))	
	return dataLabelsFinals
