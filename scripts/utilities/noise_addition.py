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
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
from collections import namedtuple
import re
import array
import math
import random
import operator
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize
import sys
import imp
import scipy.stats as stat
import time
import pickle
import sys

def noise_addition(volumePath,volumeType,file_name,sigma):
	np.random.seed(0)
	datafile=open(volumePath,"rb")
	maxval=255
	actualData=np.fromfile(datafile, np.uint8)
	noise=np.random.normal(0,sigma,len(actualData))
	actualData=np.int32(actualData)	
	noise=np.int32(noise)
	actualData_new=np.asarray(actualData)
	actualData_new=actualData_new+noise
	for i in range(len(actualData_new)):
		if(actualData_new[i]<0):
			actualData_new[i]=0
		if(actualData_new[i]>maxval):
			actualData_new[i]=maxval
	sigma=str(sigma)
	sigma=sigma.replace(".","_")
	ext="_"+str(sigma)+".raw"
	volumePath="../../data/NRRD/phantom_data/"+file_name+ext
	actualData_new=np.uint8(actualData_new)
	datafile=open(volumePath,"wb")
	datafile.write(bytearray(actualData_new))
	datafile.close()
	return file_name+ext

arguments=sys.argv
fileName="../../data/NRRD/"+str(arguments[1])+".nhdr"
fileNameFile=open(fileName,"r")
filePathData=[]
for line in fileNameFile:
    filePathData.append(line)
fileNameFile.close()
volumeType=filePathData[1]
volumePath = fileName.replace("\n", "");
volumePath = fileName.replace(".nhdr", "")
volumePath=volumePath+".raw"
base=os.path.normpath(volumePath)
os.makedirs("../../data/NRRD/phantom_data",exist_ok=True)
sigma=[0.25,0.45,0.85,1.65,3.25]

for i in sigma:
	newVolumeName=noise_addition(volumePath,volumeType,arguments[1],i)
	filePathData[6]='data file: ./'+newVolumeName+'\n'
	newVolumeName=newVolumeName.replace(".raw",".nhdr")
	f=open("../../data/NRRD/phantom_data/"+newVolumeName,"w")
	for i in filePathData:
		f.write(i)
	f.close()

