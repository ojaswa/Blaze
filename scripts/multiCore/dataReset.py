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

actualData=[]
mx=0

def dataValueReset(i):
	r=[]
	r.append(i)
	r.append(mx-actualData[i])
	return r


def dataReset(actualDataip,mxip):	
	global actualData
	actualData=actualDataip
	global mx
	mx=mxip
	dataLabelsFinals=[0 for i in range(len(actualData))]
	ptemparray=[]
	pool=Pool(processes=3)
	inputs=range(int(len(actualData)))
	t=pool.map(dataValueReset,inputs)
	for i in range(len(t)):
		a=t[i][0]
		b=t[i][1]
		dataLabelsFinals[a]=b
	return dataLabelsFinals	


database=os.path.normpath(os.getcwd() + os.sep + os.pardir)
database=os.path.join(database,"Volume data")
database=os.path.join(database,"NRRD")
print(database)

dataheadFileName="calix_resizehalf.nhdr"
datarawFileName="calix_resizehalf.raw"

datahead=open(os.path.join(database,dataheadFileName),"r")
datafile=open(os.path.join(database,datarawFileName),"rb")
actualData=np.fromfile(datafile, np.uint16)
mxa=256**2
resetdata=dataReset(actualData,mxa)
filename="reset_"+dataheadFileName
wfile=open(filename,'w')
wfile.write(resetdata)

