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

BCD=imp.load_source('bhattacharyaDistUni',"./scripts/utilities/bhattacharyaDist.py")

def distanceMatrix(stdData,meanhl):
	# print(hlData)
	# print(stdData)
	# print(len(meanhl))
	distMat=[[0 for i in range(len(meanhl))] for j in range(len(meanhl))]
	for i in range(0,len(meanhl),2):
		for j in range(0,len(meanhl),2):
			if(i==j):
				distMat[i][j]=float("inf")
				distMat[i+1][j]=float("inf")
				distMat[i][j+1]=float("inf")
				distMat[i+1][j+1]=float("inf")
	# print(distMat)
	for i in range(len(meanhl)):
		for j in range(i,len(meanhl)):
			if(distMat[i][j]==float("inf")):
				continue
			else:
				# continue
				distMat[i][j]=float(str.format("{0:.4f}",BCD.bhattacharyaDistUni(meanhl[i],meanhl[j],stdData[i],stdData[j],[1])))			
				# distMat[j][i]=distMat[i][j]
	return distMat