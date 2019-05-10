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
from scipy.stats import multivariate_normal

BCD=imp.load_source('bhattacharyaDistUni',"./scripts/utilities/bhattacharyaDist.py")

# ---------- calclating the seperability metric ---------

def seprability(graph,threshs,coordData,labels):
	sep_mean_w=[0 for i in range(len(graph))]
	sep_mean_nw=[0 for i in range(len(graph))]
	sep_mean_wfisc=[0 for i in range(len(graph))]
	sep_mean_nwfisc=[0 for i in range(len(graph))]
	sep_eig=[0 for i in range(len(graph))]
	for i in range(len(graph)):
		# print(len(graph),i)
		total=0
		total2=0
		max_w=0
		max_nw=0
		min_w=float("inf")
		min_nw=float("inf")
		temp_w=[]
		temp_nw=[]
		weights=0
		total3=0
		total4=0
		total5=0

		adjacency_matrix=np.zeros((graph[i][0],graph[i][0]))

		# --------------- For calculating distance of all materials from all ---------------	
		# for j in range(graph[i][0]):
		# 	for k in range(j,graph[i][0],1):
		# 		total=total+(graph[i][4][j]+graph[i][4][k])*BCD.bhattacharyaDistUni(graph[i][1][j],graph[i][1][k],graph[i][3][j],graph[i][3][k],[1,1])
		# 		total2=total2+BCD.bhattacharyaDistUni(graph[i][1][j],graph[i][1][k],graph[i][3][j],graph[i][3][k],[1,1])
		# 		temp_w.append((graph[i][4][j]+graph[i][4][k])*BCD.bhattacharyaDistUni(graph[i][1][j],graph[i][1][k],graph[i][3][j],graph[i][3][k],[1,1]))
		# 		temp_nw.append(BCD.bhattacharyaDistUni(graph[i][1][j],graph[i][1][k],graph[i][3][j],graph[i][3][k],[1,1]))
		# sep_mean_w[i]=(total)/((len(graph[i][1]))**2)
		# sep_mean_nw[i]=(total2)/(len(graph[i][1])**2)

		# ---------------- For calculating distances along the edges only --------------

		for j in range(len(graph[i][2])):
			bcdval=BCD.bhattacharyaDistUni(graph[i][1][graph[i][2][j][0]],graph[i][1][graph[i][2][j][1]],graph[i][3][graph[i][2][j][0]],graph[i][3][graph[i][2][j][1]],[1,1])
			weights+=graph[i][4][graph[i][2][j][0]]+graph[i][4][graph[i][2][j][1]]
			total=total+(graph[i][4][graph[i][2][j][0]]+graph[i][4][graph[i][2][j][1]])*bcdval
			total2=total2+bcdval
			total3+=graph[i][3][graph[i][2][j][0]]
			adjacency_matrix[graph[i][2][j][0],graph[i][2][j][1]]=bcdval
	

		total_mean=total/weights
		pdfs=[]
		for j in range(graph[i][0]):
			if(graph[i][3][j]==0):
				std1=1e-4
			else:
				std1=graph[i][3][j]
			pdfs.append(multivariate_normal(mean=graph[i][1][j],cov=std1))
		
		# for k in range(coordData.shape[1]):
		# 	id1=2*labels[k]
		# 	id2=2*labels[k]+1
		# 	# print(coordData.shape)
		# 	# print(id1," ",id2," ",coordData[0,k]," ",coordData[1,k])
		# 	# print(pdfs[graph[i][6][id1]].pdf(coordData[0,k]))
		# 	total4+=np.log(1+pdfs[graph[i][6][id1]].pdf(coordData[0,k]))
		# 	total4+=np.log(1+pdfs[graph[i][6][id2]].pdf(coordData[1,k]))
			# print(total4)


		# for j in range(len(graph[i][2])):
		# 	total4+=(graph[i][4][graph[i][2][j][0]]+graph[i][4][graph[i][2][j][1]])*abs(BCD.bhattacharyaDistUni(graph[i][1][graph[i][2][j][0]],graph[i][1][graph[i][2][j][1]],graph[i][3][graph[i][2][j][0]],graph[i][3][graph[i][2][j][1]],[1,1])-total_mean)
		# 	total5+=abs(BCD.bhattacharyaDistUni(graph[i][1][graph[i][2][j][0]],graph[i][1][graph[i][2][j][1]],graph[i][3][graph[i][2][j][0]],graph[i][3][graph[i][2][j][1]],[1,1])-total_mean)

		# sep_mean_wfisc[i]=total4/total3
		# sep_mean_nwfisc[i]=total5/total3

		# eigvals=np.linalg.eigvals(adjacency_matrix)
		# sep_eig[i]=np.sum(eigvals)/threshs[i]	
		# print(total4)

		sep_mean_w[i]=(total)/((len(graph[i][2])))#*total3)#*(len(graph[i][1])))
		sep_mean_nw[i]=(total2)/(len(graph[i][2]))#*total3)#*(len(graph[i][1])))
		sep_eig[i]=(total*total4)/(len(graph[i][2]))
		
	xaxis=[i+1 for i in range(len(graph))]
	
	return xaxis,sep_mean_w,sep_mean_nw,sep_eig
