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
from matplotlib.colors import LogNorm
import matplotlib.cm as cm
import matplotlib.pyplot as plt
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

# ---------------- function to subsample the LH space ----------------- 

def mean_and_std(actualData,data_array,c2,labels,volumeType):
    
    random.seed(0)
    bins=[[[] for i in range(256)] for j in range(256)]
    idxs=range(len(data_array[0]))
    maxH=np.max(data_array[0])
    minH=np.min(data_array[0])
    maxL=np.max(data_array[1])
    minL=np.min(data_array[1])
    
    # ------------- binning -------------------
    for i in range(len(data_array[0])):
        t=[]
        if(data_array[0][i]<data_array[1][i] and actualData[i]>0):
            t=[data_array[0][i],data_array[1][i]]
            if(volumeType=="1"):
                bins[int(data_array[0][i])][int(data_array[1][i])].append(t)
            else:
                bins[int((data_array[0][i]-minH)*255/(maxH-minH))][int((data_array[1][i]-minL)*255/(maxL-minL))].append(t)
    flag=[]
    labelwiseactualdatal=[]
    labelwiseactualdatah=[]
    for i in range(256*256):
        flag.append(0)
    meanhl=[]
    stdhl=[]
    num_clus=[]
    alloc_x=0
    alloc_y=0
    unique_labels=np.unique(labels)
    for i in range(len(unique_labels)):
        if(unique_labels[i]==-1):
            continue
        num_clus.append(unique_labels[i])
        
    for i in range(len(num_clus)):
        labelwiseactualdatah.append([])
        labelwiseactualdatal.append([])
    for i in range(len(num_clus)):
        mymembers=labels==num_clus[i]
        for l in c2[mymembers]:
            temph=[]
            templ=[]
            if(volumeType=="1"):
                alloc_x = int(l[0])
                alloc_y = int(l[1])
            else:
                alloc_x=int((l[0]-minH)*255/(maxH-minH))
                alloc_y=int((l[1]-minL)*255/(maxL-minL))
            if(flag[alloc_x+alloc_y*256])==0:
                for j in range(len(bins[alloc_x][alloc_y])):
                    temph.append(bins[alloc_x][alloc_y][j][1])
                    templ.append(bins[alloc_x][alloc_y][j][0])
                labelwiseactualdatah[i].extend(temph)
                labelwiseactualdatal[i].extend(templ)
                flag[alloc_x+alloc_y*256]=1

    for i in range(len(labelwiseactualdatah)):
        meanhl.append(np.mean(labelwiseactualdatal[i]))
        meanhl.append(np.mean(labelwiseactualdatah[i]))
        if(np.std(labelwiseactualdatal[i])==0):
            stdhl.append(0.1)
        else:
            stdhl.append(np.std(labelwiseactualdatal[i]))
        if(np.std(labelwiseactualdatah[i])==0):
            stdhl.append(0.1)
        else:
            stdhl.append(np.std(labelwiseactualdatah[i]))
    return meanhl,stdhl
