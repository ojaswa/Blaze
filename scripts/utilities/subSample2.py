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
from matplotlib.colors import LogNorm
import sys
import imp
from multiprocessing import Pool

def subSampling2(int_data,occ_data):

    #------------Sampling 100000 points uniformly-----------------
    num_points=int(0.08*len(int_data))
    if(num_points>100000):
        num_points=100000
    random.seed(0)
    coordData=[]
    for i in range(len(int_data)):
        coordData.append([int_data[i],occ_data[i]])
    sample_points=random.sample(range(len(coordData)),num_points)
    cData=[coordData[i] for i in sample_points]
    maxH=np.max(np.asarray(int_data))
    minH=np.min(np.asarray(int_data))
    maxL=np.max(np.asarray(occ_data))
    minL=np.min(np.asarray(occ_data))
    subSampleFrac=num_points/len(int_data)

    #-------------Code for plotting the graph of sub-sampled data------------

    x=[i[0] for i in cData]
    y=[i[1] for i in cData]
    # plt.hist2d(x,y,bins=64,norm=LogNorm())
    plt.xlim(0,maxH)
    plt.ylim(0,maxL)
    # plt.show()

    #------------Return the sub sampled data-------------------
    print("before subsampling = ",len(int_data))
    print("After Subsampling = ",len(cData)," sub sampling fraction = ",subSampleFrac)   
    return cData,x,y,num_points
