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
from matplotlib.colors import LogNorm
from sklearn.cluster import MeanShift, estimate_bandwidth
import math
import random
import operator
import sys
import imp
from multiprocessing import Pool

# ---------------- function to subsample the LH space ----------------- 

def subSampling(actualData,data_array,volumeType):
    
    #----------------Uniformly sampling 100000 points from the data------------------

    random.seed(0)
    coordData=[]
    for i in range(len(data_array[0])):
        if(data_array[0][i]<data_array[1][i]):
            coordData.append([data_array[0][i],data_array[1][i]])
    num_points=int(0.1*len(coordData))
    if(num_points>100000):
        num_points=100000
    sample_points=random.sample(range(len(coordData)),num_points)
    cData=[coordData[i] for i in sample_points]
    maxH=np.max(data_array[0])
    minH=np.min(data_array[0])
    maxL=np.max(data_array[1])
    minL=np.min(data_array[1])
    subSampleFrac=num_points/len(data_array[0])

    #-------------Code for plotting the graph of sub-sampled data------------

    x=[i[0] for i in cData]
    y=[i[1] for i in cData]
    plt.hist2d(x,y,bins=64,norm=LogNorm())
    plt.xlim(0,maxH)
    plt.ylim(0,maxL)
    # plt.show()

    #-------------Return the sampled data long with maximum L and maximum H values---------

    print("Data length beofre sub-sampling = ",len(data_array[0]))
    print("Data length after sub sampling = ",len(cData), " sub sampling fraction = ",subSampleFrac) 
    return cData,maxH,maxL,minH,minL,num_points
