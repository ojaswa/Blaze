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
from sklearn.cluster import KMeans
from sklearn.cluster import DBSCAN
from sklearn import datasets
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import MeanShift, estimate_bandwidth
from itertools import cycle
from sklearn.cluster import SpectralClustering
from sklearn.feature_extraction import image
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import AffinityPropagation
from sklearn.neighbors import kneighbors_graph
import math
import random

base=os.path.normpath(os.getcwd() + os.sep + os.pardir)
base=os.path.join(base,"Code")
base=os.path.join(base,"LH computation")
base=os.path.join(base,"data")
print(base)

headFileName="CThead_LH.hdr"
rawFileName="CThead_LH.bin"

head=open(os.path.join(base,headFileName),"r")
file=open(os.path.join(base,rawFileName),"rb")

hAr=[]
bins=[[[] for i in range(256)] for j in range(256)]
binX=[]
binY=[]
byteSize=2 #has to be manually set as per the use case
limiter=256**byteSize
divider=256**(byteSize-1)
granularity=1
subSampleFrac=0.003

for line in head:
    hAr.append(line)

length=int(hAr[0].split(" ")[2])
width=int(hAr[0].split(" ")[3])
height=int(hAr[0].split(" ")[4])

data_array = np.fromfile(file, np.uint16).reshape((-1, 2)).T
totalPts=len(data_array[0])
print(totalPts)
print(data_array)
maxH=max(data_array[0])
minH=min(data_array[0])
maxL=max(data_array[1])
minL=min(data_array[1])
divH=((maxH)/(divider-1))*granularity
divH=divH+1
divL=((maxL)/(divider-1))*granularity
divL=divL+1
print(divH)
print(divL)
# sum=0
for i in range(len(data_array[0])):
	if(data_array[0][i]<data_array[1][i]):
		t=[]
		t.append(data_array[0][i])
		t.append(data_array[1][i])
		bins[int(data_array[0][i]/divH)][int(data_array[1][i]/divL)].append(t)
		# sum=sum+1

coordData=[]    

# print(sum)
    
for j in range(256):
    for i in range(256):
        if(len(bins[i][j])>0):
            nmIdxs=int(math.floor(len(bins[i][j])*subSampleFrac))
            sampleIdxs=random.sample(range(0,len(bins[i][j])-1),nmIdxs)
            for k in sampleIdxs:
            	coordData.append(bins[i][j][k])
            	binX.append(bins[i][j][k][0])
            	binY.append(bins[i][j][k][1])	

print(len(coordData))            		
plot.figure(figsize=(9,9))

def plotData(binX,binY):			
	plot.subplot(331)
	plot.title("Sub-Sampled Data")
	plot.plot(binX,binY,'o',markersize=0.5)
	plot.xlabel("L")
	plot.ylabel("H")
	
	return

def kMeans(coordData):
	coordData = StandardScaler().fit_transform(coordData)
	k_means = KMeans(n_clusters=4)
	k_means.fit(coordData)
	labels=k_means.labels_
	n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
	unique_labels = set(labels)
	colors = plot.cm.Spectral(np.linspace(0, 1, len(unique_labels)))	
	plot.subplot(332)
	for i in range(len(coordData)):
		k=labels[i]
		if(labels[i]==-1):
			col=colors[-1]
		else:
			col=colors[k]	
		plot.plot(coordData[i][0], coordData[i][1], 'o', color=col,markersize=0.5)	 

	plot.title('Kmeans clusters: %d' % n_clusters)
	return

def dbScan(coordData):
	coordData = StandardScaler().fit_transform(coordData)
	# print(coordData)
	db = DBSCAN(eps=0.25, min_samples=5).fit(coordData)
	print(db)
	core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
	core_samples_mask[db.core_sample_indices_] = True
	labels = db.labels_
	# print(labels)	
	plot.subplot(333)
	n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
	unique_labels = set(labels)
	colors = plot.cm.Spectral(np.linspace(0, 1, len(unique_labels)))
	for k, col in zip(unique_labels, colors):
	    if k == -1:
	        col = 'k'

	    class_member_mask = (labels == k)
	    xy = coordData[class_member_mask & core_samples_mask]
	    plot.plot(xy[:, 0], xy[:, 1], 'o', color=col,markersize=0.5)

	    xy = coordData[class_member_mask & ~core_samples_mask]
	    plot.plot(xy[:, 0], xy[:, 1], 'o', color=col,markersize=0.5)

	plot.title('dbScan clusters: %d' % n_clusters_)
	
	return

def meanShift(coordData):
	coordData = StandardScaler().fit_transform(coordData)	
	bandwidth = estimate_bandwidth(coordData, quantile=0.2)
	ms = MeanShift(bandwidth=bandwidth,bin_seeding=True)
	ms.fit(coordData)
	labels = ms.labels_
	cluster_centers = ms.cluster_centers_
	labels_unique = np.unique(labels)
	n_clusters_ = len(labels_unique)
	plot.subplot(334)
	plot.title("Mean Shift, cluster=%d" % n_clusters_)
	# colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
	colors = plot.cm.Spectral(np.linspace(0, 1, len(labels_unique)))	
	print(colors)
	for i in range(len(coordData)):
		plot.plot(coordData[i][0],coordData[i][1],'o',color=colors[labels[i]],markersize=0.5)
	   
	return

def specClus(coordData):
	coordData = StandardScaler().fit_transform(coordData)
	specCluster = SpectralClustering(eigen_solver='arpack', affinity="nearest_neighbors")
	specCluster.fit(coordData)
	labels=specCluster.labels_
	print(len(labels))
	labelSet=np.unique(labels)
	print(len(labelSet))
	plot.subplot(335)
	plot.title("Spectral Clustering, cluster=%d" % len(labelSet))
	colors = plot.cm.Spectral(np.linspace(0, 1, len(labelSet)))
	for i in range(len(coordData)):
		plot.plot(coordData[i][0],coordData[i][1],'o',color=colors[labels[i]],markersize=0.5)
	return

def agglomAverage(coordData):
	coordData = StandardScaler().fit_transform(coordData)
	connectivity = kneighbors_graph(coordData, n_neighbors=10, include_self=False)
	connectivity = 0.5 * (connectivity + connectivity.T)
	average_linkage = AgglomerativeClustering(linkage="average", affinity="cityblock", connectivity=connectivity)
	average_linkage.fit(coordData)
	labels=average_linkage.labels_
	print(len(labels))
	labelSet=np.unique(labels)
	print(len(labelSet))
	plot.subplot(336)
	plot.title("agglomerative Cluster cluster=%d" %len(labelSet))
	colors = plot.cm.Spectral(np.linspace(0, 1, len(labelSet)))
	for i in range(len(coordData)):
		plot.plot(coordData[i][0],coordData[i][1],'o',color=colors[labels[i]],markersize=0.5)
	return	
	 

def hierarchicalWithConnectivity(coordData):
	coordData = StandardScaler().fit_transform(coordData)
	connectivity = kneighbors_graph(coordData, n_neighbors=10, include_self=False)
	connectivity = 0.5 * (connectivity + connectivity.T)
	ward = AgglomerativeClustering(linkage="ward", connectivity=connectivity)
	ward.fit(coordData)
	labels=ward.labels_
	print(len(labels))
	labelSet=np.unique(labels)
	print(len(labelSet))
	plot.subplot(337)
	plot.title("hierarchical with connectivity cluster=%d" %len(labelSet))
	colors = plot.cm.Spectral(np.linspace(0, 1, len(labelSet)))
	for i in range(len(coordData)):
		plot.plot(coordData[i][0],coordData[i][1],'o',color=colors[labels[i]],markersize=0.5)
	 
	return

def hierarchicalWithoutConnectivity(coordData):
	coordData = StandardScaler().fit_transform(coordData)
	ward = AgglomerativeClustering(linkage="ward")
	ward.fit(coordData)
	labels=ward.labels_
	print(len(labels))
	labelSet=np.unique(labels)
	print(len(labelSet))
	plot.subplot(338)
	plot.title("hierarchical without connectivity cluster=%d" %len(labelSet))
	colors = plot.cm.Spectral(np.linspace(0, 1, len(labelSet)))
	for i in range(len(coordData)):
		plot.plot(coordData[i][0],coordData[i][1],'o',color=colors[labels[i]],markersize=0.5)
	return

def affinityPropogation(coordData):
	coordData = StandardScaler().fit_transform(coordData)
	affinity_propagation =AffinityPropagation(damping=.9, preference=-200)
	affinity_propagation.fit(coordData)
	labels=affinity_propagation.labels_
	labelSet=np.unique(labels)
	print(len(labelSet))
	plot.subplot(339)
	plot.title("affinity Propogation cluster=%d" %len(labelSet))
	colors = plot.cm.Spectral(np.linspace(0, 1, len(labelSet)))
	for i in range(len(coordData)):
		plot.plot(coordData[i][0],coordData[i][1],'o',color=colors[labels[i]],markersize=0.5) 
	return


kMeans(coordData)
dbScan(coordData)
meanShift(coordData)
specClus(coordData)
agglomAverage(coordData)
hierarchicalWithConnectivity(coordData)
hierarchicalWithoutConnectivity(coordData)
affinityPropogation(coordData)
plotData(binX,binY)	
plot.show()