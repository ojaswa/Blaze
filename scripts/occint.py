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
from sklearn.cluster import DBSCAN
from sklearn import datasets
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import estimate_bandwidth
import math
from matplotlib.colors import LogNorm
import random
import operator
import sys
import hdbscan
import imp
import scipy.stats as stat
import time
from prettytable import PrettyTable
import pickle
from mean_shift_ import MeanShift
from itertools import cycle
from sklearn.cluster import KMeans
from scipy.spatial import ConvexHull
from scipy import interpolate
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import RandomForestClassifier

def wmean(x):
	num=0
	denom=0
	for a in x:
		num+=a[0]*a[1]
		denom+=a[1]
	return (num/denom)

def random_classify(cenetrwise_alloc,int_data,mat_data_reduce):
	train=[]
	label=[]
	for i in range(len(cenetrwise_alloc)):
		for j in range(len(cenetrwise_alloc[i])):
			train.append(cenetrwise_alloc[i][j][0])
			label.append(i)
	train=np.reshape(train,(-1,1))
	clf = RandomForestClassifier(n_jobs=-1, random_state=0)
	clf.fit(train, label)
	int_data=np.reshape(int_data,(-1,1))
	label=clf.predict(int_data)

	per=[0 for i in range(len(cenetrwise_alloc))]
	for i in range(len(int_data)):
		per[label[i]]+=1

	for i in range(len(per)):
		print("Material ",i,": mean = ",mat_data_reduce[i][0],", composition = ",per[i]/len(int_data)*100)

def occint_cluster():

	# ------------------- Setting up data ------------------

	t0=time.time()

	intn=open("./process/intensity_data.txt",'r')
	occl=open("./process/occlusion_data.txt",'r')
	materials=open("./process/materials.txt",'r')
	int_data=[]
	occ_data=[]
	mat_data=[]
	dim=[]
	dim=intn.readline().strip("\\n").split()
	dim=[int(i) for i in dim]
	print("input data dimensions = ",dim)
	x=materials.readline()
	mat_t=x.strip("\\n")
	mats=int(mat_t)
	x=materials.readline()
	tempdiv=x.strip("\\n")
	divider=float(tempdiv)
	kkii=0
	for lines in materials:
		mat_data.append(lines.split(','))
		mat_data[-1][0]=float(mat_data[-1][0])
		mat_data[-1][1]=float(mat_data[-1][1])
		mat_data[-1][2]=float(mat_data[-1][2])
		mat_data[-1].append(kkii)
		kkii+=1
	materials.close()
	edges=open("./process/material_graph.txt","r")
	vesizes_str=edges.readline().strip("\\n")
	vesizes_str=vesizes_str.split(',')
	vesizes=[int(i) for i in vesizes_str]
	gedges=[]
	for lines in edges:
		gedges.append(lines.split(","))
		gedges[-1][0]=float(gedges[-1][0])
		gedges[-1][1]=float(gedges[-1][1])
	edges.close()	
	for lines in intn:
		t=lines.strip("\\n")
		int_data.append(float(t))

	min_i=min(int_data)
	max_i=max(int_data)
	count=0
	for lines in occl:
		if(lines=="nan\n"):
			occ_data.append(0)
			count=count+1
			continue
		occ_data.append(float(lines.strip("\\n")))

	mxi=max(int_data)
	mni=min(int_data)
	mxo=max(occ_data)
	mno=min(occ_data)

	t1=time.time()
	
	output_log2=open("./process/outlog2.txt","w")
	print("max intensity= ",mxi," max occlusion= ",mxo," min intensity= ",mni," min occlusion= ",mno)
	print("Occlusion data loading time = ",t1-t0)
	output_log2.write("Occlusion data loading time = "+str(t1-t0)+"\n")
	
	# ------------- Linear Sub-sampling of data---------------------------

	t0=time.time()
		
	SS=imp.load_source('subSampling2',"./scripts/utilities/subSample2.py")
	coordData,ss_int,ss_occ,subSize=SS.subSampling2(int_data,occ_data)
		
	t1=time.time()
		
	print("Occlusion Intensity space sub sampling time = ",t1-t0)
	output_log2.write("Occlusion Intensity space sub sampling time = "+str(t1-t0)+"\n")

	m=0
	mat_data_reduce=mat_data

	# --------------- Mapping possible materials to clusters -----------	
	#cluster to material mapping based on in which material rectangle each point lies

	cenetrwise_alloc=[[] for i in range(len(mat_data_reduce))]
	for i in range(len(mat_data_reduce)):
		mu=mat_data_reduce[i][0]
		std=mat_data_reduce[i][1]
		upper=mu+2*std
		lower=mu-2*std
		for j in range(len(coordData)):
			xc=coordData[j][0]
			yc=coordData[j][1]
			# span=mat_data_reduce[j][2]/2
			if((xc<=upper and xc>lower)):
				cenetrwise_alloc[i].append(coordData[j])

	# ---------------------- Classifying the volume using random forest ------------
	random_classify(cenetrwise_alloc,int_data,mat_data_reduce)

	# -------------------- Colour and Alpha assignment -----------------
	#Color and alpha assignment based on occlusion value which is interpolated and scaled.

	
	mat_data_reduce2=mat_data_reduce
	colours_rgb=cm.plasma(np.linspace(0.1,0.9,len(mat_data_reduce2)))
	# random.seed(42)
	# indices=[]
	# indices=random.sample(range(0,len(colours_rgb)),len(mat_data_reduce2))
	alphas=[]
	occ_mat_map=[]
	size=0
	for i in range(len(mat_data_reduce2)):
		alpha_temp=[]
		y=[j[1] for j in cenetrwise_alloc[i]]
		a=np.asarray(y)
		if(a.size==0):
			occ_mat_map.append(1)
			alphas.append(int(1))
			continue
		size=size+a.size
		occ_std=np.std(a)
		occ_mean=np.mean(a)
		alpha_temp2=occ_mean
		alpha_temp.append([occ_mean,a.size])
		occ_interpolated=np.interp(np.mean(alpha_temp2),[mno,mxo],[0,0.15])
		occ_interpolated=int(round(occ_interpolated*255))
		occ_mat_map.append(wmean(alpha_temp))
		alphas.append(int(occ_interpolated))
	colours=[]
	# for i in indices:
	for i in range(len(colours_rgb)):
		temp=[]
		# colours_rgb[i][3]=alphas[i]	
		for j in range(len(colours_rgb[i])):
			t=int(round(colours_rgb[i][j]*256))
			if(t==256):
				t-=1
			temp.append(t)
		colours.append(temp)

	# ------------------ Background values ----------------

	gamma=0.7
	bkg_sel=[]
	material_sizes=[]
	for i in range(len(mat_data_reduce2)):
		# size=0
		# for j in range(len(cluster_sizes)):
		# 	if(clus_mat_map[j]==i):
		# 		size+=cluster_sizes[j]
		y=[j[1] for j in cenetrwise_alloc[i]]
		a=np.asarray(y)
		material_sizes.append(a.size)	
		bkg_sel.append((1/(a.size+1))*(gamma*occ_mat_map[i]+(1-gamma)*mat_data_reduce2[i][0]))	

	print("final number of materials= ",len(mat_data_reduce2))
	final_reduced_materials=PrettyTable()
	final_reduced_materials.add_column("Intensities",(np.mat(mat_data_reduce2)[:,0]).tolist())
	final_reduced_materials.add_column("Standard Deviation",(np.mat(mat_data_reduce2)[:,1]).tolist())
	final_reduced_materials.add_column("R colors",(np.mat(colours)[:,0]).tolist())
	final_reduced_materials.add_column("G colors",(np.mat(colours)[:,1]).tolist())
	final_reduced_materials.add_column("B colors",(np.mat(colours)[:,2]).tolist())
	final_reduced_materials.add_column("Alpha",alphas)
	final_reduced_materials.add_column("Bkg val",bkg_sel)
	final_reduced_materials.add_column("mat size",material_sizes)
	print(final_reduced_materials)
	output_log2.close()

	# ----------------- Material data output -----------------
	#Writing data about material graph and color assignment to appropriate files

	edgesredtemp=[]	
	gedges_red2=gedges
	for i in range(len(gedges_red2)):
		e1=e2=-1
		for j in range(len(mat_data_reduce2)):
			if(gedges_red2[i][0]==mat_data_reduce2[j][3]):
				e1=j
			if(gedges_red2[i][1]==mat_data_reduce2[j][3]):
				e2=j
		if(e1>=0 and e2>=0):
			edgesredtemp.append([e1,e2])
	gedges_red2=edgesredtemp					

	mat_red_file=open("./process/materials_reduced.txt","w")
	mat_str=str(len(mat_data_reduce2))+"\n"
	for i in range(len(mat_data_reduce2)):
		mat_str+=str(mat_data_reduce2[i][0])+","+str(mat_data_reduce2[i][1])+","+str(colours[i][0])+","+str(colours[i][1])+","+str(colours[i][2])+","+str(alphas[i])+","+str(bkg_sel[i])+"\n"
	
	mat_str+=str(np.min(ss_occ))+"\n"+str(np.max(ss_occ))+"\n"	
	mat_red_file.write(mat_str)
	mat_red_file.close()

	gedges_red_file=open("./process/edges_reduced.txt","w")
	gedges_str=str(len(mat_data_reduce2))+","+str(len(gedges_red2))+"\n"
	for i in range(len(gedges_red2)):
		gedges_str+=str(gedges_red2[i][0])+","+str(gedges_red2[i][1])+"\n"
	gedges_red_file.write(gedges_str)
	gedges_red_file.close()	

	# #undoing all the changes that were made to material data and edge data
	# mat_data_reduce2=mat_data_reduce
	# gedges_red2=gedges
