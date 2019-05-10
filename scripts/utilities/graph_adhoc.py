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
import math
import random
import operator
import sys
import imp

BCD=imp.load_source('bhattacharyaDistUni',"./scripts/utilities/bhattacharyaDist.py")

def change_val(i,j,flag_merg,dist_holder,dm):
	flag_merg[i][j] = 1
	dist_holder[j] = [dm[i][j],i]
	return flag_merg,dist_holder

## of j is even oe = 1; else oe = -1	
def change_rem_noparent(i,j,flag_merg,dist_holder,dm):
	flag_merg[dist_holder[j][1]][j]=0
	flag_merg[i][j] = 1
	dist_holder[j] = [dm[i][j],i]
	return flag_merg,dist_holder

def w_mean(means,sizes):
	num_sum = 0 
	denom_sum = 0
	for i in range(len(means)):
		num_sum = num_sum + means[i]*sizes[i]
		denom_sum = denom_sum + sizes[i]
	return float(num_sum)/float(denom_sum)	

def w_std(stds,sizes):
	num_sum = 0 
	denom_sum = 0
	var = np.square(stds)
	sum_vars = np.sum(var)
	return sum_vars**(0.5)

def change_and_swap_parent(i,j,flag_merg,dist_holder,oe,dm):
	flag_merg[dist_holder[i][1]][j] = 1
	# flag_merg[i][j] = 1
	dist_holder[j] = [dm[i][j] , dist_holder[i][1]]
	flag_merg[i][j+oe] = 0
	dist_holder[j+oe] = [-1,j]
	return flag_merg,dist_holder 

def change_and_swap_noparent(i,j,flag_merg,dist_holder,oe,dm):
	flag_merg[i][j]=1
	flag_merg[i][j+oe] = 0
	dist_holder[j] = [dm[i][j],i]
	dist_holder[j+oe] = [-1,j]
	return flag_merg,dist_holder


def graph_finder(dm,distThreshold,meanhl,stdhl,classSize,max_int):
	toBeMerged= [[] for i in range(len(dm))]
	flag_merg = np.zeros((len(dm) , len(dm[0])))
	dist_holder = [[-1,i] for i in range(len(dm))]
	dm=np.array(dm)
	numIter = 1 ## variable
	changed = True
	k=0

	## finding the merging vertices

	while(changed and k<numIter):
		changed = False

		for i in range(dm.shape[0]):
			for j in range(i,dm.shape[1]):
				if(dm[i][j]==float("inf" ) or flag_merg[i,j]==1):
					continue
				elif(dm[i][j]<distThreshold):
					if(j%2==0):
						oe = 1
					else:
						oe = -1	

					if(dist_holder[j][0]==-1):
						if(flag_merg[i][j+oe]==0 and flag_merg[dist_holder[i][1]][j+oe]==0):
							if(dist_holder[i][0]==-1):	
								flag_merg,dist_holder = change_val(i,j,flag_merg,dist_holder,dm)
							else:
								flag_merg,dist_holder = change_val(dist_holder[i][1],j,flag_merg,dist_holder,dm)
								flag_merg[i][j]=1
							changed=True	
						else:
							if(dm[i][j+oe]<dm[i][j]):
								continue
							else:
								if(dist_holder[i][0]==-1):
									flag_merg,dist_holder = change_and_swap_noparent(i,j,flag_merg,dist_holder,oe,dm)
								else:
									flag_merg,dist_holder = change_and_swap_parent(i,j,flag_merg,dist_holder,oe,dm)
								changed=True	
					else:
						if(dist_holder[j][0]<=dm[i][j]):
							continue
						else:
							if(dist_holder[i][0]==-1):
								if(flag_merg[i][j+oe]==0):
									flag_merg,dist_holder = change_rem_noparent(i,j,flag_merg,dist_holder,dm)
									changed=True	
								elif(dm[i][j+oe]<dm[i][j]):
									continue	
								else:
									flag_merg,dist_holder = change_and_swap_noparent(i,j,flag_merg,dist_holder,oe,dm)
									changed=True		
							else:
								if(flag_merg[dist_holder[i][1]][j+oe] == 0 and flag_merg[i][j+oe]==0):
									flag_merg,dist_holder = change_rem_noparent(dist_holder[i][1],j,flag_merg,dist_holder,dm)
									changed=True	
								elif(dm[i][j+oe]<dm[i][j]):
									continue
								else:
									flag_merg,dist_holder = change_and_swap_parent(i,j,flag_merg,dist_holder,oe,dm)
									changed=True	
				else:
					continue

		k=k+1

	## merging step	

	num_vertices = 0
	v_map = []
	for i in range(len(dist_holder)):
		if(dist_holder[i][0]==-1):
			v_map.append([num_vertices,i])
			num_vertices=num_vertices+1
	
	vertice_set = [[] for i in range(num_vertices)]
	vertice_set_dev = [[] for i in range(num_vertices)]
	vertice_set_size = [[] for i in range(num_vertices)]
	vertice_mark = [0 for i in range(len(meanhl))]
	edge_set = set()
	for i in range(len(dist_holder)):
		for j in v_map:
			if(i==j[1] or dist_holder[i][1] == j[1]):
				vertice_set[j[0]].append(meanhl[i])
				vertice_set_dev[j[0]].append(stdhl[i])
				vertice_set_size[j[0]].append(classSize[i])
				vertice_mark[i] = j[0]

	# print("distThreshold = ",distThreshold)
	# print("meanhl \n",meanhl)
	# print("stdhl \n", stdhl)
	# print("dist holder\n",dist_holder)
	# print("verice mark \n",vertice_mark)			

	for i in range(0,len(meanhl),2):
		# if([vertice_mark[i],vertice_mark[i+1]] not in edge_set):
		edge_set.add((vertice_mark[i],vertice_mark[i+1]))


	vertices_combined = []
	vertices_combined_dev = []
	vertices_combined_size = []
	vertices_combined_span = []
	vertices_combined_index = []

	for i in range(num_vertices):
		if(len(vertice_set[i])>0):
			vertices_combined.append(w_mean(vertice_set[i],vertice_set_size[i]))
			vertices_combined_dev.append(w_std(vertice_set_dev[i],vertice_set_size[i]))
			std_temp = vertices_combined_dev[-1]
			span = std_temp*6
			vertices_combined_span.append(span)
			vertices_combined_index.append(i)
			vertices_combined_size.append(np.sum(vertice_set_size[i]))
	# print("vertice set before\n",vertices_combined)	
	# print("edge set before\n",edge_set)			

	## cross check with the current vertices to see if some mach or not (just a sanity check)
	check_mapping = -1*np.ones((num_vertices,1))
	for i in range(num_vertices):
		boundi_p = vertices_combined[i] + 2*vertices_combined_dev[i]
		boundi_n = vertices_combined[i] - 2*vertices_combined_dev[i]
		for j in range(i+1,num_vertices):
			boundj_p = vertices_combined[j] + 2*vertices_combined_dev[j]
			boundj_n = vertices_combined[j] - 2*vertices_combined_dev[j]
			if(check_mapping[i] == -1):
				if(boundj_n<=0 and boundi_n<=0):
					check_mapping[j] = i
				elif(boundj_p>=max_int and boundi_p>=max_int):
					check_mapping[j] = i
				elif(boundi_n <= vertices_combined[j] and vertices_combined[i]>=vertices_combined[j]):
					check_mapping[j] = i
				elif(boundi_p >= vertices_combined[j] and vertices_combined[i]<=vertices_combined[j]):
					check_mapping[j] = i
				else:
					continue
			else:
				if(boundj_n<=0 and boundi_n<=0):
					check_mapping[j] = check_mapping[i] 
				elif(boundj_p>=max_int and boundi_p>=max_int):
					check_mapping[j] = check_mapping[i]
				elif(boundi_n <= vertices_combined[j] and vertices_combined[i]>=vertices_combined[j]):
					check_mapping[j] = check_mapping[i]
				elif(boundi_p >= vertices_combined[j] and vertices_combined[i]<=vertices_combined[j]):
					check_mapping[j] = check_mapping[i]
				else:
					continue				
								
	vertice_backmap = np.zeros((num_vertices,1))
	num_vertices_af = len(np.where(check_mapping==-1)[0])
	vert_idxs_af = np.where(check_mapping==-1)[0]

	k=0
	for i in range(num_vertices):
		if(i in vert_idxs_af):
			verts_of_i = np.where(check_mapping == i)[0]
			if(len(verts_of_i)==0):
				vertice_backmap[i] = k
				k=k+1
				continue
			vertices_combined[i] = vertices_combined[i] * vertices_combined_size[i]	
			vertices_combined_dev[i] = vertices_combined_dev[i]**2
			for j in range(len(verts_of_i)):
				vertices_combined[i] = vertices_combined[i]+ vertices_combined[verts_of_i[j]] * vertices_combined_size[verts_of_i[j]]
				vertices_combined_dev[i] = vertices_combined_dev[i] + vertices_combined_dev[verts_of_i[j]]**2
				vertices_combined_size[i] = vertices_combined_size[i] + vertices_combined_size[verts_of_i[j]]
				vertice_backmap[i] = k
				vertice_backmap[verts_of_i[j]]=k
			vertices_combined[i] = vertices_combined[i]/vertices_combined_size[i]	
			vertices_combined_span[i] = vertices_combined_dev[i]*6
			k=k+1		

	idxs_notneed = np.where(check_mapping!=-1)[0][::-1]
	for i in idxs_notneed:
		del(vertices_combined[i])
		del(vertices_combined_size[i])
		del(vertices_combined_dev[i])
		del(vertices_combined_span[i])

	edge_set = list(edge_set)
	edge_set_new = set()
	for e in edge_set:
		edge_set_new.add((int(vertice_backmap[e[0]][0]),int(vertice_backmap[e[1]][0])))
	edge_set = edge_set_new
	edge_set = list(edge_set)
	for e in edge_set:
		e = list(e)
	for i in range(len(vertice_mark)):
		vertice_mark[i] = int(list(vertice_backmap[vertice_mark[i]])[0])		
	
	# print("vertice cobined\n",vertices_combined)
	# print("edge set\n",edge_set)

	## fraph energy calculation
	inter_energy=0
	for i in range(len(edge_set)):
		m1 = edge_set[i][0]
		m2 = edge_set[i][1]
		inter_energy+=(vertices_combined_size[m1]+vertices_combined_size[m2])*BCD.bhattacharyaDistUni(vertices_combined[m1],vertices_combined[m2],vertices_combined_dev[m1],vertices_combined_dev[m2],[1])

	inter_energy=inter_energy/len(edge_set)

	intra_energy=0
	for i in range(len(meanhl)):
		m1 = [meanhl[i],stdhl[i]]
		idx2 = int(vertice_backmap[vertice_mark[i]])
		m2 = [vertices_combined[idx2],vertices_combined_dev[idx2]]
		intra_energy+=(classSize[i]+classSize[idx2])*BCD.bhattacharyaDistUni(m1[0],m2[0],m1[1],m2[1],[1])
	intra_energy=intra_energy/(len(meanhl))	
	
	graph_energy=inter_energy/intra_energy

	return vertices_combined_size,vertices_combined,vertices_combined_dev,vertices_combined_span,edge_set,vertice_mark,graph_energy	



											






