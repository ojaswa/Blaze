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
import matplotlib.cm as cm
import colorsys
# def remove_background(): # based on lowest intensity value
	

def define_colors(materials):	
	clus_n=len(materials)
	colours_hsv=cm.hsv(np.linspace(0, 1, clus_n))
	colours_rgb=[]
	for i in range(len(colours_hsv)):
		colours_rgb.append(colorsys.hsv_to_rgb(colours_hsv[i][0],colours_hsv[i][1],colours_hsv[i][2]))
	return colours_rgb

def define_alpha(occlusions,x):	
	minlim=np.min(occlusions)
	maxlim=np.max(occlusions)
	occ_interpolated=np.interp(x,[minlim,maxlim],[0,1])
	return occ_interpolated

def colalpha(materials,occlusions,centers,alloc):
	cols=define_colors(materials)
	count=0
	# alpha=define_alpha(occlusions,x)
	clalp=[]
	min_idx=-1
	min_val=float('inf')
	for i in range(len(materials)):
		if(materials[i][0]<min_val):
			min_val=materials[i][0]
			min_idx=i
	# color_alloc=[]
	# for i in range(len(centers)):

	print(cols)
	for i in range(len(centers)):
		print(define_alpha(occlusions,centers[i][1]))

	# for i in range(len(centers)):
	# 	print(i)
	# 	if(alloc[i]==materials[min_idx][0]):
	# 		clalp.append([cols[alloc[i]],0])
	# 	else:	
	# 		count+=1
	# 		clalp.append([cols[(np.where(materials[:,0]==alloc[i])[0])[0]],define_alpha(occlusions,centers[i][1])])

	print("colour alpha mappings ...")
	print(clalp)	
	print(count)
	# call all functions here
	# some post processing 
	# return them	
