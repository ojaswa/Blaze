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
import os
from collections import namedtuple
import re
import array
import math

base=os.path.normpath(os.getcwd() + os.sep + os.pardir)
base=os.path.join(base,"data")
base=os.path.join(base,"NRRD")

headFileName="occlusion_test.nhdr"
rawFileName="occlusion_test.raw"
head=open(os.path.join(base,headFileName),"r")
file=open(os.path.join(base,rawFileName),"rb")

actualData=np.fromfile(file, np.uint16)
maxi=np.max(actualData)
mini=np.min(actualData)
span=(maxi*2)/255

noise_ratio=[0.1,0.15,0.2,0.25,0.3]

for i in range(len(noise_ratio)):
	print(noise_ratio[i])
	temp=np.copy(actualData)
	for j in range(temp.shape[0]):
		if(np.random.random()<noise_ratio[i]):
			noise=np.random.normal(0,span)
			temp[j]+=noise

	tempfilename="../data/NRRD/testing/"+str(noise_ratio[i])+"_"+rawFileName
	tempnhdrfilename="../data/NRRD/testing/"+str(noise_ratio[i])+"_"+headFileName
	os.makedirs(os.path.dirname(tempfilename), exist_ok=True)
	tempfile=open(tempfilename,"wb")
	os.makedirs(os.path.dirname(tempnhdrfilename), exist_ok=True)
	tempnhdrfile=open(tempnhdrfilename,"w")
	tempnhdrfilestr=""
	ln=0
	head=open(os.path.join(base,headFileName),"r")
	for line in head:
		if(ln==6):
			tempnhdrfilestr+=str("data file: ./"+str(noise_ratio[i])+"_"+rawFileName+"\n")#.write("data file: ./"+str(noise_ratio[i])+"_"+rawFileName+"\n")
		else:
			tempnhdrfilestr+=line#.write(line)	
		ln+=1	
	print(tempnhdrfilestr)
	tempnhdrfile.write(tempnhdrfilestr)
	temp.tofile(tempfile)
	tempfile.close()
	tempnhdrfile.close()
	
head.close()
file.close()		



