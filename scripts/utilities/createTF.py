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

def retAlpha(x,meanhl,stdhl,discardValue):
    key=[(0,0),(0.03,0),(1,0)]
    intensity=meanhl[x]
    istd=stdhl[x]
    if(intensity<=discardValue):
        print("not a valid intensity point value")
        return '0'  
         
    key.append((intensity,0.6))
    bn=10
    
    if(intensity-istd>0.03):
        iminus=intensity-istd
        key.append((iminus,0.6))
    else:
        iminus=0.03
        key.append((iminus,0.6))
    if(intensity+istd<1):
        iplus=intensity+istd
        key.append((iplus,0.6))
    else:
        iplus=1   
        key.append((iplus,0.6))         

    if(iminus-0.01<=0.03):
        nminus=0.03
    else:
        nminus=iminus-0.01
    # nminus=(iminus+nminus)/2    
    if(iplus+0.01>1):
        nplus=1
    else:
        nplus=iplus+0.01    
    # nplus=(nplus+iplus)/2 
    key.append((nplus,0))
    key.append((nminus,0))
    # print("intensity="+ str(intensity) " iplus=" + str(iplus)+ " iminus="+str(iminus)+" nplus="+str(nplus)+" nminus="+str(nminus)+"\n")
    # for i in range(1,bn+1):
    #     a=((iminus-(i*(iminus-nminus)/bn)),(1-i*(1/float(bn))))
    #     b=((iplus+(i*(nplus-iplus)/bn)),(1-i*(1/float(bn))))
    #     key.append(a) 
    #     key.append(b)

    key.sort(key=operator.itemgetter(0))
    # print(key)
    keyval=([ a for a,b in key ], [ b for a,b in key ])
    # print(keyval)
    ky=keyval[0]
    val=keyval[1]    
    out=""
    for i in range(len(key)):
        out=out+"\t\t<Node Key=\""+str(ky[i])+"\" Value=\""+str(val[i])+"\"/>\n"
    # print(out)    
    return out
        
def retColor(meanhl):
    color = plot.cm.Spectral(np.linspace(0.1, 0.9, len(meanhl)))
    value=[]
    for i in color:
        tmp='#'
        tmp=tmp+str(hex(int(i[0]*255))[2:])
        tmp=tmp+str(hex(int(i[1]*255))[2:])
        tmp=tmp+str(hex(int(i[2]*255))[2:])
        value.append(tmp)
    out=""
    for i in range(len(meanhl)):
        out=out+"\t\t<Node Key=\""+str(meanhl[i])+"\" Value=\""+str(value[i])+"\"/>\n"
    # print(out)    
    return out    

def writeTfFile(dataheadFileName,x,meanhl,stdhl,discardValue):
    filename="test_"+dataheadFileName+"_"+str(x)+".tf1"
    wfile=open(filename,'w')
    alpha=retAlpha(x,meanhl,stdhl,discardValue)
    if(alpha=='0'):
        return
    colors=retColor(meanhl)
    toWrite="<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<TransferFunction1D>\n\t<Alpha>\n"+alpha+"\t</Alpha>\n\t<Color>\n"+colors+"\t</Color>\n</TransferFunction1D>"
    # print(toWrite)
    wfile.write(toWrite)
    return