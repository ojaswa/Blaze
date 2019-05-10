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
import math

# --------------- for calculating bhattacharya distance using the integral formula on a ---------
# ------------ 2D expanse matrix of the clusters given their Sigma matrix -------------

def bhattacharyaDistMul(xMu,yMu,xSigmadet,xSigmainv,ySigmadet,ySigmainv,coordData):
    sum=0.0
    for ztemp in coordData:
        p=(1/(2*math.pi*(xSigmadet**(1/2))))*np.exp((-1/2)*(ztmp-xMu).transpose()*xSigmainv*(ztmp-xMu))
        q=(1/(2*math.pi*(ySigmadet**(1/2))))*np.exp((-1/2)*(ztmp-yMu).transpose()*ySigmainv*(ztmp-yMu))
        r=(p*q)**(1/2)
        sum=sum+r
    ans=-1*(np.log(sum))
    return ans    

# ------------------ for calculating bhattacharya distance using the integral formula ------------
# ------ computationally expensive O(N) --------- 

# def bhattacharyaDistUni(xMu,yMu,xStd,yStd,coordData):
#     sum=0.0
#     for k in range(len(coordData)):
#         p=(1/((2*math.pi)**(1/2)*(xStd**(1/2))))*np.exp((-1/2)*(((coordData[k][0]-xMu)**2)/(xStd**2)))
#         q=(1/((2*math.pi)**(1/2)*(yStd**(1/2))))*np.exp((-1/2)*(((coordData[k][0]-yMu)**2)/(yStd**2)))
#         r=(p*q)**(1/2)
#         sum=sum+r
#         p2=(1/((2*math.pi)**(1/2)*(xStd**(1/2))))*np.exp((-1/2)*(((coordData[k][1]-xMu)**2)/(xStd**2)))
#         q2=(1/((2*math.pi)**(1/2)*(yStd**(1/2))))*np.exp((-1/2)*(((coordData[k][1]-yMu)**2)/(yStd**2)))
#         r2=(p2*q2)**(1/2)
#         sum=sum+r2
#     print(sum) 
#     ans=round(-1*(np.log(sum)),2)
#     return ans    

# ------------------ for calculating bhattacharya distance using simple O(1) formula that is the closed form ------------
# ------ not computationally expensive --------- 

def bhattacharyaDistUni(xMu,yMu,xStd,yStd,coordData):
    sum=0.0
    pa=0
    pb=0
    if(xStd>0 and yStd>0 ):#and (xMu>0 or yMu>0)):
        pa=(0.25)*np.log((0.25)*(((xStd**2)/(yStd**2))+((yStd**2)/(xStd**2))+2))
        pb=(0.25)*(((xMu-yMu)**2)/((xStd**2)+(yStd**2)))

    # -------------- for cases where standard deviation or mean is 0 for any or both materials ----------

    else:
        if(yStd==0 or xStd==0):
            if(yStd==0 and xStd>0 ):
                ## we consider the distance is independent of yStd
                pa=0#(0.25)*np.log((0.25)*(((xStd**2)/(yStd**2))+((yStd**2)/(xStd**2))+2))
                pb=(0.25)*(((xMu-yMu)**2)/((xStd**2)+(yStd**2)))
            elif(yStd==0 and xStd==0):
                ## considering simply squared distance between two points
                pa=0
                pb=0.25*(abs(xMu-yMu)**2)
            elif(yStd>0 and xStd==0):
                ## we consider the distance is independent of xStd
                pa=0
                pb=(0.25)*(((xMu-yMu)**2)/((xStd**2)+(yStd**2)))
        elif(xMu==0 and yMu==0):
            pa=0
            pb=0       
    sum=pa+pb
    return sum

# def bhattacharya_bound(xMu,yMu,xStd,yStd):
    
# def kullback_liebler()
    

