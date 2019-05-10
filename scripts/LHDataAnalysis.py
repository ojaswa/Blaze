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
import random
import operator
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize
import sys
import imp
# from scipy import stats
import scipy.stats as stat
import time
from prettytable import PrettyTable
import hdbscan
from itertools import cycle
import pickle
from mean_shift_ import MeanShift
import sys
from sklearn.cluster import KMeans
from sklearn.naive_bayes import GaussianNB
import shutil

t0=time.time()

# ------------ Dynamically importing scripts ---------

SSample=imp.load_source('subSampling',"./scripts/utilities/subSample.py")
DM=imp.load_source('distanceMatrix',"./scripts/utilities/distanceMatrix.py")
GSE=imp.load_source('multiCoreEdge',"./scripts/multiCore/gaussianSampling.py")
GSM=imp.load_source('multiCoreMaterial',"./scripts/multiCore/gaussianSampling.py")
BCDU=imp.load_source('bhattacharyaDistUni',"./scripts/utilities/bhattacharyaDist.py")
SPY=imp.load_source('seprability',"./scripts/utilities/seprtability.py")
GAC=imp.load_source('graph_finder',"./scripts/utilities/graph_adhoc.py")
MSD=imp.load_source('mean_and_std',"./scripts/mean_and_std.py")
print("Running material graph detection script ...")

# --------------------- Debugging toggle ------------------------
# 0 = don't skip; 1 = skip

arguments=sys.argv
mcs_value=int(arguments[1]) #Set minimum cluster size for HDBScan manually, -1 to use default parameters	
ms_value=int(arguments[2]) #Set minimum sample size fro HDBScan manually, -1 to use default parameters
restore_volume=int(arguments[3]) #Restoration of volume from results directory, 1 implies perform restore 
run_only_algo=int(arguments[4])

# --------------- opening LH data and volume intensity files ----------

base=os.path.normpath(os.getcwd())
base=os.path.join(base,"process")
headFileName="tmp_LH.hdr"
rawFileName="tmp_LH.bin"
head=open(os.path.join(base,headFileName),"r")
file=open(os.path.join(base,rawFileName),"rb")
fileNameFile=open(os.path.join(base,"path.txt"),"r")
filePathData=[]
for line in fileNameFile:
    filePathData.append(line)
volumePath=filePathData[0]
volumeType=filePathData[1]
volumePath = volumePath.replace("\n", "");
datahead=open(os.path.normpath(volumePath),"r")
volumePath = volumePath.replace(".nhdr", "")
volumePath = volumePath.replace(".nrrd", "")
volumePath=volumePath+".raw"   
datafile=open(volumePath,"rb")

# ------------------------ logging file -------------------------

output_log=open("./process/outlog1.txt","w")

# ------------- Reading LH data and intensity files -------------

datahAr=[]
for line in datahead:
    datahAr.append(line)
datalength=int(datahAr[4].split(" ")[1])
datawidth=int(datahAr[4].split(" ")[2])
dataheight=int(datahAr[4].split(" ")[3])
byteSize=2
if(volumeType=="1"):
    actualData=np.fromfile(datafile, np.uint8)
    byteSize=1
else:
    actualData=np.fromfile(datafile, np.uint16)
    byteSize=2
hAr=[]
limiter=256**byteSize
divider=256**(byteSize-1)
granularity=1
for line in head:
    hAr.append(line)
length=int(hAr[0].split(" ")[2])
width=int(hAr[0].split(" ")[3])
height=int(hAr[0].split(" ")[4])
if(volumeType=="1"):
    data_array = np.fromfile(file, np.uint8).reshape((-1, 2)).T
else:
    data_array = np.fromfile(file, np.uint16).reshape((-1, 2)).T   
totalPts=data_array.shape[1]
print("total LH points = ",totalPts)

t1=time.time()

print("Data Loading time = ",t1-t0)
output_log.write("Data Loading time = "+ str(t1-t0) +"\n")

# ----------------Clustering Algorithm---------------------------------------
#Option between meanshift and HDBScan

if(restore_volume==1):
    base_new=os.path.normpath(os.getcwd())
    results=os.path.join(base_new,"results")
    process=os.path.join(base_new,"process")
    for item in os.listdir(results):
        s=os.path.join(process,item)
        d=os.path.join(results,item)
        shutil.copy2(d,s)
	
if(restore_volume==0): 
    
    t0=time.time()

    coordData,maxH,maxL,minH,minL,subSize=SSample.subSampling(actualData,data_array,volumeType)
    
    t1=time.time()
    
    print("Sub sampling time = ",t1-t0)
    output_log.write("Sub sampling time = "+ str(t1-t0) +"\n")  

    def meanShift(coData,quantileValue,min_bin):
        transform=StandardScaler().fit(coData)
        cData = transform.transform(coData)
        c2Data= transform.inverse_transform(cData)
        print("Quantile Value = ",quantileValue, "min bin size = ",min_bin_freq)
        bandwidth = estimate_bandwidth(c2Data, quantile=quantileValue,n_jobs=-1)
        ms = MeanShift(bandwidth=bandwidth,bin_seeding=True,n_jobs=-1,min_bin_freq=min_bin,cluster_all=True)
        print(ms)
        output_log.write(str(ms))
        ms.fit(c2Data)
        labels = ms.labels_
        centers=ms.cluster_centers_ 
        return labels,c2Data

    def hdb(coData,subSize,mcs_value,ms_value):
        transform = StandardScaler().fit(coData)
        cData = transform.transform(coData)
        c2 = transform.inverse_transform(cData)
        if(mcs_value==-1):
            mcs=int(subSize*0.015)
        else:
            mcs=mcs_value
        if(ms_value==-1):
            ms=int(subSize*0.025)
        else:
            ms=ms_value
        clusterer=hdbscan.HDBSCAN(min_cluster_size=mcs,min_samples=ms,core_dist_n_jobs=1,cluster_selection_method='leaf',alpha=1.0,allow_single_cluster=True)
        labels = clusterer.fit_predict(c2)
        print(clusterer)
        labels_unique = np.unique(labels)
        n_clusters_ = len(labels_unique)
        colors = cycle('grmykgrmykbgrmykgrmyk')
        for k, col in zip(range(n_clusters_), colors):
            my_members = labels == k
            plt.plot(c2[my_members, 0], c2[my_members, 1], col + '.',alpha=0.1)
        plt.title('cluster_graph.py')
        # plt.show()
        return labels,c2

    print("running algorithm...")
    
    t0=time.time()
    
    usehdb=True
    if(usehdb):
        labels,c2=hdb(coordData,subSize,mcs_value,ms_value)
        labels_unique = np.unique(labels)
        n_clusters_ = len(labels_unique)
    else:
        quantileValue=float(0.01)
        min_bin_freq=3
        labels,c2=meanShift(coordData,quantileValue,min_bin_freq)
        labels_unique = np.unique(labels)
        n_clusters_ = len(labels_unique)
    
    t1=time.time()
    
    print("Algorithm running time = ",t1-t0," clusters = ",np.unique(labels))
    output_log.write("Algorithm running time = "+ str(t1-t0) +"\n") 
    unique_labels=np.unique(labels)

    labelsn=unique_labels.shape[0]
    labelwisedata=[]
    labelwisel=[]
    labelwiseh=[]
    
    t0=time.time()
    
    accThr=0.0001*len(coordData)
    coordData=np.matrix(coordData).T
    outliers_L=((coordData.T[np.where(labels==-1)]).T[1]).tolist()
    outliers_H=((coordData.T[np.where(labels==-1)]).T[0]).tolist()
    for i in range(labelsn):
        indexes=np.where(labels==i)
        if(indexes[0].shape[0]==0 or indexes[0].shape[0]<=accThr):
            continue
        labelwisedata.append(coordData.T[indexes])
        labelwiseh.append(((coordData.T[indexes]).T[0]).tolist())
        labelwisel.append(((coordData.T[indexes]).T[1]).tolist())
    
    #-------------------Calculating mean and standard deviation of clusters-------------------------------------------
    
    mx=max(maxH,maxL)
    meanhl=[]
    stdhl=[]
    labelwisemeanh=[]
    labelwisemeanl=[]
    xs=[]
    classSize=[]

    for i in range(len(labelwisedata)):
        if((np.std(labelwiseh[i])>=0) and (np.std(labelwisel[i])>=0)): #and len(labelwiseh[i])>=accThr and len(labelwisel[i])>=accThr):
            labelwisemeanh.append(np.mean(labelwiseh[i]))
            labelwisemeanl.append(np.mean(labelwisel[i]))
            if(np.std(labelwiseh[i])==0):
                stdhl.append(0.1)
            else:
                stdhl.append(np.std(labelwiseh[i]))
            meanhl.append(np.mean(labelwiseh[i])) 
            if(np.std(labelwisel[i])==0):
                stdhl.append(0.1)
            else:
                stdhl.append(np.std(labelwisel[i]))
            meanhl.append(np.mean(labelwisel[i]))     
            classSize.append(len(labelwiseh[i][0]))
            classSize.append(len(labelwisel[i][0]))
            xs.append(i)
# ---------------------------- Dumping the data to appropriate files ----------------------
    
    pkldmpfilew=open("./process/hlssclus.p","wb")
    pickle.dump(meanhl,pkldmpfilew)
    pickle.dump(stdhl,pkldmpfilew)
    pickle.dump(classSize,pkldmpfilew)  
    pickle.dump(coordData,pkldmpfilew)
    pickle.dump(labels,pkldmpfilew)   
    pkldmpfilew.close()   

# ---------------------------- Running material deduction algorithm ----------------------

if(run_only_algo==0):

    pkldmpfiler=open("./process/hlssclus.p","rb")
    meanhl=pickle.load(pkldmpfiler) 
    stdhl=pickle.load(pkldmpfiler)
    classSize=pickle.load(pkldmpfiler)
    print("class size total =", np.sum(classSize))
    coordData=pickle.load(pkldmpfiler)
    labels=pickle.load(pkldmpfiler)
    pkldmpfiler.close()
    print("Distance matrix length ",len(meanhl))
    dm=DM.distanceMatrix(stdhl,meanhl)

    edgeDist=[]
    for i in range(0,len(meanhl),2):
        edgeDist.append(BCDU.bhattacharyaDistUni(meanhl[i],meanhl[i+1],stdhl[i],stdhl[i+1],coordData))
    distThreshold = min(edgeDist)
    maxDist = max(edgeDist)
    if(distThreshold==0):
        distThreshold=1.0
    print("base Distance Threshold = ",distThreshold)
    distThrArray=[]
    denoiseThrArray=[]
    
    for i in range(-100,100,1):
        distThrArray.append(distThreshold*(10**(i/20)))
    
    # for i in range(9):
    #     denoiseThrArray.append(distThreshold+(i+1)*(maxDist - distThreshold)/50)    
        
    max_int=np.max(actualData)
    graphs=[]
    sep_mean_w=[]
    l=0
    # for den in denoiseThrArray:

    #     meanhl_denoised=[]
    #     stdhl_denoised=[]
    #     classsize_denoised=[]
    #     distanceMatrix_denoised=[]

    #     for i in range(len(edgeDist)):
    #         if(edgeDist[i]<den):
    #             continue
    #         else:
    #             meanhl_denoised.append(meanhl[2*i])
    #             meanhl_denoised.append(meanhl[2*i+1])
    #             stdhl_denoised.append(stdhl[2*i])
    #             stdhl_denoised.append(stdhl[2*i+1])
    #             classsize_denoised.append(classSize[2*i])
    #             classsize_denoised.append(classSize[2*i+1])

        # print(len(meanhl_denoised) )
        # print(len(stdhl_denoised) )
        # print(len(classsize_denoised) )

        # distanceMatrix_denoised = DM.distanceMatrix(stdhl_denoised,meanhl_denoised)        

    for dt in distThrArray:
        l=l+1
        # verticeClassSize,verticeSetAlpha,verticeSetAlphaDeviation,verticeSetSpan,edgeSetAlpha,verticeMarkup,graphEnergy=GAC.graph_finder(distanceMatrix_denoised ,dt,meanhl_denoised ,stdhl_denoised ,classsize_denoised ,max_int)
        verticeClassSize,verticeSetAlpha,verticeSetAlphaDeviation,verticeSetSpan,edgeSetAlpha,verticeMarkup,graphEnergy=GAC.graph_finder(dm ,dt,meanhl ,stdhl ,classSize ,max_int)
        z=[]
        if(len(edgeSetAlpha)==0 or (len(verticeSetAlpha)-1>len(edgeSetAlpha))):
            continue
        z.append(len(verticeSetAlpha))
        z.append(verticeSetAlpha)
        z.append(edgeSetAlpha)
        z.append(verticeSetAlphaDeviation)
        z.append(verticeClassSize)
        z.append(verticeSetSpan)
        z.append(verticeMarkup)
        z.append(dt)
        z.append(dt)
        sep_mean_w.append(graphEnergy)
        graphs.append(z)

    mingr=float("inf")
    mingridx=0

    grx=[]
    gry=[]

    for g in range(len(graphs)):
        if(graphs[g][0]<mingr):
            mingr=graphs[g][0]
            mingridx=g
        grx.append(g+1)
        gry.append(graphs[g][0])

    max=-1
    # usewsep=0
    # if(usewsep==0):
    for i in range(len(sep_mean_w)):
        if(sep_mean_w[i]>=max):
            max=sep_mean_w[i]
            maxGraphSepIdx=i    

    pklgrphfw=open("./process/graphsf.p","wb")
    pickle.dump(graphs,pklgrphfw) 
    pickle.dump(sep_mean_w,pklgrphfw)
    # pickle.dump(sep_mean_nw,pklgrphfw)
    pklgrphfw.close()               

    print("Selected graph: vertices= "+str(graphs[maxGraphSepIdx][0])+" edges= "+str(len(graphs[maxGraphSepIdx][2]))+" theta= "+str(graphs[maxGraphSepIdx][-1]))
    verticeSetAlpha=graphs[maxGraphSepIdx][1]
    edgeSetAlpha=graphs[maxGraphSepIdx][2]
    verticeSetAlphaDeviation=graphs[maxGraphSepIdx][3]
    verticeClassSize=graphs[maxGraphSepIdx][4]
    verticeSetSpanAlpha=graphs[maxGraphSepIdx][5]

    t1=time.time()
    print("time to run deduction algorithm (including DM matrix, deduction)",t1-t0)
    output_log.write("time to run deduction algorithm (including DM matrix, deduction) = "+ str(t1-t0) +"\n")  
    print("writing to disk ...")

    material_file_out=open(os.path.join(base,"materials.txt"),"w")
    material_out_file_str=str(len(verticeSetAlpha))+"\n"+str(divider)+"\n"
    for i in range(len(verticeSetAlpha)):
        material_out_file_str=material_out_file_str+str(verticeSetAlpha[i])+","+str(verticeSetAlphaDeviation[i])+","+str(verticeSetSpanAlpha[i])+"\n"
    material_file_out.write(material_out_file_str)
    material_file_out.close()

    edge_file_out=open(os.path.join(base,"material_graph.txt"),"w")
    edge_file_out_str=str(len(verticeSetAlpha))+","+str(len(edgeSetAlpha))+"\n"
    for i in range(len(edgeSetAlpha)):
        edge_file_out_str=edge_file_out_str+str(edgeSetAlpha[i][0])+","+str(edgeSetAlpha[i][1])+"\n"
    edge_file_out.write(edge_file_out_str)
    edge_file_out.close()

    lh_file_out=open(os.path.join(base,"lh_graph.txt"),"w")
    lh_file_out_str=str(len(meanhl))+"\n"
    for i in range(0,len(meanhl),2):
        lh_file_out_str=lh_file_out_str+str(meanhl[i])+","+str(stdhl[i])+","+str(meanhl[i+1])+","+str(stdhl[i+1])+"\n"
    lh_file_out.write(lh_file_out_str)
    lh_file_out.close()
    print("written to disk.")
    print("number of materials = ",len(verticeSetAlpha))  
    print("Final materials after first material deduction run ... ")
    final_table_part1=PrettyTable()
    final_table_part1.add_column("Intensity",verticeSetAlpha)
    final_table_part1.add_column("standard deiation",verticeSetAlphaDeviation)
    final_table_part1.add_column("class size",verticeClassSize)
    final_table_part1.add_column("class span",verticeSetSpanAlpha)
    print(final_table_part1)

    esepth=open("./process/esepth.txt","w")
    esepth.write(str(len(graphs)) +"\n")
    for i in range(len(graphs)):
        esepth.write(str(graphs[i][7])+","+str(graphs[i][8])+","+str(graphs[i][0])+","+str(len(graphs[i][2]))+","+str(sep_mean_w[i])+ "\n")
        # ,"+str(sep_mean_nw[i])+","+str(sep_mean_wfisc[i])+"\n")
    esepth.close() 

    graphs_all=open("./process/graphs.txt","w")
    graphs_all.write(str(len(graphs))+"\n")
    for i in range(len(graphs)):
        graphs_all.write(str(graphs[i][0])+","+str(graphs[i][6])+"\n")
        for j in range(graphs[i][0]):
            if(j==0):
                graphs_all.write(str(graphs[i][1][j]))
            else:
                graphs_all.write(","+str(graphs[i][1][j])) 
        graphs_all.write("\n")        
        for j in range(graphs[i][0]):
            if(j==0):
                graphs_all.write(str(graphs[i][3][j]))
            else:
                graphs_all.write(","+str(graphs[i][3][j])) 
        graphs_all.write("\n")              
        for j in range(graphs[i][0]):
            if(j==0):
                graphs_all.write(str(graphs[i][4][j]))
            else:
                graphs_all.write(","+str(graphs[i][4][j])) 
        graphs_all.write("\n")              
        for j in range(graphs[i][0]):
            if(j==0):
                graphs_all.write(str(graphs[i][5][j]))
            else:
                graphs_all.write(","+str(graphs[i][5][j]))                 
        graphs_all.write("\n")      
        for j in range(len(graphs[i][2])):
            if(j==0):
                graphs_all.write(str(graphs[i][2][j]))
            else:
                graphs_all.write(","+str(graphs[i][2][j])) 
        graphs_all.write("\n") 
    graphs_all.close()                  

output_log.close()
OIC=imp.load_source('occint_cluster',"./scripts/occint.py")
OIC.occint_cluster()

if(run_only_algo==0):
    print("Copying the results to the process folder")
    base_new=os.path.normpath(os.getcwd())
    os.makedirs("results",exist_ok=True)
    results=os.path.join(base_new,"results")
    process=os.path.join(base_new,"process")
    for item in os.listdir(process):
        s=os.path.join(process,item)
        d=os.path.join(results,item)
        shutil.copy2(s,d)
