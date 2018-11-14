#=============================================================================
# Project : MSIhunter0.0.1
# Py Name: loadConfig
# Author : Peng Jia
# Date : 18-10-24
# Email : pengjia@stu.xjtu.edu.cn 
# Description : 'load bed file and microsatellites'
#=============================================================================
from src.globalVar import *
import pandas as pd
def loadbed():
    print()
def loadMicrosatellite(pathMicosatellite):
    MicroSatellite=getMicroSatellite()
    # windowSize=100000
    print("[MSIHunter INFO] Loading Microsatellite file from " + pathMicosatellite+" ...")
    dfMicroSatellite = pd.read_csv(pathMicosatellite,index_col=0)
    for index,row in dfMicroSatellite.iterrows():
        chrID = row["chr"]
        if chrID not in MicroSatellite :MicroSatellite[chrID]={}
        pos=row["pos"]
        MicroSatellite[chrID][pos]=[row["motif"],row["motifLen"],row["repeatTimes"],row["prefix"],row["suffix"],row["threshold"]]

    # dfMicroSatellite=pd.read_table(pathMicosatellite)

    # for index,row in dfMicroSatellite.iterrows():
    #     chrID=row["chr"]
    #     if chrID not in MicroSatellite :MicroSatellite[chrID]={}
    #     pos=row["pos"]
    #     MicroSatellite[chrID][pos]=[row["motif"],row["motifLen"],row["repeatTimes"],row["prefix"],row["suffix"]]
    #     """
    #     0: motif
    #     1: motifLen
    #     2: repeatTimes
    #     3: prefix
    #     4: suffix
    #     """
    setMicroSatellite(MicroSatellite)
    print("[MSIHunter INFO] Loading Microsatellite successfully")
if __name__ == "__main__":
    print()
