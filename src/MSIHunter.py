#=============================================================================
# Project : MSIhunter0.0.1
# Py Name: MSIHunter
# Author : Peng Jia
# Date : 18-10-23
# Email : pengjia@stu.xjtu.edu.cn 
# Description : 'main of the project'
#=============================================================================
import pysam
import os
import collections
import pandas as pd
import numpy as np
from src.bam2dis import bam2dis
from src.globalVar import *
chrs=["chr"+str(i) for i in range(1,23)]+["chrX","chrY"]
def optInit():
    """
    argument procress
    """
    import argparse
    parser = argparse.ArgumentParser(description='MSIHunter: a Microsatellite Instability(MSI) detection tools')
    parser.add_argument('-i','--input_bam',required=True,type=str,nargs=1,default="CRC",help="input bam file")
    parser.add_argument('-o', '--output_prefix', required=True, type=str, nargs=1, default="NA",
                        help="prefix of the output")
    parser.add_argument('-m','--Microsatellite',required=True,type=str,nargs=1,default="NA",help="path of the Microsatellite")
    parser.add_argument('-bed', '--bed_region', type=str, nargs=1, default=".NA", help="only procress read in this bed region")
    parser.add_argument('-t','--tumor_type',type=str,nargs=1,default="CRC",help="tumor type of the case(CRC,UCEC,STAD)")
    # parser.add_argument("-v1", "--verbose", help="increase output verbosity",
    #                     action="store_true")
    # parser.add_argument("")
    args=parser.parse_args()
    arguments={}
    arguments["inputBam"]=args.input_bam[0]
    arguments["outPre"] = args.output_prefix[0]
    arguments["inputBed"] = args.bed_region[0]
    arguments["tumorType"]=args.tumor_type[0]
    arguments["Microsatellite"] = args.Microsatellite[0]
    # arguments["v"] = args.verbose
    # print("[MSIHunter I] input bam is ")
    print("[MSIHunter INFO] Initializing...  Successfully!")
    return arguments
    # print(arguments)
    # print(args.tumor_type)
    # print(args.input_bam)
def loadbed():
    print()
def loadMicrosatellite(pathMicosatellite):
    windowSize=100000
    print("[MSIHunter INFO] Loading Microsatellite file from " + pathMicosatellite+" ...")
    dfMicroSatellite=pd.read_table(pathMicosatellite)
    for index,row in dfMicroSatellite.iterrows():
        chrID=row["chr"]
        if chrID not in MicroSatellite :MicroSatellite[chrID]={}
        pos=row["pos"]
        MicroSatellite[chrID][pos]=[row["motif"],row["motifLen"],row["repeatTimes"],row["prefix"],row["suffix"]]
        """
        0: motif
        1: motifLen
        2: repeatTimes
        3: prefix
        4: suffix
        """
    # print(MicroSatellite[chrID].keys())
    print("[MSIHunter INFO] Loading Microsatellite successfully")

def main():
    import os
    argument=optInit() #argument procress
    # global MicroSatellite,Distribution
    if argument["inputBed"]!=".":
        loadbed()
    if os.path.exists(argument["Microsatellite"]):
        loadMicrosatellite(argument["Microsatellite"])
    else:
        print('[MSIHunter ERROR] Fail to load Microsatellite from "'+argument["Microsatellite"]+' "')
        return 2
    if os.path.exists(argument["inputBam"]):
        bam2dis(argument["inputBam"])
    else:
        print(argument)
if __name__ == "__main__":

    main()


