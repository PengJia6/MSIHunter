#=============================================================================
# Project : MSIhunter0.0.1
# Py Name: MSIHunter
# Author : Peng Jia
# Date : 18-10-23
# Email : pengjia@stu.xjtu.edu.cn 
# Description : 'main of the project'
#=============================================================================
import sys
import datetime
sys.path.append(sys.path[0][:-3])
sys.path.append(sys.path[0][:-3]+"data")
from src.argument import optInit
from src.bam2dis import bam2dis
from src.loadConfig import loadMicrosatellite,loadbed
from src.globalVar import *
def main():
    import os
    # global Arguments
    golInit()
    optInit() #argument procress
    Arguments=getArguments()
    if not os.path.exists(Arguments["outPre"]):
        os.makedirs(Arguments["outPre"])
    else:
        if os.path.isfile(Arguments["outPre"]):
            os.renames(Arguments["outPre"],Arguments["outPre"]+"_backup")
            os.makedirs(Arguments["outPre"])
    for filesuffix in [".dis","",".pro"]:
        file=open(Arguments["outPreF"]+filesuffix,"w")
        file.close()
    # if Arguments["inputBed"]!="NA":
    #     loadbed()
    if os.path.exists(Arguments["Microsatellite"]):
        loadMicrosatellite(Arguments["Microsatellite"])
    else:
        print('[MSIHunter ERROR] Fail to load Microsatellite from "'+Arguments["Microsatellite"]+'"')
        return -1
    if os.path.exists(Arguments["inputBam"]):
      bam2dis(Arguments["inputBam"])
    else:
        print("[Error:**]: Not such a bam file in "+Arguments["inputBam"] + " ,Please give a valid bam file!")
if __name__ == "__main__":
    # startTime=datetime.datetime.now()
    # print(startTime)
    main()
    # endTime = datetime.datetime.now()
    # print(endTime)
    # print(endTime-startTime)


