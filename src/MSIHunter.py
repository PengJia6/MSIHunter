#=============================================================================
# Project : MSIhunter0.0.1
# Py Name: MSIHunter
# Author : Peng Jia
# Date : 18-10-23
# Email : pengjia@stu.xjtu.edu.cn 
# Description : 'main of the project'
#=============================================================================
import sys
sys.path.append(sys.path[0][:-3])
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
    if Arguments["inputBed"]!="NA":
        loadbed()
    if os.path.exists(Arguments["Microsatellite"]):
        loadMicrosatellite(Arguments["Microsatellite"])
    else:
        print('[MSIHunter ERROR] Fail to load Microsatellite from "'+Arguments["Microsatellite"]+' "')
        return 2
    if os.path.exists(Arguments["inputBam"]):
      bam2dis(Arguments["inputBam"])
    else:
        print(len(Arguments))
if __name__ == "__main__":
    main()
    # print()


