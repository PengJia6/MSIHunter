#=============================================================================
# Project : MSIhunter0.0.1
# Py Name: bam2dis
# Author : Peng Jia
# Date : 18-10-24
# Email : pengjia@stu.xjtu.edu.cn 
# Description : 'conver bam file to a distribution in microsatellites regions'
#=============================================================================
import pysam
from src.globalVar import *
def bam2dis(bamPath):
    global MicroSatellite
    print("[MSIHunter INFO] Loading bam file from "+bamPath+" ...")
    bam = pysam.AlignmentFile(bamPath, "rb")
    if not bam.check_index():
        print("[Error:]: Not index file for " + bamPath + " ,Please give a valid bam index!")
        return 0
    for chrId in MicroSatellite:
        for posStart,info in MicroSatellite[chrId].items():
            motifLen=int(info[1])
            motif=info[0]
            repeatTimes=int(info[2])
            prefix=info[3]
            suffix=info[4]
            posEnd=posStart+motifLen*repeatTimes
            queryStrat=posStart-5
            queryEnd=posEnd+5
            print(queryStrat,queryEnd)

            for alignment in bam.fetch(chrId, queryStrat,
                                       queryEnd):  # (refName,start,end): read which at least has a base between  start+1 and end-1
                print(alignment)

    print("[MSIHunter INFO] Loading bam file successfully")
if __name__ == "__main__":
    print()
