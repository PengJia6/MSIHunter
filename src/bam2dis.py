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
def getRepeatTimes(alignment,motif,motifLen,prefix,suffix):
    Arguments = getArguments()
    if alignment.mapping_quality<Arguments["minimum_mapping_quality"]:
        return -1
    readString=alignment.query
    prefixState=readString.find(prefix)
    if prefixState<0:return -1
    suffixState=readString.rfind(suffix)
    if suffixState<0:return -3
    if prefixState+5>=suffixState:return -2
    repeat = readString[prefixState + 5: suffixState]
    while prefixState>=0:
        count=0
        start=prefixState+5
        while start==readString.find(motif, start):
            count+=1
            start =readString.find(motif, start)+motifLen
        if (motifLen==1 and count>=1) or (motifLen>1 and count>=1):
            if start==readString.find(suffix, start):
                # print(count, "    ", prefix,motif, suffix,repeat)
                return count
        prefixState=readString.find(prefix,prefixState+1)
    return -4
    # while prefixState>=0:
    #     count=0
    #     start0=prefixState+5
    #     start1=start0
    #     while start0==readString.find(motif, start1):
    #         count+=1
    #         start1 =readString.find(motif, start1)
    #         start1+=motifLen
    #         start0=start1
    #     if (motifLen==1 and count>=1) or (motifLen>1 and count>=1):
    #         start1=start0
    #         if start1==readString.find(suffix, start0):
    #             # print(conut, "    ", prefix,motif, suffix,repeat)
    #             return count
    #     prefixState=readString.find(prefix,prefixState+1)
    # return -4
def writeDistribution():
    Distribution = getDistribution()
    MicroSatellite = getMicroSatellite()
    Arguments = getArguments()
    ShiftProbability=getShiftProbability()
    file=open(Arguments["outPreF"]+".dis",'a')
    filePro=open(Arguments["outPreF"]+".pro",'a')
    for chrId in Distribution:
        for pos,dis in Distribution[chrId].items():
            # print(chrId,pos,dis)
            info=MicroSatellite[chrId][pos]
            # file.write(chrId+"\t"+str(pos)+"\t"+info[3].lower()+info[0]+str(info[2])+info[4].lower()+"\n"+
            #            "\t".join([str(key)+":"+str(dis[key]) for key in sorted(list(dis.keys()))])+"\n")
            file.write("\t".join([chrId,str(pos),info[3],info[0],str(info[2]),info[4]])+"\n"+
                       " ".join([str(key)+":"+str(dis[key]) for key in sorted(list(dis.keys()))])+"\n")
            filePro.write("\t".join([chrId,str(pos),info[3],info[0],str(info[2]),info[4],str(ShiftProbability[chrId][pos][0]),str(ShiftProbability[chrId][pos][1])])+"\n")
    file.close()
    filePro.close()
def calcuShiftProbability(disDict,refRepeatTimes):
    insShfit=0;delShfit=0;normal=0
    # print(refRepeatTimes)
    # print(disDict)
    for rpt in disDict:
        if rpt-refRepeatTimes>0:
            insShfit=insShfit+(rpt-refRepeatTimes)
            normal=normal+rpt
        else:
            delShfit=delShfit+(refRepeatTimes-rpt)
            normal = normal + rpt
    # print(insShfit,delShfit,normal)
    return round(delShfit/(insShfit+delShfit+normal),4),round(insShfit/(insShfit+delShfit+normal),4)

def bam2dis(bamPath):
    print("[MSIHunter INFO] Procressing bam file from " + bamPath + " ...")
    Distribution=getDistribution()
    MicroSatellite=getMicroSatellite()
    Arguments=getArguments()
    bam = pysam.AlignmentFile(bamPath, "rb")
    ShiftProbability=getShiftProbability()
    if not bam.check_index():
        print("[Error:**]: Not index file for " + bamPath + " ,Please give a valid bam index!")
        return 0
    for chrId in MicroSatellite:
        setDistribution({})
        Distribution[chrId]={}
        ShiftProbability[chrId]={}
        for posStart,info in MicroSatellite[chrId].items():
            motifLen=int(info[1])
            motif=info[0]
            repeatTimes=int(info[2])
            prefix=info[3]
            suffix=info[4]
            posEnd=posStart+motifLen*repeatTimes
            queryStrat=posStart-5
            queryEnd=posEnd+5
            # print(queryStrat,queryEnd)
            alignmentList=[]
            # print()
            for alignment in bam.fetch(chrId, queryStrat,
                                       queryEnd):  # (refName,start,end): read which at least has a base between  start+1 and end-1
                # print(alignment)
                alignmentList.append(alignment)
            if len(alignmentList) <Arguments["minimum_support_reads"]:continue
            repeatTimesDict={}
            for alignment in alignmentList:
                if alignment.is_unmapped:continue
                thisRepeatTimes=getRepeatTimes(alignment,motif,motifLen,prefix,suffix)
                if thisRepeatTimes<0:continue
                if thisRepeatTimes not in repeatTimesDict:repeatTimesDict[thisRepeatTimes]=0
                repeatTimesDict[thisRepeatTimes] += 1
            if sum(qqw.values())<Arguments["minimum_support_reads"]:continue
            else:
                Distribution[chrId][posStart]=repeatTimesDict
                ShiftProbability[chrId][posStart]=calcuShiftProbability(repeatTimesDict,MicroSatellite[chrId][posStart][2])
            if len(Distribution[chrId])>100:
                setDistribution(Distribution);
                setShiftProbability(ShiftProbability)
                writeDistribution();
                Distribution[chrId]={};
                setShiftProbability[chrId]={}
        if len(Distribution[chrId])!=0:
            setDistribution(Distribution);
            setShiftProbability(ShiftProbability)
            writeDistribution();
            Distribution[chrId] = {};
            ShiftProbability[chrId] = {}
    print("[MSIHunter INFO] Loading bam file successfully")
if __name__ == "__main__":
    # print(MicroSatellite)
    print()
