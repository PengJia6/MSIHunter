#=============================================================================
# Project : MSIhunter0.0.1
# Py Name: SelectDiscriminativeMS
# Author : PengJia
# Date : 18-11-20
# Email : pengjia@stu.xjtu.edu.cn 
# Description : 'Select discriminative microsatellites for each tumor for each tumor type'
#=============================================================================
import argparse
import pandas as pd
import numpy as np
import os
import pysam
from scipy.stats import *
from multiprocessing import Pool
global dfMicroSatellite, TrainInfo, arguments,trainDict
dfMicroSatellite=pd.DataFrame()
trainDict={}
def argumentProcress():
    """
    argument procress
    """
    global arguments
    parser = argparse.ArgumentParser(description='Select discriminative microsatellites for each type of tumor! ')
    parser.add_argument('-train','--train_set_for_sites_selection',required=True,type=str,nargs=1,help="train set information [required]")
    parser.add_argument('-w', '--workspace', required=True, type=str, nargs=1, default=["./workspace"],
                        help="prefix of the output [required]")
    parser.add_argument('-m','--Microsatellite',required=True,type=str,nargs=1,default=["NA"],help="path of the Microsatellite [required]")
    parser.add_argument('-t', '--threads', type=int, nargs=1, default=[4],
                        help="Number of additional threads to use [default:4]")
    parser.add_argument('-q', '--minimum_mapping_quality', type=int, nargs=1, default=[20],
                        help="minimum mapping quality of read [default:20]")
    parser.add_argument('-s', '--minimum_support_reads', type=int, nargs=1, default=[20],
                        help="minimum support reads of an available microsatellite [default:20]")
    args=parser.parse_args()
    arguments={}
    # print(args.input_bam)
    arguments["train"]=args.train_set_for_sites_selection[0]
    arguments["workspace"] = args.workspace[0] if args.workspace[0][-1]=="/" else args.workspace[0]+"/"
    arguments["Microsatellite"] = args.Microsatellite[0]
    arguments["threads"] = args.threads[0]
    arguments["minimum_support_reads"] = args.minimum_support_reads[0]
    arguments["minimum_mapping_quality"] = args.minimum_mapping_quality[0]
    ErrorStat=False
    if os.path.isfile(arguments["train"]):
        print("[INFO] The train set is: "+arguments["train"])
    else:
        print('[ERROR] The train set "'+arguments["train"]+'" is not exist, please check again')
        ErrorStat=True

    if os.path.isfile(arguments["Microsatellite"]):
        print("[INFO] The Microsatellites file  is: " + arguments["Microsatellite"])
    else:
        print('[ERROR] The Microsatellites file "' + arguments["Microsatellite"] + '" is not exist, please check again')
        ErrorStat = True

    if os.path.exists(arguments["workspace"]):
        print('[ERROR] The workspace is still exist! in case of overwrite files in this workspace, please give a new work space' )
        ErrorStat = True
        if ErrorStat :return False
    else:
        if ErrorStat: return False
        os.mkdir(arguments["workspace"])
        os.mkdir(arguments["workspace"]+"dis/")
        print("[INFO] The workspace is : " + arguments["workspace"])
    return True
def loadMicroSatellite():
    """

    :return:
    """
    global dfMicroSatellite,arguments
    print("[INFO] Loading Microsatellite file from " + arguments["Microsatellite"] + " ...")
    dfMicroSatellite = pd.read_csv(arguments["Microsatellite"],index_col=0)
    columns=["chr","pos","motif","motifLen","repeatTimes","prefix","suffix"]
    errorItem=[]
    for item in columns:
        if item not in dfMicroSatellite.columns:
            errorItem.append(item)
    if len(errorItem)>0:
        print("[ERROR] The item "+",".join(errorItem)+" of the Micorsatellites file is not exist! Please check again!")
        return False
    return True
def procressTrainSet():
    """

    :return:
    """
    global TrainInfo,arguments,trainDict
    print("[INFO] Loading train file from " + arguments["train"] + " ...")
    TrainInfo=pd.read_csv(arguments["train"],index_col=0,dtype="str")
    columns=["case","bamPath","cancerType","MSI"]
    errorItem = []
    for item in columns:
        if item not in TrainInfo.columns:
            errorItem.append(item)
    if len(errorItem) > 0:
        print("[ERROR] The item " + ",".join(
            errorItem) + " of the Micorsatellites file is not exist! Please check again!")
        return False
    if TrainInfo.isnull().any().any():
        print("[ERROR] The train set information is incomplete,please check again ")
        return False
    # add some judgement for
    cancerTypeList=set(TrainInfo["cancerType"])
    for cancer in cancerTypeList:
        trainDict[cancer]={}
        data=TrainInfo[TrainInfo["cancerType"]==cancer]
        # MSInum=len()
        trainDict[cancer]["MSI_P"] = list(data[data["MSI"]=="MSI-H"].index)
        trainDict[cancer]["MSI_N"] = list(data[data["MSI"] != "MSI-H"].index)
        if len(trainDict[cancer]["MSI_P"])<2:
            print("[ERROR] Too few MSI-H sample in "+cancer)
            return False
        if len(trainDict[cancer]["MSI_N"])<2:
            print("[ERROR] Too few MSS/MSI-H sample in " + cancer)
            return False
        # train_list[cancer]["MSI-H"]=list()
    # print(cancerTypeList)
    return True

def getRepeatTimes(alignment,motif,motifLen,prefix,suffix):
    """
    :param alignment:
    :param motif:
    :param motifLen:
    :param prefix:
    :param suffix:
    :return:
    """
    global arguments
    if alignment.mapping_quality<arguments["minimum_mapping_quality"]:
        return -1
    readString=alignment.query
    prefixState=readString.find(prefix)
    if prefixState<0:return -1
    suffixState=readString.rfind(suffix)
    if suffixState<0:return -3
    if prefixState+5>=suffixState:return -2
    # repeat = readString[prefixState + 5: suffixState]
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
def calcuShiftProbability(disDict,refRepeatTimes):
    """
    :param disDict:
    :param refRepeatTimes:
    :return:
    """
    insShfit=0;delShfit=0;normal=0
    # print(refRepeatTimes)
    # print(disDict)
    for rpt in disDict:
        if rpt-refRepeatTimes>0:
            insShfit=insShfit+(rpt-refRepeatTimes)*disDict[rpt]
            normal=normal+rpt*disDict[rpt]-(rpt-refRepeatTimes)*disDict[rpt]
        else:
            delShfit=delShfit+(refRepeatTimes-rpt)*disDict[rpt]
            normal = normal + rpt*disDict[rpt]-(refRepeatTimes-rpt)*disDict[rpt]
    # print(insShfit,delShfit,normal)
    return round(delShfit/(insShfit+delShfit+normal),4),round(insShfit/(insShfit+delShfit+normal),4)
def procressOneBam(caseid):
    """
    :param caseid:
    :return:
    """
    global TrainInfo,arguments,dfMicroSatellite
    try:
        infoTrain=TrainInfo.loc[caseid,:]
        lociRes=pd.DataFrame()
        # print(infoTrain)
        print("[INFO] procressing "+str(caseid)+"...")
        file=open(arguments["workspace"]+"dis/"+str(caseid)+"_"+infoTrain["case"]+".dis","w")
        file.close()
        bam=pysam.AlignmentFile(infoTrain["bamPath"], "rb")
        Distribution={}
        curentMSNum=0
        for id,info in dfMicroSatellite.iterrows():
            chrId=info["chr"]
            posStart=info["pos"]
            motifLen = int(info["motifLen"])
            motif = info["motif"]
            repeatTimes = int(info["repeatTimes"])
            prefix = info['prefix']
            suffix = info['suffix']
            posEnd = posStart + motifLen * repeatTimes
            queryStrat = posStart - 5
            queryEnd = posEnd + 5
            alignmentList = []
            for alignment in bam.fetch(chrId, queryStrat,
                                       queryEnd):  # (refName,start,end): read which at least has a base between  start+1 and end-1
                alignmentList.append(alignment)
            if len(alignmentList) < arguments["minimum_support_reads"]: continue
            repeatTimesDict = {}
            for alignment in alignmentList:
                if alignment.is_unmapped: continue
                thisRepeatTimes = getRepeatTimes(alignment, motif, motifLen, prefix, suffix)
                if thisRepeatTimes < 0: continue
                if thisRepeatTimes not in repeatTimesDict: repeatTimesDict[thisRepeatTimes] = 0
                repeatTimesDict[thisRepeatTimes] += 1
            if sum(repeatTimesDict.values()) < arguments["minimum_support_reads"]:
                continue
            else:
                curentMSNum+=1
                Distribution[id] = repeatTimesDict
                lociRes.loc[id,"p"]= calcuShiftProbability(repeatTimesDict,repeatTimes)[0]
                # print(list(repeatTimesDict.keys())[np.argmax(np.array(repeatTimesDict.values()))])
                lociRes.loc[id,"rpl"]=list(repeatTimesDict.keys())[np.argmax(np.array(repeatTimesDict.values()))]
            if curentMSNum%1000==0:
                file = open(arguments["workspace"] + "dis/" + str(caseid) + "_" + infoTrain["case"] + ".dis", "a")
                for id, dis in Distribution.items():
                    file.write(str(id) + "\n" +
                               " ".join([str(key) + ":" + str(dis[key]) for key in sorted(list(dis.keys()))]) + "\n")
                file.close()
                Distribution={}
        if len(Distribution)>0:
            file = open(arguments["workspace"] + "dis/" + str(caseid) + "_" + infoTrain["case"] + ".dis", "a")
            for id, dis in Distribution.items():
                file.write(str(id) + "\n" +
                           " ".join([str(key) + ":" + str(dis[key]) for key in sorted(list(dis.keys()))]) + "\n")
            file.close()
        lociRes.index.name="id"
        lociRes.to_csv(arguments["workspace"] + "dis/" + str(caseid) + "_" + infoTrain["case"] + ".pro")
        return True
    except:
        print("[INFO] procressing " + str(caseid) + " ERROR!!")
        return False
def moreCoverage():
    global arguments,TrainInfo,trainDict,dfMicroSatellite
    try:
        for cancer in trainDict:
            LocisMatrix = pd.DataFrame()
            indelMatrix=pd.DataFrame()
            data=TrainInfo[TrainInfo["cancerType"]==cancer]
            caseNum=len(data)
            for caseid,infoTrain in data.iterrows():
                lociRes=pd.read_csv(arguments["workspace"] + "dis/" + str(caseid) + "_" + infoTrain["case"] + ".pro",index_col=0)
                for loci,info in lociRes.iterrows():
                    LocisMatrix.loc[caseid,loci]=info["p"]
                    indelMatrix.loc[caseid,loci]=info["rpl"]
            LocisMatrix.to_csv(arguments["workspace"] + "originalLocisMatrix" + cancer + ".csv")
            indelMatrix.to_csv(arguments["workspace"] + "originalRepeatTimesMatrix" + cancer + ".csv")
            fewerCoverageLociList=[]
            for loci in LocisMatrix.columns:
                if len(LocisMatrix[loci].dropna())*3<caseNum:
                    fewerCoverageLociList.append(loci)
            LocisMatrix.drop(fewerCoverageLociList,axis=1,inplace=True)
            indelMatrix.drop(fewerCoverageLociList, axis=1, inplace=True)
            mutationParameterHigh=[]
            for loci in indelMatrix.columns:
                if abs(np.mean(indelMatrix[loci].dropna())-dfMicroSatellite.loc[loci,"repeatTimes"])>0.1 or np.std(indelMatrix[loci].dropna())>1:
                    mutationParameterHigh.append(loci)
            LocisMatrix.drop(mutationParameterHigh, axis=1, inplace=True)
            LocisMatrix.to_csv(arguments["workspace"]+"finalLocisMatrix"+cancer+".csv")
        return True
    except:
        print("[ERROR] error when processing moreCoverage()")
        return False


def moreCorrelation():
    global arguments, TrainInfo, trainDict
    try:
        for cancer in trainDict:
            LocisMatrix=pd.read_csv(arguments["workspace"]+"finalLocisMatrix"+cancer+".csv",index_col=0)
            MSIPMatrix = LocisMatrix.loc[trainDict[cancer]["MSI_P"], :]
            MSINMatrix = LocisMatrix.loc[trainDict[cancer]["MSI_N"], :]
            dfSignificantP = pd.DataFrame()
            for loci in LocisMatrix.columns:
                msiNP = np.array(MSINMatrix[loci])
                msiPP = np.array(MSIPMatrix[loci])
                if msiNP.shape[0] < 2 or msiPP.shape[0] < 2: continue
                ranksumsP = ranksums(msiPP, msiNP)[1]
                if ranksumsP < 0.5 and msiPP.mean() >= msiNP.mean():
                    # print(ranksumsP)
                    series = dfMicroSatellite.loc[loci]
                    nu = msiNP.mean()
                    # series["threshold"] = round(nu + 3 * msiNP.std(), 4)
                    # print(series)
                    dfSignificantP = dfSignificantP.append(series)
                    dfSignificantP.loc[loci,"threshold"] = round(nu + 3 * msiNP.std(), 4)
                    # print(dfSignificantP)
            dfSignificantP.to_csv(arguments["workspace"] + "MS_" + cancer + ".csv")
            return True
    except:
        print("[ERROR] error when processing moreCorrelation()")
        return False


def main():
    """
    :return:
    """
    global TrainInfo,arguments
    if not argumentProcress():
        return 1
    if not loadMicroSatellite():
        return 2
    if not procressTrainSet():
        return 3
    pool=Pool(arguments["threads"])
    pool.map(procressOneBam,list(TrainInfo.index))
    pool.close()
    if not  moreCoverage():
        return 4
    if not moreCorrelation():
        return 5
    return 0
if __name__ == "__main__":
    global arguments
    sta=main()
    if sta>1:
        os.system("rm -r "+arguments["workspace"])