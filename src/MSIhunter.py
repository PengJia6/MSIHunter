# =============================================================================
# Project : MSIhunter0.0.1
# Py Name: MSIhunter
# Author : Peng Jia
# Date : 18-11-22
# Email : pengjia@stu.xjtu.edu.cn 
# Description : 'MSI detection using only tumor sequencing data'
# =============================================================================
import argparse
import os
import pysam
import pandas as pd
import numpy as np
from multiprocessing import Pool

global arguments, dfMicroSatellites, caseInfo


def argumentProcress():
    """
    argument procress
    """
    global arguments
    parser = argparse.ArgumentParser(description='MSIHunter: a Microsatellite Instability(MSI)'
                                                 ' detection tools using only tumor sequencing data!\n'
                                                 'You can test multiple sample one time in this tools')
    parser.add_argument('-i', '--input_configure', required=True, type=str, nargs=1,
                        help="The path of input configure files [required]")
    parser.add_argument('-o', '--workspace', required=True, type=str, nargs=1, default=["./workspace"],
                        help="prefix of the output [required]")
    parser.add_argument('-mc', '--microsatellites_configure', required=True, type=str, nargs=1, default=["NA"],
                        help="path of the microsatellites configure files [required]")
    parser.add_argument('-t', '--threads', type=int, nargs=1, default=[4],
                        help="mumber of additional threads to use [default:4]")
    parser.add_argument('-q', '--minimum_mapping_quality', type=int, nargs=1, default=[20],
                        help="minimum mapping quality of read [default:20]")
    parser.add_argument('-s', '--minimum_support_reads', type=int, nargs=1, default=[20],
                        help="minimum support reads of an available microsatellite [default:20]")
    args = parser.parse_args()
    arguments = {}
    arguments["input"] = args.input_configure_file[0]
    arguments["workspace"] = args.workspace[0] if args.workspace[0][-1] == "/" else args.workspace[0] + "/"
    arguments["Microsatellite"] = args.microsatellites_configure[0]
    arguments["threads"] = args.threads[0]
    arguments["minimum_support_reads"] = args.minimum_support_reads[0]
    arguments["minimum_mapping_quality"] = args.minimum_mapping_quality[0]
    ErrorStat = False
    if os.path.isfile(arguments["input"]):
        print("[INFO] The train set is : " + arguments["input"])
    else:
        print('[ERROR] The train set "' + arguments["input"] + '" is not exist, please check again')
        ErrorStat = True
    if os.path.isfile(arguments["Microsatellite"]):
        print("[INFO] The Microsatellites file  is : " + arguments["Microsatellite"])
    else:
        print('[ERROR] The Microsatellites file "'
              + arguments["Microsatellite"]
              + '" is not exist, please check again')
        ErrorStat = True
    if os.path.exists(arguments["workspace"]):
        print(
            '[ERROR] The workspace is still exist! in case of overwrite files in this workspace, '
            'please give a new work space')
        ErrorStat = True
        if ErrorStat: return False
    else:
        if ErrorStat: return False
        os.mkdir(arguments["workspace"])
        os.mkdir(arguments["workspace"] + "detailInfo/")
        print("[INFO] The workspace is : " + arguments["workspace"])
    return True


def loadMicroSatellite():
    """
    :return:
    """
    global dfMicroSatellites, arguments
    dfMicroSatellites = {}
    pdMSConfigure = pd.read_csv(arguments["Microsatellite"], index_col=0)
    for cancer, info in pdMSConfigure.iterrows():
        if not os.path.isfile(info["Discriminative_Microsatellites_path"]):
            print("[ERROR] This file " + info[
                "Discriminative_Microsatellites_path"] + "is not exist! Please check it again!")
            return False
        dfMicroSatellites[cancer] = pd.read_csv(info["Discriminative_Microsatellites_path"], index_col=0)
        print("[INFO] Loading Microsatellite file for " + cancer + " from " + info[
            "Discriminative_Microsatellites_path"] + " ...")

        columns = ["chr", "pos", "motif", "motifLen", "repeatTimes", "prefix", "suffix"]
        errorItem = []
        for item in columns:
            if item not in dfMicroSatellites[cancer].columns:
                errorItem.append(item)
        if len(errorItem) > 0:
            print("[ERROR] The item " + ",".join(
                errorItem) + " of the Micorsatellites file is not exist! Please check again!")
            return False
        if dfMicroSatellites[cancer].isnull().any().any():
            print("[ERROR] Loading Microsatellite file for " + cancer + " from " + info[
                "Discriminative_Microsatellites_path"] + " is incomplete, please check it again!")
            return False
    return True


def getRepeatTimes(alignment, motif, motifLen, prefix, suffix):
    """
    :param alignment:
    :param motif:
    :param motifLen:
    :param prefix:
    :param suffix:
    :return:
    """
    global arguments
    if alignment.mapping_quality < arguments["minimum_mapping_quality"]:
        return -1
    readString = alignment.query
    prefixState = readString.find(prefix)
    if prefixState < 0: return -1
    suffixState = readString.rfind(suffix)
    if suffixState < 0: return -3
    if prefixState + 5 >= suffixState: return -2
    while prefixState >= 0:
        count = 0
        start = prefixState + 5
        while start == readString.find(motif, start):
            count += 1
            start = readString.find(motif, start) + motifLen
        if (motifLen == 1 and count >= 1) or (motifLen > 1 and count >= 1):
            if start == readString.find(suffix, start):
                # print(count, "    ", prefix,motif, suffix,repeat)
                return count
        prefixState = readString.find(prefix, prefixState + 1)
    return -4


def calcuShiftProbability(disDict, refRepeatTimes):
    """
    :param disDict:
    :param refRepeatTimes:
    :return:
    """
    insShfit = 0;
    delShfit = 0;
    normal = 0
    for rpt in disDict:
        if rpt - refRepeatTimes > 0:
            insShfit = insShfit + (rpt - refRepeatTimes) * disDict[rpt]
            normal = normal + rpt * disDict[rpt] - (rpt - refRepeatTimes) * disDict[rpt]
        else:
            delShfit = delShfit + (refRepeatTimes - rpt) * disDict[rpt]
            normal = normal + rpt * disDict[rpt] - (refRepeatTimes - rpt) * disDict[rpt]
    return round(delShfit / (insShfit + delShfit + normal), 4), round(insShfit / (insShfit + delShfit + normal), 4)


def procressInputConfigure():
    """
    :return:
    """
    global arguments, dfMicroSatellites, caseInfo, Result

    print("[INFO] Loading the input configure file from " + arguments["input"] + " ...")
    caseInfo = pd.read_csv(arguments["input"], index_col=0, dtype="str")
    columns = ["bamPath", "cancerType"]
    errorItem = []
    for item in columns:
        if item not in caseInfo.columns:
            errorItem.append(item)
    if len(errorItem) > 0:
        print("[ERROR] The item " + ",".join(
            errorItem) + " of the input configure file is not exist! Please check again!")
        return False
    if caseInfo.isnull().any().any():
        print("[ERROR] The information of your input configure file is incomplete,please check again ")
        return False

    cancerTypeList = set(caseInfo["cancerType"])
    for cancer in cancerTypeList:

        if cancer not in dfMicroSatellites:
            print("[ERROR] There is no microsatellite file for " + cancer + ",please check again ")
    for id, info in caseInfo.iterrows():
        if not os.path.isfile(info['bamPath']):
            print("[ERROR] This file " + info["bamPath"] + "is not exist! Please check it again!")
            return False
    Result = caseInfo
    return True


def procressOneBam(caseid):
    """
    :param caseid:
    :return:
    """
    global caseInfo, arguments, dfMicroSatellites
    try:
        infoTest = caseInfo.loc[caseid, :]
        lociRes = pd.DataFrame()
        # print(infoTrain)
        print("[INFO] procressing " + str(caseid) + "...")
        os.mkdir(arguments["workspace"] + "detailInfo/" + str(caseid))
        file = open(arguments["workspace"] + "detailInfo/" + str(caseid) + "/" + str(caseid) + ".dis", "w")
        file.close()
        cancer = infoTest.loc["cancerType"]
        bam = pysam.AlignmentFile(infoTest["bamPath"], "rb")
        Distribution = {}
        curentMSNum = 0
        for id, info in dfMicroSatellites[cancer].iterrows():
            chrId = info["chr"]
            posStart = info["pos"]
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
                                       queryEnd):
                # (refName,start,end): read which at least has a base between  start+1 and end-1
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
                curentMSNum += 1
                Distribution[id] = repeatTimesDict
                series = dfMicroSatellites[cancer].loc[id]
                series["threshold"] = round(series["threshold"], 4)
                series["p"] = round(calcuShiftProbability(repeatTimesDict, repeatTimes)[0], 4)
                lociRes = lociRes.append(series)
            if curentMSNum % 1000 == 0:
                file = open(arguments["workspace"] + "detailInfo/" + str(caseid) + "/" + str(caseid) + ".dis", "a")
                for id, dis in Distribution.items():
                    file.write(str(id) + "\n" +
                               " ".join([str(key) + ":" + str(dis[key]) for key in sorted(list(dis.keys()))]) + "\n")
                file.close()
                Distribution = {}
        if len(Distribution) > 0:
            file = open(arguments["workspace"] + "detailInfo/" + str(caseid) + "/" + str(caseid) + ".dis", "a")
            for id, dis in Distribution.items():
                file.write(str(id) + "\n" +
                           " ".join([str(key) + ":" + str(dis[key]) for key in sorted(list(dis.keys()))]) + "\n")
            file.close()
        lociRes = lociRes[
            ["chr", "pos", "motif", "motifLen", "repeatTimes", "prefix", "suffix", "threshold", "p"]]
        lociRes.index.name = "id"
        fileRes = open(arguments["workspace"] + "detailInfo/" + str(caseid) + "/" + str(caseid), "w")
        totalNum = len(lociRes)
        InstableNum = len(lociRes[lociRes["p"] > lociRes["threshold"]])
        InstablePercentage = InstableNum / totalNum * 100

        # Result.loc[caseid,"percentageOfInstableMicosatellites"]=round(InstablePercentage,4)
        fileRes.write(
            "CaseID: " + str(caseid) + "\n"
                                       "TotalNumberOfMicrosatellites: " + str(totalNum) + "\n"
                                                                                          "NumberOfInstableMicrosatellites: " + str(
                InstableNum) + "\n"
                               "PercentageOfInstableMicrosatellite: " + str(round(InstablePercentage, 2)) + " %\n"
        )
        fileRes.close()
        lociRes.to_csv(arguments["workspace"] + "detailInfo/" + str(caseid) + "/" + str(caseid) + ".pro")
        return True
    except:
        print("[INFO] procressing " + str(caseid) + " ERROR!!")
        return False


def resultOut():
    global caseInfo, arguments
    try:
        for caseid in caseInfo.index:
            caseInfo.loc[caseid, "MSIScore(%)"] = \
                open(arguments["workspace"] + "detailInfo/" + str(caseid) + "/" + str(caseid)).readlines()[-1].split(
                    ": ")[
                    -1][:-2]
        caseInfo.to_csv(arguments["workspace"] + "Result.csv")
        return True
    except:
        print("[ERROR] error when processing resultOut()")
        return False


def main():
    global caseInfo, Result, arguments
    if not argumentProcress():
        return 1
    if not loadMicroSatellite():
        return 2
    if not procressInputConfigure():
        return 3
    pool = Pool(arguments["threads"])
    pool.map(procressOneBam, list(caseInfo.index))
    pool.close()
    if not resultOut():
        return 4
    return 0


if __name__ == "__main__":
    if main() > 1:
        os.system("rm -r " + arguments["workspace"])
    print()
