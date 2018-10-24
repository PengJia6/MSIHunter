#=============================================================================
# Project : MSIhunter0.0.1
# Py Name: MSIHunter
# Author : 0.0.1
# Date : 18-10-23
# Email : pengjia@stu.xjtu.edu.cn 
# Description : 'main of the project'
#=============================================================================
import pysam
import os
import pandas as pd

MicroSatellite = pd.DataFrame()
Distribution = {}
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
    print("[MSIHunter INFO] Loading Microsatellite file from " + pathMicosatellite+" ...")
    MicroSatellite=pd.read_table(pathMicosatellite,index_col=0)

    print(MicroSatellite)
    print("[MSIHunter INFO] Loading Microsatellite successfully")
def bam2dis(bam):
    print("[MSIHunter INFO] Loading bam file from "+bam+" ...")
    print("[MSIHunter INFO] Loading bam file successfully")







    print()
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


