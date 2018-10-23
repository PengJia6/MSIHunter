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
import pandas
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
    parser.add_argument("-v1", "--verbose", help="increase output verbosity",
                        action="store_true")
    # parser.add_argument("")
    args=parser.parse_args()
    arguments={}
    arguments["inputBam"]=args.input_bam[0]
    arguments["outPre"] = args.output_prefix[0]
    arguments["inputBed"] = args.bed_region[0]
    arguments["tumorType"]=args.tumor_type[0]
    arguments["v"] = args.verbose
    # print("[MSIHunter I] input bam is ")
    return arguments
    # print(arguments)
    # print(args.tumor_type)
    # print(args.input_bam)
def loadbed():
    print()
def bam2dis(bam):


    print()
def main():
    import os

    argument=optInit() #argument procress
    if argument["inputBed"]!=".":
        loadbed()
    if os.path.exists(argument["inputBam"]):
        bam2dis(argument["inputBam"])
    else:
        print(argument)
if __name__ == "__main__":

    main()
