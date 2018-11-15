#=============================================================================
# Project : MSIhunter0.0.1
# Py Name: argument
# Author : 
# Date : 18-10-24
# Email : pengjia@stu.xjtu.edu.cn 
# Description : ''
#=============================================================================
import argparse
from src.globalVar import *
def optInit():
    """
    argument procress
    """
    print("[MSIHunter INFO] Initializing...  ")
    parser = argparse.ArgumentParser(description='MSIHunter: a Microsatellite Instability(MSI) detection tools')
    parser.add_argument('-i','--input_bam',required=True,type=str,nargs=1,default=["NA"],help="input bam file")
    parser.add_argument('-o', '--output_prefix', required=True, type=str, nargs=1, default=["NA"],
                        help="prefix of the output")
    parser.add_argument('-m','--Microsatellite',required=True,type=str,nargs=1,default=["NA"],help="path of the Microsatellite")
    parser.add_argument('-c', '--case_id', type=str, nargs=1, default=["NA"], help="case id")
    parser.add_argument('-q', '--minimum_mapping_quality', type=int, nargs=1, default=[10],
                        help="minimum mapping quality of read [default:10]")
    parser.add_argument('-s', '--minimum_support_reads', type=int, nargs=1, default=[10],
                        help="minimum support reads of a microsatellites [default:10]")
    parser.add_argument('-t','--tumor_type',type=str,nargs=1,default=["CRC"],help="tumor type of the case(CRC,UCEC,STAD)")
    # parser.add_argument("-v1", "--verbose", help="increase output verbosity",
    #                     action="store_true")

    args=parser.parse_args()
    arguments={}
    # print(args.input_bam)
    arguments["inputBam"]=args.input_bam[0]
    arguments["outPre"] = args.output_prefix[0] if args.output_prefix[0][-1]=="/" else args.output_prefix[0]+"/"
    arguments["CaseID"] = args.input_bam[0].split("/")[-1] if args.case_id[0]=="NA" else args.case_id[0]
    arguments["tumorType"]=args.tumor_type[0]
    arguments["Microsatellite"] = args.Microsatellite[0]
    arguments["minimum_support_reads"] = args.minimum_support_reads[0]
    arguments["minimum_mapping_quality"] = args.minimum_mapping_quality[0]
    arguments["outPreF"]=arguments["outPre"]+arguments["outPre"].split("/")[-2]
    print("[MSIHunter INFO] Initialization Successfully!")
    setArguments(arguments)
if __name__ == "__main__":
    print()
