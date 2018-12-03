# =============================================================================
# Project : MSIhunter0.0.1
# Py Name: ScanMicrosatellites
# Author :
# Date : 18-12-2
# Email : pengjia@stu.xjtu.edu.cn
# Description : ''
# =============================================================================
import argparse
import os
import datetime
import numpy as np

global arguments, chrstr


def argumentProcress():
    """
    argument procress
    """
    global arguments
    parser = argparse.ArgumentParser(description='Scan Microsatellites from reference genome')
    parser.add_argument('-r', '--reference', required=True, type=str, nargs=1,
                        help="path of reference file [required]")
    parser.add_argument('-m', '--microsatellites', required=True, type=str, nargs=1, default=["./workspace"],
                        help="path of output microsatellites [required]")
    parser.add_argument('--maximal_repeat_unit_size', default=[5], type=int, nargs=1,
                        help=" maximal repeat unit size [default=5]")
    parser.add_argument('--context_length', default=[5], type=int, nargs=1,
                        help=" size of prefix and suffix in output [default=5]")
    parser.add_argument('--ranges_of_repeat_times', default=["1-1:10;2-5:5"], type=str, nargs=1,
                        help="ranges_of_repeat_times [1-1:10;2-5:5]")
    parser.add_argument('-t', '--threads', type=int, nargs=1, default=[4],
                        help="Number of additional threads to use [default:4]")
    # parser.add_argument('-q', '--minimum_mapping_quality', type=int, nargs=1, default=[20],
    #                     help="minimum mapping quality of read [default:20]")
    # parser.add_argument('-s', '--minimum_support_reads', type=int, nargs=1, default=[20],
    #                     help="minimum support reads of an available microsatellite [default:20]")
    args = parser.parse_args()
    arguments = {}
    # print(args.input_bam)
    arguments["reference"] = args.reference[0]
    arguments["microsatellites"] = args.microsatellites[0]
    arguments["maximal_repeat_unit_size"] = args.maximal_repeat_unit_size[0]
    arguments["threads"] = args.threads[0]
    arguments["context_length"] = args.context_length[0]
    arguments["ranges_of_repeat_times"] = {}
    for i in args.ranges_of_repeat_times[0].split(";"):
        # print(i)
        unitRange, repeatRange = i.split(":")
        unitStart, unitEnd = tuple(map(int, unitRange.split("-")))
        repeatStart=int(repeatRange)
        # print(unitStart,unitEnd,"  ",repeatStart, repeatEnd)
        for ur in range(unitStart, unitEnd + 1):
            arguments["ranges_of_repeat_times"][ur] = {}
            arguments["ranges_of_repeat_times"][ur]= repeatStart

    for i in range(1, arguments["maximal_repeat_unit_size"] + 1):
        if i not in arguments["ranges_of_repeat_times"]:
            ErrorStat = True
            print("[ERROR] You need to give the range of repeat times when the repeat unit length is", i)

    if os.path.isfile(arguments["reference"]):
        print("[INFO] The reference is: " + arguments["reference"])
        ErrorStat = False
    else:
        print('[ERROR] The reference "' + arguments["reference"] + '" is not exist, please check again')

    if os.path.exists(arguments["microsatellites"]):
        print('[ERROR] ' + arguments["microsatellites"] + ' is still exist! in case of overwrite files in this file,'
                                                          ' please give a new path')
        ErrorStat = True
    else:
        print("[INFO] The output microsatellites file is: " + arguments["microsatellites"])
    if ErrorStat: return False
    file = open(arguments["microsatellites"], "w")
    file.write("ID,chr,pos,motif,motifLen,repeatTimes,prefix,suffix\n")
    file.close()
    return True


def scanGenome():
    global arguments
    alphabet = {"A": "", "T": "", "C": "", "G": ""}
    def procressOneContig(chrome):
        content = arguments["context_length"]
        resBuf = ""
        chrLen = len(chrstr) - 100
        if chrLen < 20: return
        print("[INFO] Processing " + chrome + " ...")
        i = 0
        while i < chrLen:
            tmp2 = chrstr[i]
            # print(tmp2)

            if tmp2 not in alphabet:
                # print(i)
                i += 1
                continue
            i0 = i
            repeatLength = 0
            while chrstr[i] == tmp2:
                repeatLength += 1
                i += 1
            if repeatLength >= arguments["ranges_of_repeat_times"][1]:
                # print(chrome, i0, repeatLength, tmp2,chrstr[i0-5:i0],chrstr[i:i+5])
                resBuf+=",".join(
                    [chrome + "_" + str(i0), chrome, str(i0), tmp2, "1", str(repeatLength),
                     chrstr[i0 - content:i0], chrstr[i:i + content] + "\n"])
            else:
                FlagRepeat = False
                FlagN = False
                for k in range(2, arguments["maximal_repeat_unit_size"] + 1):
                    i = i0
                    m = 0
                    while m < k:
                        if chrstr[i] not in alphabet:
                            FlagN = True
                            break
                        m += 1
                        i += 1
                    if FlagN:
                        i0 = i - 1
                        break
                    s00 = chrstr[i - k:i]
                    repeatLength = 0
                    i = i - k
                    while chrstr[i:i + k] == s00:
                        repeatLength += 1
                        i += k
                    if repeatLength >= arguments["ranges_of_repeat_times"][k]:
                        resBuf +=",".join([chrome + "_" + str(i0), chrome, str(i0), s00, str(k), str(repeatLength),
                                                chrstr[i0 - content:i0], chrstr[i:i + content] + "\n"])
                        FlagRepeat = True
                        break
                if FlagRepeat:
                    i = i0 + k * repeatLength
                else:
                    i = i0 + 1
        file = open(arguments["microsatellites"], "a")
        for item in resBuf:
            file.write(item)
        file.close()
        # resBuf=[]

    with open(arguments["reference"]) as ref:
        chrstr = ""
        for line in ref:
            if line[0] == ">":
                if chrstr != "":
                    chrstr = chrstr.upper()
                    chrstr = chrstr + "N" * 100
                    procressOneContig(chrome)
                    chrome = line[:-1].split(" ")[0][1:]
                    chrstr = ""
                    continue
                else:
                    chrome = line[:-1].split(" ")[0][1:]
                    chrstr = ""
                    continue
            else:
                chrstr = chrstr + line[:-1]
        chrstr = chrstr.upper()
        chrstr = chrstr + "N" * 100
        procressOneContig(chrome)
def main():
    global arguments
    if not argumentProcress():
        return 1
    scanGenome()


if __name__ == "__main__":
    startTime = datetime.datetime.now()
    print(startTime)
    main()
    endTime = datetime.datetime.now()
    print(endTime)
    print(endTime - startTime)
