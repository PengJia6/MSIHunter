#=============================================================================
# Project : MSIhunter0.0.1
# Py Name: globalVar
# Author : pengjia
# Date : 18-10-24
# Email : pengjia@stu.xjtu.edu.cn 
# Description : ''
#=============================================================================\
# class GlobalVar:
#     db_handle = None
#     mq_client = None
#     def set_db_handle(db):
#       GlobalVar.db_handle = db
#     def get_db_handle():
#       return GlobalVar.db_handle
#     def set_mq_client(mq_cli):
#       GlobalVar.mq_client = mq_cli
#     def get_mq_client():
#       return GlobalVar.mq_client
global MicroSatellite, Distribution, Arguments
def golInit():
    global MicroSatellite, Distribution, Arguments,ShiftProbability
    MicroSatellite = {}
    Distribution = {}
    Arguments = {}
    ShiftProbability = {}


def setMicroSatellite( value):
    global Arguments
    MicroSatellite = value
def getMicroSatellite():
    global Arguments
    return MicroSatellite
def setArguments( value):
    global Arguments
    Arguments = value
def getArguments():
    global Arguments
    return Arguments
def setDistribution( value):
    global Distribution
    Distribution = value
def getDistribution():
    global Distribution
    return Distribution
def setShiftProbability( value):
    global ShiftProbability
    ShiftProbability = value
def getShiftProbability():
    global ShiftProbability
    return ShiftProbability


if __name__ == "__main__":
    print()
