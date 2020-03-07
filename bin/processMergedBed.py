#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
import sys
import os
import math


def main():
    print_help()
    tel = int(sys.argv[2])
    ins = int(sys.argv[3])
    stepN = int(math.ceil((tel - ins) / 5.0)) + 1
    dict_tel = {}
    for l in open(sys.argv[1], "r"): # open merged bed
        line = l.strip().split()
        count = line[3].split("|")
        transRec = line[4].split("|")
        transN = transRec[0].split(",")[0]
        transStrand = transRec[0].split(",")[3]
        transCor = []
        transSpliced = []
        for i in transRec: transCor.append(map(int,i.split(",")[1:3]))
        for i in transRec: transSpliced.append(i.split(",")[-2])
        readSt = map(int,line[6].split("|"))
        readEd = map(int,line[7].split("|"))
        # skip if supporting read less than 2
        if len(count) < 2:
            if transSpliced[0] == "1":
                tsp = str(count[0])
            else:
                tsp = "0"
            print "\t".join([line[0],line[1],line[2],str(count[0]),"0",line[5],transN,transStrand,str(transCor[0][0]),str(transCor[0][1]),line[4],tsp,tsp])
        # else, calculate supporting read counts in each sliding window
        else:
            finalStep = 0
            finalCount = 0
            for i in range(stepN):
                if transStrand == "-": tst=i*5+1; ted=i*5+ins # forward slide
                else: tst=tel-i*5-ins+1; ted=tel-i*5 # reverse slide
                tc = 0 # calculate read number in window
                for j in range(len(transCor)):
                    if transCor[j][0] >= tst and transCor[j][1] <= ted:
                        tc += float(count[j])
                if tc > finalCount: # refresh final window when read number is higher
                    finalStep = i
                    finalCount = tc
            # anchor the final window
            if transStrand == "-":
                finalSt=finalStep*5+1; finalEd=min(finalStep*5+ins,tel)
            else:
                finalSt=max(tel-finalStep*5-ins+1,1); finalEd=tel-finalStep*5
            tpCount = 0
            fpCount = 0
            start = {}
            end = {}
            stn = {}
            edn = {}
            finalTrans = []
            teSt = 100000
            teEd = 0
            if finalStep == 0:
                teSt=finalSt
                teEd=finalEd
            for j in range(len(transCor)):
                if transCor[j][0] >= finalSt and transCor[j][1] <= finalEd:
                    tpCount += float(count[j])
                    finalTrans.append(transRec[j])
                    # anchor genome start and end coordinates
                    st = readSt[j]
                    ed = readEd[j]
                    if transSpliced[j] == "1":
                        tsp = float(count[j])
                    else:
                        tsp = 0
                    if st not in start:
                        start[st] = st
                        stn[st] = tsp
                    else:
                        start[st] = start[st] - (st-start[st])*2*float(count[j]) - 2*float(count[j])
                        stn[st] += tsp
                    if ed not in end:
                        end[ed] = ed
                        edn[ed] = tsp
                    else:
                        end[ed] = end[ed] + (end[ed]-ed)*2*float(count[j]) + 2*float(count[j])
                        edn[ed] += tsp
                    # anchor transposon start and end coordinates if finalStep is not 1
                    if finalStep != 0:
                        teSt = min(teSt, transCor[j][0])
                        teEd = max(teEd, transCor[j][1])
                else:
                    fpCount += float(count[j])
            st = min(start, key=start.get)
            ed = max(end, key=end.get)
            to = "|".join(finalTrans)
            print "\t".join([line[0],str(st),str(ed),str(tpCount),str(fpCount),line[5],transN,transStrand,str(teSt),str(teEd),to,str(stn[st]),str(edn[ed])])

# --------functions--------
def print_help():
    if len(sys.argv) < 2:
        print "usage:"
        print "%s in.merged.bed transposon_size insert_size"%sys.argv[0]
        sys.exit(0)

# --------process--------
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        bb.fun_print_error("user interrupted, abort!")
        sys.exit(0)
