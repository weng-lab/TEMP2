#! /usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
import math

def main():
    print_help()
    # read TE size
    tel = {}
    for l in open(sys.argv[5]):
        li = l.strip().split()
        tel[li[0]] = int(li[1])
    tpr = {}
    for l in open(sys.argv[6]):
        li = l.strip().split()
            tpr[l.strip()] = [li[0], int(li[1]), int(li[2]), li[5]]
    # read bdgs
    dict_all_s_bdg = read_bdg(sys.argv[1], tel)
    dict_all_a_bdg = read_bdg(sys.argv[2], tel)
    dict_sglt_s_bdg = read_bdg(sys.argv[3], tel)
    dict_sglt_a_bdg = read_bdg(sys.argv[4], tel)
    # predict candidate TP region
    for rg in tpr:
        if tpr[rg][3] == "+":
            lamda = (sum(dict_sglt_a_bdg[te]) + 1) * (sum(dict_all_a_bdg[te][tstart:tend]) + 1) / (sum(dict_all_a_bdg[te]) + 1)

    for te in tel: # for each TE
        for i in range(int((tel[te]-1)/100)): # for each 100bp window
            tstart = i*100; tend = min((i+1)*100, tel[te])
            lamda = (sum(dict_sglt_a_bdg[te]) + 1) * (sum(dict_all_a_bdg[te][tstart:tend]) + 1) / (sum(dict_all_a_bdg[te]) + 1)
            x = sum(dict_sglt_a_bdg[te][tstart:tend])
            pv.append([te, tstart, tend, x, lamda, "-"])
            tstart = max(0, tel[te]-(i+1)*100); tend = tel[te]-i*100
            lamda = (sum(dict_sglt_s_bdg[te]) + 1) * (sum(dict_all_s_bdg[te][tstart:tend]) + 1) / (sum(dict_all_s_bdg[te]) + 1)
            x = sum(dict_sglt_s_bdg[te][tstart:tend])
            pv.append([te, tstart, tend, x, lamda, "+"])
    for i in pv:
        sys.stdout.write("\t".join(map(str, i)) + "\n")


# --------functions--------
def read_bdg(path, tel):
    d = {}
    for te in tel: d[te] = [0] * tel[te]
    for l in open(path, "r"):
        li = l.strip().split()
        for i in range(int(li[1]),int(li[2])):
            try:
                d[li[0]][i] = float(li[3])
            except KeyError:
                continue
    return d


def print_help():
    if len(sys.argv) < 5:
        print sys.argv[0] + " overall.sense.bdg overall.anti.bdg sinleton.sense.bdg singleton.anti.bdg te.size TPregion"
        print "output a bed file contains observed reads and lamda in 4th and 5th colums"
        sys.exit(0)

# --------process--------
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        bb.fun_print_error("user interrupted, abort!")
        sys.exit(0)
