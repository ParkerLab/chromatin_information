#!/usr/bin/env python2

from __future__ import print_function
import sys

infile = sys.argv[1]

with open(infile, 'r') as data:
    for line in data.readlines():
        lineData = line.strip().split("\t")
        chrom = lineData[0]
        try:
            st = int(lineData[1])
            en = int(lineData[2])
        except ValueError:
            st = int(float(lineData[1]))
            en = int(float(lineData[2]))
        paired = lineData[3]
        if en - st > 8:
            en, st = en - 4, st + 4
        lineData = [chrom, st, en, paired]
        print(*lineData, sep = "\t", file = sys.stdout)
