#!/usr/bin/env python3

# mat2jaspar.py
# Converts pwms into jaspar format

## Input
# ID
# 455 46 146 353
# 174 25 759 42
# 0 941 0 59
# 5 0 988 7

## Output
# > ID
# A [455 174 0 5 ]
# C [46 25 941 0 ]
# G [146 759 0 988 ]
# T [353 42 59 7 ]

import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("infile", help = "Input matrix file")
parser.add_argument("id", help = "Motif ID used for output")
parser.parse_args()

infile   = args.infile
motif_id = args.id

with open(infile, 'r') as my_matrix:
    A_str = []
    C_str = []
    G_str = []
    T_str = []
    non_prob = re.compile(r'[A-Za-z]')
    for line in my_matrix.readlines():
        if non_prob.search is None:
            line = line.strip().split(" ")
            line = [ int(x) for x in line ]
            A_str.append(line[0])
            C_str.append(line[1])
            G_str.append(line[2])
            T_str.append(line[3])