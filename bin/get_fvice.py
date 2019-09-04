#!/usr/bin/env python3

import pathlib
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir")
    parser.add_argument("output_file")
    return parser.parse_args()

args = parse_arguments()

files = list(pathlib.Path(args.input_dir).glob("*.out"))
files = sorted(files)

with open(args.output_file, "w") as outfile:
    header_was_printed = False
    for f in files:
        motif = f.stem
        with open(str(f), "r") as infile:
            header = infile.readline().strip()
            line = infile.readline().strip()
            line = "\t".join([motif, line])
            if not header_was_printed:
                header = "\t".join(["motif", header])
                print(header, file=outfile)
                print(line, file=outfile)
                header_was_printed = True
            else:
                print(line, file=outfile)