"""Simple program to efficently filter a count table in O(n)"""
import argparse

## Set up global variables as empty to begin with, so we throw an
## error if we're missing any

## Set up argument parsing
parser = argparse.ArgumentParser()
parser.add_argument('-l', '--labels')
parser.add_argument('-f', '--file')
parser.add_argument('-o', '--out')
args = parser.parse_args()

LABELS = args.labels
FILE = args.file
OUTFILE = args.out

if LABELS == None or FILE == None or OUTFILE == None:
    print("Missing arguments. Try --help")
    quit()

## Load in our file and read in the keys as a dictionary, so we have
## constant-time behavior
DIC = {}
with open(LABELS) as f:
    for line in f:
        curr_line = line.split()
        DIC[curr_line[3]] = True

## Iterate over the file and append our new key to the end of the file
with open(OUTFILE, 'w+') as out:
    with open(FILE) as f:
        out.write(f.readline()) # Preserve the header
        for line in f:
            curr_line = line.split()
            out_line = line.strip()
            try:
                newitem = DIC[curr_line[0]]  # 3rd line is ID
            except KeyError:
                continue
            out.write(out_line + "\n")
