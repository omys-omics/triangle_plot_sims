#! /usr/bin/env python3

import argparse

# Create an ArgumentParser to handle command-line arguments
parser = argparse.ArgumentParser(description='Read and copy the first 13 lines from an input file to an output file.')

# Add arguments for the input and output file paths
parser.add_argument('input_file_path', help='Path to the input file')
parser.add_argument('output_file_path', help='Path to the output file')

# Parse the command-line arguments
args = parser.parse_args()

# Open the input file and output file
with open(args.input_file_path, 'r') as input_file, open(args.output_file_path, 'w') as output_file:
    # start counter
    i = 0
    # start dictionary
    positions = {}
    for line in input_file:
      i += 1
      # write first 13 lines
      if i < 13:
        output_file.write(line)
      # get num of individuals
      if i == 13:
        fields = line.split('\t')
        inds = len(fields) - 9
        output_file.write(line)
      # get line for each position
      if i > 13:
        fields = line.split('\t')
        positions.setdefault(fields[1], '\t'.join(fields))
    
    # write invariant genotypes
    gt = "0|0\t"
    gts = gt * (inds - 1)
    gts += "0|0\n"
    
    # add invariant genotypes to new file
    for i in range(1, 2000001):
      j = str(i)
      if j not in positions.keys():
        line = f"1\t{i}\t.\tA\t.\t1000\tPASS\tMID=-1;S=0;DOM=0.5;PO=-1;GO=-1;MT=-1;AC=0;DP=1000\tGT\t"
        output_file.write(line)
        output_file.write(gts)
      else:
        line = positions[j]
        output_file.write(line)
        









