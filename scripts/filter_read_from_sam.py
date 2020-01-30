#!/usr/bin/env python3

import pysam
import sys

if __name__ == '__main__':
  if len(sys.argv) < 3:
    print('[!] Usage : ./filter_read_from_sam.py <in.bam> <out.bam>')
    exit()

  inputBamfile = pysam.AlignmentFile(sys.argv[1], "rb")
  outputBamfile = pysam.AlignmentFile(sys.argv[2], "wb", template = inputBamfile)
  K = 20

  for read in inputBamfile:
    if(read.is_secondary):
      continue

    outputRead = False

    if(read.is_unmapped):
      outputRead = True
    else:
      indelOrClippingSum = 0
      for a,b in read.cigartuples:
        if(a != 0):
          indelOrClippingSum += b;

      if(indelOrClippingSum >= K):
        outputRead = True;

    if(outputRead):
      lowBaseCount = 0
      for c in read.qual:
        if(c <= ','):
          lowBaseCount += 1

      if(lowBaseCount < 100):
        outputBamfile.write(read);

  inputBamfile.close()
  outputBamfile.close()
