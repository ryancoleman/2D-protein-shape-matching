#!/usr/bin/env python

import glob
import string
import png
import os
import project
import sys

databaseDir = 'database'
queryFileName = sys.argv[1]
ignore = None
if len(sys.argv) > 2:
  ignore = sys.argv[2] #left, right, top, or bottom can be ignored.
  #"North is not up and East is not right" -atom and his package, sadly ignored
  if ignore not in ['top', 'bottom', 'left', 'right']:
    print "usage is to ignore one half of the image, specified by " + \
        "top, left, bottom or right as the 2nd argument"
    sys.exit(1)
 
def readPngBitString(pngName, ignore=None):
  '''read a png into a bitstring. if ignore is top, bottom, left or right,
  don't read those bits'''
  bitlist = []
  reader = png.Reader(pngName)
  if ignore is None:
    for line in reader.read()[2]:
      bitlist.append(project.listToString(list(line)))
  elif ignore == 'top' or ignore == 'bottom':
    lines = []  # save all, then delete first or second half
    for line in reader.read()[2]:
      lines.append(line)
    if ignore == 'top':
      keptLines = lines[:len(lines)/2]  # keep the bottom half
    else:  # ignore the bottom or second half
      keptLines = lines[len(lines)/2:]  # keep the top half
    for line in keptLines:
      bitlist.append(project.listToString(list(line)))
  elif ignore == 'left' or ignore == 'right':
    for line in reader.read()[2]:
      listLine = list(line)
      if ignore == 'left':
        bitlist.append(project.listToString(list(line)[len(listLine)/2:]))
      else:   # ignore must be the right half, so keep left
        bitlist.append(project.listToString(list(line)[:len(listLine)/2]))
  bitstring = string.join(bitlist, '')
  return bitstring
  
querybits = readPngBitString(queryFileName, ignore)
bestTani, bestName = 0.0, None
for count, pngFileName in enumerate(
    glob.iglob(os.path.join(databaseDir, '*', '*.png'))):
  searchbits = readPngBitString(pngFileName, ignore)
  both = 0
  either = 0
  for count in xrange(len(querybits)):
    if querybits[count] == '1' and searchbits[count] == '1':
      both += 1
    if querybits[count] == '1' or searchbits[count] == '1':
      either += 1
  tanimoto = both/float(either)
  if tanimoto > bestTani:
    bestTani, bestName = tanimoto, pngFileName
print bestTani, bestName
