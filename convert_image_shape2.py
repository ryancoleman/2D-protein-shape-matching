#!/usr/bin/env python

import pdb
import project
import png
import glob
import os
import geometry
import string

size = 256  # should be enough for anybody
shapeDir = 'shapes'

pngName = os.path.join(shapeDir, 'australia.png')
reader = png.Reader(pngName)
for line in reader.read()[2]:
  #for some reason when exporting from photoshop, every other pixel is white
  for start in xrange(0, len(line), 2):
    pixel = 0
    total = line[start]
    if total > 20:
      pixel = 1
    print pixel,
  print 
