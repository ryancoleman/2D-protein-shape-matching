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

for oneFile in glob.iglob(os.path.join(shapeDir, '*.txt')):  # every file
  oneName = string.split(os.path.split(oneFile)[-1], '.')[0]
  matrix = project.make2dMapFromShape(oneFile, size)
  pngfile = open(os.path.join(shapeDir, oneName + '.png'), 'wb')
  pngwriter = png.Writer(size, size, greyscale=True, bitdepth=1)
  pngwriter.write(pngfile, matrix)
  pngfile.close()

