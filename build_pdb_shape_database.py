#!/usr/bin/env python

import pdb
import project
import png
import glob
import os
import geometry
import string

size = 256  # should be enough for anybody
databaseDir = 'database'
pdbLocation = 'pdbs'
radius = 1.8 #radius of atoms, estimate
try:
  os.mkdir(databaseDir)
except OSError:
  pass  # directory exists, fine

for onePdb in glob.iglob(os.path.join(pdbLocation, '*.pdb')):  # every PDB
  pdbD = pdb.pdbData(onePdb)
  pdbCode = string.split(os.path.split(onePdb)[-1], '.')[0]
  points = pdbD.coords
  for vector in xrange(2):
    for thetaTen in xrange(0, 62, 3):
      vecX, vecY, vecZ = 0., 0., 0.
      if vector == 0:
        vecX = 1.
      elif vector == 1:
        vecY = 1.
      normalVec = geometry.normalizeVector((vecX, vecY, vecZ))
      theta = thetaTen / 10.  
      name = string.join([str(vecX), str(vecY), str(vecZ), str(theta)], '.')
      projectedPts = project.projectPointsOnto2D(points, normalVec, theta)
      mins, maxs = project.size2dSquare(projectedPts, radius)
      matrix = project.make2dMap(projectedPts, radius, mins, maxs, size)
      pngfile = open(
          os.path.join(databaseDir, pdbCode + '.' + name + '.png'), 'wb')
      pngwriter = png.Writer(size, size, greyscale=True, bitdepth=1)
      pngwriter.write(pngfile, matrix)
      pngfile.close()



'''
  #example code to read later
  reader = png.Reader(os.path.join(databaseDir, pdbCode + '.png'))
  for line in reader.read()[2]:
    print list(line)
'''
