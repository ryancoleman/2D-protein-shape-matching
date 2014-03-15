#!/usr/bin/env python

import pdb
import project
import png
import glob
import os
import geometry
import string

size = 256
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
  normalVec = geometry.normalizeVector((1., 0., 0.))
  theta = 0.
  projectedPts = project.projectPointsOnto2D(points, normalVec, theta)
  mins, maxs = project.size2dSquare(projectedPts, radius)
  matrix = project.make2dMap(projectedPts, radius, mins, maxs, size)
  pngfile = open(os.path.join(databaseDir, pdbCode + '.png'), 'wb')
  pngwriter = png.Writer(size, size, greyscale=True, bitdepth=1)
  pngwriter.write(pngfile, matrix)
  pngfile.close()
