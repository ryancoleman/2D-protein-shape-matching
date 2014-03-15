#!/usr/bin/env python

import pdb
import project
import png
import glob
import os
import geometry

size = 256
databaseDir = 'database'
pdbLocation = 'pdbs'
radius = 3

for onePdb in glob.iglob(os.path.join(pdbLocation, '*.pdb')):  # every PDB
  pdbD = pdb.pdbData(onePdb)
  points = pdbD.coords
  normalVec = geometry.normalizeVector((1., 0., 0.))
  theta = 0.
  projectedPts = project.projectPointsOnto2D(points, normalVec, theta)
  mins, maxs = project.size2dSquare(projectedPts, radius)
  matrix = project.make2dMap(projectedPts, radius, mins, maxs, size)
  pngfile = open('test.png', 'wb')
  pngwriter = png.Writer(size, size, greyscale=True, bitdepth=1)
  pngwriter.write(pngfile, matrix)
  pngfile.close()
