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

for onePdb in glob.iglob(os.path.join(pdbLocation, '*.pdb')):  # every PDB
  pdbD = pdb.pdbData(onePdb)
  points = pdbD.getHeavyAtomXYZ()
  normalVec = geometry.normalizeVector((1., 0., 0.))
  theta = 0.
  projectedPts = project.projectPointsOnto2D(points, normalVec, theta)
  print projectedPts
  pngwriter = png.Writer(size, size, greyscale=True, bitdepth=1)
