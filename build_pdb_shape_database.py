#!/usr/bin/env python

import pdb
import project
import png

pdbD = pdb.pdbData('pdbs/966c.pdb')
points = pdbD.getHeavyAtomXYZ()
print len(points)
