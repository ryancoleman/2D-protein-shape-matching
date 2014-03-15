#!/usr/bin/env python

import math
import geometry
import random

def projectPoints(pointList, vector, theta):
  '''
  http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/
  ArbitraryAxisRotation.html
  f(x,y,z,u,v,w,theta) =
  u(ux + vy + wz )(1 - cos theta) + x cos theta + (- wy + vz) sin theta
  v(ux + vy + wz )(1 - cos theta) + ycos theta + (wx - uz) sin theta
  w(ux +  vy + wz )(1 - cos theta) + z cos theta + (- vx + uy) sin theta
  
  use normal vector UVW, rotate theta around this, use every point in pointList

  '''
  vecU, vecV, vecW = vector
  cosTheta = math.cos(theta)
  sinTheta = math.sin(theta)



points = [[1., 2., 3.], [3., 2., 3.], [-2, 0., 1.]]
vecRot = geometry.normalizeVector(
    (random.random(), random.random(), random.random()))
theta = random.random()
print projectPoints(points, vecRot, theta)
