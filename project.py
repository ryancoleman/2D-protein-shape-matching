#!/usr/bin/env python

import math
import geometry
import random

def projectPointsOnto2D(pointList, vector, theta):
  '''
  moves 3d points around, then discared z
  '''
  newPts = projectPoints(pointList, vector, theta)
  twodpts = []
  for onePt in newPts:
    twodpts.append(onePt[:2])  # discard Z
  return twodpts

def projectPoints(pointList, vector, theta):
  '''
  moves 3d points around to other 3d points

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
  newPts = []
  for point in pointList:
    ptX, ptY, ptZ = point
    uxvywz = vecU * ptX + vecV * ptY + vecW * ptZ
    newX = vecU * uxvywz * (1 - cosTheta) + ptX * cosTheta + (
         -vecW * ptY + vecV * ptZ) * sinTheta
    newY = vecV * uxvywz * (1 - cosTheta) + ptY * cosTheta + (
         -vecW * ptX + vecU * ptZ) * sinTheta
    newZ = vecW * uxvywz * (1 - cosTheta) + ptZ * cosTheta + (
         -vecV * ptX + vecU * ptY) * sinTheta
    newPts.append((newX, newY, newZ))
  return newPts

#following is for rudimentary testing
if __name__ == "__main__":
  points = [[1., 2., 3.], [3., 2., 3.], [-2, 0., 1.]]
  vecRot = geometry.normalizeVector(
      (random.random(), random.random(), random.random()))
  theta = random.random()
  print projectPoints(points, vecRot, theta)
