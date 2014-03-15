#!/usr/bin/env python

import math
import geometry
import random
import string

def listToString(aList):
  '''turns [1, 0, 0, 1] into '1001' '''
  retStr = string.join([str(item) for item in aList], '')
  return retStr

def make2dMap(projectPts, radius, mins, maxs, size):
  '''takes list of 2d points & radius. makes a grid of the size by size.
  using mins & maxs, maps those points onto the grid, coloring points within
  the radius. returns ['0110', '0010', '0100', '0000'] for example.
  which can be directly written to a 1-bit png file'''
  oneRow = [0 for count in xrange(size)]
  matrix = [oneRow[:] for count in xrange(size)]
  scale = (maxs[0] - mins[0]) / size
  radiusScaled = radius / scale
  for point in projectPts:
    locationX = (point[0] - mins[0]) / scale
    locationY = (point[1] - mins[1]) / scale
    xMin = int(math.floor(locationX - radiusScaled))
    yMin = int(math.floor(locationY - radiusScaled))
    #this makes boxes instead of circles... so lazy
    for xPt in xrange(xMin, xMin + int(radiusScaled * 2)):
      for yPt in xrange(yMin, yMin + int(radiusScaled * 2)):
        matrix[xPt][yPt] = 1
  retMatrix = []
  for row in matrix:
    convertedRow = listToString(row)
    retMatrix.append(convertedRow)
  return retMatrix

def size2dSquare(pointList, radius):
  '''takes list of 2d points. return min & max, make them square'''
  mins = [pointList[0][0] - radius, pointList[0][1] - radius]  # first point to init
  maxs = [pointList[0][0] + radius, pointList[0][1] + radius]
  for point in pointList[1:]:  # now do the rest
    for dimension in xrange(2):
      if point[dimension] - radius < mins[dimension]:
        mins[dimension] = point[dimension] - radius
      if point[dimension] + radius > maxs[dimension]:
        maxs[dimension] = point[dimension] + radius
  #okay now make this square, left and top (mins) aligned
  xSize = maxs[0] - mins[0]
  ySize = maxs[1] - mins[1]
  if ySize > xSize:
    maxs[0] = mins[0] + ySize
  else:
    maxs[1] = mins[1] + xSize
  return mins, maxs

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
