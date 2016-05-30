import numpy as np, SimPEG as simpeg, vtk, sys, os, time
import vtk.util.numpy_support as npsup

# import polydata, extraction

def findZofXYOnPolydata(points,vtkPolydata):

	# Make the cell locator
	cellLocator = vtk.vtkCellLocator()
	cellLocator.SetDataSet(vtkPolydata)
	cellLocator.BuildLocator()
	# Find the min/max of the polydata.
	lbot, ltop = np.array(vtkPolydata.GetBounds())[4::]

	# Loop over all the locations.
	intersectList = []
	try:
		for nr, loc in enumerate(points):
			# Make line
			p1 = np.hstack((loc[0:2],ltop))
			p2 = np.hstack((loc[0:2],lbot))
			# Pre define variables as in C++
			t = vtk.mutable(0)
			pIntSect = [0.0, 0.0, 0.0]
			pcoords = [0.0, 0.0, 0.0]
			sub_id = vtk.mutable(0)
			cellLocator.IntersectWithLine(p1,p2,1e-6,t,pIntSect,pcoords,sub_id)
			intersectList.append(pIntSect)
	except KeyboardInterrupt as k:
		print 'Stopped at iteration {:d} in the for loop.'.format(nr)
		raise k

	# Return the intersects
	return np.array(intersectList)