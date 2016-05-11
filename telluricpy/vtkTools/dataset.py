import numpy as np, modelTools as mT, SimPEG as simpeg, vtk, sys, os, time
import vtk.util.numpy_support as npsup

import polydata, extraction


# Function in vtk to work with dataset
def thresholdCellId2vtp(obj,ind):
    """
        Function to return polygon from a cell in a data object.
        obj has to have a 'id' cell array
    """
    thresh = vtk.vtkThreshold()
    thresh.SetInputData(obj)
    thresh.ThresholdBetween(ind-.1,ind+.1)
    # The numbers are: 1- idx, 2-port, 3-connection, 4-fieldAssociation, 5-name
    thresh.SetInputArrayToProcess(0, 0, 0, 1, "id")
    thresh.Update()
    vtpObj = extraction.vtu2vtp(thresh.GetOutput())

    return polydata.normFilter(polydata.triangulatePolyData(vtpObj))


def getCell2vtp(obj,ind):
    """
        Function gets a cell by ind and constructs a polydata from it.

    """

    # Get the cell
    cE = obj.GetCell(ind)
    # Make the polygon
    if cE.GetCellType() == 11:
        # Use a cubeSource, much faster
        cube = vtk.vtkCubeSource()
        cube.SetBounds(cE.GetBounds())
        cube.Update()
        vtpObj = cube.GetOutput()
    else:
        polygons = vtk.vtkCellArray()
        for iF in range(cE.GetNumberOfFaces()):
            f = cE.GetFace(iF)
            poly = vtk.vtkPolygon()
            poly.GetPointIds().SetNumberOfIds(f.GetNumberOfPoints())
            for nr in range(f.GetNumberOfPoints()):
                poly.GetPointIds().SetId(nr,f.GetPointId(nr))
            polygons.InsertNextCell(poly)
        # Build the polydata
        vtpObj = vtk.vtkPolyData()
        vtpObj.SetPoints(obj.GetPoints())
        vtpObj.SetPolys(polygons)

    return polydata.normFilter(polydata.triangulatePolyData(vtpObj))
