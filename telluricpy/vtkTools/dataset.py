import numpy as np, SimPEG as simpeg, vtk, sys, os, time
import vtk.util.numpy_support as npsup

import polydata, extraction

# Functions that take any vtkData object as input.


# Function in vtk to work with dataset
def thresholdCellId2vtp(vtkObj,ind):
    """
        Function to return polygon from a cell in a data object.
        vtkObj has to have a 'id' cell array
    """

    thresObj = thresFilt(vtkObj,'ind',[ind-.1,ind+.1],thType='Between')
    vtpObj = extraction.vtu2vtp(thresObj)

    return polydata.normFilter(polydata.triangulatePolyData(vtpObj))


def getCell2vtp(vtkObj,ind):
    """
        Function gets a cell by ind and constructs a polydata from it.

    """

    # Get the cell
    cE = vtkObj.GetCell(ind)
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

# Add a numpy array to a VTKobject
def addNPDataArrays(vtkObj,arrDict,arrType='Cell'):
    """ Function to add a nparray to vtkObject"""
    for nameArr,npArr in arrDict.iteritems():
        vtkArr = npsup.numpy_to_vtk(npArr,deep=1)
        vtkArr.SetName(nameArr)
        if arrType == 'Cell':
            vtkObj.GetCellData().AddArray(vtkArr)
        elif arrType == 'Point':
            vtkObj.GetPointData().AddArray(vtkArr)
        else:
            raise Exception('Not a support arrType')

def getDataArrayNames(vtkObj,arrType='Cell'):
    """Function that returns a list of all the names of cell/point data arrays."""
    l = []
    if arrType == 'Cell':
        nameList = [ vtkObj.GetCellData().GetArrayName(i) for i in range(vtkObj.GetCellData().GetNumberOfArrays())]
    elif arrType == 'Point':
        nameList = [ vtkObj.GetPointData().GetArrayName(i) for i in range(vtkObj.GetPointData().GetNumberOfArrays())]
    else:
        raise Exception('Not a support arrType')
    return nameList

def getDataArray(vtkObj,name,arrType='Cell'):
    """Function that returns the cell/point data array. """
    return npsup.vtk_to_numpy(vtkObj.GetCellData().GetArray(name))
    if arrType == 'Cell':
        return npsup.vtk_to_numpy(vtkObj.GetCellData().GetArray(name))
    elif arrType == 'Point':
        return npsup.vtk_to_numpy(vtkObj.GetPointData().GetArray(name))
    else:
        raise Exception('Not a support arrType')


def thresFilt(vtkObj,arrName,value,thType='Upper'):
    thresFilt = vtk.vtkThreshold()
    thresFilt.SetInputData(vtkObj)
    if thType in 'Upper':
        thresFilt.ThresholdByUpper(value)
    elif thType in 'Lower':
        thresFilt.ThresholdByLower(value)
    elif thType in 'Between':
        thresFilt.ThresholdBetween(value[0],value[1])
    thresFilt.AllScalarsOn()
    thresFilt.SetInputArrayToProcess(0,0,0,vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS,arrName)
    thresFilt.Update()
    return thresFilt.GetOutput()