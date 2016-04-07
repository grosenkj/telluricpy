import numpy as np, modelTools as mT, SimPEG as simpeg, vtk, sys, os, time
import vtk.util.numpy_support as npsup
# sys.path.append('/home/gudni/Dropbox/code/python/vtkTools/')
from vtkTools.polygons import convertToImplicitPolyDataDistance

# Clip and extrude the topo
def clipDataSetWithPolygon(vtkDataSet,vtkPoly,returnImpDist=False,insideOut=True,extractBounds=False):
    """
    Function to clips cells from a vtkDataSet, given polygon/s in a vtkPolyData.
    Returns a clipped cells that fall inside or outside the polygon boundary.

    """

    # Make a implicit function
    impDist = convertToImplicitPolyDataDistance(vtkPoly)
    # Reduce the data to the bounds
    if extractBounds:
        extBoundsFilt = extractDataSetByBounds(vtkDataSet,vtkPoly)
    else:
        extBoundsFilt = extractDataSetByBounds(vtkDataSet,vtkDataSet)

    if vtkDataSet.IsA('vtkPolyData'):
        clipFilt = vtk.vtkClipPolyData()
    else:
        clipFilt = vtk.vtkClipDataSet()
    clipFilt.SetInsideOut(insideOut+0)
    clipFilt.SetInputConnection(extBoundsFilt.GetOutputPort())
    clipFilt.SetClipFunction(impDist)
    clipFilt.SetOutputPointsPrecision(vtk.vtkAlgorithm.DOUBLE_PRECISION)
    clipFilt.SetMergeTolerance(0.000001)
    clipFilt.Update()
    if returnImpDist:
        return clipFilt.GetOutput(), impDist
    else:
        return clipFilt.GetOutput()

# Cut and extrude the topo
def cutDataSetWithPolygon(vtkDataSet,vtkPoly,returnImpDist=False,extractBounds=False):
    """
    Function to cuts cells from a vtkDataSet, given polygon/s in a vtkPolyData.
    Returns a 1 dimension reduction of the cut cell the are intersected by the polygon.
    (3D cell to 2D polygon; 2D polygon to lines.)

    """
    # Make a implicit function
    impDist = convertToImplicitPolyDataDistance(vtkPoly)
    # Reduce the data to the bounds
    if extractBounds:
        extBoundsFilt = extractDataSetByBounds(vtkDataSet,vtkPoly)
    else:
        extBoundsFilt = extractDataSetByBounds(vtkDataSet,vtkDataSet)
    # Define the cutter filter
    cutFilt = vtk.vtkCutter()
    cutFilt.SetInputConnection(extBoundsFilt.GetOutputPort())
    cutFilt.SetCutFunction(impDist)
    cutFilt.SetOutputPointsPrecision(vtk.vtkAlgorithm.DOUBLE_PRECISION)
    cutFilt.Update()
    if returnImpDist:
        return cutFilt.GetOutput(), impDist
    else:
        return cutFilt.GetOutput()

# Extract cells with a polygon, needs to have normals to work
def extractDataSetWithPolygon(vtkDataSet,vtkPoly,returnImpDist=False,extInside=True,extBoundaryCells=True,extractBounds=False):
    """
    Function to extract cells from a vtkDataSet, given polygon/s in a vtkPolyData.
    Returns a full cells that fall fully inside/outside and/or on the polygon boundary.

    """

    # Make a implicit function
    impDist = convertToImplicitPolyDataDistance(vtkPoly)
    # Reduce the data to the bounds
        # Reduce the data to the bounds
    if extractBounds:
        extBoundsFilt = extractDataSetByBounds(vtkDataSet,vtkPoly)
    else:
        extBoundsFilt = extractDataSetByBounds(vtkDataSet,vtkDataSet)
    # If input is a vtkPolyData
    if vtkDataSet.IsA('vtkPolyData'):
        extractFilt = vtk.vtkExtractPolyDataGeometry()
    else:
        extractFilt = vtk.vtkExtractGeometry()
    extractFilt.SetExtractInside(extInside+0)
    extractFilt.SetExtractBoundaryCells(extBoundaryCells+0)
    extractFilt.SetInputConnection(extBoundsFilt.GetOutputPort())
    extractFilt.SetImplicitFunction(impDist)
    extractFilt.Update()
    if returnImpDist:
        return extractFilt.GetOutput(), impDist
    else:
        return extractFilt.GetOutput()

def extractDataSetByBounds(fullvtkDataSet,boundvtkDataSet):
    """
    Function to extract from a vtkDataSet within bounds of another vtkDataSet.

    Returns a extrct filter

    """

    # Define a bounding box implicit
    if boundvtkDataSet.IsA('vtkDataSet'):
        impFunc = vtk.vtkBox()
        impFunc.SetBounds(boundvtkDataSet.GetBounds())
    elif boundvtkDataSet.IsA('vtkImplicitFunction'):
        impFunc = boundvtkDataSet
    # Extract all inside and boundaries
    if fullvtkDataSet.IsA('vtkPolyData'):
        extractFilt = vtk.vtkExtractPolyDataGeometry()
    else:
        extractFilt = vtk.vtkExtractGeometry()
    extractFilt.SetExtractInside(1)
    extractFilt.SetExtractBoundaryCells(1)
    extractFilt.SetInputData(fullvtkDataSet)
    extractFilt.SetImplicitFunction(impFunc)
    extractFilt.Update()
    return extractFilt


def sumImplicitFunc(impFunctions):
    impBool = vtk.vtkImplicitBoolean()
    impBool.SetOperationTypeToDifference()
    for impFunc in impFunctions:
        impBool.AddFunction(impFunc)

    return impBool

def makePolyhedronCell(vtkPolyData,returnGrid=False):
    """ Function that makes polyhedron cell from polygon. """

    # ToDo:
    # Add check that the polygon is "Waterthight"
    # Not have to write a file, rather declare the polyhedron directly.
    # Extract the information needed from the poly data
    ptsIds = vtk.vtkIdList()
    numCellFaces = vtkPolyData.GetNumberOfCells()
    ptsIds.InsertNextId(numCellFaces) # Number of cell faces
    # Add the faces
    for cF in range(numCellFaces):
        numPtsInFace = vtkPolyData.GetCell(cF).GetNumberOfPoints()
        ptsIds.InsertNextId(numPtsInFace)
        for cPF in range(numPtsInFace):
            ptsIds.InsertNextId(vtkPolyData.GetCell(cF).GetPointId(cPF))

    # Make the grid
    UnstructPolyHed = vtk.vtkUnstructuredGrid()
    UnstructPolyHed.SetPoints(vtkPolyData.GetPoints())
    UnstructPolyHed.InsertNextCell(vtk.VTK_POLYHEDRON,ptsIds)

    vtkPolyhed = UnstructPolyHed.GetCell(0)
    if returnGrid:
        return UnstructPolyHed
    else:
        return vtkPolyhed

def geometryFilt(vtkObject):
    geoFilt = vtk.vtkGeometryFilter()
    geoFilt.SetInputData(vtkObject)
    geoFilt.Update()
    return geoFilt.GetOutput()

def vtu2vtp(vtuObject):
    vtu2vtpFilt = vtk.vtkDataSetSurfaceFilter()
    vtu2vtpFilt.SetInputData(vtuObject)
    vtu2vtpFilt.Update()
    return vtu2vtpFilt.GetOutput()

def vtp2vtuPolyhedron(vtpObject):
    """
    Convert a Polydata to individual polyhedron cells in a vtu grid
    """
    # Find individual
    conFilt = vtk.vtkPolyDataConnectivityFilter()
    conFilt.SetInputData(vtpObject)
    conFilt.SetExtractionMode(5)
    conFilt.SetColorRegions(1)
    conFilt.Update()
    #
    phAppFilt = vtk.vtkAppendFilter()
    # phAppFilt.MergePointsOn()
    for nr in np.arange(conFilt.GetNumberOfExtractedRegions()):
        thresh = vtk.vtkThreshold()
        thresh.SetInputConnection(conFilt.GetOutputPort())
        thresh.ThresholdBetween(nr-.1,nr+.1)
        thresh.SetInputArrayToProcess(1, 0, 0, 0, "Regionid")
        thresh.Update()
        # Convert to a Polyhedron and add to the append filter
        phAppFilt.AddInputData(makePolyhedronCell(vtu2vtp(thresh.GetOutput()),returnGrid=True))
    phAppFilt.Update()
    # Return the grid
    return phAppFilt.GetOutput()
