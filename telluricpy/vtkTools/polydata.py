import numpy as np, SimPEG as simpeg, vtk
import vtk.util.numpy_support as npsup
import io

# Calculate normals for polygons
def normFilter(vtkPoly):
    '''
    Calcluates normals and makes the ordering of the cells consistent (positive outwards).
    '''
    polyNormFilter = vtk.vtkPolyDataNormals()
    polyNormFilter.SetInputData(vtkPoly)
    polyNormFilter.ComputeCellNormalsOn()
    polyNormFilter.ComputePointNormalsOn()
    polyNormFilter.AutoOrientNormalsOn()
    polyNormFilter.ConsistencyOn()
    polyNormFilter.NonManifoldTraversalOff()
    polyNormFilter.SetSplitting(0)
    # polyNormFilter.FlipNormalsOn()
    polyNormFilter.Update()
    # Don't need this
    # reverseFilter = vtk.vtkReverseSense()
    # reverseFilter.SetInputConnection(polyNormFilter.GetOutputPort())
    # reverseFilter.ReverseCellsOn()
    # reverseFilter.ReverseNormalsOn()
    # reverseFilter.Update()
    return polyNormFilter.GetOutput()

def extrudePolygon(polyObject,extrudeFactor,vector=None,centered=False):
    """
    Extrude polygons by a scalar factor.
    """
    extrude = vtk.vtkLinearExtrusionFilter()
    extrude.CappingOn()
    if vector is not None:
        extrude.SetExtrusionTypeToVectorExtrusion()
        extrude.SetVector(vector)
    else:
        extrude.SetExtrusionTypeToNormalExtrusion()
        polyObject = normFilter(polyObject)
    extrude.SetInputData(polyObject)
    extrude.SetScaleFactor(extrudeFactor)
    extrude.Update()
    ext = extrude.GetOutput()
    if centered and vector is not None:
        extPts = npsup.vtk_to_numpy(ext.GetPoints().GetData())
        newPts = vtk.vtkPoints()
        newNppts = extPts - (np.array(vector)*(extrudeFactor/2))
        newPts.SetData(npsup.numpy_to_vtk(newNppts,deep=1))
        ext.SetPoints(newPts)

    return ext

def calculateVolume(polyData):
    """
    Calculates the volume of a polydata cell.

    Input:
        polyData - vtkPolyData object - Contains a cell defined by triangular polygons and
            needs to be waterthight.

    Returns:
        float - The volume of the cell in the respecitve units.
    """
    # Mass prop
    mp = vtk.vtkMassProperties()
    # Set the input
    mp.SetInputData(polyData)
    return float(mp.GetVolume())

def decimatePolygon(polyObject,reduceFactor=0.5):
    '''
    '''
    deci = vtk.vtkDecimatePro()
    deci.SetInputData(polyObject)
    deci.SetTargetReduction(reduceFactor)
    deci.BoundaryVertexDeletionOff()
    deci.PreserveTopologyOn()
    deci.SetSplitting(0)
    deci.Update()
    return deci.GetOutput()

def smoothPolyData(polyObject):
    smoother = vtk.vtkSmoothPolyDataFilter()
    smoother.SetInputData(polyObject)
    smoother.SetFeatureAngle(60.0)
    smoother.SetEdgeAngle(60.0)
    smoother.SetBoundarySmoothing(1)
    smoother.SetFeatureEdgeSmoothing(1)
    smoother.SetNumberOfIterations(50)
    smoother.Update()
    return smoother.GetOutput()

def triangulatePolyData(polyData):
    triFilter = vtk.vtkTriangleFilter() # triFilter to clean the polydata
    triFilter.SetInputData(polyData)
    triFilter.Update()
    return triFilter.GetOutput()

def convertToImplicitPolyDataDistance(polyData):
    triFilter = vtk.vtkTriangleFilter() # triFilter to clean the polydata
    triFilter.SetInputData(polyData)
    triFilter.Update()
    ImpDist = vtk.vtkImplicitPolyDataDistance()
    ImpDist.SetTolerance(1e-1)
    ImpDist.SetInput(triFilter.GetOutput())
    return ImpDist

def sinWinSmoother(polyData):
    smoother = vtk.vtkWindowedSincPolyDataFilter()
    smoother.SetInputData(polyData);
    smoother.SetNumberOfIterations(15);
    smoother.BoundarySmoothingOff();
    smoother.FeatureEdgeSmoothingOff();
    smoother.SetFeatureAngle(60.0);
    smoother.SetPassBand(0.1);
    smoother.NonManifoldSmoothingOn();
    smoother.NormalizeCoordinatesOn();
    smoother.Update()
    return smoother.GetOutput()

def transformToLocalCoords(locPts,vtkPolyObj):
    """ Function that transforms points to a local system, given by the input"""
    newPts = vtk.vtkPoints()
    for i in range(vtkPolyObj.GetPoints().GetNumberOfPoints()):
        newPts.InsertNextPoint(vtkPolyObj.GetPoints().GetPoint(i)-locPts)
    newvtkPolyObj = vtk.vtkPolyData()
    newvtkPolyObj.DeepCopy(vtkPolyObj)
    newvtkPolyObj.SetPoints(newPts)
    return newvtkPolyObj

def join2Polydata(vtkPolydata1,vtkPolydata2,threshold1='lower',threshold2='upper',outThres=False,saveIntersect=False):
    """ Function finds the intersect of two polydata and returns the append of the 2.

    """

    # NOTE: Test this to deal with inaccuracies.
    # Copy the structure from the inputs.
    vtkpolyStructure1 = vtk.vtkPolyData()
    vtkpolyStructure1.CopyStructure(vtkPolydata1)
    vtkpolyStructure2 = vtk.vtkPolyData()
    vtkpolyStructure2.CopyStructure(vtkPolydata2)

    # Make triangular mesh from the polydata.
    tri1Filt = vtk.vtkTriangleFilter()
    tri1Filt.SetInputData(vtkpolyStructure1)
    tri1Filt.Update()
    tri2Filt = vtk.vtkTriangleFilter()
    tri2Filt.SetInputData(vtkpolyStructure2)
    tri2Filt.Update()

    # Use the clean polydata filter
    if False:
        tri1CleanFilt = vtk.vtkCleanPolyData()
        tri1CleanFilt.SetInputConnection(tri1Filt.GetOutputPort())
        tri1CleanFilt.Update()
        tri2CleanFilt = vtk.vtkCleanPolyData()
        tri2CleanFilt.SetInputConnection(tri2Filt.GetOutputPort())
        tri2CleanFilt.Update()

        # Try on local grid
        locPoint = np.concatenate((np.sum(np.array(tri2Filt.GetOutput().GetBounds())[0:4].reshape((2,2)),axis=1)/2,np.array([0])))
        tri1Loc = transformToLocalCoords(locPoint,tri1CleanFilt.GetOutput())
        tri2Loc = transformToLocalCoords(locPoint,tri2CleanFilt.GetOutput())
    else:
        # Try on local grid
        locPoint = np.concatenate((np.sum(np.array(tri2Filt.GetOutput().GetBounds())[0:4].reshape((2,2)),axis=1)/2,np.array([0])))
        tri1Loc = transformToLocalCoords(locPoint,tri1Filt.GetOutput())
        tri2Loc = transformToLocalCoords(locPoint,tri2Filt.GetOutput())

    # Setup the intersect.
    polyIntersectFilt = vtk.vtkIntersectionPolyDataFilter()
    polyIntersectFilt.SplitFirstOutputOn()
    polyIntersectFilt.SplitSecondOutputOn()
    polyIntersectFilt.SetInputData(0,tri1Loc)
    polyIntersectFilt.SetInputData(1,tri2Loc)
    polyIntersectFilt.Update()

    # Temp: save the outputs
    if saveIntersect:
        io.writeVTPFile('InterSect0.vtp',polyIntersectFilt.GetOutput(0))
        io.writeVTPFile('InterSect1.vtp',polyIntersectFilt.GetOutput(1))
        io.writeVTPFile('InterSect2.vtp',polyIntersectFilt.GetOutput(2))
    # To do: need to check if they intersect.
    # Calculate the distance of the intersect
    polyDist = vtk.vtkDistancePolyDataFilter()
    polyDist.SetInputConnection(0,polyIntersectFilt.GetOutputPort(1))
    polyDist.SetInputConnection(1,polyIntersectFilt.GetOutputPort(2))
    if saveIntersect:
        polyDist.Update()
        io.writeVTPFile('distPoly0.vtp',polyDist.GetOutput(0))
        io.writeVTPFile('distPoly1.vtp',polyDist.GetOutput(1))
    poly0thr = vtk.vtkThreshold()
    poly0thr.AllScalarsOn()
    poly0thr.SetInputArrayToProcess(0,0,0,vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS,'Distance')
    poly0thr.SetInputConnection(polyDist.GetOutputPort(0))
    if threshold1=='lower':
        poly0thr.ThresholdByLower(0.0)
    else:
        poly0thr.ThresholdByUpper(0.0)
    poly0surf = vtk.vtkDataSetSurfaceFilter()
    poly0surf.SetInputConnection(poly0thr.GetOutputPort())
    poly1thr = vtk.vtkThreshold()
    poly1thr.AllScalarsOn()
    poly1thr.SetInputArrayToProcess(0,0,0,vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS,'Distance')
    poly1thr.SetInputConnection(polyDist.GetOutputPort(1))
    if threshold2=='upper':
        poly1thr.ThresholdByUpper(0.0)
    else:
        poly1thr.ThresholdByLower(0.0)
    poly1surf = vtk.vtkDataSetSurfaceFilter()
    poly1surf.SetInputConnection(poly1thr.GetOutputPort())
    # Append the polydata's
    polyAppendFilt = vtk.vtkAppendPolyData()
    polyAppendFilt.SetInputConnection(poly0surf.GetOutputPort())
    polyAppendFilt.AddInputConnection(poly1surf.GetOutputPort())
    polyAppendFilt.Update()
    # Return the appended polydata
    if outThres:
        return transformToLocalCoords(-locPoint,polyAppendFilt.GetOutput()),transformToLocalCoords(-locPoint,poly0surf.GetOutput()),transformToLocalCoords(-locPoint,poly1surf.GetOutput())
    else:
        return transformToLocalCoords(-locPoint,polyAppendFilt.GetOutput())


def cleanPolyData(polyData):
    cleanFilt = vtk.vtkCleanPolyData()
    cleanFilt.SetConvertLinesToPoints(0)
    cleanFilt.SetConvertPolysToLines(0)
    cleanFilt.SetConvertStripsToPolys(0)
    cleanFilt.SetPointMerging(1)
    cleanFilt.SetInputData(polyData)
    cleanFilt.Update()
    return cleanFilt.GetOutput()

def polyDataList2polyhedVTU(polyDataList):
    """
    Function that makes polyhedron ugrid from a list of watertight polygons stored in PolyData.



    """

    # What is needed is:
    #   allPoints -vtkPoints-
    #   pointsIds -vtkIdList- for each cell
    #   hedronFaces -vtkCellArray- with all the faces.
    # Insert cells with
    # ugrid.InsertNextCell(vtk.vtkPolyHedron,
    #           nrVerts,hedronVertPtsIds,nrFaces,faceCellArray.GetPointer())
    def makeIdList(polyData,ptsShift):
        """
        Make a vtkIdList for a polyhedron from a polydata
        """
        ptsIds = vtk.vtkIdList()
        numCellFaces = polyData.GetNumberOfCells()
        ptsIds.InsertNextId(numCellFaces) # Number of cell faces
        # Add the faces
        for cF in range(numCellFaces):
            numPtsInFace = polyData.GetCell(cF).GetNumberOfPoints()
            ptsIds.InsertNextId(numPtsInFace)
            for cPF in range(numPtsInFace):
                ptsIds.InsertNextId(polyData.GetCell(cF).GetPointId(cPF)+ptsShift)
        return ptsIds
    # Start
    unstructGrid = vtk.vtkUnstructuredGrid()
    nrPts = 0
    hedronList = []
    for vtkPolyData in polyDataList:
        pts = npsup.vtk_to_numpy(vtkPolyData.GetPoints().GetData())
        nrVert = pts.shape[0]
        cellIdList = makeIdList(vtkPolyData,nrPts)
        unstructGrid.InsertNextCell(vtk.VTK_POLYHEDRON,cellIdList)

        try:
            npPts = np.concatenate((npPts,pts))
        except NameError as e:
            npPts = pts.copy()
        nrPts = nrPts + pts.shape[0]

    # Make the grid
    # Set the points
    vtkPts = vtk.vtkPoints()
    vtkPts.SetData(npsup.numpy_to_vtk(npPts,deep=1))
    unstructGrid.SetPoints(vtkPts)

    # Return
    return unstructGrid