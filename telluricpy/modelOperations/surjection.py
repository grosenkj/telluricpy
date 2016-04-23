

# Collection of functions that surject models (nd-arrays).

def getVolumemetricSurjectMatrix(vtkDataSet1,vtkDataSet2):
    """
    Function to calculate and return a sparse matrix of surjection from vtkDataSet1 to vtkDataSet2.

    Based on a volumemetric estimation of intersecting cells.

    """
    import numpy as np, SimPEG as simpeg, vtk
    import vtk.util.numpy_support as npsup

    from telluricpy import vtkTools
    sp = simpeg.sp

    # Iniate a sparse matrix
    outMat = sp.csr_matrix((vtkDataSet2.GetNumberOfCells(),vtkDataSet1.GetNumberOfCells()),dtype=float)

    # Check the type of DataSet2
    if np.all(np.array([vtkDataSet2.GetCellType(i) for i in range(vtkDataSet2.GetNumberOfCells())]) == 11):
        useBox = True
    # Extract the indieces of the object
    indArr = npsup.vtk_to_numpy(vtkDataSet2.GetCellData().GetArray('id'))
    # Want to end up with nC2xnC1 sparse matrix
    for iV in indArr:
        # Get the base cell and calc volume.
        if useBox:
            volDict = _calculateVolumeByBoxClip(vtkDataSet1,vtkDataSet2,iV)
        else:
            volDict = _calculateVolumeByBoolean(vtkDataSet1,vtkDataSet2,iV)
        # Insert the data
        for iR in volDict.iterkeys():
            outMat[iV,iR] = volDict[iR]
    # Return the matrix
    return outMat

def _calculateVolumeByBoxClip(vtkDataSet1,vtkDataSet2,iV):
    """
    Function to calculate the volumes of a cell intersecting a mesh.

    Uses a BoxClip filter to calculate the intersection, assumes the cell to
    be a voxel (=11) and aligned to the grid. This should work well for rectilinear
    and unstructured grids with all cells as voxels

    """
    import numpy as np, SimPEG as simpeg, vtk
    import vtk.util.numpy_support as npsup

    from telluricpy import vtkTools
    # Make the outDict
    outDict = {}
    # Triangulate polygon and calc normals
    baseC = vtkTools.dataset.cell2vtp(vtkDataSet2,iV)
    baseVol = vtkTools.polydata.calculateVolume(baseC)
    # Define the box clip and clip the first mesh with the cell
    cb = baseC.GetBounds()
    boxClip = vtk.vtkBoxClipDataSet()
    boxClip.SetInputData(vtkDataSet1)
    boxClip.SetBoxClip(cb[0],cb[1],cb[2],cb[3],cb[4],cb[5])
    boxClip.Update()
    # Get the intersect grid
    intCells = boxClip.GetOutput()
    idList = npsup.vtk_to_numpy(intCells.GetCellData().GetArray('id'))
    uniIDs = np.unique(idList)
    # Extract cells from the first mesh that intersect the base cell
    # Calculate the volumes of the clipped cells and insert to the matrix
    for nrCC,iR in enumerate(uniIDs):
        vol =  vtkTools.polydata.calculateVolume(vtkTools.dataset.cell2vtp(intCells,iR))
        outDict[iR] = float(vol)/float(baseVol)
    return outDict

def _calculateVolumeByBoolean(vtkDataSet1,vtkDataSet2,iV):
    """
    Function to calculate the volumes of a cell intersecting a mesh.

    Uses a boolean polydata filter to calculate the intersection,
    a general implementation but slow.

    """
    import numpy as np, SimPEG as simpeg, vtk
    import vtk.util.numpy_support as npsup

    from telluricpy import vtkTools
    # Make the outDict
    outDict = {}

    # Triangulate polygon and calc normals
    baseC = vtkTools.dataset.cell2vtp(vtkDataSet2,iV)
    baseVol = vtkTools.polydata.calculateVolume(baseC)
    # print iV, baseVol
    # Extract cells from the first mesh that intersect the base cell
    extractCells = vtkTools.extraction.extractDataSetWithPolygon(vtkDataSet1,baseC,extInside=True,extBoundaryCells=True,extractBounds=True)
    extInd = npsup.vtk_to_numpy(extractCells.GetCellData().GetArray('id'))
    # print extInd
    # Assert if there are no cells cutv
    assert extractCells.GetNumberOfCells() > 0, 'No cells in the clip, cell id {:d}'.format(iV)
    # Calculate the volumes of the clipped cells and insert to the matrix
    for nrCC,iR in enumerate(extInd):
        tempCell = vtkTools.dataset.cell2vtp(extractCells,iR)
        # Find the intersection of the 2 cells
        boolFilt = vtk.vtkBooleanOperationPolyDataFilter()
        boolFilt.SetTolerance(1e-2)
        boolFilt.SetInputData(0,tempCell)
        boolFilt.SetInputData(1,baseC)
        boolFilt.SetOperationToIntersection()
        # If they intersect, calculate the volumes
        if boolFilt.GetOutput().GetNumberOfPoints() > 0:
            cleanInt = vtkTools.polydata.cleanPolyData(boolFilt.GetOutputPort())
            del3dFilt = vtk.vtkDelaunay3D()
            del3dFilt.SetInputData(cleanInt)
            del3dFilt.Update()
            # Get the output
            intC = vtkTools.extraction.vtu2vtp(del3dFilt.GetOutput())
            intVol = vtkTools.polydata.calculateVolume(tempCell)
            # Calculate the volume
            volVal = intVol/baseVol
            # print iR, intVol, volVal
            # Insert the value
            if volVal > 0.0:
                outDict[iR] = volVal/baseVol

def getIWDSurjectMatrix(vtkDataSet1, vtkDataSet2, leafsize = 10, nrofNear = 9, eps=0, p=1.):
    """
    Function to calculate and return a sparse matrix of surjection form vtkDataSet1 to vtkDataSet2

    Input:
        vtkDataSet1 - vtk DataSet - Mesh with cells
        vtkDataSet2 - vtk DataSet - Mesh with cells

    Output:
        scipy sparce matrix - Weight inverse distance values
    """


    # Import packages
    import numpy as np, SimPEG as simpeg, vtk
    import vtk.util.numpy_support as npsup

    from scipy.spatial import cKDTree as KDTree



    # Prep the matrices
    # Get cell centers
    ds1CC = vtk.vtkCellCenters()
    ds1CC.SetInputData(vtkDataSet1)
    ds1CC.Update()
    ds1CCpts = ds1CC.GetOutput()
    ds1Arr = npsup.vtk_to_numpy(ds1CCpts.GetPoints().GetData())

    ds2CC = vtk.vtkCellCenters()
    ds2CC.SetInputData(vtkDataSet2)
    ds2CC.Update()
    ds2CCpts = ds2CC.GetOutput()
    ds2Arr = npsup.vtk_to_numpy(ds2CCpts.GetPoints().GetData())

    # Use a inverse distance weighting of the cell centers
    # Based on this stack overflow.
    # http://stackoverflow.com/questions/3104781/inverse-distance-weighted-idw-interpolation-with-python

    # Build the KDtree
    KDtree = KDTree( ds1Arr, leafsize = leafsize )
    # Calculate the distances and indexes
    distances, ixs = KDtree.query( ds2Arr, k=nrofNear, eps=eps )
    outMat = simpeg.sp.csr_matrix((vtkDataSet2.GetNumberOfCells(),vtkDataSet1.GetNumberOfCells()),dtype=float)
    for nr,(dist, ix) in enumerate(zip( distances, ixs )):
        if nrofNear == 1:
            # Only using one point
            outMat[nr,ix] = 1
        elif dist[0] < 1e-10:
            # If points are the "same"
            outMat[nr,ix[0]] = 1
        else:  # weight z s by 1/dist --
            w = 1. / dist**p
            w /= np.sum(w)
            outMat[nr,ix] = w
    # Return the matrix
    return outMat
