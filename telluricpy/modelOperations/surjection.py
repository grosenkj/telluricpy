

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

    # Check the type of DataSet2
    if np.all(np.array([vtkDataSet2.GetCellType(i) for i in range(vtkDataSet2.GetNumberOfCells())]) == 11):
        useBox = True
    # Extract the indieces of the object
    indArr = npsup.vtk_to_numpy(vtkDataSet2.GetCellData().GetArray('id'))
    # Want to end up with nC2xnC1 sparse matrix
    iL = []
    jL = []
    valL = []
    for iV in indArr:
        # Get the base cell and calc volume.
        if useBox:
            jT,valT = _calculateVolumeByBoxClip(vtkDataSet1,vtkDataSet2,iV)
        else:
            jT,valT = _calculateVolumeByBoolean(vtkDataSet1,vtkDataSet2,iV)
        # Insert the data
        iL.append(np.ones_like(jT)*iV)
        jL.append(jT)
        valL.append(valT)
    # Return the matrix
    i = np.concatenate(iL)
    j = np.concatenate(jL)
    val = np.concatenate(valL)
    # Make the
    outMat = simpeg.sp.csr_matrix((val,(i,j)),shape=(vtkDataSet2.GetNumberOfCells(),vtkDataSet1.GetNumberOfCells()),dtype=float)
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

    # Triangulate polygon and calc normals
    baseC = vtkTools.dataset.getCell2vtp(vtkDataSet2,iV)
    # Extrct the cell by the bounds first, significant speed up...
    if vtkDataSet1.IsA('vtkRectilinearGrid'):
        # Use a faster implementation of the rectgrid extraction
        extCells = _extractRectGridByBounds(vtkDataSet1,baseC)
    else:
        extCells = vtkTools.extraction.extractDataSetByBounds(vtkDataSet1,baseC)
    # Define the box clip and clip the first mesh with the cell
    cb = baseC.GetBounds()
    boxClip = vtk.vtkBoxClipDataSet()
    boxClip.SetInputConnection(extCells.GetOutputPort())
    boxClip.SetBoxClip(cb[0],cb[1],cb[2],cb[3],cb[4],cb[5])
    boxClip.Update()
    # Get the intersect grid
    intCells = boxClip.GetOutput()
    idList = npsup.vtk_to_numpy(intCells.GetCellData().GetArray('id'))
    uniIDs = np.unique(idList)
    # Extract cells from the first mesh that intersect the base cell
    # Calculate the volumes of the clipped cells and insert to the matrix
    volList =[]
    for i in range(intCells.GetNumberOfCells()):
        c = intCells.GetCell(i)
        cPts = c.GetPoints()
        volList.append(c.ComputeVolume(cPts.GetPoint(0),cPts.GetPoint(1),cPts.GetPoint(2),cPts.GetPoint(3)))
    volArr =  np.array(volList)
    # Calculate the volumes
    volCal = np.array([np.sum(volArr[idList == curId]) for curId in uniIDs])

    return uniIDs, volCal/np.sum(volCal)

def _extractRectGridByBounds(vtrObj,boundObj):
    '''
    Function that extracts cell from a rectilinear grid (vtr) using bounds.

    Should be signifacantly faster the extractBounds method.

    '''
    import numpy as np, SimPEG as simpeg, vtk
    import vtk.util.numpy_support as npsup

    bO = boundObj.GetBounds()
    xC = npsup.vtk_to_numpy(vtrObj.GetXCoordinates())
    yC = npsup.vtk_to_numpy(vtrObj.GetYCoordinates())
    zC = npsup.vtk_to_numpy(vtrObj.GetZCoordinates())
    iL = np.where(xC <= bO[0])[0][-1]
    iU = np.where(xC >= bO[1])[0][0]
    jL = np.where(yC <= bO[2])[0][-1]
    jU = np.where(yC >= bO[3])[0][0]
    kL = np.where(zC <= bO[4])[0][-1]
    kU = np.where(zC >= bO[5])[0][0]
    extRect = vtk.vtkExtractRectilinearGrid()
    extRect.SetInputData(vtrObj)
    extRect.SetVOI((iL,iU,jL,jU,kL,kU))
    extRect.Update()
    return extRect

def _calculateVolumeByBoolean(vtkDataSet1,vtkDataSet2,iV):
    """
    Function to calculate the volumes of a cell intersecting a mesh.

    Uses a boolean polydata filter to calculate the intersection,
    a general implementation but slow.

    """
    import numpy as np, SimPEG as simpeg, vtk
    import vtk.util.numpy_support as npsup

    from telluricpy import vtkTools

    # Triangulate polygon and calc normals
    baseC = vtkTools.dataset.getCell2vtp(vtkDataSet2,iV)
    baseVol = vtkTools.polydata.calculateVolume(baseC)
    # print iV, baseVol
    # Extract cells from the first mesh that intersect the base cell
    extractCells = vtkTools.extraction.extractDataSetWithPolygon(vtkDataSet1,baseC,extInside=True,extBoundaryCells=True,extractBounds=True)
    extInd = npsup.vtk_to_numpy(extractCells.GetCellData().GetArray('id'))
    # print extInd
    # Assert if there are no cells cutv
    assert extractCells.GetNumberOfCells() > 0, 'No cells in the clip, cell id {:d}'.format(iV)
    # Calculate the volumes of the clipped cells and insert to the matrix
    volL = []
    for nrCC,iR in enumerate(extInd):
        tempCell = vtkTools.dataset.thresholdCellId2vtp(extractCells,iR)
        # Find the intersection of the 2 cells
        boolFilt = vtk.vtkBooleanOperationPolyDataFilter()
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
                volL.append(volVal)
    return extInd,np.array(volL)

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
    iL = []
    jL = []
    valL = []
    for nr,(dist, ix) in enumerate(zip( distances, ixs )):
        iL.append(nr*np.ones_like(ix))
        jL.append(ix)
        if nrofNear == 1:
            # Only using one point
            valL.append(np.ones_like(ix))
        elif dist[0] < 1e-10:
            # If points are the "same"
            valL.append(np.ones_like(ix))
        else:  # weight z s by 1/dist --
            w = 1. / dist**p
            w /= np.sum(w)
            valL.append(w)
    # Return the matrix
    i = np.concatenate(iL)
    j = np.concatenate(jL)
    val = np.concatenate(valL)
    # Make the
    outMat = simpeg.sp.csr_matrix((val,(i,j)),shape=(vtkDataSet2.GetNumberOfCells(),vtkDataSet1.GetNumberOfCells()),dtype=float)
    return outMat
