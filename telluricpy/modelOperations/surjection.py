import numpy as np, SimPEG as simpeg, vtk, sys, os, time
import vtk.util.numpy_support as npsup

import vtkTools

# Collection of functions that surject models (nd-arrays).

def getVolumemetricSurjectMatrix(vtkDataSet1,vtkDataSet2):
    """
    Function to calculate and return a sparse matrix of surjection from vtkDataSet1 to vtkDataSet2.

    Based on a volumemetric estimation of intersecting cells.

    """

    sp = simpeg.sp

    # Iniate a sparse matrix
    outMat = sp.csr_matrix((vtkDataSet2.GetNumberOfCells(),vtkDataSet1.GetNumberOfCells()),dtype=float)

    # Extract the indieces of the object
    indArr = npsup.vtk_to_numpy(vtkDataSet2.GetCellData().GetArray('id'))
    # Want to end up with nC1xnC2 sparse matrix
    for iV in indArr:
        # Get the base cell and calc volume. Triangulate polygon and calc normals

        baseC = vtkTools.dataset.c_ell2vtp(vtkDataSet2,iV)
        baseVol = vtkTools.polydata.calculateVolume(baseC)
        print iV, baseVol
        # Cut the first mesh with base cell
        clipVols = vtkTools.extration.clipDataSetWithPolygon(vtkDataSet1,baseC,insideOut=True,extractBounds=True)
        clipInd = npsup.vtk_to_numpy(clipVols.GetCellData().GetArray('id'))
        print clipInd
        # Assert if there are no cells cutv
        assert clipVols.GetNumberOfCells() > 0, 'No cells in the clip, cell id {:d}'.format(iV)
        # Calculate the volumes of the clipped cells and insert to the matrix
        for nrCC,iR in enumerate(clipInd):
            tempCell = vtkTools.dataset.cell2vtp(clipVols,nrCC)
            tempVol = vtkTools.polydata.calculateCellVolume(tempCell)
            # Calculate the volume
            volVal = tempVol/baseVol
            print iR, tempVol, volVal
            # Insert the value
            outMat[iV,iR] = volVal

    # Return the matrix
    return outMat

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
    import vtk,pdb,numpy as np, os
    from vtk.util import numpy_support as npsup
    import scipy.interpolate as sciint
    import numpy as np
    from scipy.spatial import cKDTree as KDTree



    # Prep the matrices
    # Get cell centers
    ds1CC = vtk.vtkCellCenters()
    ds1CC.SetInputData(vtkDataSet1)
    ds1CC.Update()
    ds1CCpts = compGridCC.GetOutput()
    ds1Arr = npsup.vtk_to_numpy(ds1CCpts.GetPoints().GetData())

    ds2CC = vtk.vtkCellCenters()
    ds2CC.SetInputData(vtkDataSet2)
    ds2CC.Update()
    ds2CCpts = modelCC.GetOutput()
    ds2Arr = npsup.vtk_to_numpy(ds2CCpts.GetPoints().GetData())

    # Use a inverse distance weighting of the cell centers
    # Based on this stack overflow.
    # http://stackoverflow.com/questions/3104781/inverse-distance-weighted-idw-interpolation-with-python

    # Build the KDtree
    KDtree = KDTree( ds2Arr, leafsize=10 )
    # Calculate the distances and indexes
    distances, ix = KDtree.query( ds1Arr, k=nrofNear, eps=eps )
    outMat = simpeg.sp.csr_matrix((vtkDataSet2.GetNumberOfCells(),vtkDataSet1.GetNumberOfCells()),dtype=float)
    for nr,(dist, ix) in enumerate(zip( distances, ix )):
        if nrofNear == 1:
            # Only using one point
            outMat[nr,ix] = 1
        elif dist[0] < 1e-10:
            # If points are the "same"
            outMat[nr,ix[0]] = 1
        else:  # weight z s by 1/dist --
            w = 1 / dist**p
            w /= np.sum(w)
            outMat[nr,ix] = w
    # Return the matrix
    return outMat
