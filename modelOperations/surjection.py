import numpy as np, SimPEG as simpeg, vtk, sys, os, time
import vtk.util.numpy_support as npsup

import vtkTools

# Collection of functions that surject models (nd-arrays).

def volumemetricSurject(vtkData1,vtkData2):
    """
    Function to calculate and return a sparse matrix of projection from vtkData1 to vtkData2.

    Based on a volumemetric estimation of intersecting cells.

    """

    sp = simpeg.sp

    # Iniate a sparse matrix
    outMat = sp.csc_matrix((vtkData2.GetNumberOfCells(),vtkData1.GetNumberOfCells()),dtype=float)

    # Extract the indieces of the object
    indArr = npsup.vtk_to_numpy(vtkData2.GetCellData().GetArray('id'))
    # Want to end up with nC1xnC2 sparse matrix
    for iV in indArr:
        # Get the base cell and calc volume. Triangulate polygon and calc normals

        baseC = _cell2vtp(vtkData2,iV)
        baseVol = calculateCellVolume(baseC)
        print iV, baseVol
        # Cut the first mesh with base cell
        clipVols = vtkTools.extration.clipDataSetWithPolygon(vtkData1,baseC,insideOut=True,extractBounds=True)
        clipInd = npsup.vtk_to_numpy(clipVols.GetCellData().GetArray('id'))
        print clipInd
        # Assert if there are no cells cutv
        assert clipVols.GetNumberOfCells() > 0, 'No cells in the clip, cell id {:d}'.format(iV)
        # Calculate the volumes of the clipped cells and insert to the matrix
        for nrCC,iR in enumerate(clipInd):
            tempCell = _cell2vtp(clipVols,nrCC)
            tempVol = calculateCellVolume(tempCell)
            # Calculate the volume
            volVal = tempVol/baseVol
            print iR, tempVol, volVal
            # Insert the value
            outMat[iV,iR] = volVal

    # Return the matrix
    return outMat