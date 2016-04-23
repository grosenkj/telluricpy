import numpy as np, modelTools as mT, SimPEG as simpeg, vtk, sys, os, time
import vtk.util.numpy_support as npsup

import polydata, extraction


# Function in vtk to work with dataset
def cell2vtp(obj,ind):
    """
        Function to return polygon from a cell in a data object.
        obj has to have a 'id' cell array
    """
    thresh = vtk.vtkThreshold()
    thresh.SetInputData(obj)
    thresh.ThresholdBetween(ind-.1,ind+.1)
    thresh.SetInputArrayToProcess(1, 0, 0, 0, "id")
    thresh.Update()
    vtpObj = extraction.vtu2vtp(thresh.GetOutput())

    return polydata.normFilter(polydata.triangulatePolyData(vtpObj))


