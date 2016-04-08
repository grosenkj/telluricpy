import numpy as np, vtk, vtk.util.numpy_support as npsup

# Functions that make vtk object out of XYZ locations.
def makeCylinderPtsVTP(locXYZ,radius=50,height=50,res=10):
    # Load the file
    if type(locXYZ) == np.ndarray:
        loc = locXYZ
    elif type(locXYZ) == str:
        loc = np.genfromtxt(locXYZ)
    # Make append poly filter
    appPoly = vtk.vtkAppendPolyData()
    # Loop through all the locations
    for pt in loc[:,0:3]:
    	# Make the cylinters
        cyl = vtk.vtkCylinderSource()
        cyl.SetCenter(pt)
        cyl.SetRadius(radius)
        cyl.SetHeight(height)
        cyl.SetResolution(res)
        # Rotate to be vertical
        rot = vtk.vtkTransform()
        rot.Translate(-pt)
        rot.PostMultiply()
        rot.RotateX(90)
        rot.Translate(pt)
        tranFor = vtk.vtkTransformPolyDataFilter()
        tranFor.SetInputConnection(cyl.GetOutputPort())
        tranFor.SetTransform(rot)
        tranFor.Update()
        # Append
        appPoly.AddInputConnection(tranFor.GetOutputPort())
    # Update and return.
    appPoly.Update()
    return appPoly.GetOutput()