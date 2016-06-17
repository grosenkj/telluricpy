import numpy as np, SimPEG as simpeg, vtk
import vtk.util.numpy_support as npsup


# Simple write model functions.
def writeVTPFile(fileName,vtkPolyObject):
    '''Function to write vtk polydata file (vtp).'''
    polyWriter = vtk.vtkXMLPolyDataWriter()
    if float(vtk.VTK_VERSION.split('.')[0]) >=6:
        polyWriter.SetInputData(vtkPolyObject)
    else:
        polyWriter.SetInput(vtkPolyObject)
    polyWriter.SetFileName(fileName)
    polyWriter.Update()

def writeVTUFile(fileName,vtkUnstructuredGrid,compress=True):
    '''Function to write vtk unstructured grid (vtu).'''
    Writer = vtk.vtkXMLUnstructuredGridWriter()
    if float(vtk.VTK_VERSION.split('.')[0]) >=6:
        Writer.SetInputData(vtkUnstructuredGrid)
    else:
        Writer.SetInput(vtkUnstructuredGrid)
    if not compress:
        Writer.SetCompressorTypeToNone()
        Writer.SetDataModeToAscii()
    Writer.SetFileName(fileName)
    Writer.Update()

def writeVTRFile(fileName,vtkRectilinearGrid):
    '''Function to write vtk rectilinear grid (vtr).'''
    Writer = vtk.vtkXMLRectilinearGridWriter()
    if float(vtk.VTK_VERSION.split('.')[0]) >=6:
        Writer.SetInputData(vtkRectilinearGrid)
    else:
        Writer.SetInput(vtkRectilinearGrid)
    Writer.SetFileName(fileName)
    Writer.Update()

def writeVTSFile(fileName,vtkStructuredGrid):
    '''Function to write vtk structured grid (vts).'''
    Writer = vtk.vtkXMLStructuredGridWriter()
    if float(vtk.VTK_VERSION.split('.')[0]) >=6:
        Writer.SetInputData(vtkStructuredGrid)
    else:
        Writer.SetInput(vtkStructuredGrid)
    Writer.SetFileName(fileName)
    Writer.Update()

def readVTSFile(fileName):
    '''Function to read vtk structured grid (vts) and return a grid object.'''
    Reader = vtk.vtkXMLStructuredGridReader()
    Reader.SetFileName(fileName)
    Reader.Update()
    return Reader.GetOutput()

def readVTUFile(fileName):
    '''Function to read vtk structured grid (vtu) and return a grid object.'''
    Reader = vtk.vtkXMLUnstructuredGridReader()
    Reader.SetFileName(fileName)
    Reader.Update()
    return Reader.GetOutput()

def readVTRFile(fileName):
    '''Function to read vtk structured grid (vtr) and return a grid object.'''
    Reader = vtk.vtkXMLRectilinearGridReader()
    Reader.SetFileName(fileName)
    Reader.Update()
    return Reader.GetOutput()

def readVTPFile(fileName):
    '''Function to read vtk structured grid (vtp) and return a grid object.'''
    Reader = vtk.vtkXMLPolyDataReader()
    Reader.SetFileName(fileName)
    Reader.Update()
    return Reader.GetOutput()