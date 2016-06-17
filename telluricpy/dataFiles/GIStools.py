# Script to convert GIS .shp and .dbt files to VTK vtp files.
import vtk, numpy as np, vtk.util.numpy_support as npsup, sys
from telluricpy import vtkTools

def makePolyhedron(polygon,thickness=1,elevation=0,triangulate=False,returnGrid=False):
    """
    Function to make a polyhedron from a polygon.

    Inputs
        polygon - array list of vertices in clockwise order
        thickness - float, thinkess of the polyhedron
        elevation - float or array, elevation
    """

    volpolyPolydata = makeVolumePolygon(polygon,thickness=thickness,elevation=elevation,triangulate=triangulate,cap=True)

    # Make the polyhedron
    if returnGrid:
        return _makePolyhedronCell(volpolyPolyData,True)
    else:
        return _makePolyhedronCell(volpolyPolyData)

def _makePolyhedronCell(vtkPolyData,returnGrid=False):
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

def makePolygon(polygon,elevation=0,triangulate=False):
    """
    Function to make a 2D vtk polygon PolyData from points making up a polygon.

    Inputs
        polygon - array list of vertices in clockwise order
        elevation - float, elevation of the polygon
        triangulate - boolean, triangulate the surface

    Output
        polyData - returns polydata of the surface polygon
    """

    # Make a enclosed polygons
    nrBasepts = len(polygon)
    ptsArr = np.hstack((polygon,np.ones((polygon.shape[0],1))*elevation))
    ppts = vtk.vtkPoints()
    ppts.SetData(npsup.numpy_to_vtk(ptsArr,deep=1))
    # Make the cells array of polygons
    polys = vtk.vtkCellArray()
    poly = vtk.vtkPolygon()
    poly.GetPointIds().SetNumberOfIds(len(ptsArr))
    for ind in np.arange(len(ptsArr)):
        poly.GetPointIds().SetId(ind,ind)
    polys.InsertNextCell(poly)

    polyPolyData = vtk.vtkPolyData()
    polyPolyData.SetPoints(ppts)
    polyPolyData.SetPolys(polys)
    if triangulate:
        triFilt = vtk.vtkTriangleFilter()
        triFilt.SetInputData(polyPolyData)
        triFilt.Update()
        return triFilt.GetOutput()
    else:
        return polyPolyData

def makeVolumePolygon(polygon,thickness=1,elevation=0.0,triangulate=False,cap=True):
    """
    Function to make a 3D/volume polygon from a 2D/planepolygon.

    Input:
        polygon - array list of vertices in clockwise order
        thickness - float, thinkess of the polyhedron
        elevation - float or triangulated vtkPolyData, elevation of the top of the polyhedron layer

    Output
        polyData - returns polydata of the volume polygon
    """

    # Check the elevation
    if type(elevation) is float:
        vtkPoly = makePolygon(polygon,elevation=elevation,triangulate=triangulate)
    elif 'vtkobject' in str(type(elevation)):
        if elevation.IsA('vtkPolyData'):
            eB = elevation.GetBounds()
            cutPoly = makePolygon(polygon,elevation=eB[-2]-10,triangulate=triangulate)
            cutVol = vtkTools.polydata.normFilter(vtkTools.polydata.extrudePolygon(cutPoly,eB[-1]-eB[-2]+20,[0,0,1],False))
            # Cut the elevation
            vtkPoly = vtkTools.extraction.clipDataSetWithPolygon(vtkTools.polydata.normFilter(elevation),cutVol,insideOut=True,extractBounds=True)

    # Extrude the topo surface downwards by the given thicknes.
    volPoly  = vtkTools.polydata.extrudePolygon(vtkPoly,thickness,[0,0,-1],False)
    return volPoly


def shape2polyhedron(shpFile,dbfHeader=None,thickness=1.0,elevation=0.0):
    """
    Function to convert GIS .shp and dbf to VTK vtu polyhedron.


    """

    import vtk, pysal, numpy as np, vtk.util.numpy_support as npsup

    # Read the input shp file
    shp = pysal.open(shpFile,'r')
    # Read the connected dbf file
    dbfFile = shpFile.split('.shp')[0] + '.dbf'
    dbf = pysal.open(dbfFile,'r')

    # Define the thickness and elevation of the holes
    if 'vtkobject' in str(type(elevation)):
        # Make the holes larger taller then the elevation
        eB = elevation.GetBounds()
        holeThink = eB[-1]-eB[-2] + 20.
        holeElev = eB[-1]+10.
    else:
        holeThink = thickness + 20
        holeElev = elevation + 10

    # Make a vtk object.
    # Initialize the data
    shpPHAppFilt = vtk.vtkAppendFilter()
    shpPHAppFilt.MergePointsOn()
    for nr, poly in enumerate(shp):
        mainPolygon = vtkTools.polydata.normFilter(makeVolumePolygon(np.array(poly.parts[0][:-1]),thickness,elevation,triangulate=True))
        # Deal with holes
        if poly.holes[0] == []:
            mainPHGrid = _makePolyhedronCell(mainPolygon,returnGrid=True)
        else:
            holesAppendFilt = vtk.vtkAppendPolyData()
            for hole in poly.holes:
                holePoly = vtkTools.polydata.normFilter(makeVolumePolygon(np.array(hole[:-1]),holeThink,holeElev,triangulate=True))
                holesAppendFilt.AddInputData(holePoly)
            # Make holes polygon
            holesAppendFilt.Update()
            holesPolygons = holesAppendFilt.GetOutput()
            # Cut the holes
            mainCutHolesPolygon = vtkTools.polydata.join2Polydata(mainPolygon,holesPolygons,threshold1='upper',threshold2='lower')
            # Add the cut polyhedron
            mainPHGrid = _makePolyhedronCell(mainCutHolesPolygon,returnGrid=True)
        shpPHAppFilt.AddInputData(mainPHGrid)
    shpPHAppFilt.Update()
    # Extract the vtu object.
    vtuPolyhedObj = shpPHAppFilt.GetOutput()
    # np.load(fileT)
    # Check if there is data to be associated with the cells.
    if dbfHeader:
        # Find all the headers given
        for head in dbfHeader:
            # Get the data
            nparr = np.array(dbf.by_col[head])
            # Sort the type of data
            if nparr.dtype.type is np.string_:
                vtkarr = vtk.vtkStringArray()
                for ele in ['NoDef' if l=='' else l for l in nparr ]: # Add 'NoDef' if the string is '' (empty)
                    vtkarr.InsertNextValue(ele)
            elif nparr.dtype.type is np.object:
                for nr, i in enumerate(nparr):
                    if i is None:
                        nparr[nr] = np.nan
                nparr = np.array(nparr,dtype=np.float)
                vtkarr = npsup.numpy_to_vtk(nparr,deep=1)
            else:
                nparr = np.array(nparr,dtype=np.float)
                vtkarr = npsup.numpy_to_vtk(nparr,deep=1)
            vtkarr.SetName(head)

            vtuPolyhedObj.GetCellData().AddArray(vtkarr)# Return the PH grid

    return vtuPolyhedObj



def shape2vtpFile(shpFile,dbfHeader):
    '''
    Function to convert GIS .shp and .dbt files to VTK .vtp file.

    NOTE: likely won't work, since cells in a vtp object can't have holes in them.
    '''

    import vtk, pysal, numpy as np, vtk.util.numpy_support as npsup

    # Read the input shp file
    shp = pysal.open(shpFile,'r')
    # Read the connected dbf file
    dbfFile = shpFile.split('.shp')[0] + '.dbf'
    dbf = pysal.open(dbfFile,'r')


    # Transfer all the info to the vtkPolyData
    appPolyData = vtk.vtkAppendPolyData()

    # cpid = 0
    for poly in shp:
        polyData = vtk.vtkPolyData()
        polyData.Allocate(10,10)
        polyPts = vtk.vtkPoints()
        polyPts.SetData(npsup.numpy_to_vtk(np.hstack( (poly.vertices,np.zeros((poly.len,1))) )[0:-1,:],deep=1))
        # polygons = vtk.vtkCellArray()
        # polygons.Allocate(len(shp),len(shp))

        # Initiate the vtk objects
        # p = vtk.vtkPolygon()
        # p.GetPointIds().SetNumberOfIds(poly.len-1)
        vtkidList = vtk.vtkIdList()
        for nr  in np.arange(0,poly.len -1,1):
            # cpid = cpid + nr
            vtkidList.InsertNextId(nr)
            # polyPts.InsertNextPoint(np.hstack((poly.vertices[nr],[0])))
            # p.GetPointIds().InsertNextId(cpid)

        # polygons.InsertNextCell(p)
        polyData.SetPoints(polyPts)
        polyData.InsertNextCell(vtk.VTK_POLYGON,vtkidList)
        appPolyData.AddInputData(polyData)
    # polyDat = vtk.vtkPolyData()

    # polyDat.SetPolys(polygons)
    appPolyData.Update()
    pData = appPolyData.GetOutput()
    # Set the data
    for head in dbfHeader:
        nparr = np.array(dbf.by_col[head])
        if nparr.dtype.type is np.string_:
            vtkarr = vtk.vtkStringArray()
            for ele in nparr:
                vtkarr.InsertNextValue(ele)
        elif nparr.dtype.type is np.object:
            for nr, i in enumerate(nparr):
                if i is None:
                    nparr[nr] = np.nan
            nparr = np.array(nparr,dtype=np.float)
            vtkarr = npsup.numpy_to_vtk(nparr,deep=1)
        else:
            nparr = np.array(nparr,dtype=np.float)
            vtkarr = npsup.numpy_to_vtk(nparr,deep=1)
        vtkarr.SetName(head)

        pData.GetCellData().AddArray(vtkarr)


    # Return the polydata object
    return pData


if __name__ == "__main__":
    # Main program

    shpFile = '/home/Gudni/Dropbox/Work/ISOR/Hengill/GeologyData/Jardfraedikort_tAtafla.shp'
    dbfHead = ['AREA','KENNI','atafla_LET','atafla_SKY','atafla_TYP']


    vtpObj = shape2vtpFile(shpFile,dbfHead)
    vtkTools.io.writeVTPFile('Hengill_geologicMap.vtp',vtpObj)
