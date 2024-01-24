
from pyvista.utilities.helpers import vtk_points
import vtk
import pyvista as pv
import numpy as np
import math
import copy
import os
from tqdm import tqdm
#from vtkmodules.vtkIOXML import vtkXMLPolyDataWriter
import vtkmodules.numpy_interface as np2vtk
from vtkmodules.vtkCommonDataModel import (
    vtkCellArray,
    vtkPolyData,
    vtkPolyLine,
    vtkGenericCell
)
import gmsh
from VesselInterpreter import (
    #ReadPolyData,
    ExtractLine,
    radiusArrayName,
    parallelTransportNormalsArrayName
)
from itertools import chain
from vtk_read_write import (
    read_poly_data
)

def create_gmsh_mesh_carving(clFileName,folder="temp",ofile="gmsh.msh",cx=0.0,cy=-18.0,cz=-18.0,dx=20.0,dy=36.0,dz=36.0,scalingfactor=1.0,brep_module=100):
    baseCl=read_poly_data(os.path.join(folder,clFileName))
    temp = os.path.splitext(ofile)
    var = (os.path.basename(temp[0]), temp[1])
    name_ofile=var[0]

    numberOfLines=baseCl.GetNumberOfCells()
    numberOfPoints=baseCl.GetNumberOfPoints()
    print(f" Vessel segments: {numberOfLines}, Nodes: {numberOfPoints}")
    options=['t1.geo','-tol', '1.e-14','-setnumber', 'Geometry.OCCFixDegenerated', '1','-setnumber', 'Geometry.OCCFixSmallEdges', '1','-setnumber', 'Geometry.OCCFixSmallFaces', '1']
    gmsh.initialize(argv=options)
    gmsh.model.add("DFG 3D")
   # print(gmsh.model.get  r('Geometry.Tolerance'))
    newdims=[]
    
    
    #cx,cy,cz=0,-18,-18
    #dx,dy,dz=20,36,36
    
    #p1=gmsh.model.occ.addPoint(cx,cy,cz)
    #p2=gmsh.model.occ.addPoint(cx,cy+dy,cz)
    #p3=gmsh.model.occ.addPoint(cx,cy+dy,cz+dz)
    #p4=gmsh.model.occ.addPoint(cx,cy,cz+dz)
    #p5=gmsh.model.occ.addPoint(cx+dx,cy,cz)
    #p6=gmsh.model.occ.addPoint(cx+dx,cy+dy,cz)
    #p7=gmsh.model.occ.addPoint(cx+dx,cy+dy,cz+dz)
    #p8=gmsh.model.occ.addPoint(cx+dx,cy,cz+dz)
    #gmsh.model.occ.synchronize()
    #gmsh.write(f'{folder}\model_{name_ofile}_test.brep')
#
    #gmsh.model.occ.addLine(1,2, 9)
    #gmsh.model.occ.addLine(p2,p3, 10)
    #gmsh.model.occ.addLine(p3,p4, 11)
    #gmsh.model.occ.addLine(p4,p1, 12)
#
    #gmsh.model.occ.addCurveLoop([9,10,11,12], 13)
    #gmsh.model.occ.addPlaneSurface([13], 14)
#
    #gmsh.model.occ.addLine(p2,p6, 15)
    #gmsh.model.occ.addLine(p6,p5, 16)
    #gmsh.model.occ.addLine(p5,p1, 17)
#
    #gmsh.model.occ.addCurveLoop([9,15,16,17], 18)
    #gmsh.model.occ.addPlaneSurface([18], 19)
#
    #gmsh.model.occ.addLine(p3,p7, 20)
    #gmsh.model.occ.addLine(p7,p6, 21)
#
    #gmsh.model.occ.addCurveLoop([10,20,21,-15], 22)
    #gmsh.model.occ.addPlaneSurface([22], 23)
#
    #gmsh.model.occ.addLine(p4,p8, 24)
    #gmsh.model.occ.addLine(p8,p7, 25)
#
    #gmsh.model.occ.addCurveLoop([11,24,25,-20], 26)
    #gmsh.model.occ.addPlaneSurface([26], 27)
#
    #gmsh.model.occ.addLine(p5,p8, 28)
#
    #gmsh.model.occ.addCurveLoop([12,-17,28,-24], 29)
    #gmsh.model.occ.addPlaneSurface([29], 30)
#
    #outlet1=gmsh.model.occ.addCircle(20,0,15,1,1,31,zAxis=[1,0,0])
    #outlet2=gmsh.model.occ.addCircle(20,0,-15,1,1,32,zAxis=[1,0,0])
    #inlet=gmsh.model.occ.addCircle(20,0,0,1,1,33,zAxis=[1,0,0])
#
    #gmsh.model.occ.addCurveLoop([-16,-21,-25,-28], 33)
    #gmsh.model.occ.addCurveLoop([31,32,33], 34)
#
    #gmsh.model.occ.addPlaneSurface([33,34], 35)
   #
    #gmsh.write(f'{folder}\model_{name_ofile}_test.brep')

    box=gmsh.model.occ.addBox(cx,cy,cz, dx,dy,dz)
    volumes=gmsh.model.occ.getEntities(dim=3)
    gmsh.model.occ.remove(volumes)
    volumes=gmsh.model.occ.getEntities(dim=3)
    outlet1=gmsh.model.occ.addDisk(20,0,15,1,1,zAxis=[1,0,0])
    outlet2=gmsh.model.occ.addDisk(20,0,-15,1,1,zAxis=[1,0,0])
    inlet=gmsh.model.occ.addDisk(20,0,0,1,1,zAxis=[1,0,0])
    cut=gmsh.model.occ.cut([(2,2)],[(2,outlet1),(2,outlet2),(2,inlet)],removeTool= False,removeObject= True)
    gmsh.model.occ.synchronize()
    surfaces=gmsh.model.occ.getEntities(dim=2)
    shell=gmsh.model.occ.addSurfaceLoop([8,2,9,7,3,5,1,4,6])

    #shell=gmsh.model.occ.add_surface_loop(surfaces)
    gmsh.write(f'{folder}\model_{name_ofile}_test.brep')
    gmsh.model.occ.addVolume([shell])
    #gmsh.model.occ.remove([(2,2)])
    gmsh.model.occ.synchronize()

    gmsh.write(f'{folder}\model_{name_ofile}_test.brep')


    volumes=gmsh.model.occ.getEntities(dim=3)
    count = 0
    endpoints=[]
    for i in tqdm(range(0,numberOfLines)):#numberOfLineschain(range(0,136),range(138,numberOfLines))
        line=ExtractLine(i,baseCl)
        linePtsNumber=line.GetNumberOfPoints()
        newtrunks=[]
        volumes=gmsh.model.occ.getEntities(dim=3)
        try:
        #volumes=gmsh.model.occ.getEntities(dim=3) 
            #for j in range(0,linePtsNumber-1): #linePtsNumber-1
            endpoint=line.GetPoint(6)
            if endpoint not in endpoints:
                for j in range(0,7): #linePtsNumber-1
                    startPoint=line.GetPoint(j)
                    startingRadius=line.GetPointData().GetArray(radiusArrayName).GetTuple1(j)
                    endingRadius=line.GetPointData().GetArray(radiusArrayName).GetTuple1(j+1)
                    direction=line.GetPointData().GetArray(parallelTransportNormalsArrayName).GetTuple3(j)
                    if startingRadius==endingRadius:
                        cone = gmsh.model.occ.addCylinder(startPoint[0],startPoint[1],startPoint[2],direction[0], direction[1],direction[2], startingRadius)
                        newtrunks.append((3,cone))
                    else:
                        cone = gmsh.model.occ.addCone(startPoint[0],startPoint[1],startPoint[2],direction[0], direction[1],direction[2], startingRadius,    endingRadius)
                        newtrunks.append((3,cone))
                    if j<linePtsNumber-2:
                        ball = gmsh.model.occ.addSphere(startPoint[0]+direction[0],startPoint[1]+direction[1],startPoint[2]+direction[2],endingRadius)
                        newtrunks.append((3,ball))
                    #if j==0:
                    #    ball = gmsh.model.occ.addSphere(startPoint[0],startPoint[1],startPoint[2],startingRadius*1.01)
                    #    newtrunks.append((3,ball))
                #volumes=gmsh.model.occ.getEntities(dim=3)

                #pippo=gmsh.model.occ.fuse([newtrunks[0]], newtrunks)
                gmsh.model.occ.cut(volumes,newtrunks,removeTool= True)
                gmsh.model.occ.synchronize()
                #gmsh.model.occ.cut(volumes,pippo[0],removeObject= True, removeTool= True)
                endpoints.append(endpoint)
            else:
                print(f'skipping line {i}')
                
        except:
            print(f'Problem with polyline {i}')
        volumes=gmsh.model.occ.getEntities(dim=3)
        gmsh.model.occ.synchronize()
        count=count+1
        if ((count % brep_module) == 0 or i==1):
            gmsh.write(f'{folder}\model_{name_ofile}_{i}.brep')
    #cx,cy,cz=0.8734448800133734, -2.081709684995058, 1.6177133644915518
    #cx,cy,cz=0.12864125650976543,-2.0758902928179417,-1.8674509158982868
    #dx,dy,dz=2.580985600142573,3.589531105579775,3.400352085492524
    #dx,dy,dz=1, 1, 1
    #box=gmsh.model.occ.addBox(cx,cy,cz, dx,dy,dz)
    #gmsh.model.occ.cut([(3,box)],volumes, removeTool= True) 
    gmsh.model.occ.synchronize()
    gmsh.write(f'{folder}\model_{name_ofile}_fine.brep')
    print(f'Begin labeling')
    #gmsh.model.occ.remove([(2,10)])
    gmsh.model.occ.synchronize()
    surfaces = gmsh.model.occ.getEntities(dim=2)
    root_marker,   bottom_side_marker,  lateral_sides_marker = 203,  209,  213
    walls = []
    leaves = []
    laterals=[]
    #comtop =gmsh.model.occ.getCenterOfMass(2, 2)
    for surface in tqdm(surfaces):
        com = gmsh.model.occ.getCenterOfMass(surface[0], surface[1])
        surface_type= gmsh.model.getType(surface[0], surface[1])
        if surface_type == 'Plane':
            if np.allclose(com, [cx+dx,0,0]):
                top_side_marker=gmsh.model.addPhysicalGroup(surface[0], [surface[1]])
                gmsh.model.setPhysicalName(surface[0], top_side_marker, "topSide")
            elif np.allclose(com, [cx, 0, 0]):
                gmsh.model.addPhysicalGroup(surface[0], [surface[1]], bottom_side_marker)
                gmsh.model.setPhysicalName(surface[0], bottom_side_marker, "bottomSide")
            elif np.allclose(com, [0,0,0]):
                gmsh.model.addPhysicalGroup(surface[0], [surface[1]], root_marker)
                gmsh.model.setPhysicalName(surface[0], root_marker, "root")
            elif np.isclose(com[2], cz) or np.isclose(com[1], cy+dy) or np.isclose(com[2], cz+dz) or np.isclose(com[1],cy):
                laterals.append(surface[1])
            else:
                leaves.append(surface[1])
        else:
            walls.append(surface[1])
    gmsh.model.addPhysicalGroup(2, laterals, lateral_sides_marker)
    gmsh.model.setPhysicalName(2, lateral_sides_marker, "lateral")
    leaves_marker=gmsh.model.addPhysicalGroup(2, leaves)
    gmsh.model.setPhysicalName(2, leaves_marker, "leaves")
    wall_marker=gmsh.model.addPhysicalGroup(2, walls)
    gmsh.model.setPhysicalName(2, wall_marker, "Walls")
    inlet_marker=gmsh.model.addPhysicalGroup(2, [7])
    gmsh.model.setPhysicalName(2, inlet_marker, "inlet")
    outlet_marker=gmsh.model.addPhysicalGroup(2, [4,5])
    gmsh.model.setPhysicalName(2, outlet_marker, "outlet")


    volumes = gmsh.model.occ.getEntities(dim=3)
    gmsh.model.addPhysicalGroup(3, [volumes[0][1]], 600)
    gmsh.model.setPhysicalName(3, 600, "volume")
    gmsh.model.occ.synchronize()
    gmsh.write(f'{folder}\model_{name_ofile}.brep')
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 15)
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.0005)
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.75)
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(3)
    gmsh.write(os.path.join(folder,f'model_{os.path.splitext(name_ofile)[0]}.vtk'))
    print(inlet_marker,outlet_marker,wall_marker,leaves_marker)

def create_gmsh_mesh_adding(clFileName,folder="temp",ofile="gmsh.msh",scalingfactor=1):
    baseCl=read_poly_data(clFileName)
    temp = os.path.splitext(ofile)
    var = (os.path.basename(temp[0]), temp[1])
    name_ofile=var[0]
    numberOfLines=baseCl.GetNumberOfCells()
    numberOfPoints=baseCl.GetNumberOfPoints()
    print(f" Vessel segments: {numberOfLines}, Nodes: {numberOfPoints}")
    options=['t1.geo','-tol', '1.e-4','-setnumber', 'Geometry.OCCFixDegenerated', '1','-setnumber', 'Geometry.OCCFixSmallEdges', '1','-setnumber', 'Geometry.OCCFixSmallFaces', '1']
    gmsh.initialize(argv=options)
    
    gmsh.model.add("DFG 3D")
    
    newdims=[]
    line=ExtractLine(0,baseCl)
    startPoint=line.GetPoint(0)
    startingRadius=scalingfactor*line.GetPointData().GetArray(radiusArrayName).GetTuple1(0)
    endingRadius=scalingfactor*line.GetPointData().GetArray(radiusArrayName).GetTuple1(1)    
    direction=line.GetPointData().GetArray(parallelTransportNormalsArrayName).GetTuple3(0)
    if startingRadius==endingRadius:
        trunk = gmsh.model.occ.addCylinder(startPoint[0],startPoint[1],startPoint[2],direction[0], direction[1],direction[2], startingRadius)
    else:
        trunk = gmsh.model.occ.addCone(startPoint[0],startPoint[1],startPoint[2],direction[0], direction[1],direction[2], startingRadius, endingRadius)
    ball = gmsh.model.occ.addSphere(startPoint[0]+direction[0],startPoint[1]+direction[1],startPoint[2]+direction[2],endingRadius)
    badboys=[]
    
    for i in range(0, numberOfLines-1):#chain(range(381,401),range(401,numberOfLines)):#range(0, numberOfLines):#chain(range(0,136),range(138,161),range(181,numberOfLines)):#chain(range(0,136),range(138,numberOfLines))
        line=ExtractLine(i,baseCl)
        linePtsNumber=line.GetNumberOfPoints()
        newtrunks=[]
        cones=[]
        balls=[]
        print(i)

        volumes=gmsh.model.occ.getEntities(dim=3) 
        oldVolumes=gmsh.model.occ.getEntities(dim=3) 
        #try:
        #for j in range(0,linePtsNumber-1): #linePtsNumber-1
        for j in range(0,9): #linePtsNumber-1
            startPoint=line.GetPoint(j)
            startingRadius=scalingfactor*line.GetPointData().GetArray(radiusArrayName).GetTuple1(j)
            endingRadius=scalingfactor*line.GetPointData().GetArray(radiusArrayName).GetTuple1(j+1)
            direction=line.GetPointData().GetArray(parallelTransportNormalsArrayName).GetTuple3(j)
            if startingRadius==endingRadius:
                cone = gmsh.model.occ.addCylinder(startPoint[0],startPoint[1],startPoint[2],direction[0], direction[1],direction[2], startingRadius)
                cones.append((3,cone))
            else:
                cone = gmsh.model.occ.addCone(startPoint[0],startPoint[1],startPoint[2],direction[0], direction[1],direction[2], startingRadius,   endingRadius)
                cones.append((3,cone))
            if j<linePtsNumber-1:
                ball = gmsh.model.occ.addSphere(startPoint[0]+direction[0],startPoint[1]+direction[1],startPoint[2]+direction[2],endingRadius*1.01)
                balls.append((3,ball))
                test=gmsh.model.occ.fuse(cones,balls,removeObject= True, removeTool= True)
                
            #if j==0:
            #    ball = gmsh.model.occ.addSphere(startPoint[0],startPoint[1],startPoint[2],startingRadius*1.01)
            #    newtrunks.append((3,ball))
        volumes=gmsh.model.occ.getEntities(dim=3)
        [volumes.remove(item) for item in badboys if item in volumes]
        try:
            gmsh.model.occ.fuse(test[0],volumes,removeObject= True, removeTool= True)
        except:
            print("gmsh.model.occ.fuse(volumes,test[0],removeObject= True, removeTool= True")
            #gmsh.model.occ.fuse(test[0],volumes,removeObject= True, removeTool= True)
            badboys.append(test[0][0])
        volumes=gmsh.model.occ.getEntities(dim=3)
        #gmsh.model.occ.fuse([(3,trunk)], newtrunks,removeObject= False, removeTool= True)
        #except:
        #volumes=oldVolumes
        #print(f'Problem with CL {i}')
        #print(volumes,i)
        gmsh.model.occ.synchronize()

        if ((i % 2) == 0  or i==1):
            gmsh.write(f'{folder}\{name_ofile}_{i}.step')
    #cx,cy,cz=0.2816580241399551, -1.9043948378959654, -1.3874804023695808
    #dx,dy,dz=3.195877166439165,2.9768761114761113,2.5549797529430096
    #cx,cy,cz=0.22723123693223046,-1.117517697735693,-1.5175201865734897  # villi14
    #dx,dy,dz=2.0870246230892904,2.8587875728566874,3.0326460420042394    # villi14
    #cx,cy,cz=0.11962770121213051,-0.4990486324439569,-1.4965857674994025# villi15
    #dx,dy,dz=2.2865367795282356,2.5604400762919113,2.7007688261857514# villi15
    #cx,cy,cz=0,-1.2402858552276663,-1.3034999407295926# villi16
    #dx,dy,dz=2.2,2.4,3# villi16
    #cx,cy,cz=0.2, -1.75, -1.75 # villi19
    #dx,dy,dz=2.2, 3.5, 3.5 # villi19
    cx,cy,cz=8.8,2,-5.5
    dx,dy,dz=1,1,1
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 20)
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.0005)
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.18)
    gmsh.model.mesh.generate(2)
    gmsh.write(f'{folder}\model_{name_ofile}_2.vtk')
    box=gmsh.model.occ.addBox(cx,cy,cz, dx,dy,dz)
    
    gmsh.model.occ.cut([(3,box)],volumes, removeTool= True) 
    gmsh.model.occ.synchronize()
    gmsh.write(f'{folder}\model_{name_ofile}.step')
    print(f'Begin labeling')
    surfaces = gmsh.model.occ.getEntities(dim=2)
    root_marker, leaves_marker, wall_marker, bottom_side_marker, top_side_marker, lateral_sides_marker = 203, 205, 207, 209, 211, 213
    walls = []
    leaves = []
    laterals=[]
    volumes = gmsh.model.occ.getEntities(dim=3)
    gmsh.model.addPhysicalGroup(3, [2], 300)
    gmsh.model.setPhysicalName(3, 300, "volume")
    print(volumes)
    comtop =gmsh.model.occ.getCenterOfMass(2, 1)
    for surface in tqdm(surfaces):
        com = gmsh.model.occ.getCenterOfMass(surface[0], surface[1])
        surface_type= gmsh.model.getType(surface[0], surface[1])
        if surface_type == 'Plane':
            if np.allclose(com, comtop):
                gmsh.model.addPhysicalGroup(surface[0], [surface[1]], top_side_marker)
                gmsh.model.setPhysicalName(surface[0], top_side_marker, "topSide")
            elif np.allclose(com, [cx+dx, cy+dy/2, cz+dz/2]):
                gmsh.model.addPhysicalGroup(surface[0], [surface[1]], bottom_side_marker)
                gmsh.model.setPhysicalName(surface[0], bottom_side_marker, "bottomSide")
            elif np.allclose(com, [0,0,0]):
                gmsh.model.addPhysicalGroup(surface[0], [surface[1]], root_marker)
                gmsh.model.setPhysicalName(surface[0], root_marker, "root")
            elif np.isclose(com[2], cz) or np.isclose(com[1], cy+dy) or np.isclose(com[2], cz+dz) or np.isclose(com[1],cy):
                laterals.append(surface[1])
            else:
                leaves.append(surface[1])
        else:
            walls.append(surface[1])
    gmsh.model.addPhysicalGroup(2, laterals, lateral_sides_marker)
    gmsh.model.setPhysicalName(2, lateral_sides_marker, "lateral")
    gmsh.model.addPhysicalGroup(2, leaves, leaves_marker)
    gmsh.model.setPhysicalName(2, leaves_marker, "leaves")
    gmsh.model.addPhysicalGroup(2, walls, wall_marker)
    gmsh.model.setPhysicalName(2, wall_marker, "Walls")
    gmsh.model.occ.synchronize()
    gmsh.write(f'{folder}\model_{name_ofile}_2.step')
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 20)
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.0005)
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.18)
    gmsh.model.mesh.generate(3)
    gmsh.write(f'{folder}\model_{name_ofile}_3.msh')
    gmsh.write(f'{folder}\model_{name_ofile}_3.vtk')
    gmsh.finalize()
    
def boolean_example():
    gmsh.initialize()

    gmsh.model.add("boolean")

    # from http://en.wikipedia.org/wiki/Constructive_solid_geometry

    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.4)
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.4)

    R = 1.4
    Rs = R * .7
    Rt = R * 1.25

    gmsh.model.occ.addBox(-R, -R, -R, 2 * R, 2 * R, 2 * R, 1)
    gmsh.model.occ.addSphere(0, 0, 0, Rt, 2)
    gmsh.model.occ.intersect([(3, 1)], [(3, 2)], 3)
    gmsh.model.occ.addCylinder(-2 * R, 0, 0, 4 * R, 0, 0, Rs, 4)
    gmsh.model.occ.addCylinder(0, -2 * R, 0, 0, 4 * R, 0, Rs, 5)
    gmsh.model.occ.addCylinder(0, 0, -2 * R, 0, 0, 4 * R, Rs, 6)
    gmsh.model.occ.fuse([(3, 4), (3, 5)], [(3, 6)], 7)
    gmsh.model.occ.cut([(3, 3)], [(3, 7)], 8)

    gmsh.model.occ.synchronize()

    gmsh.model.mesh.generate(3)
    #gmsh.model.mesh.refine()
    #gmsh.model.mesh.setOrder(2)
    #gmsh.model.mesh.partition(4)

    gmsh.write("boolean.msh")
    gmsh.write("boolean.brep")
    gmsh.finalize()

def example_spline_extrude():
    gmsh.initialize()

    gmsh.model.add("extrude spline")
    nturns = 2 # tested ok up to 100

    npts = 100 * nturns
    r = 1.
    rd = 0.1
    h = 1. * nturns

    for i in range(npts):
      theta = i * 2. * math.pi * nturns / npts
      gmsh.model.occ.addPoint(r * math.cos(theta), r * math.sin(theta),
                              i * h / npts, i+1)

    gmsh.model.occ.addSpline(range(1, npts), 1)
    gmsh.model.occ.addWire([1], 1)

    gmsh.model.occ.addDisk(1,0,0, rd, rd, 1)

    gmsh.model.occ.addRectangle(1+2*rd,-rd,0, 2*rd,2*rd, 2, rd/5)
    gmsh.model.occ.rotate([(2, 1), (2, 2)], 0, 0, 0, 1, 0, 0, math.pi/2)

    #gmsh.model.occ.addPipe([(2, 1), (2, 2)], 1, 'DiscreteTrihedron')
    gmsh.model.occ.addPipe([(2, 1), (2, 2)], 1, 'Frenet')

    gmsh.model.occ.remove([(2, 1), (2, 2), (1, 1)])

    gmsh.model.occ.synchronize()

    gmsh.option.setNumber('Mesh.MeshSizeMin', 0.1)
    gmsh.option.setNumber('Mesh.MeshSizeMax', 0.1)
    gmsh.option.setNumber('Geometry.NumSubEdges', npts) # nicer display of curves
    gmsh.model.mesh.generate(3)
    gmsh.write("spline.msh")
    gmsh.write("spline.brep")
    gmsh.finalize()

def main():# used for debug purpose of function in this file
    #create_gmsh_mesh_adding("vtkvessel35_1st_iter.vtp",folder="step",ofile='vtkvessel35_1st_iter.vtp.msh',scalingfactor=1)
    create_gmsh_mesh_carving("vtkVilli38_1st_iter.vtp",folder="step",ofile='vtkvessel38_1st_iter.msh',scalingfactor=1)
    #create_gmsh_mesh_carving("sub_cube.vtp")
    #create_gmsh_mesh_carving("vtkVilli15.vtp",'villi15.msh')
    #boolean_example()
    #example_spline_extrude()
if __name__=='__main__':
    
    main()