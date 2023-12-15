#!/usr/bin/env python


from __future__ import print_function, absolute_import # NEED TO STAY AS TOP IMPORT
import vtk
import sys
import math
from vmtk import vtkvmtk
import vmtk.vmtksurfaceremeshing as remeshing

def ReadPolyData(filename):
   reader = vtk.vtkXMLPolyDataReader()
   reader.SetFileName(filename)
   reader.Update()
   return reader.GetOutput()

def WritePolyData(input,filename):
   writer = vtk.vtkXMLPolyDataWriter()
   writer.SetFileName(filename)
   writer.SetInputData(input)
   writer.Write()
  

def from_voronoi_to_surf(fileIn='voronoiDiagram',fileOut='reconstructedmodel.vtp'):
   #SOME COMMON VMTK DATA ARRAY NAMES
   radiusArrayName = 'MaximumInscribedSphereRadius'
   parallelTransportNormalsArrayName = 'ParallelTransportNormals'

   #OPTIONS TO SET
   interpolationHalfSize = 3
   polyBallImageSize = [300,300,300]   #size of the image for the evaluation of the polyball function

   

   VoronoiFilename       = fileIn

   surfaceFilename                 = fileOut

   Voronoi = ReadPolyData(VoronoiFilename)



   print('Reconstructing Surface from Voronoi Diagram')
   modeller = vtkvmtk.vtkvmtkPolyBallModeller()
   modeller.SetInputData(Voronoi)
   modeller.SetRadiusArrayName(radiusArrayName)
   modeller.UsePolyBallLineOff()
   modeller.SetSampleDimensions(polyBallImageSize)
   modeller.Update()
   print('MarchingCubes')
   marchingCube = vtk.vtkMarchingCubes()
   marchingCube.SetInputData(modeller.GetOutput())
   marchingCube.SetValue(0,0.0)
   marchingCube.Update()  
   envelope = marchingCube.GetOutput()
   print('smoothing')
   smoothingFilter= vtk.vtkSmoothPolyDataFilter()
   smoothingFilter.SetInputDataObject(envelope)
   smoothingFilter.SetNumberOfIterations(20)
   smoothingFilter.Update()
   envelope_smooth=smoothingFilter.GetOutput()
   print("remeshing")
   remesher = remeshing.vmtkSurfaceRemeshing()
   remesher.Surface = envelope_smooth
   remesher.TargetAreaFactor = 1.0
   remesher.NumberOfIterations = 1
   remesher.Execute()
   WritePolyData(remesher.Surface,surfaceFilename)

def test_gpt_suggestion():


   # Read the vtk surface file
   reader = vtk.vtkXMLPolyDataReader()
   reader.SetFileName("subcube_reconstructedmodel28.vtp")
   reader.Update()

   # Create a plane
   plane = vtk.vtkPlane()
   plane.SetOrigin(7.43,8.2,6.56)  # Set the origin of the plane
   plane.SetNormal(0.7,0.012,0.68)  # Set the normal vector of the plane

   # Create a cutter
   cutter = vtk.vtkCutter()
   cutter.SetInputConnection(reader.GetOutputPort())
   cutter.SetCutFunction(plane)
   
   # Create a mapper for the cut surface
   cutMapper = vtk.vtkPolyDataMapper()
   cutMapper.SetInputConnection(cutter.GetOutputPort())
   
   # Create an actor for the cut surface
   cutActor = vtk.vtkActor()
   cutActor.SetMapper(cutMapper)
   
   # Create a boolean operation filter
   booleanOperation = vtk.vtkBooleanOperationPolyDataFilter()
   booleanOperation.SetOperationToIntersection()
   booleanOperation.SetInputData(0, reader.GetOutput())
   booleanOperation.SetInputConnection(1, cutter.GetOutputPort())
   booleanOperation.Update()
   
   # Create a mapper for the inner part
   innerMapper = vtk.vtkPolyDataMapper()
   innerMapper.SetInputData(booleanOperation.GetOutput())
   
   # Create an actor for the inner part
   innerActor = vtk.vtkActor()
   innerActor.SetMapper(innerMapper)
   innerActor.GetProperty().SetColor(0, 1, 0)  # Set color to green for better visualization
   
   # Create a renderer
   renderer = vtk.vtkRenderer()
   renderer.AddActor(cutActor)
   renderer.AddActor(innerActor)
   
   # Create a render window
   renderWindow = vtk.vtkRenderWindow()
   renderWindow.AddRenderer(renderer)
   
   # Create a render window interactor
   renderWindowInteractor = vtk.vtkRenderWindowInteractor()
   renderWindowInteractor.SetRenderWindow(renderWindow)
   
   # Set up the camera
   renderer.GetActiveCamera().Azimuth(30)
   renderer.GetActiveCamera().Elevation(30)
   
   # Set up the background color
   renderer.SetBackground(1, 1, 1)
   
   # Render the scene
   renderWindow.Render()
   
   # Start the interaction
   renderWindowInteractor.Start()

def main():
   #from_voronoi_to_surf(fileIn='vor_diag_sub_cube_28.vtp',fileOut='subcube_reconstructedmodel28.vtp')
   test_gpt_suggestion()

if __name__=='__main__':
    main()

