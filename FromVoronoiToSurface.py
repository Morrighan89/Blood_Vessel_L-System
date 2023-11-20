#!/usr/bin/env python


from __future__ import print_function, absolute_import # NEED TO STAY AS TOP IMPORT
import vtk
import sys
import math
from vmtk import vtkvmtk

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
   polyBallImageSize = [1200,1200,1200]   #size of the image for the evaluation of the polyball function

   

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
   WritePolyData(envelope,surfaceFilename)

def main():
   from_voronoi_to_surf(fileIn='voronoiDiagram6.vtp',fileOut='reconstructedmodel6.vtp')

if __name__=='__main__':
    main()

