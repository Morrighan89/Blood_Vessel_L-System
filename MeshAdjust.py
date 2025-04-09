#!/usr/bin/env python


#Elena Faggiano, 11.05.2012


import sys

import vtk
from vmtk import pypes
from vmtk import vmtkscripts


#MeshAdjust = 'meshAdjust'


class meshAdjust(pypes.pypeScript):


    def __init__(self):
        pypes.pypeScript.__init__(self)
        
        self.Surface = None
        self.CellEntityIdsArrayName = 'CellEntityIds'
        self.TargetEdgeLength = 1.0
        self.ConstraintFactor = 1.0
        self.NumberOfRings = 8        
        # Output objects
        self.Mesh = None        
        
        self.SetScriptName('meshAdjust')
        self.SetScriptDoc('delineate a region on the surface and eliminate it (remeshing the hole)')
        self.SetInputMembers([['Surface','i','vtkPolyData',1,'','the input surface','vmtksurfacereader'],['CellEntityIdsArrayName','entityidsarray','str',1],['ConstraintFactor','constraint','float',1,'','amount of influence of the shape of the surface near the boundary on the shape of the cap ("smooth" method only)'],['NumberOfRings','rings','int',1,'(0,)','number of rings composing the cap ("smooth" method only)'],['TargetEdgeLength','edgelength','float',1,'(0.0,)']])
        self.SetOutputMembers([['Mesh','o','vtkPolyData',1,'','the output surface','vmtksurfacewriter']])

    def Execute(self): 
        if (self.Surface == None):
            self.PrintError('Error: no Surface.') 

        surface=vtk.vtkPolyData()
        surface=self.Surface   
        
        print("please, delineate the area to be eliminated (ctrl + left click)")
        renderWindowInteractor = vtk.vtkRenderWindowInteractor()
        renderer = vtk.vtkRenderer()
        renderWindow = vtk.vtkRenderWindow()
        renderWindow.SetInteractor(renderWindowInteractor)
        renderWindow.AddRenderer(renderer)
        
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(surface)
        mapper.ScalarVisibilityOff()
        
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetInterpolationToFlat()
        actor.GetProperty().SetRepresentationToSurface()
        actor.GetProperty().SetEdgeVisibility(1)
        actor.GetProperty().SetEdgeColor(0 , 0, 0)
        
        
        renderer.AddActor(actor)
        renderer.SetBackground(1,1,1)

        interactorStyle = vtk.vtkInteractorStyleTrackballCamera()
        renderWindowInteractor.SetInteractorStyle(interactorStyle)

        contourWidget=vtk.vtkContourWidget()
        contourWidget.SetInteractor(renderWindowInteractor)
        rep = vtk.vtkOrientedGlyphContourRepresentation()
        rep = contourWidget.GetRepresentation()
        rep.GetLinesProperty().SetColor(0, 1, 0)
        rep.GetLinesProperty().SetLineWidth(10.0)
        pointPlacer =vtk.vtkPolygonalSurfacePointPlacer()
        pointPlacer.AddProp(actor)
        pointPlacer.GetPolys().AddItem(surface)
        rep.SetPointPlacer(pointPlacer)
        polycontour = vtk.vtkPolyData()
        polycontour = rep.GetContourRepresentationAsPolyData()
        interpolator = vtk.vtkPolygonalSurfaceContourLineInterpolator()
        interpolator.GetPolys().AddItem(surface)
        rep.SetLineInterpolator(interpolator)
        renderWindow.Render()
        renderWindowInteractor.Initialize()
        contourWidget.EnabledOn()
        renderWindowInteractor.Start()

        select = vtk.vtkSelectPolyData()
        select.SetInputData(surface)
        select.SetLoop(polycontour.GetPoints())
        select.SetSelectionModeToSmallestRegion()
        select.GenerateSelectionScalarsOff()
        select.GenerateUnselectedOutputOn()
        select.Update()

        Unselected = vtk.vtkPolyData()
        Unselected = select.GetUnselectedOutput()

        capper = vmtkscripts.vmtkSurfaceCapper()
        capper.Surface = Unselected
        capper.Method = 'smooth'
        capper.ConstraintFactor = self.ConstraintFactor
        capper.NumberOfRings = self.NumberOfRings
        #capper.CellEntityIdsArrayName=self.CellEntityIdsArrayName
        capper.Execute()

        #capper.Surface.GetCellData.GetArray(self.CellEntityIdsArrayName).FillComponent(0,1)

        remeshing = vmtkscripts.vmtkSurfaceRemeshing()
        remeshing.Surface = capper.Surface
        #remeshing.CellEntityIdsArrayName = self.CellEntityIdsArrayName
        remeshing.ElementSizeMode = 'edgelength'
        remeshing.TargetEdgeLength = self.TargetEdgeLength
        remeshing.Execute()

        #cleaner = vtk.vtkCleanPolyData()
        #cleaner.SetInput(remeshing.Surface)
        #cleaner.Update()

        triangleFilter = vtk.vtkTriangleFilter()
        #triangleFilter.SetInput(cleaner.GetOutput())
        triangleFilter.SetInputData(remeshing.Surface)
        triangleFilter.Update()

        self.Mesh = triangleFilter.GetOutput()

        renderWindowInteractor = vtk.vtkRenderWindowInteractor()
        renderer = vtk.vtkRenderer()
        renderWindow = vtk.vtkRenderWindow()
        renderWindow.SetInteractor(renderWindowInteractor)
        renderWindow.AddRenderer(renderer)

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(triangleFilter.GetOutput())
        mapper.ScalarVisibilityOff()

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetInterpolationToFlat()
        actor.GetProperty().SetRepresentationToSurface()
        actor.GetProperty().SetEdgeVisibility(1)
        actor.GetProperty().SetEdgeColor(0, 0, 0)

        renderer.AddActor(actor)
        renderer.SetBackground(1, 1, 1)

        interactorStyle = vtk.vtkInteractorStyleTrackballCamera()
        renderWindowInteractor.SetInteractorStyle(interactorStyle)
        renderWindow.Render()
        renderWindowInteractor.Initialize()
        renderWindowInteractor.Start()

if __name__=='__main__':

    main = pypes.pypeMain()
    main.Arguments = sys.argv
    print(main.Arguments)
    main.Execute()
