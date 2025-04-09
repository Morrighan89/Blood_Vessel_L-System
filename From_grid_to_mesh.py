import vtk
import numpy as np
from vtkmodules.vtkCommonDataModel import vtkStaticPointLocator
from vtkmodules.vtkFiltersPoints import (
    vtkLinearKernel,
    vtkPointInterpolator
)
from vtk.util import numpy_support
import os
from vtk_read_write import (
    read_unstructured_grid,
    read_structured_grid,
    read_structuredXML_mesh,
    write_structuredXML_mesh,
    write_poly_data,
    write_unstructuredXML_mesh
)

#def load_existing_mesh(file_path):
#    reader = vtk.vtkUnstructuredGridReader()
#    reader.SetFileName(file_path)
#    reader.Update()
#    return reader.GetOutput()
#
#def read_structured_mesh(file_path):
#    reader = vtk.vtkStructuredGridReader()
#    reader.SetFileName(file_path)
#    reader.Update()
#    return reader.GetOutput()
#
#def read_structuredXML_mesh(file_path):
#    reader = vtk.vtkXMLStructuredGridReader()
#    reader.SetFileName(file_path)
#    reader.Update()
#    return reader.GetOutput()
#
def find_cell_centers(structured_mesh):
    cell_centers = vtk.vtkCellCenters()
    cell_centers.SetInputData(structured_mesh)
    cell_centers.Update()
    # Transfer cell data to cell centers
    # Enable copying of cell data to cell centers
    cell_centers.CopyArraysOn()

    return cell_centers.GetOutput()

def generate_finer_structured_grid_manual(input_grid, refinement_factor):
    # Create a structured grid
    finer_grid = vtk.vtkStructuredGrid()

    # Get dimensions of the input grid
    dims = input_grid.GetDimensions()
    sizes = input_grid.GetBounds()
    refined_dims = [round(dim * refinement_factor) for dim in dims]

    # Set the dimensions of the finer grid
    finer_grid.SetDimensions(refined_dims)

    # Create points for the finer grid
    points = vtk.vtkPoints()
    #for k in range(refined_dims[2]):
    #    for j in range(refined_dims[1]):
    #        for i in range(refined_dims[0]):
    #            x, y, z = i / refinement_factor, j / refinement_factor, k / refinement_factor
    #            points.InsertNextPoint(x, y, z)
    x,y,z = np.meshgrid(np.linspace(sizes[0], sizes[1], refined_dims[0]),
                          np.linspace(sizes[2], sizes[3], refined_dims[1]),
                          np.linspace(sizes[4], sizes[5], refined_dims[2]), indexing='ij')
    np_coord=np.column_stack((x.ravel(order='F'), y.ravel(order='F'), z.ravel(order='F')))
    vtk_data_array = numpy_support.numpy_to_vtk(num_array=np_coord,deep=True,array_type=vtk.VTK_FLOAT)
    points.SetData(vtk_data_array)

    # Set the points for the finer grid
    finer_grid.SetPoints(points)

    ## Interpolate cell data to point data
    #cell_to_point_data = vtk.vtkCellDataToPointData()
    #cell_to_point_data.SetInputData(input_grid)
    #cell_to_point_data.Update()
#
    ## Copy point data to the finer grid
    #finer_grid.GetPointData().ShallowCopy(cell_to_point_data.GetOutput().GetPointData())

    return finer_grid
def interpolate_data(grid_fine, grid_original):
    # Create the point interpolator
    centers_fine=find_cell_centers(grid_fine)
    centers_original=find_cell_centers(grid_original)
    locator = vtkStaticPointLocator()
    locator.SetDataSet(grid_original)
    locator.BuildLocator()


    linearKernel = vtkLinearKernel()
    #linearKernel.SetKernelFootprintToNClosest()
    #linearKernel.SetNumberOfPoints(9)
    linearKernel.SetKernelFootprintToRadius()
    linearKernel.SetRadius(3.5)

    interpolator = vtk.vtkPointInterpolator()
    interpolator.SetInputData(centers_fine)
    interpolator.SetSourceData(centers_original)
    interpolator.SetLocator(locator)
    interpolator.SetKernel(linearKernel)
    interpolator.SetNullPointsStrategyToClosestPoint()
    interpolator.Update()

    # Access the interpolated values at the points of the fine grid
    interpolated_data = interpolator.GetOutput().GetPointData()

    # Set the interpolated values into the fine grid
    grid_fine.GetCellData().ShallowCopy(interpolated_data)


def interpolate_gridded_data2(existing_mesh, cell_center_points):
    locator = vtkStaticPointLocator()
    locator.SetDataSet(existing_mesh)
    locator.BuildLocator()

    linearKernel = vtkLinearKernel()
    linearKernel.SetRadius(2.5)

    interpolator3 = vtkPointInterpolator()
    interpolator3.SetInputData(cell_center_points)
    interpolator3.SetSourceData(existing_mesh)
    interpolator3.SetKernel(linearKernel)
    interpolator3.SetLocator(locator)
    interpolator3.SetNullPointsStrategyToClosestPoint()
    

    interpolator3.Update()

    return interpolator3.GetOutput()



def interpolate_gridded_data(existing_mesh, structured_mesh):
    probe_filter = vtk.vtkProbeFilter()
    probe_filter.SetInputData(existing_mesh)
    probe_filter.SetSourceData(structured_mesh)
    probe_filter.PassCellArraysOn ()
    probe_filter.Update()
    return probe_filter.GetOutput()

def visualize_mesh(mesh, scalar_name='permeability'):
    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputData(mesh)
    mapper.SetScalarModeToUseCellData()
    mapper.SelectColorArray(scalar_name)
    
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    renderer = vtk.vtkRenderer()
    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)

    render_window_interactor = vtk.vtkRenderWindowInteractor()
    render_window_interactor.SetRenderWindow(render_window)

    renderer.AddActor(actor)
    renderer.SetBackground(1.0, 1.0, 1.0)  # Set background color to white

    render_window.Render()
    render_window_interactor.Start()

#def save_mesh(mesh, file_path):
#    writer = vtk.vtkXMLUnstructuredGridWriter()
#    writer.SetFileName(file_path)
#    writer.SetInputData(mesh)
#    writer.Write()
#
#def write_structuredXML_mesh(mesh,file_path):
#    writer = vtk.vtkXMLStructuredGridWriter()
#    writer.SetFileName(file_path)
#    writer.SetInputData(mesh)
#    writer.Write()
#
#
#
#def save_cell_centers(cell_centers, file_path):
#    writer = vtk.vtkXMLPolyDataWriter()
#    writer.SetFileName(file_path)
#    writer.SetInputData(cell_centers)
#    writer.Write()

def copy_cell_data(original_mesh, target_mesh, array_name):
    # Get the cell data from the original mesh
    original_cell_data = original_mesh.GetCellData()
    
    # Check if the array exists in the original cell data
    if not original_cell_data.HasArray(array_name):
        print(f"Array '{array_name}' not found in the original mesh.")
        return
    
    # Get the array from the original cell data
    original_array = original_cell_data.GetArray(array_name)
    
    # Create a new array for the target mesh
    target_array = vtk.vtkIntArray()
    target_array.SetName(array_name)
    
    # Copy the data from the original array to the target array
    target_array.DeepCopy(original_array)
    
    # Add the array to the target mesh's cell data
    target_mesh.GetCellData().AddArray(target_array)

def from_grid_to_mesh(id=38,folder='./',ref_factor=10):
    meshID=id
    existing_mesh_fname = f'model_vtkvessel{meshID}_1st_iter_mesh.vtk'
    structured_mesh_fname = f'vtkVilli{meshID}trunc_clip_grid4.vts'
    existing_mesh_path=os.path.join(folder,existing_mesh_fname)
    structured_mesh_path=os.path.join(folder,structured_mesh_fname)
    # Set the refinement factor (e.g., doubling the resolution)
    refinement_factor = ref_factor

    # Load existing mesh and structured mesh
    existing_mesh = read_unstructured_grid(existing_mesh_path)
    structured_mesh = read_structuredXML_mesh(structured_mesh_path)
    #structured_mesh = read_structured_mesh(structured_mesh_path)

    

    # Generate a finer structured grid
    finer_grid = generate_finer_structured_grid_manual(structured_mesh, refinement_factor)


    interpolate_data(finer_grid, structured_mesh)

    write_structuredXML_mesh(finer_grid,os.path.join(folder,f'{os.path.splitext(existing_mesh_fname)[0]}_inter.vts') )
    ## Find cell centers of the structured mesh
    #cell_center_points = find_cell_centers(structured_mesh)

    ## Save the cell centers
    #write_poly_data(cell_center_points, 'cell_centers.vtp')

    ## Interpolate cell data to the cell center points
    #interpolated_mesh = interpolate_gridded_data2(existing_mesh, cell_center_points)
    #write_poly_data(interpolated_mesh, 'fubar.vtp')  # Change the filename as needed
    ## Interpolate cell data to the cell center points
    interpolated_mesh = interpolate_gridded_data(existing_mesh, finer_grid)

    
    copy_cell_data(existing_mesh, interpolated_mesh, 'CellEntityIds')
    

    # Visualize the result
    visualize_mesh(interpolated_mesh,scalar_name='Permeability [mm^2]')

    # Save the interpolated mesh to a vtk.vtkXMLUnstructuredGridWriter()

    write_unstructuredXML_mesh(interpolated_mesh, os.path.join(folder,f'{os.path.splitext(existing_mesh_fname)[0]}_inter.vtu') ) # Change the filename as needed

def interpolate_grid_mesh(structured_mesh_fname,id=38,folder='./',ref_factor=10.0):
    
    structured_mesh_path=os.path.join(folder,structured_mesh_fname)
    # Set the refinement factor (e.g., doubling the resolution)
    refinement_factor = ref_factor

    # Load existing mesh and structured mesh
    
    structured_mesh = read_structuredXML_mesh(structured_mesh_path)
    #structured_mesh = read_structured_mesh(structured_mesh_path)

    

    # Generate a finer structured grid
    finer_grid = generate_finer_structured_grid_manual(structured_mesh, refinement_factor)


    interpolate_data(finer_grid, structured_mesh)

    write_structuredXML_mesh(finer_grid,os.path.join(folder,f'{os.path.splitext(structured_mesh_fname)[0]}_inter.vts') )



def main():
    from_grid_to_mesh(id=39,folder='./case_39')


if __name__ == "__main__":
    main()