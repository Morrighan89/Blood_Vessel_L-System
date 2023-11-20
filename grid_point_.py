import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import vtk
from vtk.numpy_interface import dataset_adapter as dsa


def read_vtp(file_path):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(file_path)
    reader.Update()

    polydata = reader.GetOutput()

    points = []
    for i in range(polydata.GetNumberOfPoints()):
        point = polydata.GetPoint(i)
        points.append(point)

    return points

def translate_points(points):
    min_coords = np.min(points, axis=0)
    max_coords = np.max(points, axis=0)

    translated_points = points - min_coords

    return translated_points


def transform_points(points,grid_size,epsilon):
    min_coords = np.min(points, axis=0)-epsilon
    max_coords = np.max(points, axis=0)+epsilon
    domain_size=max_coords-min_coords
    grid_spaces=(grid_size[0],grid_size[1],grid_size[2])
    translated_points = (points - min_coords)/domain_size*(grid_spaces)

    return translated_points,domain_size,min_coords




def create_3d_grid(grid_size, points):
    grid = np.zeros(grid_size, dtype=int)

    for point in points:
        x, y, z = point
        grid_x, grid_y, grid_z = int(math.ceil(x)-1), int(math.ceil(y)-1), int(math.ceil(z)-1)
        grid[grid_x, grid_y, grid_z] += 1

    return grid

def count_points_in_grid(grid):
    counts = np.bincount(grid.flatten())
    counts = counts[1:]

    return counts


def plot_density_map(grid):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Get non-zero coordinates and counts
    x, y, z = np.nonzero(grid)
    counts = grid[x, y, z]

    ax.scatter(x, y, z, c=counts, marker='o', s=counts*10, cmap='viridis')

    ax.set_xlabel('X Axis')
    ax.set_ylabel('Y Axis')
    ax.set_zlabel('Z Axis')
    ax.set_title('3D Density Map of Points')

    plt.show()

def export_vtk_structured_grid(grid,domain_size,min_coords, counts, output_file):
    # Create a VTK StructuredGrid
    structured_grid = vtk.vtkStructuredGrid()
    structured_grid.SetDimensions((grid.shape[0]+1,grid.shape[1]+1,grid.shape[2]+1))

    # Create vtkPoints and add them to the structured grid
    points = vtk.vtkPoints()
    for i in range(grid.shape[0]+1):
        for j in range(grid.shape[1]+1):
            for k in range(grid.shape[2]+1):
                
                points.InsertNextPoint(i/grid.shape[0]*domain_size[0]+min_coords[0], j/grid.shape[1]*domain_size[1]+min_coords[1], k/grid.shape[2]*domain_size[2]+min_coords[2])

    structured_grid.SetPoints(points)

    # Convert counts to 32-bit signed integers
    counts_int = counts.astype(np.int32)
    flat_grid=np.ravel(grid)

    VTK_data = dsa.numpyTovtkDataArray(flat_grid,'PointCount')
    structured_grid.GetCellData().SetScalars(VTK_data)
    # Create vtkCellData and add the count as a numeric attribute
 #   cell_data = vtk.vtkCellData()
 #   count_array = vtk.vtkIntArray()
 #   count_array.SetName('PointCount')
 #   count_array.SetNumberOfComponents(1)
 #   vtk_data = numpy_support.numpy_to_vtk(num_array=counts_int, deep=True, array_type=vtk.VTK_INT)
 #   count_array.SetArray(vtk_data, len(vtk_data), True)
#
 #  
#
 #   cell_data.AddArray(count_array)
 #   structured_grid.GetCellData().ShallowCopy(cell_data)

    # Write the StructuredGrid to a VTK file
    writer = vtk.vtkXMLStructuredGridWriter()
    writer.SetFileName(output_file)
    writer.SetInputData(structured_grid)
    writer.Write()
    




def main():
    grid_size = (20, 20, 20)
    vtp_file_path = './vtkVilli21.vtp'

    # Read points from the VTP file
    points = read_vtp(vtp_file_path)

    # Translate points to fit within the specified volume
    translated_points,domain_size,min_coords = transform_points(points,grid_size,1.e-4)

    # Create the 3D grid and count points in each grid element
    grid = create_3d_grid(grid_size, translated_points)
    counts = count_points_in_grid(grid)

    # Print the result
    for i, count in enumerate(counts):
        print(f"Grid element {i + 1}: {count} points")

    # Plot the density map
    plot_density_map(grid)

    # Export the grid as VTK StructuredGrid with the count as a numeric attribute
    output_vtk_file = 'output_grid21.vts'
    export_vtk_structured_grid(grid,domain_size,min_coords, counts, output_vtk_file)

if __name__ == "__main__":
    main()
