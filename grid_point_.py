import numpy as np
import math
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import vtk
from vtk.numpy_interface import dataset_adapter as dsa
from VesselInterpreter import CreateVoronoiDiagram


def read_vtp_file(file_path):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(file_path)
    reader.Update()

    return reader.GetOutput()

def WritePolyData(input,filename,flag_ascii=0):
   writer = vtk.vtkXMLPolyDataWriter()
   writer.SetFileName(filename)
   writer.SetInputData(input)
   if flag_ascii:
    writer.SetDataModeToAscii()
   writer.Write()

def get_point_attribute(polydata, attribute_name):
    point_data = polydata.GetPointData()
    attribute_array = point_data.GetArray(attribute_name)

    if attribute_array:
        num_points = polydata.GetNumberOfPoints()
        attribute_values=[]
        for i in range(num_points):
            attribute_value = attribute_array.GetValue(i)
            attribute_values.append(attribute_value)
            #print(f"Point {i + 1}: {attribute_name} = {attribute_value}")
        return attribute_values
    else:
        print(f"Attribute '{attribute_name}' not found.")
        return None

def get_vtp_points(polydata):

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

def points_data_in_3d_grid(grid_size, points,data):
    grid = np.empty(grid_size, dtype=object)
    for i in range(grid_size[0]):
        for j in range(grid_size[1]):
            for k in range(grid_size[2]):
                grid[i, j, k] = []
    for point,value in zip(points,data):
        x, y, z = point
        grid_x, grid_y, grid_z = int(math.ceil(x)-1), int(math.ceil(y)-1), int(math.ceil(z)-1)
        grid[grid_x, grid_y, grid_z].append(value)
    return grid

def lines_in_3d_grid(grid_size,domain_size, points,lines):
    grid = np.empty(grid_size, dtype=object)
    percentage=np.empty(grid_size, dtype=object)
    for i in range(grid_size[0]):
        for j in range(grid_size[1]):
            for k in range(grid_size[2]):
                grid[i, j, k] = []
                percentage[i, j, k] = []
    for cell_id in range(lines.GetNumberOfCells()):
        cell = lines.GetCell(cell_id)
        x1, y1, z1 =  my_get_points(points,cell,0)
        x2, y2, z2 =  my_get_points(points,cell,1)
        grid_x1, grid_y1, grid_z1 = my_get_indices(x1, y1, z1)
        grid_x2, grid_y2, grid_z2 = my_get_indices(x2, y2, z2)
        if [grid_x1, grid_y1, grid_z1]==[grid_x2, grid_y2, grid_z2]:
            grid[grid_x1, grid_y1, grid_z1].append(cell_id)
            percentage[grid_x1, grid_y1, grid_z1].append(1)
        else:
            #grid[grid_x1, grid_y1, grid_z1].append(cell_id)
            #grid[grid_x2, grid_y2, grid_z2].append(cell_id)
            distance=pairwise_distances_numpy(np.array([x1, y1, z1]),np.array([x2, y2, z2]))
            direction=np.array([x1, y1, z1])-np.array([x2, y2, z2])
            for i in range(min(grid_x1,grid_x2),max(grid_x1,grid_x2)+1):
                for j in range(min(grid_y1,grid_y2),max(grid_y1,grid_y2)+1):
                    for k in range(min(grid_z1,grid_z2),max(grid_z1,grid_z2)+1):
                        intersection1,intersection2=intersect_segment_cube([x1, y1, z1],[x2, y2, z2],[i+0.5, j+0.5, k+0.5],[0.5,0.5,0.5])
                        distance_intersections=pairwise_distances_numpy(intersection1,intersection2)
                        grid[i,j,k].append(cell_id)
                        percentage[i,j,k].append(distance_intersections/distance)
            #inersection1,intersection2=intersect_segment_cube([x1, y1, z1],[x2, y2, z2],[grid_x1+0.5, grid_y1+0.5, grid_z1+0.5],[0.5,0.5,0.5])
            #inersection2,intersection1=intersect_segment_cube([x2, y2, z2],[x1, y1, z1],[grid_x2+0.5, grid_y2+0.5, grid_z2+0.5],[0.5,0.5,0.5])
            #distance1=pairwise_distances_numpy(np.array([x1, y1, z1]),intersection1)
            #distance2=pairwise_distances_numpy(np.array([x2, y2, z2]),intersection2)
            #percentage[grid_x1, grid_y1, grid_z1].append(distance1/distance)
            #percentage[grid_x2, grid_y2, grid_z2].append(distance2/distance)

    return grid, percentage
def my_get_points(points,cell,id):
    x,y,z=points[cell.GetPointId(id)][0],points[cell.GetPointId(id)][1],points[cell.GetPointId(id)][2]
    return x,y,z
def my_get_indices(x, y, z):
    grid_x, grid_y, grid_z = int(math.ceil(x)-1), int(math.ceil(y)-1), int(math.ceil(z)-1)
    return grid_x, grid_y, grid_z
def intersect_segment_cube(p1, p2, cube_center, cube_half_size):
    # p1 and p2 are the endpoints of the line segment
    # cube_center is the center of the cube
    # cube_half_size is half the size of the cube along each dimension
    
    # Direction vector of the line segment
    direction = np.array(p2) - np.array(p1)
    distance=pairwise_distances_numpy(np.array(p1),np.array(p2))
    # Calculate the minimum and maximum t-values for each face of the cube
    t_min = np.full_like(direction, -np.inf)
    t_max = np.full_like(direction, np.inf)

    for i in range(3):
        if np.abs(direction[i]) < 1e-10:  # Ray is parallel to the face
            if p1[i] < cube_center[i] - cube_half_size[i] or p1[i] > cube_center[i] + cube_half_size[i]:
                # Ray is outside the cube along this dimension
                return None  # No intersection
        else:
            t1 = (cube_center[i] - cube_half_size[i] - p1[i]) / direction[i]
            t2 = (cube_center[i] + cube_half_size[i] - p1[i]) / direction[i]

            t_min[i] = max(t_min[i], min(t1, t2))
            t_max[i] = min(t_max[i], max(t1, t2))

    # Check if the segment intersects the cube
    if np.any(t_min > t_max):
        return None  # No intersection

    # Find the entry and exit points of the segment on the cube
    t_entry = np.max(t_min)
    t_exit = np.min(t_max)

    # Check if the segment is behind the viewer
    if t_exit < 0:
        return None  # No intersection
    if t_entry < 0:
        intersection_point2 = np.array(p1)
    else:
        intersection_point2 = np.array(p1) + t_entry * direction

    # Calculate the intersection point
    if np.linalg.norm(t_exit * direction)<distance:
        intersection_point1 = np.array(p1) + t_exit * direction
    else:
        intersection_point1 = np.array(p2)
    return intersection_point1,intersection_point2

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

def break_polyline_to_lines(polyline,attribute_name=""):
    # Create a vtkPolyData to store the polyline
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(polyline.GetPoints())
    polydata.SetLines(polyline.GetLines())
    
    if attribute_name:
        point_data = polyline.GetPointData()
        attribute_array = point_data.GetArray(attribute_name)
    # Create a vtkCellArray to store the lines
    lines = vtk.vtkCellArray()

    # Iterate over the cells in the original polyline
    for cell_id in range(polyline.GetNumberOfCells()):
        cell = polyline.GetCell(cell_id)

        # Check if the cell is a line
        if cell.GetCellType() == vtk.VTK_LINE:
            # Add line as is
            lines.InsertNextCell(cell)
        elif cell.GetCellType() == vtk.VTK_POLY_LINE:
            # Add each line segment to the new cell array
            for i in range(cell.GetNumberOfPoints() - 1):
                line = vtk.vtkLine()
                line.GetPointIds().SetId(0, cell.GetPointId(i))
                line.GetPointIds().SetId(1, cell.GetPointId(i + 1))
                lines.InsertNextCell(line)
        else:
            print("thiscell is not a line ignored")


    # Create a new polydata to store the lines
    lines_polydata = vtk.vtkPolyData()
    lines_polydata.SetPoints(polyline.GetPoints())
    lines_polydata.SetLines(lines)
    if attribute_name:
        lines_polydata.GetPointData().AddArray(attribute_array)

    return lines_polydata


def export_vtk_structured_grid(grid,attributes,domain_size,min_coords, counts, output_file):
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
    VTK_data1 = dsa.numpyTovtkDataArray(flat_grid,'PointCount')
    structured_grid.GetCellData().AddArray(VTK_data1)
    
    flat_grid=np.ravel(attributes)
    VTK_data = dsa.numpyTovtkDataArray(flat_grid,'Porosities')
    structured_grid.GetCellData().AddArray(VTK_data)
    # Create vtkCellData and add the corunt as a numeric attribute
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
def compute_hydraulic_conductivity(lines_polidata,grid_lines,percentages,grid_size,domain_size,scale_factor):
    """
    compute the tissue hydraulic onductivity for the darcy equation using the kozeny carmann formula and tissue density of the villous space
    """
    point_data = lines_polidata.GetPointData()
    attribute_array = point_data.GetArray("MaximumInscribedSphereRadius")
    attribute_np_array=vtk.util.numpy_support.vtk_to_numpy(attribute_array)
    grid_volume=np.prod(domain_size/grid_size)
    porosities= np.ones(grid_size, dtype=float)
    for i in range(grid_size[0]):
        for j in range(grid_size[1]):
            for k in range(grid_size[2]):
                if grid_lines[i][j][k]:
                    p0s_id=[]
                    p1s_id=[]
                    p0s=[]
                    p1s=[]
                    r0s=[]
                    r1s=[]
                    element_percentages=[]
                    for cell_id,percentage in zip(grid_lines[i][j][k],percentages[i][j][k]):
                        cell = lines_polidata.GetCell(cell_id)
                        p0_id=(cell.GetPointId(0))
                        p1_id=(cell.GetPointId(1))
                        p0s_id.append(cell.GetPointId(0))
                        p1s_id.append(cell.GetPointId(1))
                        p0s.append(lines_polidata.GetPoint(p0_id))
                        p1s.append(lines_polidata.GetPoint(p1_id))
                        r0s.append(attribute_np_array[p0_id])
                        r1s.append(attribute_np_array[p1_id])
                        element_percentages.append(percentage)
                    p0s_np=np.array(p0s)
                    p1s_np=np.array(p1s)
                    r0s_np=np.array(r0s)
                    r1s_np=np.array(r1s)
                    distances=pairwise_distances_numpy(p0s_np,p1s_np)
                    volumes=compute_volumes(r0s_np,r1s_np,distances)
                    occupied_volume=np.dot(volumes,np.array(element_percentages)) #the weighted sum is obtained as the dot product between the qty array and the weights array
                    porosity=1-(occupied_volume/grid_volume)
                    if porosity<0:
                        porosities[i,j,k]=0
                    else:
                        porosities[i,j,k]=porosity
                    #print("fatto?")
    print("finito!")          
    return porosities
                    


    #k=d**2/180*(1-phi)**3/phi**2
#def compute_hydraulic_conductivity2(grid_lines,scale_factor):
#    """
#    compute the tissue hydraulic conductivity for the darcy equation using the kozeny carmann formula and local porisity according to  10.1016/j.placenta.#2022.09.008
#    """
#    k=d**2/180*(phi)**3/(1-phi)**2
def compute_volumes(r0s,r1s,heights):
     """
    Calculate volumes of 3D parallel faced truncated cones using NumPy.

    Parameters:
    - r0s: vector containing the radii of face 0.
    - r1s: vector containing the radii of face 1.
    - heights vector containing the heights of the cones.
    Returns:
    - volumes: vector containing the volumes.
    """
     volumes=1/3*np.pi*heights*(r0s**2+r1s**2+r0s*r1s)
     return volumes

def pairwise_distances_numpy(matrix1, matrix2):
    """
    Calculate pairwise distances between corresponding points in two matrices using NumPy.

    Parameters:
    - matrix1: 3 by n matrix containing coordinates of the first set of points.
    - matrix2: 3 by n matrix containing coordinates of the second set of points.

    Returns:
    - distances: n by n matrix containing pairwise distances.
    """
    matrix1 = matrix1.T  # Transpose to make each column represent a point
    matrix2 = matrix2.T  # Transpose to make each column represent a point

    # Calculate squared differences for each pair of coordinates
    squared_diff = (matrix1- matrix2) ** 2

    # Sum the squared differences along the coordinate axis
    sum_squared_diff = np.sum(squared_diff, axis=0)

    # Take the square root to get the distances
    distances = np.sqrt(sum_squared_diff)

    return distances

def select_cube_lines_and_points(indices,lines_polydata,lines_in_grid,fname='sub_cube.vtp'):
    """
    Given a coplex tree like polyline structure, and the gridding containin the line indices, extracts the lines and nodes of a grid element and write it as a separtate vtk polyline
    """
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(lines_polydata.GetPoints())
    polydata.SetLines(lines_polydata.GetLines())
    
    point_data = lines_polydata.GetPointData()
    attribute_array = point_data.GetArray("MaximumInscribedSphereRadius")
    # Create a vtkCellArray to store the lines
    subcube_lines = vtk.vtkCellArray()
    subcube_points = vtk.vtkPoints()
    old_point_id_list=[]
    # Iterate over the cells in the original polyline
    for cell_id in lines_in_grid[indices[0]][indices[1]][indices[2]]:
        cell = lines_polydata.GetCell(cell_id)
        #print(cell_id)
        for i in range(cell.GetNumberOfPoints() - 1):
            subcube_line = vtk.vtkLine()
            beg_point_id=cell.GetPointId(i)
            if beg_point_id not in old_point_id_list:
                old_point_id_list.append(beg_point_id)
                new_point=polydata.GetPoint(beg_point_id)
                subcube_points.InsertNextPoint(new_point)
                new_beg_point_id=len(old_point_id_list)-1

            else:
                new_beg_point_id=old_point_id_list.index(beg_point_id)
            end_point_id=cell.GetPointId(i + 1)
            if end_point_id not in old_point_id_list:
                old_point_id_list.append(end_point_id)
                new_point=polydata.GetPoint(end_point_id)
                subcube_points.InsertNextPoint(new_point)
                new_end_point_id=len(old_point_id_list)-1

            else:
                new_end_point_id=old_point_id_list.index(end_point_id)
            subcube_line.GetPointIds().SetId(0, new_beg_point_id)
            subcube_line.GetPointIds().SetId(1, new_end_point_id)
            subcube_lines.InsertNextCell(subcube_line)
    
    newRadiusArray = vtk.vtkDoubleArray()
    newRadiusArray.SetNumberOfComponents(1)
    newRadiusArray.SetNumberOfTuples(len(old_point_id_list))
    newRadiusArray.SetName("MaximumInscribedSphereRadius")
    newRadiusArray.FillComponent(0,0.0)
    for i in range(len(old_point_id_list)):
        radius=attribute_array.GetValue(old_point_id_list[i])
        newRadiusArray.SetTuple1(i,radius)
    subcube_lines_polydata = vtk.vtkPolyData()
    subcube_lines_polydata.SetPoints(subcube_points)
    subcube_lines_polydata.SetLines(subcube_lines)
    subcube_lines_polydata.GetPointData().AddArray(newRadiusArray)
    WritePolyData(subcube_lines_polydata,fname)
    #    # Check if the cell is a line
    #    if cell.GetCellType() == vtk.VTK_LINE:
    #        # Add line as is
    #        lines.InsertNextCell(cell)
    #    elif cell.GetCellType() == vtk.VTK_POLY_LINE:
    #        # Add each line segment to the new cell array
    #        for i in range(cell.GetNumberOfPoints() - 1):
    #            line = vtk.vtkLine()
    #            line.GetPointIds().SetId(0, cell.GetPointId(i))
    #            line.GetPointIds().SetId(1, cell.GetPointId(i + 1))
    #            lines.InsertNextCell(line)
    #    else:
    #        print("thiscell is not a line ignored")
#
#
    ## Create a new polydata to store the lines
    #subcube_lines_polydata = vtk.vtkPolyData()
    #subcube_lines_polydata.SetPoints(polyline.GetPoints())
    #subcube_lines_polydata.SetLines(lines)
    #return

def main():
    grid_size = (15, 15, 15)
    vtp_file_folder = './'
    vtp_ifile_name='vtkVilli29trunc.vtp'
    vtp_ifile_path=os.path.join(vtp_file_folder,vtp_ifile_name)
    vtp_ofile_path=os.path.join(vtp_file_folder,f'{os.path.splitext(vtp_ifile_name)[0]}_ls.vtp')
    attribute_name = "MaximumInscribedSphereRadius"

    # Read points from the VTP file
    polydata=read_vtp_file(vtp_ifile_path)
    points = get_vtp_points(polydata)
    attribute_array = get_point_attribute(polydata, attribute_name)
    # Translate points to fit within the specified volume
    translated_points,domain_size,min_coords = transform_points(points,grid_size,1.e-4)

    # Create the 3D grid and count points in each grid element
    grid = create_3d_grid(grid_size, translated_points)
    radius_to_grid=points_data_in_3d_grid(grid_size, translated_points,attribute_array)
    lines_polydata=break_polyline_to_lines(polydata,attribute_name)

    WritePolyData(lines_polydata,vtp_ofile_path)
    lines_in_grid,percentages=lines_in_3d_grid(grid_size,domain_size,translated_points,lines_polydata)
    porosities=compute_hydraulic_conductivity(lines_polydata,lines_in_grid,percentages,grid_size,domain_size,1)
    indices=np.unravel_index(np.argmax(grid, axis=None), grid.shape)
    select_cube_lines_and_points(indices,lines_polydata,lines_in_grid,'sub_cube29.vtp')
    CreateVoronoiDiagram("sub_cube29.vtp",32,'vor_diag_sub_cube_29.vtp')
    counts = count_points_in_grid(grid)

    # Print the result
    for i, count in enumerate(counts):
        print(f"Grid element {i + 1}: {count} points")

    # Plot the density map
    plot_density_map(grid)

    # Export the grid as VTK StructuredGrid with the count as a numeric attribute
    #output_vtk_file = 'output_grid28.vts'
    export_vtk_structured_grid(grid,porosities,domain_size,min_coords, counts,f'{os.path.splitext(vtp_ifile_name)[0]}_grid2.vts')

if __name__ == "__main__":
    main()
