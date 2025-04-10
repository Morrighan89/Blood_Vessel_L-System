import vtk
from tqdm import tqdm
from vtk_read_write import (
    read_fluent_data,
    read_unstructuredXML_grid,

    write_unstructuredXML_mesh
)
from From_grid_to_mesh import interpolate_gridded_data
# --- 1. Load the CFD EnSight Gold mesh ---
cfd_path=r"E:\LifeSaver_sim_results\Nabovati_data\ivs_porous_v2_nab.case"
cfd_output = read_fluent_data(cfd_path)

# --- 2. Unpack the multiblock dataset ---
# Collect all blocks into one merged unstructured grid
append_filter = vtk.vtkAppendFilter()

if isinstance(cfd_output, vtk.vtkMultiBlockDataSet):
    for i in range(cfd_output.GetNumberOfBlocks()):
        block = cfd_output.GetBlock(i)
        if isinstance(block, vtk.vtkDataSet):
            append_filter.AddInputData(block)
else:
    append_filter.AddInputData(cfd_output)  # Single block case

append_filter.Update()
merged_cfd_output = append_filter.GetOutput()

# --- 3. Load the simple tetrahedral mesh (e.g., target_mesh.vtu) ---
mesh_path=r"E:\Riccardo\Git\Blood_Vessel_L-System\case_38\model_vtkvessel38_1st_iter_mesh_inter.vtu"
target_mesh = read_unstructuredXML_grid(mesh_path)

# --- Transform target mesh ---
transform = vtk.vtkTransform()
transform.Scale(0.001, 0.001, 0.001)  # Example scale

transform_filter = vtk.vtkTransformFilter()
transform_filter.SetInputData(target_mesh)
transform_filter.SetTransform(transform)
transform_filter.Update()
scaled_mesh = transform_filter.GetOutput()


## --- Use vtkProbeFilter for interpolation ---
#def progress_callback(caller, event):
#    progress = caller.GetProgress() * 100
#    print(f"Progress: {progress:.1f}%")
#
## --- Setup the probe filter as before ---
#probe = vtk.vtkProbeFilter()
#probe.SetSourceData(merged_cfd_output)
#probe.SetInputData(scaled_mesh)
#
## Attach the progress observer
#probe.AddObserver("ProgressEvent", progress_callback)
#
## Run the filter (will print progress updates)
#probe.Update()

output_mesh=interpolate_gridded_data(scaled_mesh,merged_cfd_output)

# Get the probed output


# Copy original arrays from the target mesh
target_point_data = scaled_mesh.GetPointData()
probed_point_data = output_mesh.GetPointData()

for i in range(target_point_data.GetNumberOfArrays()):
    array = target_point_data.GetArray(i)
    # Avoid overwriting same-named arrays
    if not probed_point_data.HasArray(array.GetName()):
        probed_point_data.AddArray(array)


target_cell_data = scaled_mesh.GetCellData()
probed_cell_data = output_mesh.GetCellData()

for i in range(target_cell_data.GetNumberOfArrays()):
    array = target_cell_data.GetArray(i)
    if not probed_cell_data.HasArray(array.GetName()):
        probed_cell_data.AddArray(array)


## --- 4. Set up a point locator (VTK’s spatial search structure) ---
#point_locator = vtk.vtkPointLocator()
#point_locator.SetDataSet(merged_cfd_output)
#point_locator.BuildLocator()
#
## --- 5. Create new arrays to store interpolated data ---
#with tqdm(total=merged_cfd_output.GetPointData().GetNumberOfArrays()*target_mesh.GetNumberOfPoints(),miniters=250, desc="Overall Progress") as pbar:
#    for i in range(merged_cfd_output.GetPointData().GetNumberOfArrays()):
#        array = merged_cfd_output.GetPointData().GetArray(i)
#        name = array.GetName()
#
#        interpolated_array = vtk.vtkFloatArray()
#        interpolated_array.SetName(name)
#        interpolated_array.SetNumberOfComponents(array.GetNumberOfComponents())
#        interpolated_array.SetNumberOfTuples(target_mesh.GetNumberOfPoints())
#
#        # Interpolate by nearest neighbor
#        for pt_id in range(target_mesh.GetNumberOfPoints()):
#            pt = target_mesh.GetPoint(pt_id)
#            nearest_id = point_locator.FindClosestPoint(pt)
#            interpolated_array.SetTuple(pt_id, array.GetTuple(nearest_id))
#            pbar.update(1)
#        target_mesh.GetPointData().AddArray(interpolated_array)
#
# --- 6. Write the new interpolated mesh ---
write_unstructuredXML_mesh(output_mesh,"fluent_interp.vtu")