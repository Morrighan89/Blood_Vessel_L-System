import vtk
from tqdm import tqdm
from vtk_read_write import (
    read_fluent_data,
    read_unstructured_grid,

    write_unstructuredXML_mesh
)
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
mesh_path=r"E:\Riccardo\Git\Blood_Vessel_L-System\case_39\model_vtkvessel39_1st_iter_mesh.vtk"
target_mesh = read_unstructured_grid(mesh_path)

# --- 4. Set up a point locator (VTK’s spatial search structure) ---
point_locator = vtk.vtkPointLocator()
point_locator.SetDataSet(merged_cfd_output)
point_locator.BuildLocator()

# --- 5. Create new arrays to store interpolated data ---
with tqdm(total=merged_cfd_output.GetPointData().GetNumberOfArrays()*target_mesh.GetNumberOfPoints(),miniters=250, desc="Overall Progress") as pbar:
    for i in range(merged_cfd_output.GetPointData().GetNumberOfArrays()):
        array = merged_cfd_output.GetPointData().GetArray(i)
        name = array.GetName()

        interpolated_array = vtk.vtkFloatArray()
        interpolated_array.SetName(name)
        interpolated_array.SetNumberOfComponents(array.GetNumberOfComponents())
        interpolated_array.SetNumberOfTuples(target_mesh.GetNumberOfPoints())

        # Interpolate by nearest neighbor
        for pt_id in range(target_mesh.GetNumberOfPoints()):
            pt = target_mesh.GetPoint(pt_id)
            nearest_id = point_locator.FindClosestPoint(pt)
            interpolated_array.SetTuple(pt_id, array.GetTuple(nearest_id))
            pbar.update(1)
        target_mesh.GetPointData().AddArray(interpolated_array)

# --- 6. Write the new interpolated mesh ---
write_unstructuredXML_mesh(target_mesh,"fluent_interp.vtu")