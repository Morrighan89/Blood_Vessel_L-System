import vtk
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor
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

# --- 4. Interpolation function ---
def interpolate_field(array_name):
    source_array = merged_cfd_output.GetPointData().GetArray(array_name)

    interp_array = vtk.vtkFloatArray()
    interp_array.SetName(array_name)
    interp_array.SetNumberOfComponents(source_array.GetNumberOfComponents())
    interp_array.SetNumberOfTuples(target_mesh.GetNumberOfPoints())

    # Define a worker for a chunk of point IDs
    def process_range(start, end):
        for pt_id in range(start, end):
            pt = target_mesh.GetPoint(pt_id)
            nearest_id = point_locator.FindClosestPoint(pt)
            interp_array.SetTuple(pt_id, source_array.GetTuple(nearest_id))

    # Split work into chunks
    n_points = target_mesh.GetNumberOfPoints()
    n_threads = 64
    chunk_size = (n_points + n_threads - 1) // n_threads

    with ThreadPoolExecutor(max_workers=n_threads) as executor:
        futures = []
        for i in range(n_threads):
            start = i * chunk_size
            end = min(start + chunk_size, n_points)
            futures.append(executor.submit(process_range, start, end))
        [f.result() for f in futures]  # wait for all to finish

    return interp_array

# --- 5. Run interpolation in parallel over all fields ---
field_names = [
    merged_cfd_output.GetPointData().GetArray(i).GetName()
    for i in range(merged_cfd_output.GetPointData().GetNumberOfArrays())
]

with ThreadPoolExecutor(max_workers=8) as executor:  # parallel over fields too
    results = list(executor.map(interpolate_field, field_names))

# --- 6. Add arrays to mesh ---
for array in results:
    target_mesh.GetPointData().AddArray(array)
# --- 6. Write the new interpolated mesh ---
write_unstructuredXML_mesh(target_mesh,"fluent_interp.vtu")