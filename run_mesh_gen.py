import Clipp_with_box
from grid_point_ import (
    run_villi_to_grid
)
from From_grid_to_mesh import (
    from_grid_to_mesh
)
from Build_GMSH_mesh_from_PolyLine_maintrunc import (
    create_gmsh_mesh_carving
)
import os

def main():
    id=39
    folder='case_39'
    box_boundaries=(0,20,-18,18,-18,18)
    grid_size=(9,19,19)

    ## Clipps the random villi structure to the defined bounding bos for the simulation
    #original_villi=Clipp_with_box.read_poly_data(os.path.join(folder,f'vtkVilli{id}Trunc.vtp'))
    #clipped_polyline=Clipp_with_box.clip_polyline_with_box(original_villi,box_boundaries)
    #Clipp_with_box.write_poly_data(clipped_polyline, os.path.join(folder,f'vtkVilli{id}Trunc_clip.vtp'),flag_ascii=1)
    #### Readas the clipped villi structure  and superimpose the grid to compute the spatially averaged quantities for the darcy model
    run_villi_to_grid(grid_size,folder,id)
    ### Create a tetrahedral mesh of the box including the first branches of the villi and add inletand outlet arteries
    #create_gmsh_mesh_carving(f'vtkVilli{id}_1st_iter.vtp',folder,ofile=f'vtkvessel{id}_1st_iter_mesh.msh',scalingfactor=1,brep_module=400)
    ### interpolates the properties from the grid to the tetrahedral mesh
    from_grid_to_mesh(id=id,ref_factor=10,folder=folder)

if __name__ == "__main__":
    main()