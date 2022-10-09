# L-system : Probabilistic BloodVessel 3d Generator using Python
This project allows you to define [L-systems] to define blood vessel networks in a fractal like fashion.

## Example

Even branching            |  Uneven branching
:-------------------------:|:-------------------------:
![Alt text](Vessel_Ex3.PNG?raw=true "Blood Vessel Example") |  ![Alt text](Vessel_Ex2.PNG?raw=true "Blood Vessel Example")



### More information
Basic understanding of the L-System to create fractal structure can be found in the following articles

1. [L-systems : draw nice fractals and plants (part I)](https://medium.com/@hhtun21/l-systems-draw-your-first-fractals-139ed0bfcac2)
2. [L-systems: draw a stochastic plant (part II)](https://medium.com/@hhtun21/l-systems-draw-a-stochastic-plant-ii-f322df2ea3c5)

# Installation note
To use this scripts the following packages are required

- NumPy
- imageio
- appdirs
- meshio
- scooby
- tqdm
- MatplotLib
- vtk
- PyVista

For more details refer to the single packages installation requirements.

Python **3.6.+** is **required**

We suggest to create a new environment using anaconda and proceed to install all the packages from conda forge using the example code provided


```
conda create -n myenv python=3.6
conda activate myenv
conda install -c conda-forge numpy
conda install -c conda-forge imageio
conda install -c conda-forge appdirs
conda install -c conda-forge meshio
conda install -c conda-forge scooby
conda install -c conda-forge tdqm
conda install -c conda-forge matplotlib
conda install -c conda-forge vtk
conda install -c conda-forge PyVista
conda install -c vmtk vtk itk vmtk
```

Please refer also to the individual packages installation notes

# Steps and instruction

1. (conda activate myenv) and Run
   ```python
    Vessel_Generator.py
    ```
    until desired branch is generated.

    Branching parameters can be tweaked in the script Vessel_generator.py acting on the rule weights. Rule can be modified in Lsystem_class.py
	
	Now you have a nice set of branching points and bloodvessel diameters, as well as blood vessels centerlines. If you need a surface go further.
2.  Run 
    ```python
    FromVoronoiToSurface.py
    ```
    Modify the script FromVoronoiToSurface.py cpassing the correct input filename in line 518  to use the file generated previously in Vessel_Generator.py in this line
     ```python
     VesselInterpreter.CreateVoronoiDiagram(clFileName,numberOfInterpolationPoints,ofile="voronoiDiagram.vtp")
     ```
	
3. use vmtk pypad to decimate the obtained surface
    ```
    vmtksurfacedecimation -ifile reconstructedmodel.vtp -ofile  reconstructedmodel_dec.vtp
    ```
4. compute the centerlines from the decimated surface and the list of output points provided by the script Vessel_Generator.py
    ```
    vmtksurfacereader -ifile reconstructedmodel_dec.vtp --pipe     vmtkcenterlines  -seedselector pointlist -sourcepoints 0 0 0    -targetpoints 1 1 1 2 2 2 3 3 3 -ofile reconstructedmodel05_dec_cl.vtp
    ```
5. remesh the surface with vmtk (use the centerlines to obtain a radius adaptive mesh)
    ```
    vmtkdistancetocenterlines -ifile reconstructedmodel_dec.vtp     -centerlinesfile reconstructedmodel_dec_cl.vtp -useradius 1    --pipe vmtksurfaceremeshing -elementsizemode edgelengtharray   -edgelengtharray DistanceToCenterlines -edgelengthfactor 0.3 -ofile   reconstructedmodel_re.vtp
    ```
6. Recompute nice and smooth centerlines
    ```
    vmtksurfacereader -ifile reconstructedmodel_re.vtp --pipe  vmtkcenterlines  -seedselector pointlist -sourcepoints 0 0 0     -targetpoints 1 1 1 2 2 2 3 3 3 -ofile reconstructedmodel_re_cl.vtp
    ```
7. use them to open the ends of the surface and create inlet and outlets
    ```
    vmtkendpointextractor -ifile reconstructedmodel_re_cl.vtp -radiusarray    MaximumInscribedSphereRadius -ofile reconstructedmodel_re_cl_cli.vtp  --pipe vmtkbranchclipper -ifile reconstructedmodel_re.vtp --pipe    vmtksurfaceconnectivity -cleanoutput 1 --pipe vmtksurfacewriter -ofile     reconstructedmodel_re_ct.vtp
    ```
8. generate volumic mesh
    ```
    vmtkdistancetocenterlines -ifile reconstructedmodel_re_ct.vtp     -centerlinesfile reconstructedmodel05_re_cl.vtp -useradius 1    --pipe vmtkmeshgenerator -elementsizemode edgelengtharray   -edgelengtharray DistanceToCenterlines -edgelengthfactor 0.3 -ofile   reconstructedmodel_vol.vtu
    ```
9. Inspect bounbdary
    ```
    vmtkmeshboundaryinspector -ifile reconstructedmodel_vol.vtu -entityidsarray CellEntityIds 
    ```
    

