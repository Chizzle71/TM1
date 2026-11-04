# Thermal metamaterial framework

This repo conatains all the relevant models used in the paper "Adaptive Routing and Cloaking through Thermal Metamaterial Intelligence".

# Overall description

All files are separated into the 2D and 3D folders and the folders and files within are described in the subsections below (and in much more detail in the paper).

## Geometry

This Java FEM model is to be imported into a COMSOL method within the application builder and once called will automatically ouptput the .txt files with the temperature field data. The two files in this subsection have different methods of inputting the geometry, the matrix/tensor version takes a 2D or 4D array (for the 2D or 3D model respectively) in the code itself for the theta values (defined in figure 2 of the paper) and the files version takes a single 2D space separated array .txt file or (for the 3D model) 2 .txt files per xy layer of the lattice (so for a 6 by 6 by 6 there are 2 folders as there are 2 theta typesfor each unit cell, each folder containing 6 lots of .txt files each containing a 6 by 6 space separated array). The key paramaters that can be changed in this file are in the define paramters section at the top of the files and have comments to describe them, also in the add mesh section towards the end the coarseness of the mesh can be changed in the autoMeshSize() function.

## Equivalent

This Java FEM model for the sake of running it is functionally identical to the geometry files but inside the model the complex geometry is represented by anisotropic conductivity matrices instead to save computational load.

## FDM

FDM that is a Python file that also uses the anisotropic conductivity model. There is a plot file which takes a geometry input and outputs the temperature field produced, also with changeable geometry dimensions of the unit cells, lattice size and boundary conditions. The optimiser file takes a starting geometry of usually random thetas and has the same other changeable input variables as the plot file with the addition of the thermal requirements, the output file is  either one or a set (2D or 3D) of .txt files which contain the optimised geometry.

To run the optimiser file yourself it requires a little more explanation. The thermal requirements listed at the start of the file are what are later used in the file to discriminate the performance of the geometry and therefore iterate the geometry. There are two categories of requirements, cloaks which have the goal to make the temperature in the area specified uniform, and targets which are useful for the routing application as they have the goal to reach the

## ML

Contains a Python file to run the respective models in batch mode from the terminal/powershell to obtain the sample data for training the ML models.

## Menu

Contains 2D/4D arrays for the 2D/3D models of common thermal devices or testing setups that can be used in the matrix input files for convenience.






The optimiser FDM files are the optimising code that require thermal requirements and iterate the geometry and the plot FDM files take the optimised geometry files and run them with higher divsions and convergence criteria and then plots the results - the thermal requirements/features need to be input in this file below the solver if labels want to be made on the graphs.
General recommendation for the pikl is do around 300 generations (100 at a time) (untill the score stops decreasing) at a standard deviation (sigma) of 5 and population size of 8 then decrease sigma to 1 or 2 and increase population size to 16-20 (a higher population size isn't reuired but is recommended for the 3D model due to the increased number of variables).
The optimiser files currently run at 5 division and a low convergence criteria for faster simulation times.
Change number of cores to your preference.





# Dependencies

Python 3 will be required for the FDM models and ML. COMSOL will be required for the FEM models.

To install the required python libraries for this repo:

```bash
pip install numba numpy matplotlib cma
```