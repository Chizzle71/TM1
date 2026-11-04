# Thermal metamaterials

This repo conatains all the relevant files used in the paper "Adaptive Routing and Cloaking through Thermal Metamaterial Intelligence".

## Description

All files are separated into the 2D and 3D folders and the folders/files within are described below:

- Geometry

- Equivalent

- FDM

- ML

Contains a python file to run the respective models in batch mode to obtain the sample data for training the ML models.

- Menu

All files are separated within the 2D and 3D folders, with an additional folder for ML which contains a python file to run the respective models in batch mode to obtain the sample data for training the ML models, and an additional 2D/3DMenu file which contains 2D/4D arrays of common thermal devices or testing setups that can be used in the matrix input files. Each model has 2 ways to enter the geometry of the unit cells, either through an array in the code or by reading .txt files. The Java files are to be imported into the application builder in COMSOL and made into a method then the method can be ran to obtain results.
The first model is called Geometry and is a FEM to be imported into Java.
The second model is called Equivalent and is also a FEM to be imported into Java but where the compex geometry is represented by anistropic conductivity matrics insted.
The third model is FDM and is a Python file that also uses the anisotropic conductivity model. There is a plot file which takes a geometry input and outputs the temperature field produced, also with changeable geometry dimensions of the unit cells, lattice size and boundary conditions. The optimiser file takes a starting geometry of usually random thetas and has the same other changeable input variablea as the pllot file with the addition of the thermal requirements and the output file is instead a set of .txt files with the optimised geometry.
Each model is explained in much more detail in the paper.

These are the first main files. The optimiser FDM files are the optimising code that require thermal requirements and iterate the geometry and the plot FDM files take the optimised geometry files and run them with higher divsions and convergence criteria and then plots the results - the thermal requirements/features need to be input in this file below the solver if labels want to be made on the graphs.
General recommendation for the pikl is do around 300 generations (100 at a time) (untill the score stops decreasing) at a standard deviation (sigma) of 5 and population size of 8 then decrease sigma to 1 or 2 and increase population size to 16-20 (a higher population size isn't reuired but is recommended for the 3D model due to the increased number of variables).
The optimiser files currently run at 5 division and a low convergence criteria for faster simulation times.
Change number of cores to your preference.

## Dependencies

Python 3 will be required for the FDM models and ML. COMSOL will be required for the FEM models.

To install the required python libraries for this repo:

```bash
pip install numba numpy matplotlib cma
```