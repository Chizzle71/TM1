# Thermal metamaterials

This repository conatains all the relevant files used in the paper "Adaptive Routing and Cloaking through Thermal Metamaterial Intelligence".

## Descprition

These are the first main files. The optimiser FDM files are the optimising code that require thermal requirements and iterate the geometry and the plot FDM files take the optimised geometry files and run them with higher divsions and convergence criteria and then plots the results - the thermal requirements/features need to be input in this file below the solver if labels want to be made on the graphs.
General recommendation for the pikl is do around 300 generations (100 at a time) (untill the score stops decreasing) at a standard deviation (sigma) of 5 and population size of 8 then decrease sigma to 1 or 2 and increase population size to 16-20 (a higher population size isn't reuired but is recommended for the 3D model due to the increased number of variables).
The optimiser files currently run at 5 division and a low convergence criteria for faster simulation times.
Change number of cores to your preference.

# Dependencies

```bash
pip install numba numpy matplotlib cma
```