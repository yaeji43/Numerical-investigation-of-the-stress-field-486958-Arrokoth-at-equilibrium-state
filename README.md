# Finite Element Model (FEM) - (486958) Arrokoth

This FEM package is developed to calculate the stress field of (486958) Arrokoth. 
The stress calculation is based on linear-elastic deformation. 
It is considered to have no large shape deformation, which means the body keeps its initial shape all the time while applying the loading vectors on the object. 
The general formulation of the pacakge is well described in [1]. 
This package provides the stress field of Arrokoth at the equlibrium state which will used to generate the advanced stress calculation (Figure 1 in [2]).


[1] [Yaeji Kim and Masatoshi Hirabayashi. A numerical approach using a finite element model to constrain the possible interior layout of (16) psyche. The Planetary Science Journal, 3(5):122, 2022.](https://iopscience.iop.org/article/10.3847/PSJ/ac6b39/meta) 

[2] [AGU Fall Meeting 2022, held in Chicago, IL, 12-16 December 2022, id. P26A-08.](https://baas.aas.org/pub/2022n8i410p01/release/1)


# How to run the package
1. Run the shell script <br/>
```./exe.sh```

2. Run the generated executive file <br/>
```./main```

3. Then the output file named 'Output.txt' will be generated <br/>
4. [Optional] You can use 'Visualize_static_FEM.m' to visuazlie the stress field of Arrokoth of the entire and sliced structure with the output file [Output.txt].


# Output file structure

If you successfully run the FEM package, you will get the output file <Output.txt>. This file includes the node, stress vector, density, yield stress, and mass information on all nodes.

%% +++++++++++++++++++ SIMULATION INFO +++++++++++++++++++++ <br/>
%% Node    number:= # <br/>
%% Element number:= # <br/>
%% Spin period:= 15.92 <br/>
%% Node Info ================================= <br/> <br/>
_Node data. <br/>
_Each row has an x-, y-, and z- component of the assigned node (i) - (x_i, y_i, z_i) .. <br/> <br/>
%% Stress Info ================================= <br/> <br/>
_Stress data .. <br/>
_Each row has a stress field acting on the assigned node (i) - (sigma_xx_i, sigma_yy_i, sigma_zz_i, sigma_xy_i, sigma_yz_i, sigma_zx_i) .. <br/> <br/>
%% Density Info ============================= <br/> <br/>
_Density data .. <br/>
_Each row has a density value at the assigned node - (density_i) .. <br/> <br/>
%%  Yield stress ============================= <br/> <br/>
_Yield stress data .. <br/>
_Each row has the yield stress at the assigned node (i) ..  <br/>
_The yield stress is calculated based on the drucker-prager criterion - (Yield_stress_i) ..  <br/> <br/>
%% Mass distribution ============================= <br/> <br/>
_Mass data .. <br/>
_Each row has the mass value at the assigned node (i) - (mass_i) ..  <br/>






