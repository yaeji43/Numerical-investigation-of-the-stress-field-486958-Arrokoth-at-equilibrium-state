# Static_FEM_Arrokoth

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
