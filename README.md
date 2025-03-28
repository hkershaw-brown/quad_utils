## DART Interpolation on quads

Plots from https://github.com/NCAR/DART/issues/833

For notebook plot_func this uses the fortran code quad_utils_mod.f90  
To use this, compile with f2py

``
f2py -c -m quad_utils quad_utils_mod.f90
``

IWD : Inverse Weighted Distance

Notebooks:

| ipynb         | Description | 
|---------------| ------------|
| plot_func.ipynb |  call fortran quad_utils_mod |
| plot-interps.ipynb | IWD plot numpy |
| plot-interps-great-circle.ipynb | IWD plot numpy |
| globe-plot.ipynb | location of the quad on the Earth |


There is also quad_barycentric_mod.f90, test_barycentric.f90, test_jacobian.f90 to play 
with quad_barycentric. 