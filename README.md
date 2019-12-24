# Trilinos-Sacado-Testing
Testing capabilities of Sacado mainly in combination with tensors

aim: compute e. g. Tangent

<a href="https://www.codecogs.com/eqnedit.php?latex=\overset{4}{C}&space;=&space;\frac{\partial\boldsymbol{\sigma}}{\partial\boldsymbol{\varepsilon&space;}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\overset{4}{C}&space;=&space;\frac{\partial\boldsymbol{\sigma}}{\partial\boldsymbol{\varepsilon&space;}}" title="\overset{4}{C} = \frac{\partial\boldsymbol{\sigma}}{\partial\boldsymbol{\varepsilon }}" /></a>

as a fourth order tensor on the material point level (quadrature point) based on the implementation of the stress-equation only.

You can use Sacado to compute the derivatives of functions (with or without tensors) with respect to variables (double, Tensors, ...).

The Doxygen documentation for the code can be found in the documentation.zip. Just download this zip-file unpack it and open the index.html in the folder /html in your webbrowser.

ToDo:
- Wrapper
