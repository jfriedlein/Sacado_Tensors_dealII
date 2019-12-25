# Testing Trilinos::Sacado mainly with tensors
Testing capabilities of Sacado mainly in combination with tensors

## The Aim
The aim is to compute e. g. the tangent

<a href="https://www.codecogs.com/eqnedit.php?latex=\overset{4}{C}&space;=&space;\frac{\partial\boldsymbol{\sigma}}{\partial\boldsymbol{\varepsilon&space;}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\overset{4}{C}&space;=&space;\frac{\partial\boldsymbol{\sigma}}{\partial\boldsymbol{\varepsilon&space;}}" title="\overset{4}{C} = \frac{\partial\boldsymbol{\sigma}}{\partial\boldsymbol{\varepsilon }}" /></a>

as a fourth order tensor on the material point level (quadrature point) based on the implementation of the stress-equation sigma(eps,phi) only. Similarly, we can compute the tangent in combinatin with a scalar variable, such as the scalar damage phi(eps):

<a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{A}&space;=&space;\frac{\partial\boldsymbol{\sigma}}{\partial{\phi&space;}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{A}&space;=&space;\frac{\partial\boldsymbol{\sigma}}{\partial{\phi&space;}}" title="\boldsymbol{A} = \frac{\partial\boldsymbol{\sigma}}{\partial{\phi }}" /></a>

 or 
 
 <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{B}&space;=&space;\frac{\partial{\phi}}{\partial\boldsymbol{\varepsilon&space;}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{B}&space;=&space;\frac{\partial{\phi}}{\partial\boldsymbol{\varepsilon&space;}}" title="\boldsymbol{B} = \frac{\partial{\phi}}{\partial\boldsymbol{\varepsilon }}" /></a>
 
You can use Sacado to compute general derivatives of functions (with or without tensors) with respect to variables (double, Tensors, ...).

## The Documentation
The Doxygen documentation for the code can be found here https://jfriedlein.github.io/Sacado-Testing/html/index.html. (Please ignore the folder /documentation and take a closer look at the folder /docs, if you are interested in the details.)

## ToDo:
- reference the Wrapper
- enable LaTeX equations in the documentation hosted via GitHub
- Rework the design and structure of the Wrapper
- implement more data types (e.g. Vector, nonsym tensor) and compatibility with more combinations (e.g. vector-vector, vector-double).
