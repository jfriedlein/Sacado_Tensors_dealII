# Trilinos::Sacado: Computing derivatives incorporating tensors
Examples and functions to compute derivatives of tensor-valued equations with Sacado

## The Aim
The aim is to compute e. g. the tangent

<a href="https://www.codecogs.com/eqnedit.php?latex=\overset{4}{C}&space;=&space;\frac{\partial\boldsymbol{\sigma}}{\partial\boldsymbol{\varepsilon&space;}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\overset{4}{C}&space;=&space;\frac{\partial\boldsymbol{\sigma}}{\partial\boldsymbol{\varepsilon&space;}}" title="\overset{4}{C} = \frac{\partial\boldsymbol{\sigma}}{\partial\boldsymbol{\varepsilon }}" /></a>

as a fourth order tensor on the material point level (quadrature point) based on the implementation of the stress-equation sigma=sigma(eps,phi) only. Similarly, we can compute the tangent in combination with a scalar variable, such as the scalar damage phi and compute second derivatives (see feature list).
 
The here shown code only implements functions (and finally the Wrapper) to pack the derivatives related to tensors into a nice format, pass them to Sacado to compute the derivatives and unpack the results back into tensors.

This approach might be useful when you want to compute e. g. the tangent modulus at quadrature points and keep everything in an enclosed material model function/file. An example: You developed or found a material model, lets say elasto-plasticity with some saturated hardening, and want to implement this model to compare it to another one. So you implement the equations and possibly subiterations on the material point level into your material model that gets the strain and the history variables as inputs. Now you would need to derive the tangents (e.g. d_stress/d_strain) to be able to assemble the tangent. However, the latter can take some effort especially for more complicated models or if the equations are still being developed as we speak. With the Sacado_Wrapper you implement the equations (e.g. stress equation) as you would normally do, change the data types (when you use templated functions here, that is done automatically) and set the dofs as shown in the documentation. In the end you call get_tangent() and there you have your tangent (quadratic convergence without effort). When you're satisfied with the material model, you can derive the tangent by hand and simply change the data types back and replace the get_tangent() call with your analytical tangent equation. To sum up, the Wrapper enables you keep the standard material model framework and still use the power (and beauty) of automatic differentation.

It will probably be significantly more efficient if possible to assemble the residuum and compute its derivatives as shown in, for instance, the deal.ii tutorial 33 https://www.dealii.org/current/doxygen/deal.II/step_33.html with already implemented deal.ii-functionality. Also keep in mind that by automatically computing the tangent from the residuum you don't have to derive the linearisation of the residuum by hand, which can be advantageous for complicated PDEs. When you use this Wrapper, you still have to derive the linearisation of the residuum by hand and fill it with the tangents from the material model.

## The Documentation
The Doxygen documentation for the code can be found here https://jfriedlein.github.io/Sacado-Testing/html/index.html. It shows a few examples and describes how to use the Sacado_Wrapper.

## Current features of the Sacado_Wrapper
### Legend:

<a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{\varepsilon}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{\varepsilon}" title="\boldsymbol{\varepsilon}" /></a> (2nd order symmetric tensor, here: strain tensor)

<a href="https://www.codecogs.com/eqnedit.php?latex=\varphi" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\varphi" title="\varphi" /></a> (scalar, here: global damage variable)

<a href="https://www.codecogs.com/eqnedit.php?latex=\Psi&space;=&space;\Psi(\boldsymbol{\varepsilon},\varphi)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\Psi&space;=&space;\Psi(\boldsymbol{\varepsilon},\varphi)" title="\Psi = \Psi(\boldsymbol{\varepsilon},\varphi)" /></a> (scalar, here: free energy)

<a href="https://www.codecogs.com/eqnedit.php?latex=\overset{4}{C}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\overset{4}{C}" title="\overset{4}{C}" /></a> (4th order sym. tensor, here: consistent tangent moduli)

<a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{A,&space;B,&space;E,&space;F}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{A,&space;B,&space;E,&space;F}" title="\boldsymbol{A, B, E, F}" /></a> (2nd order sym. tensors)

<a href="https://www.codecogs.com/eqnedit.php?latex=D,&space;G" target="_blank"><img src="https://latex.codecogs.com/gif.latex?D,&space;G" title="D, G" /></a> (scalar derivatives)

### Features:

- Compute derivatives of equations with respect to a single scalar.

<a href="https://www.codecogs.com/eqnedit.php?latex={D}&space;=&space;\frac{\partial{\Psi}}{\partial{\varphi&space;}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?{D}&space;=&space;\frac{\partial{\Psi}}{\partial{\varphi&space;}}" title="{D} = \frac{\partial{\Psi}}{\partial{\varphi }}" /></a>
 ; 
<a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{A}&space;=&space;\frac{\partial\boldsymbol{\sigma}}{\partial{\varphi&space;}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{A}&space;=&space;\frac{\partial\boldsymbol{\sigma}}{\partial{\varphi&space;}}" title="\boldsymbol{A} = \frac{\partial\boldsymbol{\sigma}}{\partial{\varphi }}" /></a>

- Compute tangents from equations with respect to a single SymmetricTensor<2,dim>. (Example 3B)

<a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{B}&space;=&space;\frac{\partial{\Psi}}{\partial\boldsymbol{\varepsilon&space;}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{B}&space;=&space;\frac{\partial{\Psi}}{\partial\boldsymbol{\varepsilon&space;}}" title="\boldsymbol{B} = \frac{\partial{\Psi}}{\partial\boldsymbol{\varepsilon }}" /></a>
 ; 
<a href="https://www.codecogs.com/eqnedit.php?latex=\overset{4}{C}&space;=&space;\frac{\partial\boldsymbol{\sigma}}{\partial\boldsymbol{\varepsilon&space;}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\overset{4}{C}&space;=&space;\frac{\partial\boldsymbol{\sigma}}{\partial\boldsymbol{\varepsilon&space;}}" title="\overset{4}{C} = \frac{\partial\boldsymbol{\sigma}}{\partial\boldsymbol{\varepsilon }}" /></a>

- Compute tangents from equations with respect to a single SymmetricTensor<2,dim> and a single scalar. (Example 4)

<a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{A,&space;B,&space;}&space;\overset{4}{C}&space;\textrm{&space;and&space;}&space;$D$&space;\text{&space;at&space;the&space;same&space;time}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{A,&space;B,&space;}&space;\overset{4}{C}&space;\textrm{&space;and&space;}&space;$D$&space;\text{&space;at&space;the&space;same&space;time}" title="\boldsymbol{A, B, } \overset{4}{C} \textrm{ and } $D$ \text{ at the same time}" /></a>

- Compute second derivatives of equations with respect to a single SymmetricTensor<2,dim> and a single scalar. (Example 7, 8)

<a href="https://www.codecogs.com/eqnedit.php?latex=\overset{4}{C}&space;=&space;\frac{\partial^2\Psi}{\partial\boldsymbol{\varepsilon&space;}^2}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\overset{4}{C}&space;=&space;\frac{\partial^2\Psi}{\partial\boldsymbol{\varepsilon&space;}^2}" title="\overset{4}{C} = \frac{\partial^2\Psi}{\partial\boldsymbol{\varepsilon }^2}" /></a>
;
<a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{E}&space;=&space;\frac{\partial^2\Psi}{\partial\boldsymbol{\varepsilon&space;}&space;\partial\varphi&space;}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{E}&space;=&space;\frac{\partial^2\Psi}{\partial\boldsymbol{\varepsilon&space;}&space;\partial\varphi&space;}" title="\boldsymbol{E} = \frac{\partial^2\Psi}{\partial\boldsymbol{\varepsilon } \partial\varphi }" /></a>
;
<a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{F}&space;=&space;\frac{\partial^2\Psi}{\partial\varphi&space;\partial\boldsymbol{\varepsilon&space;}&space;}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{F}&space;=&space;\frac{\partial^2\Psi}{\partial\varphi&space;\partial\boldsymbol{\varepsilon&space;}&space;}" title="\boldsymbol{F} = \frac{\partial^2\Psi}{\partial\varphi \partial\boldsymbol{\varepsilon } }" /></a>
;
<a href="https://www.codecogs.com/eqnedit.php?latex=G&space;=&space;\frac{\partial^2\Psi}{\partial\varphi^2&space;}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?G&space;=&space;\frac{\partial^2\Psi}{\partial\varphi^2&space;}" title="G = \frac{\partial^2\Psi}{\partial\varphi^2 }" /></a>

## How to start
- First you should check whether the herin described concepts and examples fit your needs. The easiest way to do this is by looking at the Doxygen documentation linked above and the list of current features.
- If you're interested in the details and want to get a look under the hood of the Sacado_Wrapper, then please consider the examples 1, 2, 3, 6 and 7
- If you right away want to use the Sacado_Wrapper, then consider the examples 3b, 4 and 8 and follow the instructions on how to get the Wrapper into your own code.

## Instructions on how to "install" the Sacado_Wrapper
1. Download the file "Sacado_Wrapper.h" and place it into your working directory, e.g. where your main code such as Sacado_examples.cc is.
2. Include the .h file with your other headers in your code via '#include "Sacado_Wrapper.h"'.
3. Copy one of the examples into your code for testing.
4. Compile, run and play around with the code.

## Remarks
- Currently it is not possible to use "Vector<fad_double>" or similar Sacado data types. You would probably have to first create instances of "Vector" for the desired data type as described here https://www.dealii.org/current/doxygen/deal.II/Instantiations.html. To avoid this you can simply use "std::vector<fad_double>".

## ToDo:
- maybe add the init to zero for the constructors of the SW data types
- add a note on the "reset_its_deriv" and "reset_other_deriv" issue
- add note on initialisation of value AND derivative via evolution equations with zero values
- check whether we can temporarily add Sacado dofs to the variable to e.g. compute some intermediate derivative and then later on delete these derivatives.
- template the Sacado data type e.g. Sacado::Fad::DFad (most flexible), Sacado::Fad::SFad (most efficient)  [according to https://docs.trilinos.org/latest-release/packages/sacado/doc/html/index.html]
- instead of using Dofs SymTensor+double+double+double, check using SymTensor+std::vector<double>

- find a more suitable name
- enable LaTeX equations in the documentation hosted via GitHub
- Rework the design and structure of the Wrapper
- implement more data types (e.g. Vector, nonsym tensor) and compatibility with more combinations (e.g. vector-vector, vector-double).
- add some text that: Every variable computed from the variables set as dofs contains the derivatives as shown in the figure. Hence, you can compute the tangents for every variable.
- add note on the efficiency/computation time
- update links in the documentation
- add link to https://arxiv.org/pdf/1811.05031.pdf
- remove todo for factor 0.5 in the beginning

## Test cases:
- ..
