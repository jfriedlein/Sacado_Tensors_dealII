/**

@{

\brief Example usage of Sacado for tensor calculus

\author jfriedlein

@tableofcontents

@mainpage Trilinos::Sacado example documentation

@section intro Introduction
The way we see and use Sacado here is as follows: \n
If you usually compute the following equation \n
\f[	c = a + b \f]
for instance with data types double as \n
\f[ 1.0 + 1.0 \rightarrow 2.0 \f]
your results is just a double number \f$ c \f$ that contains the value \f$ 2 \f$. \n

Using Sacado, on the other hand, the variable \f$ c_{fad} \f$ is now of, for example, data type
@code
	Sacado::Fad::DFad<double> c_fad;
@endcode
As a result, \f$ c_{fad} \f$ now contains not just the number \f$ 2 \f$, but also all the derivatives of \f$ c_{fad} \f$
with respect to the previously defined degrees of freedom (set via command *.diff(*)). \n
The following figure tries to visualize this for first derivatives:
\image html Sacado_data-type.png
and for second derivatives (requires different data type):
\image html Sacado_data-type_1_2_derivatives.png
Every variable that was declared as such a data type now contains besides the actual value \f$ c \f$, also its derivatives.

@subsection subsec_overview Overview
This overview shall give you a first impression what to expect from each of the examples.
The background/basics group gives you the promised look under the hood.
Whereas the application group shows you how you can use the Sacado_Wrapper to quickly compute tangents.
\image html overview_of_examples.png

If you right away want to use Sacado, then you might skip the first examples and jump to \ref Ex3B "example 3B".
There we show how to use the "Sacado_Wrapper" that does everything from \ref Ex2 "example 2" and \ref Ex3 "example 3" in just a view lines of code. This does not mean that
the here shown approach is the fastest or most efficient, it is just simple and easy to use.

Furthermore, if you, for instance, compute problems with two-fields (e.g. displacement and scalar damage) and you need
tangents with respect to both a tensor (e.g. strain tensor) and a scalar (e.g. damage variable), you can use the Sacado_Wrapper as shown in \ref Ex4 "example 4".

@subsection subsec_more_basics Some more basics
One can access the double value of \f$ c_{fad} \f$ with the Sacado command *.val():
@code
	double c_value = c_fad.val();
@endcode
The derivatives of \f$ c_{fad} \f$ can be accessed with the command *.dx():
@code
	double d_c_d_a = c_fad.dx(0);
	double c_c_d_b = c_fad.dx(1);
@endcode
The arguments of \a dx, namely 0 and 1 are the numbers corresponding to the dof that belong to \a a and \a b. More details on how to set this up
and use it, are given in \ref Ex1 "example 1".

@subsection subsec_resources Some resources/links
You can use Sacado to compute general derivatives of functions (with or without tensors) with respect to variables (double, Tensors, ...).
@todo Link is broken

- Basics on Sacado, automatic differentation and coding examples (important remarks, e.g. on point of non-differentiability): \n
	https://software.sandia.gov/SESS/past_seminars/111307_Phipps.pdf
- Basics and available autodiff libraries provided by deal.ii: \n
	https://www.dealii.org/current/doxygen/deal.II/group__auto__symb__diff.html
- Template-based generic programming: \n
	downloads.hindawi.com/journals/sp/2012/202071.pdf
- Tutorial 33 by deal.ii: introduction to Sacado and implementation to assemble the residuum and compute its derivative \n
	https://www.dealii.org/current/doxygen/deal.II/step_33.html

The here shown examples shall solely show how Sacado can be applied and give some background and a look under the hood.
The code is neither elegant nor efficient, but it works. A more user-friendly version is provided by means of the "Sacado_Wrapper". \n
@todo Check if factor 0.5 is also necessary for d_sigma / d_phi

@note This documentation and code only protocol my first steps with Sacado. They are not guaranteed to be correct neither are they verified.
Any comments, criticism, corrections, feedback, improvements, ... are appreciated and very well welcomed.

@section code The commented program

\code

\endcode

@section END The End

Hosted via GitHub according to https://goseeky.wordpress.com/2017/07/22/documentation-101-doxygen-with-github-pages/ \n
Design of the documentation inspired by the deal.ii tutorial programs and created as shown here: https://github.com/jfriedlein/Custom_Doxygen_Documentation

@}
*/
