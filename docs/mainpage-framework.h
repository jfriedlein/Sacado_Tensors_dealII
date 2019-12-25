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

Using Sacado, on the other hand, the variable \f$ c \f$ is now of, for example, data type @code Sacado::Fad::DFad<double> @endcode.
As a result, \f$ c \f$ now contains not just the number \f$ 2 \f$, but also all the derivatives of \f$ c \f$
with respect to the previously defined degrees of freedom (set via command *.diff(*)). \n
The following figure tries to visualize this:
\image html Sacado_data-type.png
@todo redo this figure

If you right away want to use Sacado, then you might skip the first examples and jump to Ex3B.
There we show how to use the "Sacado_Wrapper" that does everything from Ex2 and Ex3 in just a view lines of code.

Furthermore, if you, for instance, compute problems with two-fields (e.g. displacement and scalar damage) and you need
tangents with respect to both a tensor (e.g. strain tensor) and a scalar (e.g. damage variable), you can use the Sacado_Wrapper as shown in Ex4.

Some resources/links: \n
@todo link the Sacado and DII pages

The here shown examples shall solely show how Sacado can be applied and give some background and a look under the hood.
The code is neither elegant nor efficient, but it works. A more user-friendly version is provided by means of the "Sacado_Wrapper". \n
@todo add list of files and an overview
@todo explain how to use the Wrapper (download the file Sacado_Wrapper.h, #include, ...)

@note This documentation and code only protocol my first steps with Sacado. They are not guaranteed to be correct neither are they verified.
Any comments, criticism, corrections, feedback, improvements, ... are very well welcomed.

@section code The commented program

\code

\endcode

@section END The End

Hosted via GitHub according to https://goseeky.wordpress.com/2017/07/22/documentation-101-doxygen-with-github-pages/

@}
*/
