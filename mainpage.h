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

Some resources/links: \n
@todo link the Sacado and DII pages

The here shown examples shall solely show how Sacado can be applied and give some background and a look under the hood.
The code is neither elegant nor efficient, but it works. A more user-friendly version is available (check-out XXX) \n
@todo add the link to the wrapper

@note This documentation and code only protocol my first steps with Sacado. They are not guaranteed to be correct neither are they verified.
Any comments, criticism, corrections, feedback, improvements, ... are very well welcomed.

@section code The commented program

\code
 
/*
 * Author: jfriedlein, 2019
 * 		dsoldner, 2019
 */
 
\endcode
@section includes Include Files
The data type SymmetricTensor and some related operations, such as trace, symmetrize, deviator, ... for tensor calculus
\code
#include <deal.II/base/symmetric_tensor.h>
 
\endcode
C++ headers (some basics, standard stuff)
\code
#include <iostream>
#include <fstream>
#include <cmath>
 
\endcode
Sacado (from Trilinos, data types, operations, ...)
\code
#include <Sacado.hpp>
 
#include "Sacado_Wrapper.h"
 
\endcode
Those headers are related to data types and autodiff, but don't seem to be needed
\code
//#  include <deal.II/base/numbers.h>
//#  include <deal.II/differentiation/ad/ad_number_traits.h>
//#  include <deal.II/differentiation/ad/sacado_number_types.h>
 
\endcode
According to the basics of deal.ii-programming (see dealii.org and https://www.dealii.org/current/doxygen/deal.II/step_1.html for a start)
\code
using namespace dealii;
 
\endcode
Defining a data type for the Sacado variables (here we simply used the standard types from the deal.ii step-33 tutorial's introduction)
\code
using fad_double = Sacado::Fad::DFad<double>;	// this data type now represents a double, but also contains the derivative of this variable with respect to the defined dofs (set via command *.diff(*))
 
 
\endcode
@section Ex1 1. example: simple scalar equation
1. example: simple scalar equation from deal.ii-tutorial step-33 (see the introduction there to get a first impression, https://www.dealii.org/current/doxygen/deal.II/step_33.html)
@todo clean up the documentation of the classes

\code
void sacado_test_scalar ()
{
	std::cout << "Scalar Test:" << std::endl;
\endcode
define the variables used in the computation (inputs: a, b; output: c; auxiliaries: *) as the Sacado-data type
\code
	 fad_double a,b,c;
\endcode
initialize the input variables a and b; This (a,b) = (1,2) will be the point where the derivatives are computed.
Compare: y=x² -> (dy/dx)(\@x=1) = 2. We can only compute the derivative numerically at a certain point.
\code
	 a = 1;
	 b = 2;
 
	a.diff(0,2);  // Set a to be dof 0, in a 2-dof system.
	b.diff(1,2);  // Set b to be dof 1, in a 2-dof system.
\endcode
Our equation here is very simply. But you can use nested equations and many standard mathematical operations, such as sqrt, pow, sin, ...
\code
	c = 2*a + std::cos(a*b);
	double *derivs = &c.fastAccessDx(0); // Access the derivatives of
\endcode
Output the derivatives of c with respect to the two above defined degrees of freedom (dof)
\code
	std::cout << "Derivatives at the point (" << a << "," << b << ")" << std::endl;
	std::cout << "dc/da = " << derivs[0] << ", dc/db=" << derivs[1] << std::endl;
}
 
 
\endcode
@section Ex2 2. example: Preparation for the use of Sacado with tensors
Here we want to introduce tensors for the first time. Hence, we limit ourselves to a trivial equation relating the strain tensor \a eps
with dim x dim components with the stress tensor \a sigma. Both here used tensors are symmetric, hence we use the SymmetricTensor class and
have to keep some details in mind (see below factor 0.5 related to Voigt-Notation).
\code
/*
 * 2. example: use of tensors
 */
void sacado_test_2 ()
{
	std::cout << "Tensor Test 2:" << std::endl;
 
\endcode
First we set the dimension \a dim: 2D->dim=2; 3D->dim=3 \n This defines the "size" of the tensors and the number of dofs. Ex2 only works in 3D, whereas the following Ex3 is set up dimension-independent.
\code
	const unsigned int dim = 3;
 
\endcode
Declare our input, auxiliary and output variables as SymmetricTensors consisting of fad_doubles (instead of the standard SymmetricTensor out of doubles)
\code
	SymmetricTensor<2,dim, fad_double> sigma, eps;
 
\endcode
Init the strain tensor (the point at which the derivative shall be computed)
\code
	eps[0][0] = 1;
	eps[1][1] = 2;
	eps[2][2] = 3;
	eps[0][1] = 4;
	eps[0][2] = 5;
	eps[1][2] = 6;
 
\endcode
Now we declare the dofs. The derivative to a tensor requires all components, therefore we set the components of the strain tensor here one by one as the dofs.
Because our tensors are symmetric, we only need 6 components in 3D instead of 9 for a full second order tensor
\code
	eps[0][0].diff(0,6);
	eps[1][1].diff(1,6);
	eps[2][2].diff(2,6);
	eps[0][1].diff(3,6);
	eps[0][2].diff(4,6);
	eps[1][2].diff(5,6);
 
\endcode
The equation describing the stresses (here just a simple test case)
\code
	sigma = eps;
 
\endcode
Let's output the computed stress tensor.
\code
	std::cout << sigma << std::endl;
\endcode
The resulting values of \a sigma are fairly boring, due to our simple equation. It is the additional output generated by
this, that is interesting here: \n
output: \n
1 [ 1 0 0 0 0 0 ] 4 [ 0 0 0 1 0 0 ] 5 [ 0 0 0 0 1 0 ] 4 [ 0 0 0 1 0 0 ] 2 [ 0 1 0 0 0 0 ] 6 [ 0 0 0 0 0 1 ] 5 [ 0 0 0 0 1 0 ] 6 [ 0 0 0 0 0 1 ] 3 [ 0 0 1 0 0 0 ] \n
The numbers 1, 4, 5, 4, ... are the entries in the stress tensor \a sigma. In square brackets we see the derivatives of sigma with respect to all the dofs set previously
given in the order we defined them above. Meaning: The first entry in the square brackets corresponds to the 0-th dof set by
@code eps[0][0].diff(0,6); @endcode referring to the component (0,0) in the strain tensor \a eps.
 
Computing the derivatives for certain components of the resulting tangent modulus: \n
We now access these lists of derivatives (output above in square brackets) for one component of the stress tensor \a sigma at a time.
\code
	{
\endcode
Access the derivatives corresponding to the component (0,0) of the stress tensor \a sigma
\code
		double *derivs = &sigma[0][0].fastAccessDx(0);
\endcode
The following output will show us the same derivatives that we already saw above, just formatted differently \n
output: d_sigma[0][0]/d_eps = 1 , 0 , 0 , 0 , 0 , 0 ,
\code
		std::cout << "d_sigma[0][0]/d_eps = ";
		for ( unsigned int i=0; i<6; ++i)
			std::cout << derivs[i] << " , ";
		std::cout << std::endl;
	}
	{
\endcode
Access the derivatives corresponding to the component (1,2) of the stress tensor \a sigma
\code
		double *derivs = &sigma[1][2].fastAccessDx(0);
\endcode
output: d_sigma[1][2]/d_eps = 0 , 0 , 0 , 0 , 0 , 1 ,
\code
		std::cout << "d_sigma[1][2]/d_eps = ";
		for ( unsigned int i=0; i<6; ++i)
			std::cout << derivs[i] << " , ";
		std::cout << std::endl;
	}
}
 
 
\endcode
@section Ex3 3. example: Using a slightly more complicated stress equation
\code
void sacado_test_3 ()
{
	std::cout << "Tensor Test 3:" << std::endl;
 
	const unsigned int dim = 3;
 
\endcode
Here we also define some constant, for instance the bulk modulus \a kappa and the second Lamè parameter \a mu.
We now also define one of our constants as fad_double. By doing this we can use the normal multiplication (see below).
\code
	double kappa_param = 5;
	fad_double kappa (kappa_param);
\endcode
The second constant remains as a double just to show the difference.
\code
	double mu = 2;
 
	SymmetricTensor<2,dim, fad_double> sigma, eps;
 
\endcode
To simplify the access to the dofs we define a map that relate the components of our strain tensor to the dof-nbr
\code
    std::map<unsigned int,std::pair<unsigned int,unsigned int>> std_map_indicies;
 
\endcode
The point at which the derivative shall be computed: \n
As mentioned previously, we will implement this example for 2D and 3D, hence we once have to set up a strain tensor
and the derivatives for 3D with 6 independent components ...
\code
	if(dim==3)
	{
		eps[0][0] = 1;
		eps[1][1] = 2;
		eps[2][2] = 3;
 
		eps[0][1] = 4;
		eps[0][2] = 5;
		eps[1][2] = 6;
 
 
		eps[0][0].diff(0,6);
		eps[0][1].diff(1,6);
		eps[0][2].diff(2,6);
		eps[1][1].diff(3,6);
		eps[1][2].diff(4,6);
		eps[2][2].diff(5,6);
 
\endcode
By using the map and the following pairs, we have to set up the relation between strain components and dofs only once
and can use the map to access the entries of the list later, without possibly mixing up indices and creating errors
\code
		std::pair<unsigned int, unsigned int> tmp_pair;
		tmp_pair.first=0; tmp_pair.second=0;
		std_map_indicies[0] = tmp_pair;
 
		tmp_pair.first=0; tmp_pair.second=1;
		std_map_indicies[1] = tmp_pair;
 
		tmp_pair.first=0; tmp_pair.second=2;
		std_map_indicies[2] = tmp_pair;
 
		tmp_pair.first=1; tmp_pair.second=1;
		std_map_indicies[3] = tmp_pair;
 
		tmp_pair.first=1; tmp_pair.second=2;
		std_map_indicies[4] = tmp_pair;
 
		tmp_pair.first=2; tmp_pair.second=2;
		std_map_indicies[5] = tmp_pair;
	}
\endcode
... and once for 2D with just 3 independent components.
\code
	else if(dim==2)
	{
		eps[0][0] = 1;
		eps[1][1] = 2;
 
		eps[0][1] = 4;
 
 
		eps[0][0].diff(0,3);
		eps[0][1].diff(1,3);
		eps[1][1].diff(2,3);
 
        std::pair<unsigned int, unsigned int> tmp_pair;
        tmp_pair.first=0; tmp_pair.second=0;
        std_map_indicies[0] = tmp_pair;
 
        tmp_pair.first=0; tmp_pair.second=1;
        std_map_indicies[1] = tmp_pair;
 
        tmp_pair.first=1; tmp_pair.second=1;
        std_map_indicies[2] = tmp_pair;        
	}
	else
	{
		throw std::runtime_error("only dim==2 or dim==3 allowed");
	}
 
\endcode
Instead of calling the *.diff(*) on the components one-by-one we could also use the following for-loop, so
we also use the map to set the dofs
@code
for ( unsigned int x=0; x<((dim==2)?3:6); ++x )
{
	unsigned int i=std_map_indicies[x].first;
	unsigned int j=std_map_indicies[x].second;
	eps[i][j].diff(x,((dim==2)?3:6));
}
@endcode
 
For our slightly more complicated stress equation we need the unit and deviatoric tensors.
We can simply define them by writing the values of the already existing deal.ii functions into newly
defined SymmetricTensors build from fad_doubles.
\code
	SymmetricTensor<2,dim, fad_double> stdTensor_I (( unit_symmetric_tensor<dim,fad_double>()) );
	SymmetricTensor<4,dim, fad_double> stdTensor_Idev ( (deviator_tensor<dim,fad_double>()) );
 
\endcode
With everything set and defined, we can compute our stress \a sigma according to:
\f[ \sigma = \kappa \cdot trace(\varepsilon) \cdot \boldsymbol{I} + 2 \cdot \mu \cdot \varepsilon^{dev} \f]
Here you can see that we can directly multiply the constant and the tensors when kappa is also declared as fad_double
\code
	sigma = kappa * (trace(eps) *  stdTensor_I);
\endcode
We didn't do the same for mu to once again emphasize the difference between constants as double and as fad_double. \n
The remaining code uses a normal double constant.
\code
    SymmetricTensor<2,dim,fad_double> tmp = deviator<dim,fad_double>(symmetrize<dim,fad_double>(eps)); tmp*=(mu*2);
    sigma +=  tmp;
\endcode
The fairly cumbersome computation is caused by the way the operators are set up for tensors out of fad_doubles.
 
\code
	std::cout << "sigma=" << sigma << std::endl;
 
\endcode
Now we want to actually build our tangent modulus called \a C_Sacado that contains all the derivatives and relates
the stress tensor with the strain tensor. \n
The fourth-order tensor \a C_Sacado is our final goal, we don't have to compute anything that is related to Sacado with
this tensor, so we can finally return to our standard SymmetricTensor out of doubles. The latter is necessary to use
the tangent in the actual FE code.
\code
	SymmetricTensor<4,dim> C_Sacado;
 
\endcode
As in Ex2 we access the components of the stress tensor one by one. In order to capture all of them we sum over the
components i and j of the stress tensor.
\code
	for ( unsigned int i=0; i<dim; ++i)
		for ( unsigned int j=0; j<dim; ++j )
		{
			double *derivs = &sigma[i][j].fastAccessDx(0); // Access the derivatives of the (i,j)-th component of \a sigma
 
\endcode
To visually ensure that every stress component has in fact all 6 derivatives for 3D or 3 for 2D, we output the size:
\code
            std::cout<<"size: "<<sigma[i][j].size()<<std::endl;
 
\endcode
We loop over all the dofs. To be able to use this independent of the chosen dimension \a dim, we use a ternary operator
to decide whether we have to loop over 6 derivatives or just 3.
\code
            for(unsigned int x=0;x<((dim==2)?3:6);++x)
            {
                unsigned int k=std_map_indicies[x].first;
                unsigned int l=std_map_indicies[x].second;
 
                if(k!=l)/*Compare to Voigt notation since only SymmetricTensor instead of Tensor*/
                {
                    C_Sacado[i][j][k][l] = 0.5*derivs[x];
                    C_Sacado[i][j][l][k] = 0.5*derivs[x];
                }
                else
                    C_Sacado[i][j][k][l] = derivs[x];
            }            
 
		}
 
\endcode
After resembling the fourth-order tensor, we now have got our tangent saved in \a C_Sacado ready to be used
 
To ensure that Sacado works properly, we can compute the analytical tangent for comparison
\code
	double kappa_d = 5;
	double mu_d = 2;
\endcode
Our stress equation in this example is still simple enough to derive the tangent analytically by hand:
\f[ \overset{4}{C_{analy}} = \kappa \cdot \boldsymbol{I} \otimes \boldsymbol{I} + 2 \cdot \mu \cdot \overset{4}{I^{dev}} \f]
\code
	SymmetricTensor<4,dim> C_analy = kappa_d * outer_product(unit_symmetric_tensor<dim>(), unit_symmetric_tensor<dim>()) + 2* mu_d * deviator_tensor<dim>();
 
 
\endcode
We again define our strain tensor \a eps_d (*_d for standard double in contrast to fad_double)
\code
	SymmetricTensor<2,dim> eps_d;
 
	if(dim==3)
	{
		eps_d[0][0] = 1;
		eps_d[1][1] = 2;
		eps_d[2][2] = 3;
 
		eps_d[0][1] = 4;
		eps_d[0][2] = 5;
		eps_d[2][1] = 6;
 
	}
	else if(dim==2)
	{
		eps_d[0][0] = 1;
		eps_d[1][1] = 2;
 
		eps_d[1][0] = 4;
 
	}
	else
	{
		throw std::runtime_error("only dim==2 or dim==3 allowed");
	}
\endcode
@todo use boldsymbol for tensors

To output the stress tensor we first have to compute it. We do this here via
\f[ \sigma = \overset{4}{C_{analy}} : \varepsilon \f]
The output exactly matched the result obtained with Sacado.
@note Checking the Sacado stress tensor against an analytically computed or otherwise determined stress tensor is absolutely no way to check whether
the tangent computed via Sacado is correct. When we compute the stress tensor with Sacado and for example mix up a + and - sign, this might not matter
at all if the number that is added or subtracted is small. However, for the tangent this nasty sign can be very critical. Just keep in mind: the
tangent has 81 components and the stress tensor just 9, so how does one want to verify 81 variables by comparing 9?

\code
    std::cout << "sigma_analy: " << (C_analy*eps_d) << std::endl;
 
\endcode
That's the reason we compare all the entries in the Sacado and the analytical tensor one by one
\code
	for (unsigned int i=0; i<dim; ++i)
		for ( unsigned int j=0; j<dim; ++j)
			for ( unsigned int k=0; k<dim; ++k)
				for ( unsigned int l=0; l<dim; ++l)
					std::cout << "C_analy["<<i<<"]["<<j<<"]["<<k<<"]["<<l<<"] = " << C_analy[i][j][k][l] << " vs C_Sacado: " << C_Sacado[i][j][k][l] << std::endl;
 
 
\endcode
To simplify the comparison we compute a scalar error as the sum of the absolute differences of each component
\code
	double error_Sacado_vs_analy=0;
	for (unsigned int i=0; i<dim; ++i)
		for ( unsigned int j=0; j<dim; ++j)
			for ( unsigned int k=0; k<dim; ++k)
				for ( unsigned int l=0; l<dim; ++l)
					error_Sacado_vs_analy += std::fabs(C_Sacado[i][j][k][l] - C_analy[i][j][k][l]);
 
 
\endcode
As desired: The numerical error is zero (0 in double precision) and the tensor components are equal
\code
	std::cout << "numerical error: " << error_Sacado_vs_analy << std::endl;
}
 
 
\endcode
@section Ex3B 3B. Example: Using the wrapper for Ex3
\code
void sacado_test_3B ()
{
	std::cout << "Tensor Test 3B:" << std::endl;
    const unsigned int dim=3;
 
\endcode
The following declarations are usually input arguments. So you receive the strain tensor and the constants out of doubles.
\code
    SymmetricTensor<2,dim> eps_d;
	eps_d[0][0] = 1;
	eps_d[1][1] = 2;
	eps_d[2][2] = 3;
 
	eps_d[0][1] = 4;
	eps_d[0][2] = 5;
	eps_d[1][2] = 6;
 
	double kappa = 5;
	double mu = 2;
 
\endcode
Now we start working with Sacado: \n
When we use the index notation to compute e.g. our stress we do not need to declare our constants (here kappa, mu) as
fad_double.
 
We declare our strain tensor as the special data type Sacado_Wrapper::SymTensor from the file "Sacado_Wrapper.h"
where this data type was derived from the SymmetricTensor<2,dim,fad_double>.
\code
	 Sacado_Wrapper::SymTensor<dim> eps;
 
\endcode
Next we initialize our Sacado strain tensor with the values of the inputed double strain tensor:
\code
	 eps.init(eps_d);
 
\endcode
We define all the entries in the symmetric tensor \a eps as the dofs. So we can later derive any variable
with respect to the strain tensor \a eps.
\code
	 eps.set_dofs();
 
\endcode
Now we declare our output and auxiliary variables as Sacado-Tensors.
\code
	 SymmetricTensor<2,dim,fad_double> sigma;
 
	 SymmetricTensor<2,dim, fad_double> stdTensor_I (( unit_symmetric_tensor<dim,fad_double>()) );
 
\endcode
Our stress equation is now computed in index notation to simplfiy the use of the constants and
especially the use of the \a deviator.
\code
	  for ( unsigned int i=0; i<dim; ++i)
		for ( unsigned int j=0; j<dim; ++j )
			sigma[i][j] = kappa * trace(eps) *  stdTensor_I[i][j] + 2. * mu * deviator(eps)[i][j];
 
\endcode
Finally we declare our desired tangent as the fourth order tensor \a C_Sacado and compute the tangent via
the command \a get_tangent.
\code
	 SymmetricTensor<4,dim> C_Sacado;
	 eps.get_tangent(C_Sacado, sigma);
 
\endcode
We could again compare the herein computed tangent with the analytical tangent from Ex2, but as before
the results are fairly boring, because Sacado hits the analytical tangent exactly --- no surprise for such
simple equations.
 
And that's it. By using the Sacado_wrapper we can achieve everything from Ex2 (besides the equations)
with just four lines of code namely:
- eps.init(eps_d);    // To initialize the Sacado strain tensor
- eps.set_dofs();     // To declare the components of eps as the dofs
- eps.get_tangent(*); // To get the tangent
\code
}
 
 
\endcode
@section Ex4 4. Example: Computing derivatives with respect to a tensor and a scalar
\code
void sacado_test_4 ()
{
	std::cout << "Tensor Test 4:" << std::endl;
	const unsigned int dim=3;
 
\endcode
The following declarations are usually input arguments. So you receive the strain tensor \q eps_d,
the damage variable \a phi and the constants \a kappa and \a mu out of doubles.
\code
	SymmetricTensor<2,dim> eps_d;
	eps_d[0][0] = 1;
	eps_d[1][1] = 2;
	eps_d[2][2] = 3;
 
	eps_d[0][1] = 4;
	eps_d[0][2] = 5;
	eps_d[1][2] = 6;
 
	double phi = 0.3;
 
	double kappa = 5;
	double mu = 2;
 
\endcode
We set up our strain tensor as in Ex3B.
\code
	Sacado_Wrapper::SymTensor<dim> eps;
	eps.init(eps_d);
	eps.set_dofs();
 
\endcode
In order to also compute derivatives with respect to the scalar \a phi, we add this scalar to our list
of derivatives:
@todo CONTINUE HERE
 
\code
}
 
 
\endcode
@section Ex5 5. Example: Using a vector-valued equation
\code
void sacado_test_5 ()
{
    const unsigned int dim=3;
	std::cout << "Tensor Test 5:" << std::endl;
    Tensor<1,dim,fad_double> c;
	fad_double a,b;
    unsigned int n_dofs=2;
	a = 1; b = 2;	// at the point (a,b) = (1,2)
	a.diff(0,2);  // Set a to be dof 0, in a 2-dof system.
	b.diff(1,2);  // Set b to be dof 1, in a 2-dof system.
\endcode
c is now a vector with three components
\code
	c[0] = 2*a+3*b;
    c[1] = 4*a+5*b;
    c[2] = 6*a+7*b;
 
\endcode
Access to the derivatives works as before.
\code
    for(unsigned int i=0;i<dim;++i)
    {
        const fad_double &derivs = c[i]; // Access derivatives
        for(unsigned int j=0;j<n_dofs;++j)
        {
            std::cout << "Derivatives at the point (" << a << "," << b << ") for "
            <<i<<"th component wrt "<<j<<"th direction "<< std::endl;
            std::cout << "dc_i/dxj = " << derivs.fastAccessDx(j) << std::endl;            
        }
    }
}
 
 
/*
 * The main function just calls all the examples and puts some space between the outputs.
 */
int main ()
{
	sacado_test_scalar ();
 
	std::cout << std::endl;
 
	sacado_test_2 ();
 
	std::cout << std::endl;
 
	sacado_test_3 ();
 
    std::cout << std::endl;
 
	sacado_test_3B ();
 
    std::cout << std::endl;
 
    sacado_test_4();
 
    std::cout << std::endl;
 
    sacado_test_5();
 
}
\endcode

@section END The End

@}
*/
