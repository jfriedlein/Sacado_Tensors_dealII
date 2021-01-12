/* ---------------------------------------------------------------------
 *
 * The first coding assignment to get familiar with tensor calculus related
 * deal.II class templates. This includes:
 * 
 * - Tensor<1,dim>
 * - Tensor<2,dim>
 * - Tensor<4,dim>
 * - Vector<double>
 * - FullMatrix<double>
 * 
 * dim is a template variable which allows to vary e.g. between the two- and
 * three dimensional case. As described in the brief repetition of the 
 * essentials in C++, the respective tensor class templates allow different
 * dimensions as input, i.e. dim=1,dim=2,dim=3 for each respective rank
 *
 * ---------------------------------------------------------------------
 */


#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/table_indices.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/timer.h>

#include <iostream>
#include <vector>

// Sacado
#include <Sacado.hpp>
#include "Sacado_Wrapper.h"

// defining a data type for the Sacada variables (simply used the standard types from the step-33 tutorial's introduction)
using fad_double = Sacado::Fad::DFad<double>;

using namespace dealii;



//----------------------------------------------------

int main ()
{
    /* For now the three dimensional case is considered.
     * This information is essential for deal.II since
     * is suited for arbitrary dimensions due to its
     * template character. Use this variable for all
     * deal.II templates that will be created in the
     * sequent
     * 
     * It is always helpful to read the manual and see
     * how functions are implemented, i.e. what is the
     * return value, how is the function called or what
     * is the input.
     * 
     * Further some functions may be declared 
     * "DEPRECATED" which means they are still usable
     * but will be removed in future releases of the 
     * library -> not recommended to use those
     */
    
    
    const int dim=3;
    //------------------------------------------------
    /* START IMPLEMENTATION HERE                    */
    //------------------------------------------------
    
    
    //------------------------------------------------
    //          EX - 1
    /* Create two tensors of rank one and name them
     * u and v respectively and print them to the screen.
     * Therefore consider the available documentation
     * and manual on the deal.II website
     */
    
    {
	Sacado_Wrapper::SymTensor<dim> right_cauchy_green_sym;

    right_cauchy_green_sym[0][0].val() = 1.;
    right_cauchy_green_sym[1][1].val() = 2.;
    right_cauchy_green_sym[2][2].val() = 3.;
    right_cauchy_green_sym[0][1].val() = 0.4;
    right_cauchy_green_sym[0][2].val() = 0.5;
    right_cauchy_green_sym[1][2].val() = 0.6;

    right_cauchy_green_sym.set_dofs();
    
 	 std::vector<fad_double> eigenvalues (3);
 	 std::vector< Tensor<1,3,fad_double> > eigenvector (3);
 	 std::vector< SymmetricTensor<2,3,fad_double> > eigenbasis (3);
    
	for (unsigned int i = 0; i < 3; ++i)
	{
		eigenvalues[i] = eigenvectors(right_cauchy_green_sym)[i].first;
		eigenvector[i] = eigenvectors(right_cauchy_green_sym)[i].second;
	}
	
	for (unsigned int a = 0; a < 3; ++a)
		for (unsigned int i = 0; i < 3; ++i)
			for (unsigned int j = 0; j < 3; ++j)
				eigenbasis[a][i][j] += eigenvector[a][i] * eigenvector[a][j];
	
	SymmetricTensor<2,3,fad_double> hencky_strain_3D; // full 3D components
	 for (unsigned int a = 0; a < 3; ++a)
		for (unsigned int i = 0; i < 3; ++i)
			for (unsigned int j = 0; j < 3; ++j)
				hencky_strain_3D[i][j] += log(eigenvalues[a]) * eigenbasis[a][i][j];
	
	std::cout << hencky_strain_3D << std::endl;
    }
    
    {
	Sacado_Wrapper::SymTensor2<dim> right_cauchy_green_sym;

	SymmetricTensor<2,3> C;
	C[0][0] = 1.;
    C[1][1] = 2.;
    C[2][2] = 3.;
    C[0][1] = 0.4;
    C[0][2] = 0.5;
    C[1][2] = 0.6;

    right_cauchy_green_sym.init_set_dofs(C);
    
 	 std::vector< Sacado::Fad::DFad<DFadType> > eigenvalues (3);
 	 std::vector< Tensor<1,3,Sacado::Fad::DFad<DFadType> > > eigenvector (3);
 	 std::vector< SymmetricTensor<2,3,Sacado::Fad::DFad<DFadType> > > eigenbasis (3);
    
	for (unsigned int i = 0; i < 3; ++i)
	{
		eigenvalues[i] = eigenvectors(right_cauchy_green_sym)[i].first;
		eigenvector[i] = eigenvectors(right_cauchy_green_sym)[i].second;
	}
	
	for (unsigned int a = 0; a < 3; ++a)
		for (unsigned int i = 0; i < 3; ++i)
			for (unsigned int j = 0; j < 3; ++j)
				eigenbasis[a][i][j] += eigenvector[a][i] * eigenvector[a][j];
	
	SymmetricTensor<2,dim, Sacado::Fad::DFad<DFadType> > hencky_strain_3D;
	 for (unsigned int a = 0; a < 3; ++a)
		for (unsigned int i = 0; i < 3; ++i)
			for (unsigned int j = 0; j < 3; ++j)
				hencky_strain_3D[i][j] += log(eigenvalues[a]) * eigenbasis[a][i][j];
	
	 SymmetricTensor<4,dim> H_C;
	 right_cauchy_green_sym.get_tangent(H_C, hencky_strain_3D);
	 std::cout << "H_C=" << H_C << std::endl;
	
	 Tensor<6,dim> H_CC;
	 right_cauchy_green_sym.get_curvature(H_CC, hencky_strain_3D);
	 std::cout << "H_CC=" << H_CC << std::endl;
    }
}

//----------------------------------------------------