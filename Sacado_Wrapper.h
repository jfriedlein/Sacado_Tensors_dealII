
// @section includes Include Files
// The data type SymmetricTensor and some related operations, such as trace, symmetrize, deviator, ... for tensor calculus
#include <deal.II/base/symmetric_tensor.h>

// @todo Check whether the following three headers are needed at all
#include <iostream>
#include <fstream>
#include <cmath>

// Sacado (from Trilinos, data types, operations, ...)
#include <Sacado.hpp>

using namespace dealii;

using fad_double = Sacado::Fad::DFad<double>;	// this data type now represents a double, but also contains the derivative of this variable with respect to the defined dofs (set via command *.diff(*))

namespace Sacado_Wrapper
{
	template <int dim>
	class SymTensor: public SymmetricTensor<2,dim, fad_double>
	{
	public:
		SymTensor( )
		{
			std::pair<unsigned int, unsigned int> tmp_pair;

			switch ( dim )
			{
			case 2:
		        tmp_pair.first=0; tmp_pair.second=0;
		        std_map_indicies[0] = tmp_pair;
		        
		        tmp_pair.first=0; tmp_pair.second=1;
		        std_map_indicies[1] = tmp_pair;
		        
		        tmp_pair.first=1; tmp_pair.second=1;
		        std_map_indicies[2] = tmp_pair; 
				break;
			case 3:
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
				break;
			}
		}
		
		// To simplify the access to the dofs we define a map that relate the components of our strain tensor to the dof-nbr
	    std::map<unsigned int,std::pair<unsigned int,unsigned int>> std_map_indicies;
		
	    void init ( SymmetricTensor<2,dim> &tensor_double );
	    
		void set_dofs();
		
		void get_tangent (SymmetricTensor<4,dim> &Tangent, SymmetricTensor<2,dim, fad_double> &sigma);
	};
	
	
	/*
	 * Initialization of the SymTensor data type with the values from a normal double SymmetricTensor
	 */
	template<int dim>
	void SymTensor<dim>::init( SymmetricTensor<2,dim> &tensor_double )
	{
		 for ( unsigned int i=0; i<dim; ++i)
			for ( unsigned int j=0; j<dim; ++j )
				(*this)[i][j] = tensor_double[i][j];
	}
	
	
	/*
	 * Set the dofs as the components of the SymTensor data type
	 */
	template<int dim>
	void SymTensor<dim>::set_dofs()
	{
		// Instead of calling the *.diff(*) on the components one-by-one we could also use the following for-loop, so
		// we also use the map to set the dofs
		for ( unsigned int x=0; x<((dim==2)?3:6); ++x )
		{
			unsigned int i=std_map_indicies[x].first;
			unsigned int j=std_map_indicies[x].second;
			(*this)[i][j].diff(x,((dim==2)?3:6));
		}
	}
	
	
	/*
	 * Assemble the tangent based on the derivatives saved in the argument sigma
	 * @param Tangent The desired tangent that will be computed as d_sigma/d_this, where "this" could be your strain tensor
	 * @param sigma The input argument used to extract the derivatives
	 */
	template<int dim>
	void SymTensor<dim>::get_tangent( SymmetricTensor<4,dim> &Tangent, SymmetricTensor<2,dim, fad_double> &sigma )
	{
		for ( unsigned int i=0; i<dim; ++i)
			for ( unsigned int j=0; j<dim; ++j )
			{
				double *derivs = &sigma[i][j].fastAccessDx(0); // Access the derivatives of the (i,j)-th component of \a sigma

	            // We loop over all the dofs. To be able to use this independent of the chosen dimension \a dim, we use a ternary operator
	            // to decide whether we have to loop over 6 derivatives or just 3.
	            for(unsigned int x=0;x<((dim==2)?3:6);++x)
	            {
	                unsigned int k=std_map_indicies[x].first;
	                unsigned int l=std_map_indicies[x].second;

	                if(k!=l)/*Compare to Voigt notation since only SymmetricTensor instead of Tensor*/
	                {
	                    Tangent[i][j][k][l] = 0.5*derivs[x];
	                    Tangent[i][j][l][k] = 0.5*derivs[x];
	                }
	                else
	                	Tangent[i][j][k][l] = derivs[x];
	            }            
			}
	}
	
	
}