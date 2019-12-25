
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

		unsigned int start_index = 0;
		static const unsigned int n_dofs = ((dim==2)?3:6); // ToDo: Find a way to use eps.n_independent_components(); that gives ((dim==2)?3:6) for a SymmetricTensor

	    void init ( SymmetricTensor<2,dim> &tensor_double );

		void set_dofs( unsigned int nbr_total_dofs=n_dofs );

		void get_tangent (SymmetricTensor<4,dim> &Tangent, SymmetricTensor<2,dim, fad_double> &sigma);

		void get_tangent( SymmetricTensor<2,dim> &Tangent, fad_double &argument );
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
	void SymTensor<dim>::set_dofs( unsigned int nbr_total_dofs )
	{
		// Instead of calling the *.diff(*) on the components one-by-one we could also use the following for-loop, so
		// we also use the map to set the dofs
		for ( unsigned int x=start_index; x<(start_index+n_dofs); ++x )
		{
			unsigned int i=std_map_indicies[x].first;
			unsigned int j=std_map_indicies[x].second;
			(*this)[i][j].diff(x,nbr_total_dofs);
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
	            for(unsigned int x=start_index;x<(start_index+n_dofs);++x)
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

	template<int dim>
	void SymTensor<dim>::get_tangent( SymmetricTensor<2,dim> &Tangent, fad_double &argument )
	{
		double *derivs = &argument.fastAccessDx(0); // Access derivatives
		for(unsigned int x=start_index;x<(start_index+n_dofs);++x)
		{
			unsigned int k=std_map_indicies[x].first;
			unsigned int l=std_map_indicies[x].second;

			if(k!=l)/*Compare to Voigt notation since only SymmetricTensor instead of Tensor*/
			{
				Tangent[k][l] = 0.5*derivs[x];
				Tangent[l][k] = 0.5*derivs[x];
			}
			else
				Tangent[k][l] = derivs[x];
		}
	}



	//###########################################################################################################//


	/*
	 * The same as the class above, just for data type fad_double
	 */
	template<int dim>
	class SW_double: public fad_double
	{
	public:
//		SW_double( double &sdf );
//		:
//		(*this)(sdf)
//		{
//		}

		// Assignment operator: \n
		// According to https://stackoverflow.com/questions/31029007/c-assignment-operator-in-derived-class
		// the operator= is not derived from the base class fad_double and needs to be set explicitly:
		SW_double & operator=(double double_init) { fad_double::operator =( double_init ) ;return *this;}

		static const unsigned int n_dofs = 1; // a double is just a single number, so it represents a single dof
		unsigned int start_index = 0;

		void init ( const double &double_init );

		void set_dofs ( unsigned int nbr_total_dofs=n_dofs );

		void get_tangent (SymmetricTensor<2,dim> &Tangent, SymmetricTensor<2,dim, fad_double> &sigma);

		void get_tangent ( double &Tangent, fad_double &argument );
	};

	template<int dim>
	void SW_double<dim>::init ( const double &double_init )
	{
		(*this) = double_init;
	}

	template<int dim>
	void SW_double<dim>::set_dofs ( unsigned int nbr_total_dofs )
	{
		(*this).diff( this->start_index, nbr_total_dofs );
	}

	template<int dim>
	void SW_double<dim>::get_tangent (SymmetricTensor<2,dim> &Tangent, SymmetricTensor<2,dim, fad_double> &sigma)
	{
		// reassemble the tangent as a SECOND order tensor
		 for ( unsigned int i=0; i<dim; ++i)
			for ( unsigned int j=0; j<dim; ++j )
			{
				double *derivs = &sigma[i][j].fastAccessDx(0); // Access derivatives
				Tangent[i][j] = derivs[ this->start_index ]; // ToDo: check whether the 0.5* is necessary here too
			}
	}

	template<int dim>
	void SW_double<dim>::get_tangent ( double &Tangent, fad_double &argument )
	{
		 double *derivs = &argument.fastAccessDx(0);
		 Tangent = derivs[ this->start_index ];
	}


	//###########################################################################################################//


	template<int dim>
	class DoFs_summary
	{
	public:
		void set_dofs(SymTensor<dim> &eps, SW_double<dim> &double_arg );
	};

	template<int dim>
	void DoFs_summary<dim>::set_dofs(SymTensor<dim> &eps, SW_double<dim> &double_arg)
	{
		const unsigned int nbr_total_dofs = eps.n_dofs + double_arg.n_dofs;

		eps.start_index = 0;
		double_arg.start_index = eps.n_dofs;

		eps.set_dofs( nbr_total_dofs );
		double_arg.set_dofs( nbr_total_dofs );
	}
}







