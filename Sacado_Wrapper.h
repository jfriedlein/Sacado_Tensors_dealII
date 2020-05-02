#ifndef Sacado_Wrapper_H
#define Sacado_Wrapper_H

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
typedef Sacado::Fad::DFad<double> DFadType;
/*
 * A function to set up the map that relates the indices (e.g. first entry in vector is the (0,0) component, whereas the second is the (0,1) component
 */
template<int dim>
void get_index_map( std::map<unsigned int,std::pair<unsigned int,unsigned int>> &std_map_indicies )
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


namespace Sacado_Wrapper
{
	template <int dim>
	class SymTensor: public SymmetricTensor<2,dim, fad_double>
	{
	public:
		SymTensor( )
		{
			get_index_map<dim>( std_map_indicies );
		}

		// To simplify the access to the dofs we define a map that relate the components of our strain tensor to the dof-nbr
	    std::map<unsigned int,std::pair<unsigned int,unsigned int>> std_map_indicies;

		unsigned int start_index = 0;
		static const unsigned int n_dofs = ((dim==2)?3:6); // ToDo: Find a way to use eps.n_independent_components(); that gives ((dim==2)?3:6) for a SymmetricTensor

	    void init ( SymmetricTensor<2,dim> &tensor_double );

		void set_dofs( unsigned int nbr_total_dofs=n_dofs );

		void get_tangent (SymmetricTensor<4,dim> &Tangent, SymmetricTensor<2,dim, fad_double> &sigma);

		void get_tangent( SymmetricTensor<2,dim> &Tangent, fad_double &argument );
		
		void get_values ( SymmetricTensor<2,dim> &tensor_double );

		SymmetricTensor<2,dim> get_value ();
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

        	Tangent[k][l] = derivs[x]; // ToDo: check whether the 0.5* is necessary here too

        	// Correct the off-diagonal terms by the factor of 0.5
             if(k!=l)/*Compare to Voigt notation since only SymmetricTensor instead of Tensor*/
            	Tangent[k][l] *= 0.5; // ToDo: check whether the 0.5* is necessary here too
		}
	}


	template<int dim>
	void SymTensor<dim>::get_values ( SymmetricTensor<2,dim> &tensor_double )
	{
		for ( unsigned int i=0; i<dim; ++i)
			for ( unsigned int j=0; j<dim; ++j )
				tensor_double[i][j] = ((*this)[i][j]).val();
	}

	template<int dim>
	SymmetricTensor<2,dim> SymTensor<dim>::get_value ()
	{
		SymmetricTensor<2,dim> tmp;
		(*this).get_values(tmp);
		return tmp;
	}


	//###########################################################################################################//

		 
	template <int dim>
	class SymTensor2: public SymmetricTensor<2,dim, Sacado::Fad::DFad<DFadType> >
	{
	public:
		SymTensor2( )
		{
			get_index_map<dim>( std_map_indicies );
		}

		// To simplify the access to the dofs we define a map that relate the components of our strain tensor to the dof-nbr
		std::map<unsigned int,std::pair<unsigned int,unsigned int>> std_map_indicies;

		unsigned int start_index = 0;
		static const unsigned int n_dofs = ((dim==2)?3:6); // ToDo: Find a way to use eps.n_independent_components(); that gives ((dim==2)?3:6) for a SymmetricTensor

		void init_set_dofs ( SymmetricTensor<2,dim> &tensor_double, const unsigned int nbr_total_dofs=n_dofs );

		void get_tangent (SymmetricTensor<4,dim> &Tangent, SymmetricTensor<2,dim, Sacado::Fad::DFad<DFadType> > &sigma);

		void get_tangent( SymmetricTensor<2,dim> &Tangent, Sacado::Fad::DFad<DFadType> &argument );
		
		void get_tangent( SymmetricTensor<2,dim, Sacado::Fad::DFad<DFadType> > &Tangent, Sacado::Fad::DFad<DFadType> &argument );
		
		void get_curvature( SymmetricTensor<4,dim> &Curvature, Sacado::Fad::DFad<DFadType> &argument );

//			
//		void get_values ( SymmetricTensor<2,dim> &tensor_double );
	};

	/*
	 * Initialization of the \a SymTensor2 data type with the values from a normal double \a SymmetricTensor
	 */
	template<int dim>
	void SymTensor2<dim>::init_set_dofs( SymmetricTensor<2,dim> &tensor_double, const unsigned int nbr_total_dofs )
	{
		 for ( unsigned int x=0; x<n_dofs; ++x )
		 {
		 	unsigned int i=std_map_indicies[x].first;
		 	unsigned int j=std_map_indicies[x].second;
			((*this)[i][j]).diff( x, nbr_total_dofs);	// set up the "inner" derivatives
			((*this)[i][j]).val() = fad_double(nbr_total_dofs, x, tensor_double[i][j]); // set up the "outer" derivatives
		 }
	}
	
	
	template<int dim>
	void SymTensor2<dim>::get_tangent( SymmetricTensor<2,dim> &Tangent, Sacado::Fad::DFad<DFadType> &argument )
	{
		for ( unsigned int x=0; x<n_dofs; ++x )
		 {
			unsigned int i=std_map_indicies[x].first;
			unsigned int j=std_map_indicies[x].second;
			if ( i!=j )
				Tangent[i][j] = 0.5 * argument.dx(x).val();
			else
				Tangent[i][j] = argument.dx(x).val();
		 }
	}
	
	template<int dim>
	void SymTensor2<dim>::get_tangent( SymmetricTensor<4,dim> &Tangent, SymmetricTensor<2,dim, Sacado::Fad::DFad<DFadType> > &argument )
	{
		 for(unsigned int x=0;x<n_dofs;++x)
			for(unsigned int y=0;y<n_dofs;++y)
			{
				const unsigned int i=std_map_indicies[y].first;
				const unsigned int j=std_map_indicies[y].second;
				const unsigned int k=std_map_indicies[x].first;
				const unsigned int l=std_map_indicies[x].second;

				if(k!=l)/*Compare to Voigt notation since only SymmetricTensor instead of Tensor*/
					Tangent[i][j][k][l] = 0.5 * argument[k][l].dx(y).val();		// ToDo: check this on paper (was more like a gut feeling)
				else
					Tangent[i][j][k][l] = argument[k][l].dx(y).val();
			}
	}


	/*
	 * Compute the tangent still containing all the second derivatives
	 */
	template<int dim>
	void SymTensor2<dim>::get_tangent( SymmetricTensor<2,dim, Sacado::Fad::DFad<DFadType> > &Tangent, Sacado::Fad::DFad<DFadType> &argument )
	{
		for ( unsigned int x=0; x<n_dofs; ++x )
		 {
			unsigned int i=std_map_indicies[x].first;
			unsigned int j=std_map_indicies[x].second;
			if ( i!=j )
				Tangent[i][j] = 0.5 * argument.dx(x);
			else
				Tangent[i][j] = argument.dx(x);
		 }
	}
	
	
	template<int dim>
	void SymTensor2<dim>::get_curvature( SymmetricTensor<4,dim> &Curvature, Sacado::Fad::DFad<DFadType> &argument )
	{
		 for(unsigned int x=0;x<n_dofs;++x)
			for(unsigned int y=0;y<n_dofs;++y)
			{
				const unsigned int i=std_map_indicies[y].first;
				const unsigned int j=std_map_indicies[y].second;
				const unsigned int k=std_map_indicies[x].first;
				const unsigned int l=std_map_indicies[x].second;

				double deriv = argument.dx(x).dx(y); // Access the derivatives of the (i,j)-th component of \a sigma

				if ( k!=l && i!=j )
					Curvature[i][j][k][l] = 0.25* deriv;
				else if(k!=l)/*Compare to Voigt notation since only SymmetricTensor instead of Tensor*/
				{
					Curvature[i][j][k][l] = 0.5*deriv;
					Curvature[i][j][l][k] = 0.5*deriv;
				}
				else
					Curvature[i][j][k][l] = deriv;
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
		SW_double & operator=(fad_double fad_assignment) { fad_double::operator =( fad_assignment ) ;return *this;}

		static const unsigned int n_dofs = 1; // a double is just a single number, so it represents a single dof
		unsigned int start_index = 0;

		void init ( const double &double_init );

		void set_dofs ( unsigned int nbr_total_dofs=n_dofs );

		void reset_its_deriv ( );
		void reset_other_deriv ( const SW_double<dim> &other_SW_double );

		void get_tangent (SymmetricTensor<2,dim> &Tangent, SymmetricTensor<2,dim, fad_double> &sigma);

		void get_tangent ( double &Tangent, fad_double &argument );
		
		void get_values ( double &return_double );
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
	void SW_double<dim>::reset_its_deriv ( )
	{
		double *derivs = &(this)->fastAccessDx(0);
		derivs[this->start_index] = 1.;
	}
	template<int dim>
	void SW_double<dim>::reset_other_deriv ( const SW_double<dim> &other_SW_double )
	{
		double *derivs = &(this)->fastAccessDx(0);
		derivs[other_SW_double.start_index] = 0.;
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

	template<int dim>
	void SW_double<dim>::get_values ( double &return_double )
	{
		return_double = (*this).val();
	}

	//###########################################################################################################//

	

	/*
	 * The same as the class above, just for data type fad_double
	 */
	template<int dim>
	class SW_double2: public Sacado::Fad::DFad<DFadType>
	{
	public:
//		SW_double2( double &sdf );
//		:
//		(*this)(sdf)
//		{
//		}

		// Assignment operator: \n
		// According to https://stackoverflow.com/questions/31029007/c-assignment-operator-in-derived-class
		// the operator= is not derived from the base class fad_double and needs to be set explicitly:
		SW_double2 & operator=(double double_init) { Sacado::Fad::DFad<DFadType>::operator =( double_init ) ;return *this;}
		SW_double2 & operator=(Sacado::Fad::DFad<DFadType> fad_assignment) { Sacado::Fad::DFad<DFadType>::operator =( fad_assignment ) ;return *this;}

		static const unsigned int n_dofs = 1; // a double is just a single number, so it represents a single dof
		unsigned int start_index = 0;

		void init_set_dofs ( const double &double_init, unsigned int nbr_total_dofs=n_dofs );

		void get_tangent (double &Tangent, Sacado::Fad::DFad<DFadType> &argument);
		void get_tangent (SymmetricTensor<2,dim> &Tangent, SymmetricTensor<2,dim, Sacado::Fad::DFad<DFadType> > &argument, SymTensor2<dim> &eps);


		void get_curvature (double &Curvature, Sacado::Fad::DFad<DFadType> &argument);
		
		void get_curvature (SymmetricTensor<2,dim> &Curvature, SymmetricTensor<2,dim, Sacado::Fad::DFad<DFadType> > &argument, SymTensor2<dim> &eps );
		void get_curvature (SymmetricTensor<2,dim> &Curvature, Sacado::Fad::DFad<DFadType> &argument,                          SymTensor2<dim> &eps );

//		void get_values ( double &return_double );
	};

	template<int dim>
	void SW_double2<dim>::init_set_dofs ( const double &double_init, unsigned int nbr_total_dofs )
	{
		(*this).diff( start_index , nbr_total_dofs );
		(*this).val() = fad_double(nbr_total_dofs, start_index, double_init);
	}

	
	template<int dim>
	void SW_double2<dim>::get_tangent (double &Tangent, Sacado::Fad::DFad<DFadType> &argument)
	{
		Tangent = argument.dx(this->start_index).val();
	}
	
	
	template<int dim>
	void SW_double2<dim>::get_tangent (SymmetricTensor<2,dim> &Tangent, SymmetricTensor<2,dim, Sacado::Fad::DFad<DFadType> > &argument, SymTensor2<dim> &eps)
	{
		for(unsigned int x=0;x<eps.n_dofs;++x)
		{
			const unsigned int i=eps.std_map_indicies[x].first;
			const unsigned int j=eps.std_map_indicies[x].second;
			// ToDo: find a better way to loop over the indices than using eps as input argument
			Tangent[i][j] = argument[i][j].dx(start_index).val();
		}
	}


	template<int dim>
	void SW_double2<dim>::get_curvature (double &Curvature, Sacado::Fad::DFad<DFadType> &argument)
	{
		Curvature = argument.dx(this->start_index).dx(this->start_index);
	}

	
	/*
	 * Compute Curvature d_argument / d_phi, where argument is already the derivative with respect to d_eps
	 */
	template<int dim>
	void SW_double2<dim>::get_curvature (SymmetricTensor<2,dim> &Curvature, SymmetricTensor<2,dim, Sacado::Fad::DFad<DFadType> > &argument, SymTensor2<dim> &eps )
	{
		for(unsigned int x=0;x<eps.n_dofs;++x)
		{
			const unsigned int i=eps.std_map_indicies[x].first;
			const unsigned int j=eps.std_map_indicies[x].second;

//			if ( i!=j )
//				Curvature[i][j] = 0.5 * argument[i][j].val().dx(start_index);	// ToDo: check the factor 0.5
//			else
				Curvature[i][j] = argument[i][j].val().dx(start_index);
		}
	}
	
	/*
	 * Compute Curvature d2_argument / d_phi_d_eps
	 */
	template<int dim>
	void SW_double2<dim>::get_curvature (SymmetricTensor<2,dim> &Curvature, Sacado::Fad::DFad<DFadType> &argument, SymTensor2<dim> &eps )
	{
		for(unsigned int x=0;x<eps.n_dofs;++x)
		{
			const unsigned int i=eps.std_map_indicies[x].first;
			const unsigned int j=eps.std_map_indicies[x].second;

			if ( i!=j )
				Curvature[i][j] = 0.5 * argument.dx(start_index).dx(x);	// ToDo: check the factor 0.5
			else
				Curvature[i][j] = argument.dx(start_index).dx(x);
		}
	}
		
	
	
	//###########################################################################################################//
	
	

	template<int dim>
	class DoFs_summary
	{
	public:
		void set_dofs( SymTensor<dim> &eps, SW_double<dim> &double_arg );
		void init_set_dofs( SymTensor2<dim> &eps, SymmetricTensor<2,dim> &eps_init, SW_double2<dim> &double_arg, double &double_init );

		void get_curvature( SymmetricTensor<2,dim> &Curvature, Sacado::Fad::DFad<DFadType> &argument, SymTensor2<dim> &eps,        SW_double2<dim> &double_arg );
		void get_curvature( SymmetricTensor<2,dim> &Curvature, Sacado::Fad::DFad<DFadType> &argument, SW_double2<dim> &double_arg, SymTensor2<dim> &eps        );

		
		void set_dofs( SymTensor<dim> &eps, SW_double<dim> &double_arg1, SW_double<dim> &double_arg2, SW_double<dim> &double_arg3 );
		// e.g. for strain, gamma_p and gamma_d
		 void set_dofs( SymTensor<dim> &eps, SW_double<dim> &double_arg1, SW_double<dim> &double_arg2 );
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
	
	template<int dim>
	void DoFs_summary<dim>::init_set_dofs(SymTensor2<dim> &eps, SymmetricTensor<2,dim> &eps_init, SW_double2<dim> &double_arg, double &double_init )

	{
		const unsigned int nbr_total_dofs = eps.n_dofs + double_arg.n_dofs;

		eps.start_index = 0;
		double_arg.start_index = eps.n_dofs;
		
		eps.init_set_dofs( eps_init, nbr_total_dofs);
		double_arg.init_set_dofs( double_init, nbr_total_dofs );
	}
	
	template<int dim>
	void DoFs_summary<dim>::set_dofs(SymTensor<dim> &eps, SW_double<dim> &double_arg1, SW_double<dim> &double_arg2, SW_double<dim> &double_arg3 )
	{
		const unsigned int nbr_total_dofs = eps.n_dofs + double_arg1.n_dofs + double_arg2.n_dofs + double_arg3.n_dofs ;

		eps.start_index = 0;
		double_arg1.start_index = eps.n_dofs;
		double_arg2.start_index = eps.n_dofs + double_arg1.n_dofs;
		double_arg3.start_index = eps.n_dofs + double_arg1.n_dofs + double_arg2.n_dofs;

		eps.set_dofs( nbr_total_dofs );
		double_arg1.set_dofs( nbr_total_dofs );
		double_arg2.set_dofs( nbr_total_dofs );
		double_arg3.set_dofs( nbr_total_dofs );
	}
	
	template<int dim>
	void DoFs_summary<dim>::set_dofs(SymTensor<dim> &eps, SW_double<dim> &double_arg1, SW_double<dim> &double_arg2 )
	{
		const unsigned int nbr_total_dofs = eps.n_dofs + double_arg1.n_dofs + double_arg2.n_dofs ;

		eps.start_index = 0;
		double_arg1.start_index = eps.n_dofs;
		double_arg2.start_index = eps.n_dofs + double_arg1.n_dofs;

		eps.set_dofs( nbr_total_dofs );
		double_arg1.set_dofs( nbr_total_dofs );
		double_arg2.set_dofs( nbr_total_dofs );
	}
	
	/*
	 * 
	 * example call: DoFs_summary.get_curvature(d2_energy_d_eps_d_phi, energy, eps_fad, phi_fad)
	 */
	template<int dim>
	void DoFs_summary<dim>::get_curvature( SymmetricTensor<2,dim> &Curvature, Sacado::Fad::DFad<DFadType> &argument, SymTensor2<dim> &eps, SW_double2<dim> &double_arg )
	{
		SymmetricTensor<2,dim, Sacado::Fad::DFad<DFadType> > d_arg_d_eps;
		eps.get_tangent(d_arg_d_eps, argument);
		double_arg.get_curvature(Curvature, d_arg_d_eps, eps);
	}

	
	/*
	 * 
	 * example call: DoFs_summary.get_curvature(d2_energy_d_phi_d_eps, energy, phi_fad, eps_fad)
	 */
	template<int dim>
	void DoFs_summary<dim>::get_curvature( SymmetricTensor<2,dim> &Curvature, Sacado::Fad::DFad<DFadType> &argument, SW_double2<dim> &double_arg, SymTensor2<dim> &eps  )
	{
		double_arg.get_curvature(Curvature, argument, eps);
	}
}






#endif // Sacado_Wrapper_H

