#ifndef Sacado_Wrapper_H
#define Sacado_Wrapper_H

// @section includes Include Files
// The data type SymmetricTensor and some related operations, such as trace, symmetrize, deviator, ... for tensor calculus
#include <deal.II/base/symmetric_tensor.h>
// This header is crucial, it enables tensor operations also for fad_double data types (somehow)
#include <deal.II/differentiation/ad.h>

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

namespace SacadoQP
{
	/*
	 * Extract values from SymmetricTensor of Sacado or standard double data types
	 */
	template<int dim>
	SymmetricTensor<2,dim,double> get_value ( const SymmetricTensor<2,dim,fad_double> &symtensor )
	{
		SymmetricTensor<2,dim,double> symtensor_double;
		for ( unsigned int i=0; i<dim; ++i)
			for ( unsigned int j=i; j<dim; ++j ) // start at \a i only possible for symmetric tensors
				symtensor_double[i][j] = (symtensor[i][j]).val();

		return symtensor_double;
	}
	template<int dim>
	SymmetricTensor<2,dim,double> get_value ( const SymmetricTensor<2,dim,double> &symtensor )
	{
		return symtensor;
	}


	/*
	 * Extract values from Sacado or standard double data types
	 */
	double get_value ( const fad_double &faddouble )
	{
		return faddouble.val();
	}
	double get_value ( const double &doubleData )
	{
		return doubleData;
	}

	/*
	 * Set values for Sacado or standard double data types
	 */
	void set_value ( fad_double &faddouble, double &doubleData )
	{
		faddouble.val()=doubleData;
	}
	void set_value ( double &faddouble, double &doubleData )
	{
		faddouble = doubleData;
	}


	template<int dim, typename Number>
	class BuCa
	{
	public:
		// Use different BuCa constructors for (ten,ten), (ten,LM), (ten,LM,LM), ...
		/**
		 * Only a single symmetric tensor as AD-dofs
		 * @param variable_x
		 */
		BuCa( const SymmetricTensor<2,dim,Number> &variable_x )
		:
		start_indices(1),
		n_dofs_list(1)
		{
			get_index_map<dim>( std_map_indicies );

			start_indices[0]=0;
			n_dofs_list[0] = ((dim==2)?3:6); // ToDo: Find a way to use eps.n_independent_components(); that gives ((dim==2)?3:6) for a SymmetricTensor

			nbr_total_dofs = n_dofs_list.back();
		}

		/**
		 * One symmetric tensor and a single scalar as AD-dofs
		 * @param variable_x1 SymmetricTensor
		 * @param variable_x2 Scalar
		 */
		BuCa( const SymmetricTensor<2,dim,Number> &variable_x1, const Number &variable_x2 )
		:
		start_indices(2),
		n_dofs_list(2)
		{
			get_index_map<dim>( std_map_indicies );

			start_indices[0] = 0;
			n_dofs_list[0] = ((dim==2)?3:6); // ToDo: Find a way to use eps.n_independent_components(); that gives ((dim==2)?3:6) for a SymmetricTensor

			start_indices[1] = n_dofs_list[0];
			n_dofs_list[1] = 1;

			nbr_total_dofs = n_dofs_list[0] + n_dofs_list[1];
		}
		// Everything else could be captured by std::vector<SymmetricTensor> and std::vector<Number>

		//~BuCa(){}

		// To simplify the access to the dofs we define a map that relate the components of our strain tensor to the dof-nbr
	     std::map<unsigned int,std::pair<unsigned int,unsigned int>> std_map_indicies;

	//private:
		std::vector<unsigned int> start_indices;
		std::vector<unsigned int> n_dofs_list;
		unsigned int nbr_total_dofs;

	public:
		// These should be pure accessor functions (somehow I still cannot pass "const BuCa" as argument in the fncs below
		unsigned int start_index( const enums::enum_SacadoDoFs &SacDoFs ) const
		{
			return start_indices[SacDoFs];
		};
		unsigned int n_ADdofs( const enums::enum_SacadoDoFs &SacDoFs ) const
		{
			return n_dofs_list[SacDoFs];
		};
		unsigned int n_total_ADdofs(  ) const
		{
			return nbr_total_dofs;
		};
	};


	/*
	 * Set the components of the tensor as Sacado dofs
	 */
	/**
	 *
	 * @param variable_x Cannot be used as \a const, because we set the \a diff members
	 * @param SacDoFs
	 * @param bucket
	 */
	template<int dim>
	void set_dofs( SymmetricTensor<2,dim,fad_double> &variable_x, const enums::enum_SacadoDoFs SacDoFs, BuCa<dim,fad_double> &bucket )
	{
		// Instead of calling the *.diff(*) on the components one-by-one we could also use the following for-loop, so
		// we also use the map to set the dofs
		const unsigned int start_index = bucket.start_index(SacDoFs);
		const unsigned int end_index = start_index + bucket.n_ADdofs(SacDoFs);

        for( unsigned int x=start_index; x < end_index; ++x )
		{
			unsigned int i=bucket.std_map_indicies[x].first;
			unsigned int j=bucket.std_map_indicies[x].second;
			variable_x[i][j].diff( x, bucket.n_total_ADdofs() );
		}
	}
	template<int dim>
	void set_dofs( SymmetricTensor<2,dim,double> &variable_x, const enums::enum_SacadoDoFs SacDoFs, BuCa<dim,double> &bucket )
	{
		AssertThrow(false,ExcMessage("When using data type double, you don't need to set_dofs for an analytical tangent."));
	}

	/**
	 *
	 * @param variable_x Cannot be used as \a const, because we set the \a diff members
	 * @param SacDoFs
	 * @param bucket
	 */
	template<int dim>
	void set_dofs( fad_double &variable_x2, const enums::enum_SacadoDoFs SacDoFs, BuCa<dim,fad_double> &bucket )
	{
		// Instead of calling the *.diff(*) on the components one-by-one we could also use the following for-loop, so
		// we also use the map to set the dofs
		const unsigned int start_index = bucket.start_index(SacDoFs);

		variable_x2.diff( start_index, bucket.n_total_ADdofs() );
	}
	template<int dim>
	void set_dofs( double &variable_x2, const enums::enum_SacadoDoFs SacDoFs, BuCa<dim,double> &bucket )
	{
		AssertThrow(false,ExcMessage("When using data type double, you don't need to set_dofs for an analytical tangent."));
	}


	/*
	 * Extract Derivatives
	 */
	// @todo Currently we cannot differ between a tangent wrt to a tensor -> 4th order or wrt to a scalar -> 2th order output
	template<int dim>
	SymmetricTensor<4,dim> get_tangent ( SymmetricTensor<2,dim,fad_double> &function_f, const enums::enum_SacadoDoFs SacDoFs, BuCa<dim,fad_double> &bucket )
	{
		SymmetricTensor<4,dim> Tangent;
		const unsigned int start_index = bucket.start_index(SacDoFs);
		const unsigned int end_index = start_index + bucket.n_ADdofs(SacDoFs);

		for ( unsigned int i=0; i<dim; ++i)
			for ( unsigned int j=0; j<dim; ++j )
			{
				double *derivs = &function_f[i][j].fastAccessDx(0); // Access the derivatives of the (i,j)-th component of \a sigma
	            // We loop over all the dofs. To be able to use this independent of the chosen dimension \a dim, we use a ternary operator
	            // to decide whether we have to loop over 6 derivatives or just 3.
	            for( unsigned int x=start_index; x < end_index; ++x )
	            {
	                unsigned int k=bucket.std_map_indicies[x].first;
	                unsigned int l=bucket.std_map_indicies[x].second;

	                if(k!=l)/*Compare to Voigt notation since only SymmetricTensor instead of Tensor*/
	                {
	                    Tangent[i][j][k][l] = 0.5*derivs[x];
	                    Tangent[i][j][l][k] = 0.5*derivs[x];
	                }
	                else
	                	Tangent[i][j][k][l] = derivs[x];
	            }
			}
		return Tangent;
	}
	// @todo-optimize This shall just catch the error message, so we could call the get_tanget with the double data type -> better ideas?!
	template<int dim>
	SymmetricTensor<4,dim> get_tangent ( const SymmetricTensor<2,dim,double> &function_f, const enums::enum_SacadoDoFs SacDoFs, BuCa<dim,double> &bucket )
	{
		SymmetricTensor<4,dim> Tangent;

		AssertThrow(false,ExcMessage("When using data type double, you cannot use Sacado for the tangents."));
		return Tangent;
	}


	template<int dim>
	double get_tangent ( fad_double &function_f, const enums::enum_SacadoDoFs SacDoFs, BuCa<dim,fad_double> &bucket )
	{
		const unsigned int start_index = bucket.start_index(SacDoFs);

		double *derivs = &function_f.fastAccessDx(0); // Access the derivatives of the (i,j)-th component of \a sigma
		return derivs[start_index];
	}
	template<int dim>
	double get_tangent ( double &function_f, const enums::enum_SacadoDoFs SacDoFs, BuCa<dim,double> &bucket )
	{
		double Tangent;

		AssertThrow(false,ExcMessage("When using data type double, you cannot use Sacado for the tangents."));
		return Tangent;
	}

	/**
	 * We set the derivative of the given \a SacDoF in the argument \a variable_x to the value \a value,
	 * or to zero in case no value is input.
	 * @param variable_x
	 * @param SacDoFs
	 * @param bucket
	 * @param value
	 */
	template<int dim>
	void set_deriv ( fad_double &variable_x, const enums::enum_SacadoDoFs SacDoF, BuCa<dim,fad_double> &bucket, double value=0. )
	{
		double *derivs = &variable_x.fastAccessDx(0);
		derivs[bucket.start_index(SacDoF)] = value;
	}
	template<int dim>
	void set_deriv ( double &variable_x, const enums::enum_SacadoDoFs SacDoF, BuCa<dim,double> &bucket, double value=0. )
	{
		AssertThrow(false,ExcMessage("When using data type fad_double, you cannot use this function."));
	}
}




#endif // Sacado_Wrapper_H

