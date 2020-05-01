#ifndef Sacado_auxiliary_functions_H
#define Sacado_auxiliary_functions_H

#include <deal.II/base/symmetric_tensor.h>
#include <iostream>
#include <Sacado.hpp>
using fad_double = Sacado::Fad::DFad<double>;

using namespace dealii;

template<int dim>
SymmetricTensor<2,dim,fad_double> init_Sacado(const SymmetricTensor<2,dim,double> &SymTensor)
{
	SymmetricTensor<2,dim,fad_double> SymTensor_Sacado;
	 for ( unsigned int i=0; i<dim; ++i)
		for ( unsigned int j=i; j<dim; ++j )
			SymTensor_Sacado[i][j] = SymTensor[i][j];
	 
	 return SymTensor_Sacado;
}


template<int dim>
SymmetricTensor<2,dim> extract_value_from_Sacado(const SymmetricTensor<2,dim,fad_double> &SymTensor_Sacado)
{
	SymmetricTensor<2,dim> SymTensor;
	// .val() extracts the \a double value from the Sacado variable
	 for ( unsigned int i=0; i<dim; ++i)
		for ( unsigned int j=i; j<dim; ++j )
			SymTensor[i][j] = SymTensor_Sacado[i][j].val();
	 
	 return SymTensor;
}


template<int dim>
Vector<fad_double> reinit(unsigned int a, bool b)
{
	Vector<fad_double> abc (5);
	 return abc;
}

#endif //Sacado_auxiliary_functions_H
