
// framework based on step-1 tutorial program by deal.II
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <iostream>
#include <fstream>
#include <cmath>

// Sacado
#include <Sacado.hpp>


// did not seem to make a difference:
#  include <deal.II/base/numbers.h>
#  include <deal.II/differentiation/ad/ad_number_traits.h>
#  include <deal.II/differentiation/ad/sacado_number_types.h>



using namespace dealii;

// defining a data type for the Sacada variables (simply used the standard types from the step-33 tutorial's introduction)
using fad_double = Sacado::Fad::DFad<double>;



// 1. example: simple scalar equations from step-33
void sacado_test_scalar ()
{
	std::cout << "Scalar Test:" << std::endl;
	fad_double a,b,c;
	a = 1; b = 2;	// at the point (a,b) = (1,2)
	a.diff(0,2);  // Set a to be dof 0, in a 2-dof system.
	b.diff(1,2);  // Set b to be dof 1, in a 2-dof system.
	c = 2*a+cos(a*b);
	double *derivs = &c.fastAccessDx(0); // Access derivatives
	std::cout << "Derivatives at the point (" << a << "," << b << ")" << std::endl;
	std::cout << "dc/da = " << derivs[0] << ", dc/db=" << derivs[1] << std::endl;
}



// 2. example: use of tensors
void sacado_test_tensor1 ()
{
	std::cout << "Tensor Test 1:" << std::endl;

	const unsigned int dim = 3;

	SymmetricTensor<2,dim, fad_double> sigma, eps;

	// init the strain tensor (the point at which the derivative shall be computed)
	eps[0][0] = 1;
	eps[1][1] = 2;
	eps[2][2] = 3;
	eps[0][1] = 4;
	eps[0][2] = 5;
	eps[1][2] = 6;

	// This is already very cumbersome -> easier option?
	eps[0][0].diff(0,6);
	eps[1][1].diff(1,6);
	eps[2][2].diff(2,6);
	eps[0][1].diff(3,6);
	eps[0][2].diff(4,6);
	eps[1][2].diff(5,6);
	
	

	// The equation describing the stresses (here just a simple test case)
	sigma = eps;

	std::cout << sigma << std::endl;

	// computing the derivatives for certain components of the resulting tangent modulus
	{
	double *derivs = &sigma[0][0].fastAccessDx(0); // Access derivatives
    std::cout << "d_sigma[0][0]/d_eps = ";
    for ( unsigned int i=0; i<6; ++i)
    	std::cout << derivs[i] << " , ";
    std::cout << std::endl;
	}
	{
	double *derivs = &sigma[1][2].fastAccessDx(0); // Access derivatives
    std::cout << "d_sigma[1][2]/d_eps = ";
    for ( unsigned int i=0; i<6; ++i)
    	std::cout << derivs[i] << " , ";
    std::cout << std::endl;
	}
}



// using a slightly more complicated stress equation
void sacado_test_tensor2 ()
{
	std::cout << "Tensor Test 2:" << std::endl;

	const unsigned int dim = 3;

	double kappa = 5;
	double mu = 2;


	SymmetricTensor<2,dim, fad_double> sigma, eps;
    std::map<unsigned int,std::pair<unsigned int,unsigned int>> std_map_indicies;

	// The point at which the derivative shall be computed:
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



	

	SymmetricTensor<2,dim, fad_double> stdTensor_I (( unit_symmetric_tensor<dim,fad_double>()) );
	SymmetricTensor<4,dim, fad_double> stdTensor_Idev ( (deviator_tensor<dim,fad_double>()) );
    
	sigma =(trace(eps) *  stdTensor_I);
    sigma *= kappa;
    SymmetricTensor<2,dim,fad_double> tmp = deviator<dim,fad_double>(symmetrize<dim,fad_double>(eps)); tmp*=(mu*2);
    sigma +=  tmp;

	
	std::cout << "sigma=" << sigma << std::endl;

	SymmetricTensor<4,dim> C_Sacado;

	//test
	for ( unsigned int i=0; i<dim; ++i)
		for ( unsigned int j=0; j<dim; ++j )
		{
			double *derivs = &sigma[i][j].fastAccessDx(0); // Access derivatives
            std::cout<<"size: "<<sigma[i][j].size()<<std::endl;

            for(unsigned int x=0;x<6;++x)
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

	// compute the analytical tangent for comparison
	double kappa_d = 5;
	double mu_d = 2;
	SymmetricTensor<4,dim> C_analy;
	C_analy = kappa_d * outer_product(unit_symmetric_tensor<dim>(), unit_symmetric_tensor<dim>()) + 2* mu_d * deviator_tensor<dim>(); //identity_tensor<dim>();


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

    std::cout << "sigma_analy: " << (C_analy*eps_d) << std::endl;
	


	for (unsigned int i=0; i<dim; ++i)
		for ( unsigned int j=0; j<dim; ++j)
			for ( unsigned int k=0; k<dim; ++k)
				for ( unsigned int l=0; l<dim; ++l)
					std::cout << "C_analy["<<i<<"]["<<j<<"]["<<k<<"]["<<l<<"] = " << C_analy[i][j][k][l] << " vs C_Sacado: " << C_Sacado[i][j][k][l] << std::endl;


	// compare the analytical with the numerical tangent and compute the error
	double error_Sacado_vs_analy=0;
	for (unsigned int i=0; i<dim; ++i)
		for ( unsigned int j=0; j<dim; ++j)
			for ( unsigned int k=0; k<dim; ++k)
				for ( unsigned int l=0; l<dim; ++l)
				{
					error_Sacado_vs_analy += std::fabs(C_Sacado[i][j][k][l] - C_analy[i][j][k][l]);
				}
					

	std::cout << "numerical error: " << error_Sacado_vs_analy << std::endl;
}



// using a slightly more complicated stress equation
// in version 2B some more stuff is tested and extended based on 2
void sacado_test_tensor2B ()
{
	std::cout << "Tensor Test 2B:" << std::endl;

	const unsigned int dim = 3;

	double kappa_param = 5;
	fad_double kappa (kappa_param);
	double mu = 2;


	SymmetricTensor<2,dim, fad_double> sigma, eps;
    std::map<unsigned int,std::pair<unsigned int,unsigned int>> std_map_indicies;

	// The point at which the derivative shall be computed:
	if(dim==3)
	{
		eps[0][0] = 1;
		eps[1][1] = 2;
		eps[2][2] = 3;

		eps[0][1] = 4;
		eps[0][2] = 5;
		eps[1][2] = 6;

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
	else if(dim==2)
	{
		eps[0][0] = 1;
		eps[1][1] = 2;

		eps[0][1] = 4;

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


    for ( unsigned int x=0; x<((dim==2)?3:6); ++x )
	{
		unsigned int i=std_map_indicies[x].first;
		unsigned int j=std_map_indicies[x].second;
		eps[i][j].diff(x,((dim==2)?3:6));
	}


	SymmetricTensor<2,dim, fad_double> stdTensor_I (( unit_symmetric_tensor<dim,fad_double>()) );
	SymmetricTensor<4,dim, fad_double> stdTensor_Idev ( (deviator_tensor<dim,fad_double>()) );

	sigma = kappa * (trace(eps) *  stdTensor_I);
    SymmetricTensor<2,dim,fad_double> tmp = deviator<dim,fad_double>(symmetrize<dim,fad_double>(eps)); tmp*=(mu*2);
    sigma +=  tmp;


	std::cout << "sigma=" << sigma << std::endl;

	SymmetricTensor<4,dim> C_Sacado;

	// reassemble the tangent as a fourth order tensor
	 for ( unsigned int i=0; i<dim; ++i)
		for ( unsigned int j=0; j<dim; ++j )
		{
			double *derivs = &sigma[i][j].fastAccessDx(0); // Access derivatives
            //std::cout<<"size: "<<sigma[i][j].size()<<std::endl;

            for(unsigned int x=0;x<6;++x)
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

	// compute the analytical tangent for comparison
	double kappa_d (kappa_param);
	double mu_d = 2;
	SymmetricTensor<4,dim> C_analy;
	C_analy = kappa_d * outer_product(unit_symmetric_tensor<dim>(), unit_symmetric_tensor<dim>()) + 2* mu_d * deviator_tensor<dim>(); //identity_tensor<dim>();


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

    std::cout << "sigma_analy: " << (C_analy*eps_d) << std::endl;



	for (unsigned int i=0; i<dim; ++i)
		for ( unsigned int j=0; j<dim; ++j)
			for ( unsigned int k=0; k<dim; ++k)
				for ( unsigned int l=0; l<dim; ++l)
					std::cout << "C_analy["<<i<<"]["<<j<<"]["<<k<<"]["<<l<<"] = " << C_analy[i][j][k][l] << " vs C_Sacado: " << C_Sacado[i][j][k][l] << std::endl;


	// compare the analytical with the numerical tangent and compute the error
	double error_Sacado_vs_analy=0;
	for (unsigned int i=0; i<dim; ++i)
		for ( unsigned int j=0; j<dim; ++j)
			for ( unsigned int k=0; k<dim; ++k)
				for ( unsigned int l=0; l<dim; ++l)
				{
					error_Sacado_vs_analy += std::fabs(C_Sacado[i][j][k][l] - C_analy[i][j][k][l]);
				}


	std::cout << "numerical error: " << error_Sacado_vs_analy << std::endl;
}


void sacado_test_tensor3 ()
{
    const unsigned int dim=3;
	std::cout << "Tensor Test 3:" << std::endl;
    Tensor<1,dim,fad_double> c;
	fad_double a,b;
    unsigned int n_dofs=2;
	a = 1; b = 2;	// at the point (a,b) = (1,2)
	a.diff(0,2);  // Set a to be dof 0, in a 2-dof system.
	b.diff(1,2);  // Set b to be dof 1, in a 2-dof system.
	c[0] = 2*a+3*b;
    c[1] = 4*a+5*b;
    c[2] = 6*a+7*b;
    
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


//debugger for the poor:
//#define lout(name) loutf(#name, (name))
//void loutf(char *name, double value) {
//    std::cout << name << " = " << value << std::endl;
//}



int main ()
{
	sacado_test_scalar ();

	std::cout << std::endl;

	sacado_test_tensor1 ();

	std::cout << std::endl;

	sacado_test_tensor2 ();

    std::cout << std::endl;

    sacado_test_tensor2B ();

    std::cout << std::endl;

    sacado_test_tensor3();

}
