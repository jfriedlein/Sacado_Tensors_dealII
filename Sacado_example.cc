
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


#include "../MaterialModel.h"


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


// using the pure damage stress equation
// ToDo: e.g. use TimerOutput to compare the needed computation time for this Sacado approach vs using the analytical tangent
void sacado_test_tensor4 ()
{
	std::cout << "Tensor Test 4:" << std::endl;

	const unsigned int dim = 3;

	SymmetricTensor<2,dim, fad_double> sigma, eps;
    std::map<unsigned int,std::pair<unsigned int,unsigned int>> std_map_indicies;

    {
	// The point at which the derivative shall be computed:
	if(dim==3)
	{
		eps[0][0] = 1e-3;
		eps[1][1] = 2e-3;
		eps[2][2] = 3e-3;

		eps[0][1] = 4e-3;
		eps[0][2] = 5e-3;
		eps[1][2] = 6e-3;

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
		eps[0][0] = 1e-3;
		eps[1][1] = 2e-3;

		eps[0][1] = 4e-3;

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
    }


	// Setting the material parameters
	double kappa_param = 10e3;
	fad_double kappa ( kappa_param );
	fad_double mu = 7.5e3;
	fad_double eta_1 = 2.;
	fad_double eta_2 = 0.4;
	fad_double eta_3 = 10;
	fad_double eta_4 = 0.6;
	fad_double d_n = 0.1;	// from the last load step (history)
	fad_double gamma_d_n = 5.49126223844152; //5.4913636605639; // those values would come from e.g. the last iteration ...
	fad_double q_min=0.05;
	fad_double d_n_k = 5.59126223844152; //5.5913636605639; // ... while this accuracy is sufficient
	fad_double K = 2e3;
	fad_double beta_inf = 10;
	fad_double alpha_n = 0;

	// declare some auxiliary variables
	 SymmetricTensor<2,dim,fad_double> eps_dev_e, sigma_dev_t, n;
	 fad_double f_1_k, f_2_k,
				d, gamma_d_k, T11, r_d_k,
				J00, J10, J11, J01, DET,
				M11, M10, M01, M00,
				Q_gamma_d, Q_phi_d, Phi_d, q_k;

	 std::vector<fad_double> eta_i (5);
	 eta_i[1]=eta_1;
	 eta_i[2]=eta_2;
	 eta_i[3]=eta_3;
	 eta_i[4]=eta_4;

	 std::vector<fad_double> f_i (5);
	 std::vector<fad_double> f_diff_i (5);
	 std::vector<fad_double> f_diffdiff_i (5);

	// ######################################################################################################################################################

	// compute the stress tensor based on the given equations
	 for ( unsigned int i=1; i<5; ++i)	// run over eta 1...4 (ignore first entry i=0)
	 {
		 fad_double tmp = std::exp( - eta_i[i] * d_n_k );
		 f_i[i]          =                       tmp;
		 f_diff_i[i]     = - eta_i[i] *          tmp;
		 f_diffdiff_i[i] = eta_i[i] * eta_i[i] * tmp;
	 }


	 for ( unsigned int i=0; i<dim; ++i)
		for ( unsigned int j=0; j<dim; ++j )
		{
			sigma_dev_t[i][j] = 2 * mu * deviator(eps)[i][j];
			eps_dev_e[i][j] = deviator(eps)[i][j];

		}

	// n is only computable after eps_dev_e is known completely, hence the separate for-loop
	 for ( unsigned int i=0; i<dim; ++i)
		for ( unsigned int j=0; j<dim; ++j )
			n[i][j] = eps_dev_e[i][j] / eps_dev_e.norm();


	q_k = - 0.5 * f_diff_i[1] * kappa * trace(eps) * trace(eps)
		  - f_diff_i[2] * mu * eps_dev_e.norm() * eps_dev_e.norm();

	//lout (kappa_param);

	Phi_d = q_k - q_min * std::pow(1-f_i[3], 2./3.);
	fad_double Q_phi_p = -1;
	Q_phi_d = Phi_d / (std::sqrt(Phi_d*Phi_d + gamma_d_n*gamma_d_n));

	fad_double Q_gamma_p = 0;
	Q_gamma_d = gamma_d_n / (std::sqrt(Phi_d*Phi_d + gamma_d_n*gamma_d_n));

	M00 = - 2 * mu * f_i[2] / ( f_i[4]*f_i[4] )
		  - ( 2./3. * K * ( 1 - K/beta_inf * alpha_n ))
			/ ( std::pow( (1 + std::sqrt(2./3.) * K/beta_inf * 0 ) ,2) );

	M01 = ( f_diff_i[2]*f_i[4] - f_i[2]*f_diff_i[4] )/(f_i[4]*f_i[4]) * sigma_dev_t.norm() - 0;

	// compute the double contraction of eps_dev_e and n
	 fad_double eps_dev_e_n = 0;
	 for ( unsigned int i=0; i<dim; ++i)
		for ( unsigned int j=0; j<dim; ++j )
			eps_dev_e_n += eps_dev_e[i][j] * n[i][j];

	M10 = 2 * mu * f_diff_i[2] / f_i[4] * ( eps_dev_e_n );

	M11 = - f_diffdiff_i[1] * 0.5 * kappa * trace(eps) * trace(eps)
		  - f_diffdiff_i[2] * mu  * eps_dev_e.norm() * eps_dev_e.norm()
		  + 2./3. * q_min * f_diff_i[3] / (std::pow(1-f_i[3], 1./3.));

	J00 = (Q_phi_p+1)*M00 + Q_gamma_p - 1;
	J01 = (Q_phi_p+1)*M01;
	J10 = (Q_phi_d+1)*M10;
	J11 = (Q_phi_d+1)*M11 + Q_gamma_d - 1;
	DET = J00*J11-J10*J01;
	T11 =  J00 / DET;

	r_d_k = std::sqrt(Phi_d*Phi_d + gamma_d_n*gamma_d_n) + Phi_d - gamma_d_n;

	gamma_d_k = gamma_d_n - T11 * r_d_k;
	d = d_n + gamma_d_k;

	f_1_k = std::exp(-eta_1 * d);
	f_2_k = std::exp(-eta_2 * d);

	SymmetricTensor<2,dim, fad_double> stdTensor_I (( unit_symmetric_tensor<dim,fad_double>()) );
	SymmetricTensor<4,dim, fad_double> stdTensor_Idev ( (deviator_tensor<dim,fad_double>()) );

	// computing the stress
	for ( unsigned int i=0; i<dim; ++i)
		for ( unsigned int j=0; j<dim; ++j )
		{
			sigma[i][j] = f_1_k * kappa * trace(eps) * stdTensor_I[i][j] + 2 * f_2_k * mu * eps_dev_e[i][j];
		}

	SymmetricTensor<4,dim> C_Sacado;

	{
	// reassemble the tangent as a fourth order tensor
	 for ( unsigned int i=0; i<dim; ++i)
		for ( unsigned int j=0; j<dim; ++j )
		{
			double *derivs = &sigma[i][j].fastAccessDx(0); // Access derivatives

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
	}

	// ######################################################################################################################################################

	// compute the analytical tangent
	enum enum_plasti_dmg
	{
		plasti = 0,
		dmg =    1
	};

	double kappa_d = 10e3;
	double mu_d = 7.5e3;
	double d_n_d = 0.1;
	double gamma_d_n_d = 0.2;
	double q_min_d = 0.05;
	Tensor<1,2> f_gammas_k;
	f_gammas_k[dmg] = gamma_d_n_d;
	Vector<double> eta_i_d (5);
	eta_i_d[1] = 2.0;
	eta_i_d[2] = 0.4;
	eta_i_d[3] = 10;
	eta_i_d[4] = 0.6;
	SymmetricTensor<2,dim> eps_p_n_d; // is zero

	SymmetricTensor<2,dim> eps_d;

		if(dim==3)
	{
		eps_d[0][0] = 1e-3;
		eps_d[1][1] = 2e-3;
		eps_d[2][2] = 3e-3;

		eps_d[0][1] = 4e-3;
		eps_d[0][2] = 5e-3;
		eps_d[2][1] = 6e-3;

	}
	else if(dim==2)
	{
		eps_d[0][0] = 1e-3;
		eps_d[1][1] = 2e-3;

		eps_d[1][0] = 4e-3;

	}
	else
	{
		throw std::runtime_error("only dim==2 or dim==3 allowed");
	}

	//SymmetricTensor<4,dim> C_analy = compute_C_puredmg ( eps_d, kappa_d, mu_d, d_n_d, f_gammas_k, q_min_d, eta_i_d );

	SymmetricTensor<2,dim> sigma_n1;
	double stress_vM, some_Data;

	MaterialModel<dim> material( mu_d, 5e3, 0.2,
									 1e6, 2000,
									 eta_i_d, q_min_d, 10,
									 3.0, 10.0 );

	double alpha_n_d = 0;
	bool eval_QP = false;

	SymmetricTensor<4,dim> C_analy = material.get_StiffnessTensor_C_local_damage(
										eps_d, alpha_n_d, d_n_d,
										eps_p_n_d, sigma_n1,
										stress_vM, eval_QP, some_Data );

    //std::cout << "sigma_analy: " << (C_analy*eps_d) << std::endl;


	for (unsigned int i=0; i<dim; ++i)
		for ( unsigned int j=0; j<dim; ++j)
			for ( unsigned int k=0; k<dim; ++k)
				for ( unsigned int l=0; l<dim; ++l)
					std::cout << "C_analy["<<i<<"]["<<j<<"]["<<k<<"]["<<l<<"] = " << C_analy[i][j][k][l]
							  << " vs C_Sacado: " << C_Sacado[i][j][k][l] << std::endl;


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
	std::cout << " Remark: error decreases to around 1e-5 for a lower FBR-tolerance" << std::endl;
}



// using the complete damage stress equation
// ToDo: e.g. use TimerOutput to compare the needed computation time for this Sacado approach vs using the analytical tangent
void sacado_test_tensor5 ()
{
	std::cout << "Tensor Test 5:" << std::endl;

	const unsigned int dim = 3;

	SymmetricTensor<2,dim, fad_double> sigma, eps_n1;
	fad_double phi_n1 = .5;
    std::map<unsigned int,std::pair<unsigned int,unsigned int>> std_map_indicies;

    {
	// The point at which the derivative shall be computed:
	if(dim==3)
	{
		eps_n1[0][0] = 1e-3;
		eps_n1[1][1] = 2e-3;
		eps_n1[2][2] = 3e-3;

		eps_n1[0][1] = 4e-3;
		eps_n1[0][2] = 5e-3;
		eps_n1[1][2] = 6e-3;

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
		eps_n1[0][0] = 1e-3;
		eps_n1[1][1] = 2e-3;

		eps_n1[0][1] = 4e-3;

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
		eps_n1[i][j].diff(x,((dim==2)?3:6));
	}
    }


	// Setting the material parameters
	double kappa_param = 10e3;
	fad_double kappa ( kappa_param );
	fad_double mu = 7.5e3;
	fad_double eta_1 = 2.;
	fad_double eta_2 = 0.4;
	fad_double eta_3 = 10;
	fad_double eta_4 = 0.6;
	fad_double d_n = 0.1;	// from the last load step (history)
	fad_double gamma_d_n = 5.49126223844152; //5.4913636605639; // those values would come from e.g. the last iteration ...
	fad_double q_min=0.05;
	fad_double d_n_k = 5.59126223844152; //5.5913636605639; // ... while this accuracy is sufficient
	fad_double K = 2e3;
	fad_double beta_inf = 10;
	fad_double alpha_n = 0;//0.5; // #c
	fad_double alpha_k = 0;//0.7; // #c
	fad_double gamma_p_n = 0;//0.2; // #c
	fad_double beta_d = 0; //10;
	fad_double yield_stress = 1e6;//40;

	SymmetricTensor<2,dim,fad_double> eps_p_n, eps_p_k; // ToDo: init this with something

	// declare some auxiliary variables
	 SymmetricTensor<2,dim,fad_double> eps_dev_e, sigma_dev_t, sigma_dev_k, n, eps_dev_e_k, eps_p_k1;
	 fad_double f_1_k1, f_2_k1, f_4_k1,
				d_n1,
				gamma_p_k1, gamma_d_k1,
				T00, T01, T10, T11,
				r_p_k, r_d_k,
				J00, J10, J11, J01, DET,
				M11, M10, M01, M00,
				Q_gamma_p, Q_gamma_d, Q_Phi_p, Q_Phi_d, Phi_p, Phi_d, R_k, q_k;

	 std::vector<fad_double> eta_i (5);
	 eta_i[1]=eta_1;
	 eta_i[2]=eta_2;
	 eta_i[3]=eta_3;
	 eta_i[4]=eta_4;

	 std::vector<fad_double> f_i (5);
	 std::vector<fad_double> f_diff_i (5);
	 std::vector<fad_double> f_diffdiff_i (5);

	// ######################################################################################################################################################

	// compute the stress tensor based on the given equations
	 for ( unsigned int i=1; i<5; ++i)	// run over eta 1...4 (ignore first entry i=0)
	 {
		 fad_double tmp = std::exp( - eta_i[i] * d_n_k );
		 f_i[i]          =                       tmp;
		 f_diff_i[i]     = - eta_i[i] *          tmp;
		 f_diffdiff_i[i] = eta_i[i] * eta_i[i] * tmp;
	 }


	 for ( unsigned int i=0; i<dim; ++i)
		for ( unsigned int j=0; j<dim; ++j )
		{
			sigma_dev_t[i][j] = 2 * mu * ( deviator(eps_n1)[i][j] - eps_p_n[i][j] );
			eps_dev_e[i][j] = deviator(eps_n1)[i][j] - eps_p_n[i][j];
			sigma_dev_k[i][j] = 2 * mu * f_i[2] * ( deviator(eps_n1)[i][j] - eps_p_k[i][j] );
		}

	// n is only computable after eps_dev_e is known completely, hence the separate for-loop
	 for ( unsigned int i=0; i<dim; ++i)
		for ( unsigned int j=0; j<dim; ++j )
			n[i][j] = eps_dev_e[i][j] / eps_dev_e.norm();

	fad_double sqrtR_2_3 = std::sqrt( 2./3. );
	fad_double F2over3 = 2./3.;


	q_k = - 0.5 * f_diff_i[1] * kappa * trace(eps_n1) * trace(eps_n1)
		  + beta_d * ( phi_n1 - d_n_k )
		  - f_diff_i[2] * mu * eps_dev_e.norm() * eps_dev_e.norm();

	R_k = - K * alpha_k;

	Phi_p = 1./f_i[4] * sigma_dev_k.norm()
			- sqrtR_2_3 * ( yield_stress - R_k );
	Phi_d = q_k - q_min * std::pow(1-f_i[3], 2./3.);


	fad_double Q_denom_p = std::sqrt(Phi_p*Phi_p + gamma_p_n*gamma_p_n);
	Q_Phi_p = Phi_p / Q_denom_p;
	Q_gamma_p = gamma_p_n / Q_denom_p;

	fad_double Q_denom_d = std::sqrt(Phi_d*Phi_d + gamma_d_n*gamma_d_n);
	Q_Phi_d = Phi_d / Q_denom_d;
	Q_gamma_d = gamma_d_n / Q_denom_d;


	M00 = - 2 * mu * f_i[2] / ( f_i[4]*f_i[4] )
		  - ( F2over3 * K * ( 1 - K/beta_inf * alpha_n ))
			/ ( std::pow( (1 + sqrtR_2_3 * K/beta_inf * gamma_p_n ) ,2) );

	M01 = ( f_diff_i[2]*f_i[4] - f_i[2]*f_diff_i[4] )
		   / (f_i[4]*f_i[4]) * sigma_dev_t.norm()
		  - 2 * mu * ( f_diff_i[2] * f_i[4] - 2 * f_i[2] * f_diff_i[4])
		   / std::pow(f_i[4],3) * gamma_p_n;

	// compute the double contraction of eps_dev_e and n
	 fad_double eps_dev_e_DC_n = 0;
	 for ( unsigned int i=0; i<dim; ++i)
		for ( unsigned int j=0; j<dim; ++j )
			eps_dev_e_DC_n += eps_dev_e[i][j] * n[i][j];

	M10 = 2 * mu * f_diff_i[2] / f_i[4] * ( eps_dev_e_DC_n );

	M11 = - f_diffdiff_i[1] * 0.5 * kappa * trace(eps_n1) * trace(eps_n1)
		  - f_diffdiff_i[2] * mu  * eps_dev_e.norm() * eps_dev_e.norm()
		  - f_diff_i[2] * 2 * mu * f_diff_i[4] / ( f_i[4]*f_i[4] ) * gamma_p_n * eps_dev_e_DC_n
		  - beta_d
		  + F2over3 * q_min * f_diff_i[3] / (std::pow(1-f_i[3], 1./3.));

	J00 = (Q_Phi_p+1)*M00 + Q_gamma_p - 1;
	J01 = (Q_Phi_p+1)*M01;
	J10 = (Q_Phi_d+1)*M10;
	J11 = (Q_Phi_d+1)*M11 + Q_gamma_d - 1;
	DET = J00*J11 - J10*J01;

	T00 =  J11 / DET;
	T01 = -J01 / DET;
	T10 = -J10 / DET;
	T11 =  J00 / DET;

	r_p_k = Q_denom_p + Phi_p - gamma_p_n;
	r_d_k = Q_denom_d + Phi_d - gamma_d_n;

	gamma_p_k1 = gamma_p_n - T00 * r_p_k - T01 * r_d_k;
	gamma_d_k1 = gamma_d_n - T10 * r_p_k - T11 * r_d_k;

	d_n1 = d_n + gamma_d_k1;

	f_1_k1 = std::exp(-eta_i[1] * d_n1);
	f_2_k1 = std::exp(-eta_i[2] * d_n1);
	f_4_k1 = std::exp(-eta_i[4] * d_n1);

	for ( unsigned int i=0; i<dim; ++i)
		for ( unsigned int j=0; j<dim; ++j )
		{
			eps_p_k1[i][j] = eps_p_n[i][j] + 1./f_4_k1 * gamma_p_k1 * n[i][j];
			eps_dev_e_k[i][j] = deviator(eps_n1)[i][j] - eps_p_k1[i][j];
		}

	SymmetricTensor<2,dim, fad_double> stdTensor_I (( unit_symmetric_tensor<dim,fad_double>()) );
	SymmetricTensor<4,dim, fad_double> stdTensor_Idev ( (deviator_tensor<dim,fad_double>()) );

	// computing the stress
	for ( unsigned int i=0; i<dim; ++i)
		for ( unsigned int j=0; j<dim; ++j )
		{
			sigma[i][j] = f_1_k1 * kappa * trace(eps_n1) * stdTensor_I[i][j] + 2 * f_2_k1 * mu * eps_dev_e_k[i][j];
		}

	SymmetricTensor<4,dim> C_Sacado;

	{
	// reassemble the tangent as a fourth order tensor
	 for ( unsigned int i=0; i<dim; ++i)
		for ( unsigned int j=0; j<dim; ++j )
		{
			double *derivs = &sigma[i][j].fastAccessDx(0); // Access derivatives

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
	}

	// ######################################################################################################################################################

	// compute the analytical tangent
	enum enum_plasti_dmg
	{
		plasti = 0,
		dmg =    1
	};

	double kappa_d = 10e3;
	double mu_d = 7.5e3;
	double d_n_d = 0.1;
	double gamma_d_n_d = 0.2;
	double q_min_d = 0.05;
	Tensor<1,2> f_gammas_k;
	f_gammas_k[dmg] = gamma_d_n_d;
	Vector<double> eta_i_d (5);
	eta_i_d[1] = 2.0;
	eta_i_d[2] = 0.4;
	eta_i_d[3] = 10;
	eta_i_d[4] = 0.6;
	SymmetricTensor<2,dim> eps_p_n_d; // is zero

	SymmetricTensor<2,dim> eps_d;

		if(dim==3)
	{
		eps_d[0][0] = 1e-3;
		eps_d[1][1] = 2e-3;
		eps_d[2][2] = 3e-3;

		eps_d[0][1] = 4e-3;
		eps_d[0][2] = 5e-3;
		eps_d[2][1] = 6e-3;

	}
	else if(dim==2)
	{
		eps_d[0][0] = 1e-3;
		eps_d[1][1] = 2e-3;

		eps_d[1][0] = 4e-3;

	}
	else
	{
		throw std::runtime_error("only dim==2 or dim==3 allowed");
	}

	//SymmetricTensor<4,dim> C_analy = compute_C_puredmg ( eps_d, kappa_d, mu_d, d_n_d, f_gammas_k, q_min_d, eta_i_d );

	SymmetricTensor<2,dim> sigma_n1;
	double stress_vM, some_Data;

	MaterialModel<dim> material( mu_d, 5e3, 0.2,
									 1e6, 2000,
									 eta_i_d, q_min_d, 10,
									 3.0, 10.0 );

	double alpha_n_d = 0;
	bool eval_QP = false;

	SymmetricTensor<4,dim> C_analy = material.get_StiffnessTensor_C_local_damage(
										eps_d, alpha_n_d, d_n_d,
										eps_p_n_d, sigma_n1,
										stress_vM, eval_QP, some_Data );

    //std::cout << "sigma_analy: " << (C_analy*eps_d) << std::endl;


	for (unsigned int i=0; i<dim; ++i)
		for ( unsigned int j=0; j<dim; ++j)
			for ( unsigned int k=0; k<dim; ++k)
				for ( unsigned int l=0; l<dim; ++l)
					std::cout << "C_analy["<<i<<"]["<<j<<"]["<<k<<"]["<<l<<"] = " << C_analy[i][j][k][l]
							  << " vs C_Sacado: " << C_Sacado[i][j][k][l] << std::endl;


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



void call_to_MatMod()
{
	const unsigned int dim = 3;

	SymmetricTensor<2,dim> eps_n1;

	fad_double phi_n1 = .5;


	// The point at which the derivative shall be computed:
	if(dim==3)
	{
		eps_n1[0][0] = 1e-3;
		eps_n1[1][1] = 2e-3;
		eps_n1[2][2] = 3e-3;

		eps_n1[0][1] = 4e-3;
		eps_n1[0][2] = 5e-3;
		eps_n1[1][2] = 6e-3;
	}

	// Setting the material parameters
	double kappa_param = 10e3;
	double kappa ( kappa_param );
	double mu = 7.5e3;
	double eta_1 = 2.;
	double eta_2 = 0.4;
	double eta_3 = 10;
	double eta_4 = 0.6;
	Vector<double> eta_i (5);
	eta_i[1] = eta_1;
	eta_i[2] = eta_2;
	eta_i[3] = eta_3;
	eta_i[4] = eta_4;

	double d_n = 0.1;	// from the last load step (history)
	double gamma_d_n = 5.49126223844152; //5.4913636605639; // those values would come from e.g. the last iteration ...
	double q_min=0.05;
	double d_n_k = 5.59126223844152; //5.5913636605639; // ... while this accuracy is sufficient
	double K = 2e3;
	double beta_inf = 10;
	double alpha_n = 0;//0.5; // #c
	double alpha_k = 0;//0.7; // #c
	double gamma_p_n = 0;//0.2; // #c
	double beta_d = 0; //10;
	double yield_stress = 1e6;//40;

	SymmetricTensor<2,dim> sigma_n1, eps_p_n;
	double stress_vM, some_Data;

	MaterialModel<dim> material( mu, 5e3, 0.2,
									yield_stress, K,
									 eta_i, q_min, beta_inf,
									 3.0, beta_d );

	double alpha_n_d = 0;
	bool eval_QP = false;

	SymmetricTensor<4,dim> C_Sacado = material.get_StiffnessTensor_C_local_damage(
										eps_n1, alpha_n, d_n,
										eps_p_n, sigma_n1,
										stress_vM, eval_QP, some_Data );

	for (unsigned int i=0; i<dim; ++i)
		for ( unsigned int j=0; j<dim; ++j)
			for ( unsigned int k=0; k<dim; ++k)
				for ( unsigned int l=0; l<dim; ++l)
					std::cout << "C["<<i<<"]["<<j<<"]["<<k<<"]["<<l<<"] = " << C_Sacado[i][j][k][l] << std::endl;
}



int main ()
{
//	sacado_test_scalar ();
//
//	std::cout << std::endl;
//
//	sacado_test_tensor1 ();
//
//	std::cout << std::endl;

//	sacado_test_tensor2 ();
//
//    std::cout << std::endl;
//
//    sacado_test_tensor2B ();

//    std::cout << std::endl;
//
//    sacado_test_tensor3();

//    std::cout << std::endl;
//
//    sacado_test_tensor4();
//
//    std::cout << std::endl;
//
//    sacado_test_tensor5();

	  call_to_MatMod();
}
