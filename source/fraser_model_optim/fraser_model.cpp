#include <vector>
#include <tuple>
#include <iostream>
#include <cmath>
#include <nlopt.hpp>
#include <Eigen/Dense>
#include <boost/math/special_functions/binomial.hpp>

double Phi(const double &x, const double &alpha, const double &beta, const
	double &k, const int n)
{
	double mean = beta/pow(n,alpha);
	return pow(k/(k+x*mean),k);
}

/* Returns a (n+1)*(n+1) matrix such that F(i,j) = F_i^{n,j}. Columns correspond
to different values of s0, rows are the m values. As such, F is upper 
triangular, and F(i,j) for j>i is undefined. */

Eigen::MatrixXd solve_triangular(const std::vector<double> &v,
	const int n)
	{
		Eigen::MatrixXd result(n+1,n+1);
		for (int s0=0; s0 <= n; ++s0)
		{
			Eigen::MatrixXd A(s0+1,s0+1);
			Eigen::VectorXd B(s0+1);
			for(int i=0; i<=s0; ++i)
			{
				B(i) = boost::math::binomial_coefficient<double>(s0,i);
				for (int j = 0; j <= i; ++j)
				{
					A(i,j) = boost::math::binomial_coefficient<double>
					(s0-j,i-j)
					/(Phi(s0-i,v[6],v[5],v[4],n)*(pow(v[0],s0-i)));
				}
			}
			Eigen::VectorXd sol = A.triangularView<Eigen::Lower>().solve(B);
			sol.conservativeResize(n+1);
			result.col(s0)=sol;
		}

		for (int i = 0; i <= n; ++i)
		{
			for (int j = 0; j < i; ++j)
			{
				result(i,j)=0;
			}
		}
 		return result;
	}

/* Returns the value of T_{(m,n)}, the expected frequency of households of size 
n with m infected individuals. Arguments are a tuple (n,m), the running 
parameter vector and the precomputed matrix containing the F_m^{n,s0}. */
double expected_prob(const std::tuple<int,int> &sizes, const
	std::vector<double> &parameters, const Eigen::MatrixXd &distributions)
	{
	const double n = std::get<0>(sizes); // n \ge 1
	const double m = std::get<1>(sizes); // 0 \le m \le s0

	double pasx = parameters[1];
	double ppr = parameters[2];

	double sum=0.;

	for(int t=0; t<=(n-m); ++t)
	{
		for(int r=0; r<=(n-m-t); ++r)
		{
			for (int l = 0; l<= (n-m-r); ++l)
			{
				double prob1 = boost::math::binomial_coefficient<double>
				(m+t,t)*pow(pasx,t)*pow(1-pasx,m);
				double prob2 = boost::math::binomial_coefficient<double>
				(n,r)*pow(ppr,r)*pow(1-ppr,n-r);
				double expected_prob1 = distributions(m+t,n-r-l);
				double expected_prob2 = distributions(l,n-r);
				sum += prob1*prob2*expected_prob1*expected_prob2;
			}
		}
	}
	return sum;
}	

double Dev(const std::vector<double> &v, std::vector<double> &grad, void* 
	my_func_data)
	{
		/* Contains the observed final sizes k_{(m,n)}. Be careful : 
		households(i,j) contains k_{(i,j+1)}. Also, households(i,j)=0 
		if i>j. */
		const Eigen::MatrixXd *households =
			static_cast<Eigen::MatrixXd *>(my_func_data); 

		int max_n = (*households).cols();

		/*	Computes the expected final size distributions F(n; s0,m) for n \ge
		1, m\le s0 \le n. Output is a vector whose (n+1)-th element contains the
		matrix of the F(n;s0,m) for a fixed n. */
		Eigen::MatrixXd F = solve_triangular(v, max_n); 


		double sum = 0.;



		return sum; 
	}

int main(int argc, char const *argv[])
{
	std::cout << "Début du programme" << std::endl << std::flush;

	std::vector<Eigen::MatrixXd> households; /* Contiendra les valeurs de 
	tailles finales */
	
	std::cout << "Définition de l'objet" << std::endl << std::flush;

	nlopt::opt opt(nlopt::LN_BOBYQA, 7);

	void* f_data;

	std::cout << "Définition de l'objectif" << std::endl << std::flush;

	opt.set_min_objective(Dev,f_data); // Définit Dev comme fonction à minimiser

	/* Les paramètres sont dans l'ordre (cf Fraser et al. 2011) :

	- La proba d'échappement Q
	- La proba d'infection asymptomatique et infectieuse p_asx
	- La proba d'infection asymptomatique et non-infectieuse p_r
	- La proba d'échappement à une infection antérieure Q_prior
	- Hétérogénéité des infectivités k
	- Facteur de taille du domicile beta
	- Exposant de la taille du domicile alpha

	*/

	std::cout << "Définition des bornes" << std::endl;

	std::vector<double> lb(7,0.);
	std::vector<double> ub { 1.,1.,1.,1.,HUGE_VAL,HUGE_VAL,HUGE_VAL };
	opt.set_lower_bounds(lb);
	opt.set_upper_bounds(ub);

	std::cout << "Estimation initiale" << std::endl;

	std::vector<double> init_guess { 0.5,0.5,0.5,0.5,1,1,1 } ;

	std::cout << "Test" << std::endl <<std::flush;

	// opt.optimize(init_guess);

	std::cout << "Valeurs : " << opt.last_optimum_value() << std::endl; 

	return 0;
}