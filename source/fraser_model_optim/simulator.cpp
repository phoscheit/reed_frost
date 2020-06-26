#include <vector>
#include <random>
#include <cmath>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>

std::random_device rd;
std::mt19937 rng(rd());	

int reed_frost(int n, int s0, double Q, double beta, double alpha, double k)
{
	std::binomial_distribution<int> outside(s0,1-Q);

	double gamma_mean = beta/pow(n,alpha);
	double gamma_shape = k;
	std::gamma_distribution<double> hazard(gamma_shape,
		gamma_mean/gamma_shape);

	int infected = outside(rd);
	int total_infected = infected;
	int susceptibles = s0-total_infected;

	while(infected!=0)
	{
		double cum_hazard=0.;
		for (int i = 0; i < infected; ++i)
		{
			double h = hazard(rd);
			cum_hazard += h;
		}
		double trans_prob = std::exp(-cum_hazard);
		std::binomial_distribution<int> next_gen(susceptibles,1-trans_prob);
		infected = next_gen(rd);
		total_infected += infected;
		susceptibles -= infected;
	}
	return total_infected;
}

int single_household(int n,std::vector<double> &params)
{

	/* Les paramètres sont dans l'ordre (cf Fraser et al. 2011) :

	- La proba d'échappement Q
	- La proba d'infection asymptomatique et infectieuse p_asx
	- La proba d'infection asymptomatique et non-infectieuse p_r
	- La proba d'échappement à une infection antérieure Q_prior
	- Hétérogénéité des infectivités k
	- Facteur de taille du domicile beta
	- Exposant de la taille du domicile alpha

	*/

	std::binomial_distribution<int> protection(n,params[2]);

	int prot = protection(rd);
	int susceptibles = n-prot;

	int prev_infected = reed_frost(n,susceptibles,params[3],params[5],params
		[6],params[4]);

	susceptibles -= prev_infected;

	int new_infected = reed_frost(n,susceptibles,params[0],params[5],params
		[6],params[4]);

	std::binomial_distribution<int> asympt(new_infected,params[1]);
	int asymptomatics = asympt(rd);
	int symptomatics = new_infected - asymptomatics;
	return symptomatics;
}


int main(int argc, char const *argv[])
{

	std::vector<double> test_params { 0.8,0,0,0.88,0.94,0.37,0.35 } ;
	std::vector<double> params=test_params;

	std::vector<int> test_numbers {244,1126,1342,1378,1017,680,434,
	  284,117,74,25,23,5,3,0,0,1 };

	int max_n = test_numbers.size();

	Eigen::MatrixXi table(max_n+1,max_n);

	for (int n = 0; n < max_n; ++n)
	{
		std::vector<int> symptomatics;
		int num_households = test_numbers[n];
		for (int i = 0; i < num_households; ++i)
		{
			int infected = single_household(n+1,params);
			symptomatics.push_back(infected);
		}

		for (int m = 0; m <= n+1; ++m)
		{
			int occurrences = std::count(symptomatics.begin(),
				symptomatics.end(),m);
			table(m,n) = occurrences;
		}
	}

	std::cout << table << std::endl;

	std::ofstream outfile("sim_matrix.txt");
	outfile << table;

	return 0;
}