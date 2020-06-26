#include <vector>
#include <tuple>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <nlopt.hpp>
#include <Eigen/Dense>
#include <boost/math/special_functions/binomial.hpp>


#define MAXBUFSIZE  ((int) 1e6)

unsigned int function_calls =0;

std::vector< std::pair<int,std::vector<double> > > read_prior(const char*
	filename)
	{
		std::ifstream infile;
		infile.open(filename);

		std::vector< std::pair<int,std::vector<double> > > acc;

		while(!infile.eof())
		{
			std::string line;
			std::getline(infile,line);

			std::stringstream stream(line);
			int size_household;
			stream >> size_household;

			std::vector<double> priors;
			double probab;
			while(stream >> probab)
			{
				priors.push_back(probab);
			}
			acc.push_back(std::make_pair(size_household,priors));
		}
		return acc;
	}

Eigen::MatrixXd readMatrix(const char *filename)
    {
    int cols = 0, rows = 0;
    double buff[MAXBUFSIZE];

    // Read numbers from file into buffer.
    std::ifstream infile;
    infile.open(filename);
    while (! infile.eof())
        {
        std::string line;
        std::getline(infile, line);

        int temp_cols = 0;
        std::stringstream stream(line);
        while(! stream.eof())
            stream >> buff[cols*rows+temp_cols++];

        if (temp_cols == 0)
            continue;

        if (cols == 0)
            cols = temp_cols;

        rows++;
        }

    infile.close();

    // rows--;

    // Populate matrix with numbers.
    Eigen::MatrixXd result(rows,cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result(i,j) = buff[ cols*i+j ];

    return result;
    }

double Phi(const double &x, const double &alpha, const double &beta, const
	double &k, const int n)
{
	double mean = beta/pow(n,alpha);
	return pow(k/(k+x*mean),k);
}

/* Returns a (n+1)*(n+1) matrix such that F(i,j) = F_i^{n,j}. Columns correspond
to different values of s0, rows are the m values. As such, F is upper 
triangular, and F(i,j) for i>j is undefined. */

std::pair<Eigen::MatrixXd,Eigen::MatrixXd> solve_triangular(const
	std::vector<double> &v, const int n)
	{
		Eigen::MatrixXd result1(n+1,n+1);
		Eigen::MatrixXd result2(n+1,n+1);		
		for (int s0=0; s0 <= n; ++s0)
		{
			Eigen::MatrixXd A1(s0+1,s0+1);
			Eigen::MatrixXd A2(s0+1,s0+1);			
			Eigen::VectorXd B(s0+1);
			for(int i=0; i<=s0; ++i)
			{
				B(i) = boost::math::binomial_coefficient<double>(s0,i);
				for (int j = 0; j <= i; ++j)
				{
					A1(i,j) = boost::math::binomial_coefficient<double>
					(s0-j,i-j)
					/(pow(Phi(s0-i,v[6],v[5],v[4],n),j)*(pow(v[0],s0-i)));
					A2(i,j) = boost::math::binomial_coefficient<double>
					(s0-j,i-j)
					/(pow(Phi(s0-i,v[6],v[5],v[4],n),j)*(pow(v[3],s0-i)));
				}
			}
			Eigen::VectorXd sol1 = A1.triangularView<Eigen::Lower>().solve(B);
			sol1.conservativeResize(n+1);
			result1.col(s0)=sol1;
			Eigen::VectorXd sol2 = A2.triangularView<Eigen::Lower>().solve(B);
			sol2.conservativeResize(n+1);
			result2.col(s0)=sol2;			
		}

		for (int i = 0; i <= n; ++i)
		{
			for (int j = 0; j < i; ++j)
			{
				result1(i,j)=1;
				result2(i,j)=0;
			}
		}
 		return std::pair<Eigen::MatrixXd,Eigen::MatrixXd> {result1, result2};
	}

/* Returns the value of T_{(m,n)}, the expected frequency of households of size 
n with m infected individuals. Arguments are a tuple (n,m), the running 
parameter vector and the precomputed matrix containing the F_m^{n,s0}. */
double expected_prob(const std::tuple<int,int> &sizes, const
	std::vector<double> &parameters, const
	std::pair<Eigen::MatrixXd,Eigen::MatrixXd> &distributions)
	{
	const double n = std::get<0>(sizes); // n \ge 1
	const double m = std::get<1>(sizes); // 0 \le m \le s0

	double pasx = parameters[1];
	double ppr = parameters[2];

	Eigen::MatrixXd dis_current = distributions.first;
	Eigen::MatrixXd dis_prior = distributions.second;

	double sum=0.;

	for(int t=0; t<=(n-m); ++t)
	{
		for(int r=0; r<=(n-m-t); ++r)
		{
			for (int l = 0; l<= (n-m-r-t); ++l)
			{
				double prob1 = boost::math::binomial_coefficient<double>
				(m+t,t)*pow(pasx,t)*pow(1-pasx,m);
				double prob2 = boost::math::binomial_coefficient<double>
				(n,r)*pow(ppr,r)*pow(1-ppr,n-r);
				if(m+t>n-r-l) std::cout << "n=" << n << " m=" << m << " t=" << t
					<< " r=" << r << " l=" << l << std::endl;
				double expected_prob1 = dis_current(m+t,n-r-l);
				double expected_prob2 = dis_prior(l,n-r);				
				sum += prob1*prob2*expected_prob1*expected_prob2;
			}
		}
	}
	return sum;
}	

double Dev(const std::vector<double> &v, std::vector<double> &grad, void* 
	my_func_data)
	{
		if(++function_calls%1000 ==0)
		{
			std::cout << "Function calls: " << function_calls ;
			for(auto values:v)
			{
				std::cout << " Valeurs : " << values;
			}
			std::cout << std::endl;
		}
		/* Contains the observed final sizes k_{(m,n)}. Be careful : 
		households(i,j) contains k_{(i,j+1)}. Also, households(i,j)=0 
		if i>j. */
		const Eigen::MatrixXd *households =
			static_cast<Eigen::MatrixXd *>(my_func_data); 

		int max_n = (*households).cols();

		std::vector< std::pair<Eigen::MatrixXd,Eigen::MatrixXd > >
		distribution_matrices;

		/*	Computes the expected final size distributions F(n; s0,m) for 1
		\le n \le n_max and m\le s0 \le n. Output is a vector whose (n+1)-th
		element contains the matrix of the F(n;s0,m) for a fixed n. */

		for (int n = 1; n <= max_n; ++n)
		{
			std::pair<Eigen::MatrixXd,Eigen::MatrixXd> temp = solve_triangular
			(v,n);
			distribution_matrices.push_back(temp);
		}

		double sum = 0.;

		auto colSums = (*households).colwise().sum();

		for (int n = 1; n <= max_n; ++n)
		{
			double col_sum = colSums(n-1);
			for(int m=0; m<=n; ++m)
			{
				double observed = (*households)(m,n-1);
				if(observed==0.) continue;
				std::tuple<int,int> obs_tuple {n,m};
				double expected = expected_prob(obs_tuple, v,
					distribution_matrices[n-1]);
				sum += observed*(std::log10(observed/col_sum)-
					std::log10(expected));
			}
		}
		if(function_calls%1000 ==0) std::cout << sum << std::endl;
		return sum; 
	}

int main(int argc, char const *argv[])
{
	auto priors = read_prior("test_prior.txt");

	std::vector<Eigen::MatrixXd> households; /* Contiendra les valeurs de 
	tailles finales */

	nlopt::opt opt(nlopt::LN_BOBYQA, 7);

	void* f_data;

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

	std::vector<double> lb { 0.,0.,0.,0.,0.,0.,0. };
	std::vector<double> ub { 1.,1.,1.,1.,5.,5.,5. };
	opt.set_lower_bounds(lb);
	opt.set_upper_bounds(ub);

	std::vector<std::string> var_names
		{"Q","p_asx","p_pr","Q_prior","k","beta","alpha"};
	std::vector<double> init_guess { 0.5,0.5,0.5,0.5,1.,1.,1. } ;

	Eigen::MatrixXd test_matrix = readMatrix("test_matrix.txt");

	void* my_func_data = static_cast<void *>(&test_matrix);

	std::cout << "Optimisation en cours..." << std::endl;

	opt.set_min_objective(Dev,my_func_data);
	double minf =0.;

	try
	{
		nlopt::result result = opt.optimize(init_guess,minf);
	}	catch(const nlopt::roundoff_limited& e)
	{
		for(int i=0; i<7; ++i)
		{
			std::cout << var_names[i] << " : " << init_guess[i] << std::endl;
		}
		std::cout << "Vraisemblance : " << minf << std::endl;
		std::cout << "Nombre d'appels à Dev : " << function_calls << std::endl;
	}



	return 0;
}