#include "Condition.h"
#include "Compute.h"




int main(int argc, char *argv[])
{
	std::cout.unsetf ( std::ios::floatfield );                // floatfield not set
	std::cout.precision(3);

	string input_str = argv[1];
	if (input_str.compare("ga") == 0)
	{
		compute_mode = ga;
		cout << "\nCompute mode is Gaussian approximation.\n" << endl;
	}
	else
		if (input_str.compare("tr") == 0)
		{
			compute_mode = tr;
			cout << "\nCompute mode is truncation approximation.\n" << endl;
		}
		else
		{
			cout << "\nBad compute mode. The first argument should be either tr or ga\n" << endl;
			exit(0);
		}

	initial_rho = atof(argv[2]);
	g = atof(argv[3]);
	alpha = atof(argv[4]);
	int T_size = (argc - 5)/2;
	float T_list[50];
	long int equilibrium_time[50];
	for (int i = 0; i < T_size; i++)
	{
		T_list[i] = atof(argv[2*i+5]);
		equilibrium_time[i] = atoi(argv[2*i+6]);
	}
	E = v / 2;

	stringstream address;
	address.str("");
	if (compute_mode == tr)
		address << "tr-";
	else
		address << "ga-";
	if (T_list[0] > 0.5)
		address << "cooling-";
	else
		address << "heating-";
	address << "rho=" << initial_rho << "-alpha=" << alpha << "-v=" << v << "-Dr=" << T << "Lx=" << Lx << "Ly=" << Ly << "Nx=" << Nx << "Ny=" << Ny << ".dat";
	ofstream outfile;
	outfile.open(address.str().c_str());

	Init(initial_rho, time(NULL));
//	Write_To_File();

	int interval = default_interval;
	for (int i = 0; i < T_size; i++)
	{
		T = T_list[i];
//		if (fabs(T-0.2) < 0.00001)
//			K = D = 5;
//		if (compute_mode == tr)
//			D = (v*v) / (16*T) + K;
//		else
//			D = K;
		cout << "\nrho=" << initial_rho << "\t g=" << g << "\t appha=" << alpha << "\t T=" << T_list[i] << "\t K=" << K << "\t D=" << D << endl;
		interval = default_interval;
//		if (fabs(T-2.36) < 0.00001)
//			interval = 128;
//		if (fabs(T-0.2) < 0.00001)
//			interval = 256;
//		else
//			interval = default_interval;
//		if (fabs(T-2.38) < 0.00001)
//			interval = 128;
		cout << "interval:\t" << interval << endl;
		dt = 1.0 / interval;
		if ((fabs(T-0.2) < 0.00001) && (T_list[0] > 0.5))
		{
//			Equilibrium(equilibrium_time[i], 4*interval);
			K = D = 5;
			Equilibrium(equilibrium_time[i], 512*interval);
			K = D = 0.125;
		}
		Equilibrium(equilibrium_time[i], 512*interval);
		long double polarization,err;
		Compute(equilibrium_time[i]/16, 512*interval, interval, polarization, err);
		outfile << T << "\t" << polarization << "\t" << err << endl;
	}
	outfile.close();
	Destroy_FFT_Plans();
}

