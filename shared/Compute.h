#ifndef _COMPUTE_
#define _COMPUTE_

#include "Condition.h"
#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>      // std::stringstream

double Apply_Boundary_Condition()
{
	clock_t start_time, end_time;
	start_time = clock();

//	Hard_Wall(layer);
//	Soft_Wall(layer);
//	Ghost_Slip_Boundary_Condition();
//	Ghost_No_Slip_Boundary_Condition();
//	No_Slip_Boundary_Condition();
//	Slip_Boundary_Condition();
//	Tangential_Boundary_Condition();
//	Periodic_Boundary_Condition();
//	Vortex_Boundary_Condition();

//	Forward_Fourier_Transform();
	end_time = clock();
	double duration = round(time_digits*(end_time - start_time) / CLOCKS_PER_SEC) / time_digits;
	return (duration);
}

void Thresholding()
{
	for (int x = 0; x < Nx; x++)
		for (int y = 0; y < Ny; y++)
		{
			int id = index(x,y);
			if (rho.value[id] <= epsilone)
			{
				Wx.value[id] = 0;
				Wy.value[id] = 0;
			}
		}
}

void Init(double input_rho, int seed)
{
	gsl_rng_env_setup();
	gsl_T = gsl_rng_default;
	gsl_r = gsl_rng_alloc (gsl_T);
	gsl_rng_set(gsl_r, seed);

	Find_k();

	rho = input_rho;

//	Apply_Boundary_Condition();

//	Square_Vortex_Initial_Condition(input_rho, 1.0, M_PI/6.0,0.01);
//	Circular_Vortex_Initial_Condition(input_rho);
//	Horizontal_Band_Initial_Condition(input_rho, input_rho-0.2, 20, 40);
//	Vertical_Band_Initial_Condition(input_rho, input_rho-0.2, 40, 40);
//	Moving_Clump_Initial_Condition(input_rho, 0.8, 60, 60);
//	Random_Initial_Condition_With_Hole(input_rho, 0.5*input_rho,Nx_box/2);
	Polar_Initial_Condition(input_rho);
//	Random_Initial_Condition(input_rho);
//	Add_Noise(0.001);
//	Empty_Wall_Nearby(layer);
	ifstream in_data("snapshot-69600.dat");
	File_Initial_Condition(Nx_box, Ny_box, in_data);
//	Wx *= 0.05; Wy *= 0.05;

//	Apply_Boundary_Condition();
	Forward_Fourier_Transform();
	Cut(fWx, kcx, kcy);
	Cut(fWy, kcx, kcy);
	Cut(frho, kcx, kcy);
//	Backward_Fourier_Transform();
//	Compute_P();
//	Empty_Wall_Nearby(layer);
//	Forward_Fourier_Transform();
	Find_Evolution_Matrix();
}

double Compute_Nonlinear_Changes_GA()
{
	clock_t start_time, end_time;
	start_time = clock();

//	Non linear terms
	P2 = (Px*Px) + (Py*Py);

// Deriving term
	temp1 = (0.5*(1-alpha)*g)*rho;
	temp1 += (-0.5*(1-alpha)*g)*rho*(P2*P2);
	nonlinear_x = temp1*Wx;
	nonlinear_y = temp1*Wy;

// Polarity diffusion
	temp1 = (0.0625*(1-alpha)*g)*(rho);
	nonlinear_x += temp1*LWx;
	nonlinear_y += temp1*LWy;
//The removal of the following is correct when W << 1 (close to the transition). This is like bertin that W ~ epsilon, grad ~ epsilon. Then this term becomes epsilon^6.
	temp1 = (0.0625*(1-alpha)*g)*rho;
	temp2 = Px*LWx + Py*LWy;
	nonlinear_x += temp1*(P2 * (P2*LWx - 2.0*(Px*temp2)) );
	nonlinear_y += temp1*(P2 * (P2*LWy - 2.0*(Py*temp2)) );

// Polarity spread in higher order of nabla (nabla^4)
//	temp1 = ((1-alpha)*g / 384.0)*(rho);
//	nonlinear_x += temp1*L2Wx;
//	nonlinear_y += temp1*L2Wy;
//	temp1 = ((1-alpha)*g / 384.0)*(rho);
//	temp2 = Px*L2Wx + Py*L2Wy;
//	nonlinear_x += temp1*(P2 * (P2*L2Wx - 2.0*(Px*temp2)) );
//	nonlinear_y += temp1*(P2 * (P2*L2Wy - 2.0*(Py*temp2)) );

// Repulsion
	if (fabs(alpha) > 1e-10)
	{
		temp1 = (-alpha*g/6.0)*rho;
		nonlinear_x += temp1*drhodx;
		nonlinear_y += temp1*drhody;

		temp1 = (alpha*g/6.0)*(P2*rho);
		temp2 = Px*Px - Py*Py;
		temp3 = 2.0*(Px*Py);
		nonlinear_x += temp1*(temp3*drhody + temp2*drhodx);
		nonlinear_y += temp1*(temp3*drhodx - temp2*drhody);
	}

// Advection
	double a = 1; // Removing grad \rho P^4
	temp1 = (0.5*v*a)*(P2*P2);
	nonlinear_x += temp1*drhodx;
	nonlinear_y += temp1*drhody;

	temp1 = (2*v*a)*(P2*rho);
	nonlinear_x += temp1*(Px*dPxdx + Py*dPydx);
	nonlinear_y += temp1*(Px*dPxdy + Py*dPydy);

	double b = 1; // This is another pressure like term with the form of grad \rho P^2
	temp1 = (-v*b)*(P2*(Px*drhodx + Py*drhody) + 2.0*rho*(Px*Px*dPxdx + Py*Px*dPxdy + Px*Py*dPydx + Py*Py*dPydy));
	nonlinear_x += temp1*Px;
	nonlinear_y += temp1*Py;

	temp1 = (-v)*(rho*P2);
	temp2 = dPxdx + dPydy;
	nonlinear_x += temp1*(Px*temp2 + Px*dPxdx + Py*dPxdy);
	nonlinear_y += temp1*(Py*temp2 + Px*dPydx + Py*dPydy);

	Add_Noise(sqrt(dt)*T*1.0);

	fftw_execute(nonlinear_x_fwd);
	fftw_execute(nonlinear_y_fwd);
	fftw_execute(nonlinear_x_2_fwd);
	fftw_execute(nonlinear_y_2_fwd);
	Cut(fnonlinear_x, kcx, kcy);
	Cut(fnonlinear_y, kcx, kcy);
	Cut(fnonlinear_x_2, kcx, kcy);
	Cut(fnonlinear_y_2, kcx, kcy);

	fnonlinear_x = fnonlinear_x_2 + fnonlinear_x;
	fnonlinear_y = fnonlinear_y_2 + fnonlinear_y;

	end_time = clock();
	double duration = round(time_digits*(end_time - start_time) / CLOCKS_PER_SEC) / time_digits;
	return (duration);
}

double Compute_Nonlinear_Changes_Truncation()
{
	clock_t start_time, end_time;
	start_time = clock();

// Non linear terms

// Deriving terms
	temp1 = ((1-alpha)*g/2)*rho;
	temp1 += ((-(1-alpha)*(1-alpha)*g*g/(8*T))*(Wx*Wx + Wy*Wy));
	nonlinear_x = temp1*Wx;
	nonlinear_y = temp1*Wy;

// Repulsion
	if (fabs(alpha) > 1e-10)
	{
		temp1 = (((1-alpha)*alpha*g*g)/(12*T))*(Wx*drhodx + Wy*drhody);
		temp1 += (alpha*g*v/(24*T))*(Lrho);
		temp1 +=  (-alpha*alpha*g*g/(72*T))*(drhodx*drhodx + drhody*drhody);
		nonlinear_x += temp1*Wx;
		nonlinear_y += temp1*Wy;

		temp1 = (-alpha*g/6)*rho;
		nonlinear_x += (temp1*drhodx);
		nonlinear_y += (temp1*drhody);

		nonlinear_x += ((alpha*v*g)/(48*T))*(drhodx*dWxdx + drhody*dWxdy - 3.0*drhody*dWydx + 3.0*drhodx*dWydy);
		nonlinear_y += ((alpha*v*g)/(48*T))*(drhody*dWydy + drhodx*dWydx - 3.0*drhodx*dWxdy + 3.0*drhody*dWxdx);
	}


// Momentum flux spread.
	temp1 = (0.0625*(1-alpha)*g)*rho;
	nonlinear_x += temp1*(LWx);
	nonlinear_y += temp1*(LWy);
// Momentum flux spread. High epsilon order
//	temp1 -= ((1 - alpha)*(1 - alpha)*g*g/(32*T))*(Wx*(LWx) + Wy*(LWy));
//	nonlinear_x += temp1*Wx;
//	nonlinear_y += temp1*Wy;

// Advection With \grad W^2;
	nonlinear_x -= ((v*(1-alpha)*g)/(16*T))*(3.0*Wx*dWxdx + 3.0*Wy*dWxdy + 5.0*Wx*dWydy - 5.0*Wy*dWydx);
	nonlinear_y -= ((v*(1-alpha)*g)/(16*T))*(3.0*Wy*dWydy + 3.0*Wx*dWydx + 5.0*Wy*dWxdx - 5.0*Wx*dWxdy);
//	Without \grad W^2;
//	nonlinear_x -= ((v*(1-alpha)*g)/(16*T))*(8.0*Wx*dWxdx + 3.0*Wy*dWxdy + 5.0*Wx*dWydy);
//	nonlinear_y -= ((v*(1-alpha)*g)/(16*T))*(8.0*Wy*dWydy + 3.0*Wx*dWydx + 5.0*Wy*dWxdx);

	Add_Noise(sqrt(dt)*T*1.0);

	fftw_execute(nonlinear_x_fwd);
	fftw_execute(nonlinear_y_fwd);
	fftw_execute(nonlinear_x_2_fwd);
	fftw_execute(nonlinear_y_2_fwd);
	Cut(fnonlinear_x, kcx, kcy);
	Cut(fnonlinear_y, kcx, kcy);
	Cut(fnonlinear_x_2, kcx, kcy);
	Cut(fnonlinear_y_2, kcx, kcy);

	fnonlinear_x = fnonlinear_x_2 + fnonlinear_x;
	fnonlinear_y = fnonlinear_y_2 + fnonlinear_y;

	end_time = clock();
	double duration = round(time_digits*(end_time - start_time) / CLOCKS_PER_SEC) / time_digits;
	return (duration);
}

double Add_Non_Linear_Changes()
{
	fWx +=  (mu_n00*fnonlinear_x + mu_n01*fnonlinear_y);
	fWy +=  (mu_n10*fnonlinear_x + mu_n11*fnonlinear_y);
	frho += (mu_n20*fnonlinear_x + mu_n21*fnonlinear_y);
}

double Add_Linear_Changes()
{
	Cut(fWx, kcx, kcy);
	Cut(fWy, kcx, kcy);
	Cut(frho, kcx, kcy);

	temp_fWx = fWx;
	temp_fWy = fWy;
	temp_frho = frho;

	fWx =  (mu_l00*temp_fWx + mu_l01*temp_fWy + mu_l02*temp_frho);
	fWy =  (mu_l10*temp_fWx + mu_l11*temp_fWy + mu_l12*temp_frho);
	frho = (mu_l20*temp_fWx + mu_l21*temp_fWy + mu_l22*temp_frho);

	Cut(fWx, kcx, kcy);
	Cut(fWy, kcx, kcy);
	Cut(frho, kcx, kcy);
}

double Compute_Changes()
{
	clock_t start_time, end_time;
	start_time = clock();

	if (compute_mode == ga)
		Compute_Nonlinear_Changes_GA();
	else
		Compute_Nonlinear_Changes_Truncation();
// Linear implicit stepping
	Add_Non_Linear_Changes();
//	fAdd_Noise(sqrt(dt)*T*1.0);
	Add_Linear_Changes();

	end_time = clock();
	double duration = round(time_digits*(end_time - start_time) / CLOCKS_PER_SEC) / time_digits;
	return (duration);
}


void One_Step(double t[5])
{
	t[1] = Backward_Fourier_Transform();
//	t[5] = Apply_Boundary_Condition();
//	Thresholding();
	t[3] = Compute_Derivatives_Shared();
	if (compute_mode == ga)
	{
		t[2] = Compute_P();
		t[3] += Compute_Derivatives_GA();
	}
	else
	{
		t[3] += Compute_Derivatives_Truncation();
	}

	t[4] = Compute_Changes();
}

void Write_To_File(bool saving = false)
{
	stringstream address;
	address.str("");
	if (saving)
		address << "rho=" << initial_rho << "-alpha=" << alpha << "-v=" << v << "-T=" << T << ".dat";
	else
		address << "data.dat";
	ofstream data_file(address.str().c_str());

	data_file << "#x y Wx Wy rho curl" << endl;
	for (int i = 0; i < Nx_box; i++)
	{
		for (int j = 0; j < Ny_box; j++)
			data_file << (1.0*i / (Nx_box-1)) - 0.5 << "\t" << (1.0*j / (Ny_box-1)) - 0.5 << "\t" << Wx_2.value[i*Ny+j] << "\t" << Wy_2.value[i*Ny+j] << "\t" << rho_2.value[i*Ny+j] << "\t" << endl;
		data_file << endl;
	}

	data_file.close();
}

void Save_Snapshot(double& polarization)
{
	static long int frame = 0;
	stringstream address;
	address.str("");
	address << "snapshot/snapshot-";
	if (frame < 10)
		address << "0000" << frame;
	else
		if (frame < 100)
			address << "000" << frame;
		else
			if (frame < 1000)
				address << "00" << frame;
			else
				if (frame < 10000)
					address << "0" << frame;
			else
				address << frame;
	address << ".dat";
	ofstream snapshot_file(address.str().c_str());
	snapshot_file << "# Dr\t" << T << "\tK\t" << K << "\t" << frame << "\t" << lapsed_time << "\tpolarization\t" << polarization << endl;

	for (int i = 0; i < Nx_box; i++)
	{
		for (int j = 0; j < Ny_box; j++)
			snapshot_file << (1.0*i / (Nx_box-1)) - 0.5 << "\t" << (1.0*j / (Ny_box-1)) - 0.5 << "\t" << Wx_2.value[i*Ny+j] / rho_2.value[i*Ny+j] << "\t" << Wy_2.value[i*Ny+j] / rho_2.value[i*Ny+j] << "\t" << rho_2.value[i*Ny+j] << "\t" << endl;
		snapshot_file << endl;
	}

	snapshot_file.close();
	frame++;
}

double Rho_Mean()
{
	double result = rho_2.Get_Average();
//	result = creal(frho.value[index(0,0)]);
//	result /= Nx*Ny;
	return(result);
}


void Equilibrium(long int integration_duration, int snapshot_period)
{
	clock_t start_time, end_time;
	start_time = clock();
	double t[5];

	Find_Evolution_Matrix();
	long int number_of_step = round(integration_duration / dt);
	for (long int i = 0; i < number_of_step; i++)
	{
			One_Step(t);
			if (i % snapshot_period == 0)
			{
				lapsed_time += snapshot_period*dt;
				end_time = clock();
				Write_To_File(false);
				double polarization = sqrt((creal(fWx.value[0])*creal(fWx.value[0]) + creal(fWy.value[0])*creal(fWy.value[0]))/(N*N)) / initial_rho;
				Save_Snapshot(polarization);
//				cout <<  "step: " << i << " fraction:\t" << round(i*1000.0 / number_of_step)/10 << "%\t" << round((100.0*(end_time - start_time)) / CLOCKS_PER_SEC)/100 << "\t" << t[1] << "\t" << t[2] << "\t" << t[3] << "\t" << t[4] << "\t" << t[5] << "\tp: " << polarization << "\r" << flush;
				cout <<  "step: " << i << " fraction:\t" << round(i*1000.0 / number_of_step)/10 << "%\t" << round((100.0*(end_time - start_time)) / CLOCKS_PER_SEC)/100 << "\tp: " << polarization << "\r" << flush;
				start_time = clock();
				if (i % (snapshot_period*10) == 0)
				{
					Truncate_Field(fWx);
					Truncate_Field(fWy);
					Truncate_Field(frho);
				}
			}
	}
	end_time = clock();
	float tt = (float) (end_time - start_time) / CLOCKS_PER_SEC;
	double mean_rho = Rho_Mean();
	cout << "\nFinal average density is " << mean_rho << "\tin " << tt << " s" << "\t total clocks:" << end_time - start_time << "\n" << endl;
	if (mean_rho != mean_rho)
	{
		cout << "Nan values!." << endl;
		cout << "Exiting the program" << endl;
		exit(0);
	}
//	Write_To_File(true);
}


void Compute(long int integration_duration, int snapshot_period, int saving_period, long double& mean_p, long double& err)
{
	clock_t start_time, end_time;
	start_time = clock();
	double t[5];

	mean_p = 0;
	long double mean_p2 = 2;

	Find_Evolution_Matrix();
	long int number_of_step = round(integration_duration / dt);
	for (long int i = 0; i < number_of_step; i++)
	{
			One_Step(t);
			if (i % snapshot_period == 0)
			{
				lapsed_time += snapshot_period*dt;
				end_time = clock();
				Write_To_File(false);
				double polarization = sqrt((creal(fWx.value[0])*creal(fWx.value[0]) + creal(fWy.value[0])*creal(fWy.value[0]))/(N*N)) / initial_rho;
				Save_Snapshot(polarization);
				cout <<  "step: " << i << " fraction:\t" << round(i*1000.0 / number_of_step)/10 << "%\t" << round((100.0*(end_time - start_time)) / CLOCKS_PER_SEC)/100 << "\tp: " << polarization << "\r" << flush;
				start_time = clock();
				if (i % (snapshot_period*10) == 0)
				{
					Truncate_Field(fWx);
					Truncate_Field(fWy);
					Truncate_Field(frho);
				}
			}
			if (i % saving_period == 0)
			{
				double polarization = sqrt((creal(fWx.value[0])*creal(fWx.value[0]) + creal(fWy.value[0])*creal(fWy.value[0]))/(N*N)) / initial_rho;
				mean_p += polarization;
				mean_p2 += polarization*polarization;
			}
	}

	int n = number_of_step / saving_period;
	mean_p /= n;
	mean_p2 /= n;
	err = sqrt(mean_p2 - mean_p*mean_p);
	err /= sqrt(n/10);

	end_time = clock();
	float tt = (float) (end_time - start_time) / CLOCKS_PER_SEC;
	double mean_rho = Rho_Mean();
	cout << "\nFinal average density is " << mean_rho << "\tin " << tt << " s" << "\t total clocks:" << end_time - start_time << "\n" << endl;

	Write_To_File();
}


#endif
