#ifndef _FFT_
#define _FFT_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <complex.h>
#include <fftw3.h>
#include "CField.h"
#include "Parameters.h"

void Find_k()
{
	for (int i = 0; i < Nkx; i++)
	{
			if (i < (Nkx / 2))
				kx[i] = (M_PI*i/Lx);
			else
				kx[i] = (M_PI*(i - Nkx)/Lx);
	}
	for (int j = 0; j < Nky; j++)
		ky[j] = (M_PI*(j) / Ly);
}

void Cut(Field<fftw_complex, Nkx, Nky>& f, int ikcx, int ikcy)
{
//	int ikc = 2*Nky/3;
//	ikc = Nky - 1;
//	ikc = Nkx;

	for (int i = ikcx; i <= (Nkx - ikcx); i++)
		for (int j = 0; j < Nky; j++)
			f.value[i*Nky+j] = 0 + 0*I;

	for (int i = 0; i < Nkx; i++)
		for (int j = ikcy; j < Nky; j++)
			f.value[i*Nky+j] = 0 + 0*I;
}

void Derivative_x(const Field<fftw_complex, Nkx, Nky>& f, Field<fftw_complex, Nkx, Nky>& dfdx)
{
		for (int i = 0; i < Nkx; i++)
			for (int j = 0; j < Nky; j++)
					dfdx.value[i*Nky+j] = I*kx[i]*f.value[i*Nky+j];
		for (int j = 0; j < Nky; j++)
			dfdx.value[(Nkx/2)*Nky+j] = 0;
}

void Derivative_y(const Field<fftw_complex, Nkx, Nky>& f, Field<fftw_complex, Nkx, Nky>& dfdy)
{
		for (int i = 0; i < Nkx; i++)
			for (int j = 0; j < Nky; j++)
				dfdy.value[i*Nky+j] = I*ky[j]*f.value[i*Nky+j];
}

void Derivative_x2(const Field<fftw_complex, Nkx, Nky>& f, Field<fftw_complex, Nkx, Nky>& dfdx)
{
		for (int i = 0; i < Nkx; i++)
			for (int j = 0; j < Nky; j++)
					dfdx.value[i*Nky+j] = -(kx[i]*kx[i])*f.value[i*Nky+j];
		for (int j = 0; j < Nky; j++)
			dfdx.value[(Nkx/2)*Nky+j] = 0;
}

void Derivative_y2(const Field<fftw_complex, Nkx, Nky>& f, Field<fftw_complex, Nkx, Nky>& dfdy)
{
		for (int i = 0; i < Nkx; i++)
			for (int j = 0; j < Nky; j++)
				dfdy.value[i*Nky+j] = -(ky[j]*ky[j])*f.value[i*Nky+j];
}

void Laplacian(const Field<fftw_complex, Nkx, Nky>& f, Field<fftw_complex, Nkx, Nky>& Lf)
{
		for (int i = 0; i < Nkx; i++)
			for (int j = 0; j < Nky; j++)
				if (i != Nkx/2)
					Lf.value[i*Nky+j] = -(kx[i]*kx[i] + ky[j]*ky[j])*f.value[i*Nky+j];
				else
					Lf.value[i*Nky+j] = -(ky[j]*ky[j])*f.value[i*Nky+j];
}

void Laplacian_Square(const Field<fftw_complex, Nkx, Nky>& f, Field<fftw_complex, Nkx, Nky>& Lf)
{
		for (int i = 0; i < Nkx; i++)
			for (int j = 0; j < Nky; j++)
				if (i != Nkx/2)
					Lf.value[i*Nky+j] = (kx[i]*kx[i] + ky[j]*ky[j])*(kx[i]*kx[i] + ky[j]*ky[j])*f.value[i*Nky+j];
				else
					Lf.value[i*Nky+j] = (ky[j]*ky[j])*(ky[j]*ky[j])*f.value[i*Nky+j];
}

void Truncate_Field(const Field<fftw_complex, Nkx, Nky>& f)
{
	for (int i = 0; i < Nkx; i++)
	{
		for (int j = 0; j < Nky; j++)
		{
			if (cabs(f.value[i*Nky+j]) < 1e-25)
				f.value[i*Nky+j] = 0;
		}
	}
}

/////////////////////////////// Definition of our fields

Field<Real, Nx, Ny> Wx, Wy, rho, Wx_2, Wy_2, rho_2;
Field<Real, Nx, Ny> drhodx, drhody;
Field<Real, Nx, Ny> dWxdx, dWxdy, dWydx, dWydy; // dWdr is for Bertin
Field<Real, Nx, Ny> nonlinear_x, nonlinear_y, nonlinear_x_2, nonlinear_y_2;
Field<Real, Nx, Ny> temp1, temp2, temp3, Px, Py, Px_2, Py_2, P2, dPxdx, dPxdy, dPydx, dPydy, laplacian_Wx, laplacian_Wy; // this is a temporarily field for computation of nonlinear terms. All p's are for GA
Field<Real, Nx, Ny> Mean_Px, Mean_Py, Mean_rho;
Field<fftw_complex, Nkx, Nky> fWx,fWy,frho, fPx, fPy; // fPx and fPy are for GA
Field<fftw_complex, Nkx, Nky> fdrhodx, fdrhody, fdWxdx, fdWxdy, fdWydx, fdWydy, fdPxdx, fdPxdy, fdPydx, fdPydy; // Again p's are for GA
Field<fftw_complex, Nkx, Nky> fnonlinear_x, fnonlinear_y, fnonlinear_x_2, fnonlinear_y_2;
Field<fftw_complex, Nkx, Nky> ftemp, dftempxdx, dftempxdy, dftempydx, dftempydy, temp_frho, temp_fWx, temp_fWy, temp_fPx, temp_fPy; // These temporarly fields is required for 1) Backward fourier transformation to avoid data loss 2) for stepping of the fields. The temp_fPx, and temp_fPy are for GA

Field<fftw_complex, Nkx, Nky> mu_l00, mu_l01, mu_l02, mu_l10, mu_l11, mu_l12, mu_l20, mu_l21, mu_l22, mu_n00, mu_n01, mu_n02, mu_n10, mu_n11, mu_n12, mu_n20, mu_n21, mu_n22;

Field<Real, Nx, Ny> L2Wx, L2Wy, LWx, LWy, Lrho; // Laplacian and Higher order derivatives (Laplacian^2)
Field<fftw_complex, Nkx, Nky> fL2Wx, fL2Wy, fLWx, fLWy, fLrho; // Laplacian and Higher order derivatives (Laplacian^2)

/////////////////////////////// Fourier Plan deffination

fftw_plan rho_fwd = fftw_plan_dft_r2c_2d(Nx, Ny, rho.value, frho.value,  algorithm);
fftw_plan Wx_fwd = fftw_plan_dft_r2c_2d(Nx, Ny, Wx.value, fWx.value,  algorithm);
fftw_plan Wy_fwd = fftw_plan_dft_r2c_2d(Nx, Ny, Wy.value, fWy.value,  algorithm);
fftw_plan Px_fwd = fftw_plan_dft_r2c_2d(Nx, Ny, Px.value, fPx.value,  algorithm);
fftw_plan Py_fwd = fftw_plan_dft_r2c_2d(Nx, Ny, Py.value, fPy.value,  algorithm);

fftw_plan nonlinear_x_2_fwd = fftw_plan_dft_r2c_2d(Nx, Ny, nonlinear_x_2.value, fnonlinear_x_2.value,  algorithm);
fftw_plan nonlinear_y_2_fwd = fftw_plan_dft_r2c_2d(Nx, Ny, nonlinear_y_2.value, fnonlinear_y_2.value,  algorithm);
fftw_plan nonlinear_x_fwd = fftw_plan_dft_r2c_2d(Nx, Ny, nonlinear_x.value, fnonlinear_x.value,  algorithm);
fftw_plan nonlinear_y_fwd = fftw_plan_dft_r2c_2d(Nx, Ny, nonlinear_y.value, fnonlinear_y.value,  algorithm);

fftw_plan rho_2_bwd = fftw_plan_dft_c2r_2d(Nx, Ny, temp_frho.value, rho_2.value,  algorithm);
fftw_plan Wx_2_bwd = fftw_plan_dft_c2r_2d(Nx, Ny, temp_fWx.value, Wx_2.value,  algorithm);
fftw_plan Wy_2_bwd = fftw_plan_dft_c2r_2d(Nx, Ny, temp_fWy.value, Wy_2.value,  algorithm);
fftw_plan Px_2_bwd = fftw_plan_dft_c2r_2d(Nx, Ny, temp_fPx.value, Px_2.value,  algorithm); // GA
fftw_plan Py_2_bwd = fftw_plan_dft_c2r_2d(Nx, Ny, temp_fPy.value, Py_2.value,  algorithm); // GA

fftw_plan rho_bwd = fftw_plan_dft_c2r_2d(Nx, Ny, temp_frho.value, rho.value,  algorithm);
fftw_plan Wx_bwd = fftw_plan_dft_c2r_2d(Nx, Ny, temp_fWx.value, Wx.value,  algorithm);
fftw_plan Wy_bwd = fftw_plan_dft_c2r_2d(Nx, Ny, temp_fWy.value, Wy.value,  algorithm);
fftw_plan Px_bwd = fftw_plan_dft_c2r_2d(Nx, Ny, temp_fPx.value, Px.value,  algorithm); // GA
fftw_plan Py_bwd = fftw_plan_dft_c2r_2d(Nx, Ny, temp_fPy.value, Py.value,  algorithm); // GA

fftw_plan drhodx_bwd = fftw_plan_dft_c2r_2d(Nx, Ny, fdrhodx.value, drhodx.value,  algorithm);
fftw_plan drhody_bwd = fftw_plan_dft_c2r_2d(Nx, Ny, fdrhody.value, drhody.value,  algorithm);

fftw_plan dWxdx_bwd = fftw_plan_dft_c2r_2d(Nx, Ny, fdWxdx.value, dWxdx.value,  algorithm); // Bertin
fftw_plan dWxdy_bwd = fftw_plan_dft_c2r_2d(Nx, Ny, fdWxdy.value, dWxdy.value,  algorithm); // Bertin
fftw_plan dWydx_bwd = fftw_plan_dft_c2r_2d(Nx, Ny, fdWydx.value, dWydx.value,  algorithm); // Bertin
fftw_plan dWydy_bwd = fftw_plan_dft_c2r_2d(Nx, Ny, fdWydy.value, dWydy.value,  algorithm); // Bertin

fftw_plan dPxdx_bwd = fftw_plan_dft_c2r_2d(Nx, Ny, fdPxdx.value, dPxdx.value,  algorithm); // GA
fftw_plan dPxdy_bwd = fftw_plan_dft_c2r_2d(Nx, Ny, fdPxdy.value, dPxdy.value,  algorithm); // GA
fftw_plan dPydx_bwd = fftw_plan_dft_c2r_2d(Nx, Ny, fdPydx.value, dPydx.value,  algorithm); // GA
fftw_plan dPydy_bwd = fftw_plan_dft_c2r_2d(Nx, Ny, fdPydy.value, dPydy.value,  algorithm); // GA

fftw_plan LWx_bwd = fftw_plan_dft_c2r_2d(Nx, Ny, fLWx.value, LWx.value,  algorithm); // Laplacian
fftw_plan LWy_bwd = fftw_plan_dft_c2r_2d(Nx, Ny, fLWy.value, LWy.value,  algorithm); // Laplacian
fftw_plan L2Wx_bwd = fftw_plan_dft_c2r_2d(Nx, Ny, fL2Wx.value, L2Wx.value,  algorithm); // Laplacian^2
fftw_plan L2Wy_bwd = fftw_plan_dft_c2r_2d(Nx, Ny, fL2Wy.value, L2Wy.value,  algorithm); // Laplacian^2

fftw_plan Lrho_bwd = fftw_plan_dft_c2r_2d(Nx, Ny, fLrho.value, Lrho.value,  algorithm); // Laplacian


void Destroy_FFT_Plans()
{
	fftw_destroy_plan(Wx_fwd);
	fftw_destroy_plan(Wy_fwd);
}

double Compute_P()
{
	clock_t start_time, end_time;
	start_time = clock();

	for (int x = 0; x < Nx; x++)
		for (int y = 0; y < Ny; y++)
		{
			int id = index(x,y);
			if (rho.value[id] > epsilone)
			{
				Px.value[id] = Wx_2.value[id] / rho_2.value[id];
				Py.value[id] = Wy_2.value[id] / rho_2.value[id];
				double pi = sqrt(Px.value[id]*Px.value[id] + Py.value[id]*Py.value[id]);
				double wi = sqrt(Wx_2.value[id]*Wx_2.value[id] + Wy_2.value[id]*Wy_2.value[id]);
				if (pi > 1)
				{
//					Px.value[id] /= pi;
//					Py.value[id] /= pi;
//					pi = sqrt(Px.value[id]*Px.value[id] + Py.value[id]*Py.value[id]);
					cout << "Danger!: Pi = " << pi << "\tw = " << wi << "\trho = " << rho_2.value[id] << endl;
				}
			}
			else
			{
				Px.value[id] = Wx_2.value[id] / rho_2.value[id];
				Py.value[id] = Wy_2.value[id] / rho_2.value[id];
				static int count = 0;
				if (count % 100000 == 0)
				{
					cout << endl << count << " Warning Low rho\t" << rho.value[id] << "\t" << rho_2.value[id] << "\t" << Px.value[id] << "\t" << Py.value[id] << endl;
//					exit(0);
				}
				count++;
				if (count > 300000)
					exit(0);
				Px.value[id] = 0;
				Py.value[id] = 0;
			}
		}

	fftw_execute(Px_fwd);
	fftw_execute(Py_fwd);
	temp_fPx = fPx;
	temp_fPy = fPy;
	Cut(temp_fPx, 2*(Nkx/2+1)/3, 2*Nky/3);
	Cut(temp_fPy, 2*(Nkx/2+1)/3, 2*Nky/3);
	fftw_execute(Px_2_bwd);
	fftw_execute(Py_2_bwd);

	temp_fPx = fPx;
	temp_fPy = fPy;
	Cut(temp_fPx, kcx, kcy);
	Cut(temp_fPy, kcx, kcy);
	fftw_execute(Px_bwd);
	fftw_execute(Py_bwd);

	Px_2 /= N;
	Py_2 /= N;

	Px /= N;
	Py /= N;

	end_time = clock();
	double duration = round(time_digits*(end_time - start_time) / CLOCKS_PER_SEC) / time_digits;
	return (duration);
}

void Forward_Fourier_Transform()
{
	fftw_execute(Wx_fwd);
	fftw_execute(Wy_fwd);
	fftw_execute(rho_fwd);
}

double Backward_Fourier_Transform()
{
	clock_t start_time, end_time;
	start_time = clock();

	temp_fWx = fWx;
	temp_fWy = fWy;
	temp_frho = frho;

	Cut(temp_fWx, 2*(Nkx/2+1)/3, 2*Nky/3);
	Cut(temp_fWy, 2*(Nkx/2+1)/3, 2*Nky/3);
	Cut(temp_frho, 2*(Nkx/2+1)/3, 2*Nky/3);
	fftw_execute(Wx_2_bwd);
	fftw_execute(Wy_2_bwd);
	fftw_execute(rho_2_bwd);

	temp_fWx = fWx;
	temp_fWy = fWy;
	temp_frho = frho;

	Cut(temp_fWx, kcx, kcy);
	Cut(temp_fWy, kcx, kcy);
	Cut(temp_frho, kcx, kcy);
	fftw_execute(Wx_bwd);
	fftw_execute(Wy_bwd);
	fftw_execute(rho_bwd);

	Wx_2 /= N;
	Wy_2 /= N;
	rho_2 /= N;

	Wx /= N;
	Wy /= N;
	rho /= N;

	end_time = clock();
	double duration = round(time_digits*(end_time - start_time) / CLOCKS_PER_SEC) / time_digits;
	return (duration);
}

double Compute_Derivatives_Shared()
{
	clock_t start_time, end_time;
	start_time = clock();

// Density gradient
	Derivative_x(frho,fdrhodx);
	Derivative_y(frho,fdrhody);
	Cut(fdrhodx, kcx, kcy);
	Cut(fdrhody, kcx, kcy);
	fftw_execute(drhodx_bwd);
	fftw_execute(drhody_bwd);
	drhodx /= N;
	drhody /= N;

// Ws Laplacian
	Laplacian(fWx,fLWx);
	Laplacian(fWy,fLWy);
	Cut(fLWx, kcx, kcy);
	Cut(fLWy, kcx, kcy);
	fftw_execute(LWx_bwd);
	fftw_execute(LWy_bwd);
	LWx /= N;
	LWy /= N;

	end_time = clock();
	double duration = round(time_digits*(end_time - start_time) / CLOCKS_PER_SEC) / time_digits;

	return (duration);
}

// Derivatives required by GA only
double Compute_Derivatives_GA()
{
	clock_t start_time, end_time;
	start_time = clock();

// For higher nabla expansion
//	Laplacian_Square(fWx,fL2Wx);
//	Laplacian_Square(fWy,fL2Wy);
//	Cut(fL2Wx, kcx, kcy);
//	Cut(fL2Wy, kcx, kcy);
//	fftw_execute(L2Wx_bwd);
//	fftw_execute(L2Wy_bwd);
//	L2Wx /= N;
//	L2Wy /= N;

// Ps gradients
	Derivative_x(fPx,fdPxdx);
	Derivative_y(fPx,fdPxdy);
	Derivative_x(fPy,fdPydx);
	Derivative_y(fPy,fdPydy);
	Cut(fdPxdx, kcx, kcy);
	Cut(fdPxdy, kcx, kcy);
	Cut(fdPydx, kcx, kcy);
	Cut(fdPydy, kcx, kcy);
	fftw_execute(dPxdx_bwd);
	fftw_execute(dPxdy_bwd);
	fftw_execute(dPydx_bwd);
	fftw_execute(dPydy_bwd);
	dPxdx /= N;
	dPxdy /= N;
	dPydx /= N;
	dPydy /= N;

	end_time = clock();
	double duration = round(time_digits*(end_time - start_time) / CLOCKS_PER_SEC) / time_digits;

	return (duration);
}

double Compute_Derivatives_Truncation()
{
	clock_t start_time, end_time;
	start_time = clock();

// Density laplacian
	Laplacian(frho,fLrho);
	Cut(fLrho, kcx, kcy);
	fftw_execute(Lrho_bwd);
	Lrho /= N;

// Ws gradients
	Derivative_x(fWx,fdWxdx);
	Derivative_y(fWx,fdWxdy);
	Derivative_x(fWy,fdWydx);
	Derivative_y(fWy,fdWydy);
	Cut(fdWxdx, kcx, kcy);
	Cut(fdWxdy, kcx, kcy);
	Cut(fdWydx, kcx, kcy);
	Cut(fdWydy, kcx, kcy);
	fftw_execute(dWxdx_bwd);
	fftw_execute(dWxdy_bwd);
	fftw_execute(dWydx_bwd);
	fftw_execute(dWydy_bwd);
	dWxdx /= N;
	dWxdy /= N;
	dWydx /= N;
	dWydy /= N;

	end_time = clock();
	double duration = round(time_digits*(end_time - start_time) / CLOCKS_PER_SEC) / time_digits;

	return (duration);
}

#endif
