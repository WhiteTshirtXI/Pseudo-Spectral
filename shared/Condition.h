#ifndef _CONDITION_
#define _CONDITION_

#include "FFT.h"

void Add_Noise(double amplitude)
{
	for (int x = 0; x < Nx_box; x++)
		for (int y = 0; y < Ny_box; y++)
		{
			if (rho.value[index(x,y)] > epsilone*initial_rho)
			{
				double Wsquare = Wx.value[index(x,y)]*Wx.value[index(x,y)] + Wy.value[index(x,y)]*Wy.value[index(x,y)]; // Noise of GA
				double P2 = Wsquare / (rho.value[index(x,y)]*rho.value[index(x,y)]);
//				double P4 = P2*P2;
				double temp_amplitude = rho.value[index(x,y)]*(1 - P2)*amplitude; // Noise of GA
//				double temp_amplitude = initial_rho*amplitude;
				temp1.value[index(x,y)] = temp_amplitude*gsl_ran_ugaussian(gsl_r);
				temp2.value[index(x,y)] = temp_amplitude*gsl_ran_ugaussian(gsl_r);
			}
		}
	nonlinear_x_2 = temp1;
	nonlinear_y_2 = temp2;
}

void fAdd_Noise(double amplitude)
{
	for (int i = 0; i < Nkx; i++)
		for (int j = 0; j < Nky; j++)
		{
			temp_fWx.value[i*Nky+j] = amplitude*(gsl_ran_ugaussian(gsl_r)+I*gsl_ran_ugaussian(gsl_r));
			temp_fWy.value[i*Nky+j] = amplitude*(gsl_ran_ugaussian(gsl_r)+I*gsl_ran_ugaussian(gsl_r));
		}
	Cut(temp_fWx, kcx, kcy);
	Cut(temp_fWy, kcx, kcy);
	fWx += temp_fWx;
	fWy += temp_fWy;
}



void Periodic_Boundary_Condition()
{
}

void Ghost_Slip_Boundary_Condition()
{
	for (int x = 1; x < Nx_box-1; x++)
		for (int y = 1; y < Ny_box-1; y++)
		{
			rho.value[index(x+Nx_box-1,y)] = rho.value[index(Nx_box-1-x,y)];
			rho.value[index(x,y+Ny_box-1)] = rho.value[index(x,Ny_box-1-y)];
			rho.value[index(x+Nx_box-1,y+Ny_box-1)] = rho.value[index(Nx_box-1-x,Ny_box-1-y)];

			Wx.value[index(x+Nx_box-1,y)] = -Wx.value[index(Nx_box-1-x,y)];
			Wx.value[index(x,y+Ny_box-1)] = Wx.value[index(x,Ny_box-1-y)];
			Wx.value[index(x+Nx_box-1,y+Ny_box-1)] = -Wx.value[index(Nx_box-1-x,Ny_box-1-y)];

			Wy.value[index(x+Nx_box-1,y)] = Wy.value[index(Nx_box-1-x,y)];
			Wy.value[index(x,y+Ny_box-1)] = -Wy.value[index(x,Ny_box-1-y)];
			Wy.value[index(x+Nx_box-1,y+Ny_box-1)] = -Wy.value[index(Nx_box-1-x,Ny_box-1-y)];
		}
	for (int x = 0; x < Nx_box; x++)
	{
		Wy.value[index(x,0)] = 0;
		Wy.value[index(x,Ny_box-1)] = 0;
	}
	for (int y = 0; y < Ny_box; y++)
	{
		Wx.value[index(0,y)] = 0;
		Wx.value[index(Nx_box-1,y)] = 0;
	}
}

void Ghost_No_Slip_Boundary_Condition()
{
	for (int x = 1; x < Nx_box-1; x++)
		for (int y = 1; y < Ny_box-1; y++)
		{
			rho.value[index(x+Nx_box-1,y)] = rho.value[index(Nx_box-1-x,y)];
			rho.value[index(x,y+Ny_box-1)] = rho.value[index(x,Ny_box-1-y)];
			rho.value[index(x+Nx_box-1,y+Ny_box-1)] = rho.value[index(Nx_box-1-x,Ny_box-1-y)];

			Wx.value[index(x+Nx_box-1,y)] = -Wx.value[index(Nx_box-1-x,y)];
			Wx.value[index(x,y+Ny_box-1)] = -Wx.value[index(x,Ny_box-1-y)];
			Wx.value[index(x+Nx_box-1,y+Ny_box-1)] = Wx.value[index(Nx_box-1-x,Ny_box-1-y)];

			Wy.value[index(x+Nx_box-1,y)] = -Wy.value[index(Nx_box-1-x,y)];
			Wy.value[index(x,y+Ny_box-1)] = -Wy.value[index(x,Ny_box-1-y)];
			Wy.value[index(x+Nx_box-1,y+Ny_box-1)] = Wy.value[index(Nx_box-1-x,Ny_box-1-y)];
		}
}

void Hard_Wall(int depth)
{
	for (int x = 0; x < Nx_box; x++)
	{
		for (int d = 0; d < depth; d++)
		{
			int k = index(x,d);
			Wx.value[k] = 0;
			Wy.value[k] = rho.value[k];
			k = index(x,Ny_box-1-d);
			Wx.value[k] = 0;
			Wy.value[k] = -rho.value[k];
		}
	}

	for (int y = 0; y < Ny_box; y++)
	{
		for (int d = 0; d < depth; d++)
		{
			int k = index(d,y);
			Wy.value[k] = 0;
			Wx.value[k] = rho.value[k];
			k = index(Nx_box-1-d,y);
			Wy.value[k] = 0;
			Wx.value[k] = -rho.value[k];
		}
	}

	for (int dx = 0; dx < depth; dx++)
		for (int dy = 0; dy < depth; dy++)
		{
			double amplitude = sqrt(dy*dy+dx*dx+0.1);
			int k = index(dx,dy);
			Wx.value[k] = rho.value[k]*dy/amplitude;
			Wy.value[k] = rho.value[k]*dx/amplitude;
			Wx.value[k] = 0;
			Wy.value[k] = 0;
			k = index(Nx_box-1-dx,dy);
			Wx.value[k] = -rho.value[k]*dy/amplitude;
			Wy.value[k] = rho.value[k]*dx/amplitude;
			Wx.value[k] = 0;
			Wy.value[k] = 0;
			k = index(dx,Ny_box-1-dy);
			Wx.value[k] = rho.value[k]*dy/amplitude;
			Wy.value[k] = -rho.value[k]*dx/amplitude;
			Wx.value[k] = 0;
			Wy.value[k] = 0;
			k = index(Nx_box-1-dx,Ny_box-1-dy);
			Wx.value[k] = -rho.value[k]*dy/amplitude;
			Wy.value[k] = -rho.value[k]*dx/amplitude;
			Wx.value[k] = 0;
			Wy.value[k] = 0;
		}
}

//w^2 n - w.n w
void Soft_Wall(int depth)
{
	double r = 0.3*g;
	int k;
	for (int x = 0; x < Nx_box; x++)
	{
		for (int d = 0; d < depth; d++)
		{
			k = index(x,d);
			Wx.value[k] -= r*dt*(Wy.value[k]*Wx.value[k]);
			Wy.value[k] += r*dt*(Wx.value[k]*Wx.value[k]);

			k = index(x,Ny_box - 1 - d);
			Wx.value[k] += r*dt*(Wy.value[k]*Wx.value[k]);
			Wy.value[k] -= r*dt*(Wx.value[k]*Wx.value[k]);
		}		
	}

	for (int y = 0; y < Ny_box; y++)
	{
		for (int d = 0; d < depth; d++)
		{
			k = index(d,y);
			Wx.value[k] += r*dt*(Wy.value[k]*Wy.value[k]);
			Wy.value[k] -= r*dt*(Wx.value[k]*Wy.value[k]);

			k = index(Nx_box - 1 - d,y);
			Wx.value[k] -= r*dt*(Wy.value[k]*Wy.value[k]);
			Wy.value[k] += r*dt*(Wx.value[k]*Wy.value[k]);
		}		
	}
}

void Empty_Wall_Nearby(int depth)
{
	for (int x = 0; x < Nx_box; x++)
	{
		for (int d = 0; d < depth; d++)
		{
			rho.value[index(x,d)] = 0;
			rho.value[index(x,Ny_box-1-d)] = 0;
			rho.value[index(x,d)] = 0;
			rho.value[index(x,Ny_box-1-d)] = 0;
		}
	}

	for (int y = 0; y < Ny_box; y++)
	{
		for (int d = 0; d < depth; d++)
		{
			rho.value[index(d,y)] = 0;
			rho.value[index(Nx_box-1-d,y)] = 0;
			rho.value[index(d,y)] = 0;
			rho.value[index(Nx_box-1-d,y)] = 0;
		}
	}
}

void No_Slip_Boundary_Condition()
{
	for (int x = 0; x < Nx_box; x++)
	{
		Wx.value[index(x,0)] = 0;
		Wy.value[index(x,0)] = 0;
		Wx.value[index(x,Ny_box-1)] = 0;
		Wy.value[index(x,Ny_box-1)] = 0;

//		rho.value[index(x,0)] = 0;
//		rho.value[index(x,Ny-1)] = 0;
	}

	for (int y = 0; y < Ny_box; y++)
	{
		Wx.value[index(0,y)] = 0;
		Wy.value[index(0,y)] = 0;
		Wx.value[index(Nx_box-1,y)] = 0;
		Wy.value[index(Nx_box-1,y)] = 0;

//		rho.value[index(0,y)] = 0;
//		rho.value[index(Nx-1,y)] = 0;
	}

		Wx.value[index(1,1)] = 0;
		Wy.value[index(1,1)] = 0;

		Wx.value[index(Nx_box-2,1)] = 0;
		Wy.value[index(Nx_box-2,1)] = 0;

		Wx.value[index(1,Ny_box-2)] = 0;
		Wy.value[index(1,Ny_box-2)] = 0;

		Wx.value[index(Nx_box-2,Ny_box-2)] = 0;
		Wy.value[index(Nx_box-2,Ny_box-2)] = 0;
}

void Slip_Boundary_Condition()
{
	double factor = 0;
	for (int x = 0; x < Nx_box; x++)
	{
//		if (Wy.value[index(x,0)] < 0)
			Wy.value[index(x,0)] = 0;
//		if (Wy.value[index(x,Ny_box-1)] > 0)
			Wy.value[index(x,Ny_box-1)] = 0;
	}
	for (int y = 0; y < Ny_box; y++)
	{
//		if (Wx.value[index(0,y)] < 0)
			Wx.value[index(0,y)] = 0;
//		if (Wx.value[Nx-1,y)] > 0)
			Wx.value[index(Nx-1,y)] = 0;
	}
}

void Tangential_Boundary_Condition()
{
	double ldx,ldy,lux,luy,rdx,rdy,rux,ruy;

	if (Wx.value[index(1,0)] < 0)
	{
		double projection = (Wy.value[index(1,0)] - Wx.value[index(1,0)]) / 2;
		Wy.value[index(1,0)] = projection;
		Wx.value[index(1,0)] = -projection;
	}
	else
		Wy.value[index(1,0)] = 0;

	if (Wy.value[index(0,1)] < 0)
	{
		double projection = (Wx.value[index(0,1)] - Wy.value[index(0,1)]) / 2;
		Wx.value[index(0,1)] = projection;
		Wy.value[index(0,1)] = -projection;
	}
	else
		Wx.value[index(0,1)] = 0;


	if (Wx.value[index(Nx_box-2,0)] > 0)
	{
		double projection = (Wx.value[index(Nx_box-2,0)] + Wy.value[index(Nx_box-2,0)]) / 2;
		Wx.value[index(Nx_box-2,0)] = projection;
		Wy.value[index(Nx_box-2,0)] = projection;
	}
	else
		Wy.value[index(Nx_box-2,0)] = 0;

	if (Wy.value[index(Nx_box-1,1)] < 0)
	{
		double projection = (Wx.value[index(Nx_box-1,1)] + Wy.value[index(Nx_box-1,1)]) / 2;
		Wx.value[index(Nx_box-1,1)] = projection;
		Wy.value[index(Nx_box-1,1)] = projection;
	}
	else
		Wx.value[index(Nx_box-1,1)] = 0;

	if (Wy.value[index(Nx_box-1,Ny_box-2)] > 0)
	{
		double projection = (-Wx.value[index(Nx_box-1,Ny_box-2)] + Wy.value[index(Nx_box-1,Ny_box-2)]) / 2;
		Wx.value[index(Nx_box-1,Ny_box-2)] = -projection;
		Wy.value[index(Nx_box-1,Ny_box-2)] = projection;
	}
	else
		Wx.value[index(Nx_box-1,Ny_box-2)] = 0;

	if (Wx.value[index(Nx_box-2,Ny_box-1)] > 0)
	{
		double projection = (-Wx.value[index(Nx_box-2,Ny_box-1)] + Wy.value[index(Nx_box-2,Ny_box-1)]) / 2;
		Wx.value[index(Nx_box-2,Ny_box-1)] = -projection;
		Wy.value[index(Nx_box-2,Ny_box-1)] = projection;
	}
	else
		Wy.value[index(Nx_box-2,Ny_box-1)] = 0;

	if (Wx.value[index(1,Ny_box-1)] < 0)
	{
		double projection = (Wx.value[index(1,Ny_box-1)] - Wy.value[index(1,Ny_box-1)]) / 2;
		Wx.value[index(1,Ny_box-1)] = projection;
		Wy.value[index(1,Ny_box-1)] = -projection;
	}
	else
		Wy.value[index(1,Ny_box-1)] = 0;

	if (Wy.value[index(0,Ny_box-2)] > 0)
	{
		double projection = (Wx.value[index(0,Ny_box-2)] - Wy.value[index(0,Ny_box-2)]) / 2;
		Wx.value[index(1,Ny_box-1)] = projection;
		Wy.value[index(1,Ny_box-1)] = -projection;
	}
	else
		Wx.value[index(0,Ny_box-2)] = 0;

	double factor = 0;
	for (int x = 2; x < (Nx_box - 2); x++)
	{
		if (Wy.value[index(x,0)] < 0)
			Wy.value[index(x,0)] *= factor;
		if (Wy.value[index(x,Ny-1)] > 0)
			Wy.value[index(x,Ny-1)] *= factor;

//			Wy.value[index(x,0)] *= factor;
//			Wy.value[index(x,Ny_box-1)] *= factor;
	}
	for (int y = 2; y < (Ny_box - 2); y++)
	{
		if (Wx.value[index(0,y)] < 0)
			Wx.value[index(0,y)] *= factor;
		if (Wx.value[index(Nx-1,y)] > 0)
			Wx.value[index(Nx-1,y)] *= factor;

//			Wx.value[index(0,y)] *= factor;
//			Wx.value[index(Nx_box-1,y)] *= factor;
//		Wy.value[index(0,y)] = rho.value[index(0,y)];
//		Wy.value[Ny-1,y)] = -rho.value[Ny-1,y)];
	}

	Wx.value[index(0,0)] = Wx.value[index(Nx_box-1,0)] = Wx.value[index(Nx_box-1,Ny_box-1)] = Wx.value[index(0,Ny_box-1)] = 0;
	Wy.value[index(0,0)] = Wy.value[index(Nx_box-1,0)] = Wy.value[index(Nx_box-1,Ny_box-1)] = Wy.value[index(0,Ny_box-1)] = 0;

//	Wx.value[index(0,0)] = (ldx - ldy)/2;
//	Wy.value[index(0,0)] = (-ldx + ldy)/2;

//	Wx.value[index(Nx_box-1,Ny_box-1)] = (rux - ruy)/2;
//	Wy.value[index(Nx_box-1,Ny_box-1)] = (-rux + ruy)/2;

//	Wx.value[index(Nx_box-1,0)] = (rdx + rdy)/2;
//	Wy.value[index(Nx_box-1,0)] = (rdx + rdy)/2;

//	Wx.value[index(0,Ny_box-1)] = (lux + luy)/2;
//	Wy.value[index(0,Ny_box-1)] = (lux + luy)/2;
}

void Vortex_Boundary_Condition()
{
	double amplitude = 0.2;
	for (int x = 0; x < Nx_box; x++)
	{
		Wx.value[index(x,0)] = amplitude;
		Wx.value[index(x,Ny_box-1)] = -amplitude;
		Wy.value[index(x,0)] = Wy.value[index(x,Ny_box-1)] = 0;
	}
	for (int y = 0; y < Ny_box; y++)
	{
		Wy.value[index(Nx_box-1,y)] = amplitude;
		Wy.value[index(0,y)] = -amplitude;
		Wx.value[index(Nx_box-1,y)] = Wx.value[index(0,y)] = 0;
	}

	Wx.value[index(0,0)] = Wx.value[index(Nx_box-1,0)] = Wx.value[index(Nx_box-1,Ny_box-1)] = Wx.value[index(0,Ny_box-1)] = 0;
	Wy.value[index(0,0)] = Wy.value[index(Nx_box-1,0)] = Wy.value[index(Nx_box-1,Ny_box-1)] = Wy.value[index(0,Ny_box-1)] = 0;

	Wx.value[index(0,1)] = amplitude / sqrt(2.0);
	Wy.value[index(0,1)] = -amplitude / sqrt(2.0);
	Wx.value[index(Nx_box-2,0)] = amplitude / sqrt(2.0);
	Wy.value[index(Nx_box-2,0)] = amplitude / sqrt(2.0);
	Wx.value[index(Nx_box-1,Ny_box-2)] = -amplitude / sqrt(2.0);
	Wy.value[index(Nx_box-1,Ny_box-2)] = amplitude / sqrt(2.0);
	Wx.value[index(1,Ny_box-1)] = -amplitude / sqrt(2.0);
	Wy.value[index(1,Ny_box-1)] = -amplitude / sqrt(2.0);

}

void Reflective_Boundary_Condition()
{
	double factor = -1;
	for (int x = 0; x < Nx_box; x++)
	{
		if (Wy.value[index(x,0)] < 0)
			Wy.value[index(x,0)] *= factor;
		if (Wy.value[index(x,Ny_box-1)] > 0)
			Wy.value[index(x,Ny_box-1)] *= factor;

//		Wx.value[index(x,0] = -rho.value[index(x,0];
//		Wx.value[index(x,Ny-1] = rho.value[index(x,Ny-1];
	}
	for (int y = 0; y < Ny_box; y++)
	{
		if (Wx.value[index(0,y)] < 0)
			Wx.value[index(0,y)] *= factor;
		if (Wx.value[index(Nx_box-1,y)] > 0)
			Wx.value[index(Nx_box-1,y)] *= factor;

//		Wy.value[index(0,y)] = rho.value[index(0,y)];
//		Wy.value[Ny-1,y] = -rho.value[Ny-1,y];
	}
}


void Random_Initial_Condition(double input_rho)
{
	rho = input_rho;
	Wx = 0;
	Wy = 0;

	for (int id = 0; id < Nx*Ny; id++)
	{
		double r = gsl_ran_flat(gsl_r,-M_PI, M_PI);
		Wx.value[id] = rho.value[id]*cos(r);
		Wy.value[id] = rho.value[id]*sin(r);
	}
	Wx *= 0.1;
	Wy *= 0.1;
}

void Random_Initial_Condition_With_Hole(double input_rho,double low_rho, int width)
{
	rho = ((input_rho*Nx_box*Ny_box - low_rho*(width-1)*(width-1)) / (Nx_box*Ny_box - width*width));
//	rho = input_rho;
	Wx = 0;
	Wy = 0;

	for (int x = (Nx_box - width)/2; x < (Nx_box + width)/2; x++)
		for (int y = (Ny_box - width)/2; y < (Ny_box + width)/2; y++)
		{
			rho.value[index(x,y)] = low_rho;
			Wx.value[index(x,y)] = 0;
			Wy.value[index(x,y)] = 0;
		}
}

void Polar_Initial_Condition(double input_rho)
{
	rho = input_rho;
	Wx = input_rho;
	Wy = 0;
}

void Horizontal_Band_Initial_Condition(double rho_high, double rho_low, double width_ordered, double width_unordered)
{
	rho = rho_low;
	Wx = 0;
	Wy = 0;
	int int_width_ordered = (int) (Nx*(0.5*width_ordered/Lx));
	int int_width_unordered = (int) (Nx*(0.5*width_unordered/Lx));
	int l = int_width_unordered + int_width_ordered;
	int number = Nx / l;
	for (int i = 0; i < number; i++)
	{
		for (int x = int_width_unordered + i*l; x < (l+i*l); x++)
			for (int y = 0; y < Ny; y++)
			{
				rho.value[index(x,y)] = rho_high;
				Wx.value[index(x,y)] = (rho_high - rho_low);
			}
	}
}

void Vertical_Band_Initial_Condition(double rho_high, double rho_low, double width_ordered, double width_unordered)
{
	rho = rho_low;
	Wx = 0;
	Wy = 0;
	int int_width_ordered = (int) (Ny*(0.5*width_ordered/Ly));
	int int_width_unordered = (int) (Ny*(0.5*width_unordered/Ly));
	int l = int_width_unordered + int_width_ordered;
	int number = Ny / l;
	for (int i = 0; i < number; i++)
	{
		for (int x = 0; x < Nx; x++)
			for (int y = int_width_unordered + i*l; y < (l+i*l); y++)
			{
				rho.value[index(x,y)] = rho_high;
				Wy.value[index(x,y)] = (rho_high - rho_low);
			}
	}
}

void Moving_Clump_Initial_Condition(double rho_high, double rho_low, double width_clump, double height_clump)
{
	rho = rho_low;
	Wx = 0;
	Wy = 0;
	int int_width_clump = (int) (Nx*(0.5*width_clump/Lx));
	int int_height_clump = (int) (Nx*(0.5*height_clump/Ly));
	int l = int_width_clump + int_width_clump;
	int number = Nx / l;
	for (int x = (Nx - int_width_clump) / 2; x < (Nx + int_width_clump) / 2; x++)
				for (int y = (Ny - int_height_clump) / 2; y < (Ny + int_height_clump) / 2; y++)
				{
					rho.value[index(x,y)] = rho_high;
					Wx.value[index(x,y)] = (rho_high - rho_low);
				}
}

void Circular_Vortex_Initial_Condition(double input_rho)
{
	rho = input_rho;
	double W_magnitude = input_rho;
	for (int x = layer; x < Nx_box - layer; x++)
		for (int y = layer; y < Ny_box - layer; y++)
		{
			Wx.value[index(x,y)] = (y - 0.5*(Ny_box-layer-1))/(Ny_box - layer);
			Wy.value[index(x,y)] = -(x - 0.5*(Nx_box-layer-1))/(Nx_box - layer);
			double norm = sqrt(Wx.value[index(x,y)]*Wx.value[index(x,y)] + Wy.value[index(x,y)]*Wy.value[index(x,y)]) / W_magnitude;
			if (norm != 0)
			{
				Wx.value[index(x,y)] /= norm;
				Wy.value[index(x,y)] /= norm;
			}
		}
//	int width = 10;
//	for (int x = Nx_box/2 - width; x < Nx_box/2 + width; x++)
//		for (int y = Ny_box/2 - width; y < Ny_box/2 + width; y++)
//		{
//			rho.value[index(x,y)] = 0.1;
//		}
}

void Square_Vortex_Initial_Condition(double input_rho, double W_magnitude, double angle,double fraction)
{
	rho = input_rho;

	double xr, yr;

	double h = 2*Ly_box*cos(angle)*sin(angle);
	double l = 2*Ly_box*cos(angle)*cos(angle);

	cout << h << "\t\n" << l << endl;

	for (int x = 0; x < Nx_box; x++)
		for (int y = 0; y < Ny_box; y++)
		{
			xr = (2.0*x*Lx_box / (Nx_box-1)) - Lx_box;
			yr = (2.0*y*Ly_box / (Ny_box-1)) - Ly_box;

//			cout << xr << "\t" << yr << endl;

			if ((xr + Lx_box) < h)
			{
				if ((yr + Ly_box) < l)
				{
					if (((xr + Lx_box) / (yr + Ly_box)) < tan(angle))
					{
						rho.value[index(x,y)] = input_rho / ((1 + yr + Ly_box)*tan(angle));
						Wx.value[index(x,y)] = 0.3*(1 + xr / Lx_box);
						Wy.value[index(x,y)] = 1;

						rho.value[index(Nx_box - 1 - x,Ny_box - 1 - y)] = rho.value[index(x,y)];
						Wx.value[index(Nx_box - 1 - x,Ny_box - 1 - y)] = -Wx.value[index(x,y)];
						Wy.value[index(Nx_box - 1 - x,Ny_box - 1 - y)] = -Wy.value[index(x,y)];

						rho.value[index(y,Ny_box - 1 - x)] = rho.value[index(x,y)];
						Wx.value[index(y,Ny_box - 1 - x)] = Wy.value[index(x,y)];
						Wy.value[index(y,Ny_box - 1 - x)] = -Wx.value[index(x,y)];

						rho.value[index(Nx_box - 1 - y,x)] = rho.value[index(x,y)];
						Wx.value[index(Nx_box - 1 - y,x)] = -Wy.value[index(x,y)];
						Wy.value[index(Nx_box - 1 - y,x)] = Wx.value[index(x,y)];
					}
				}
				else
				{
					if ((yr + Ly_box - l) < ((h - (xr + Lx_box))*tan(angle)))
					{
						rho.value[index(x,y)] = input_rho / h;
						Wx.value[index(x,y)] = 0.3*(1 + xr / Lx_box);
						Wy.value[index(x,y)] = 1;

						rho.value[index(Nx_box - 1 - x,Ny_box - 1 - y)] = rho.value[index(x,y)];
						Wx.value[index(Nx_box - 1 - x,Ny_box - 1 - y)] = -Wx.value[index(x,y)];
						Wy.value[index(Nx_box - 1 - x,Ny_box - 1 - y)] = -Wy.value[index(x,y)];

						rho.value[index(y,Ny_box - 1 - x)] = rho.value[index(x,y)];
						Wx.value[index(y,Ny_box - 1 - x)] = Wy.value[index(x,y)];
						Wy.value[index(y,Ny_box - 1 - x)] = -Wx.value[index(x,y)];

						rho.value[index(Nx_box - 1 - y,x)] = rho.value[index(x,y)];
						Wx.value[index(Nx_box - 1 - y,x)] = -Wy.value[index(x,y)];
						Wy.value[index(Nx_box - 1 - y,x)] = Wx.value[index(x,y)];
					}
				}
			}
		}

	for (int x = 0; x < Nx_box; x++)
		for (int y = 0; y < Ny_box; y++)
		{
			double norm = sqrt(Wx.value[index(x,y)]*Wx.value[index(x,y)] + Wy.value[index(x,y)]*Wy.value[index(x,y)]);

			if (norm > 0)
			{
				Wx.value[index(x,y)] /= norm;
				Wy.value[index(x,y)] /= norm;
				Wx.value[index(x,y)] *= W_magnitude;
				Wy.value[index(x,y)] *= W_magnitude;
			}
			else
				rho.value[index(x,y)] = 0.6;
		}

	Wx.Add_Noise(fraction, gsl_r);
	Wy.Add_Noise(fraction, gsl_r);

}

void File_Initial_Condition(const int nx, const int ny, ifstream& data_file)
{
	double input_rho[nx][ny];
	double input_Wx[nx][ny];
	double input_Wy[nx][ny];
	double r;
	string line;
	getline(data_file, line);
	for (int x = 0; x < nx; x++)
	{
		for (int y = 0; y < ny; y++)
		{
			data_file >> r;
			data_file >> r;
			data_file >> input_Wx[x][y];
			data_file >> input_Wy[x][y];
			data_file >> input_rho[x][y];
			rho.value[index(x,y)] = input_rho[x][y];
			Wx.value[index(x,y)] = input_Wx[x][y]*input_rho[x][y];
			Wy.value[index(x,y)] = input_Wy[x][y]*input_rho[x][y];
		}
		getline(data_file, line);
	}
/*	for (int x = 0; x < Nx_box; x++)*/
/*	{*/
/*		for (int y = 0; y < Ny_box; y++)*/
/*		{*/
/*			double xr = (x*nx*1.0) / Nx_box;*/
/*			double yr = (y*ny*1.0) / Ny_box;*/
/*			int x_base = (int) xr;*/
/*			int y_base = (int) yr;*/
/*			rho.value[index(x,y)] = input_rho[x_base][y_base];*/
/*			Wx.value[index(x,y)] = input_Wx[x_base][y_base];*/
/*			Wy.value[index(x,y)] = input_Wy[x_base][y_base];*/
/*		}*/
/*	}*/
}


void Find_Evolution_Matrix()
{
	E = v / 2;

	for (int i = 0; i < Nkx; i++)
	{
		for (int j = 0; j < Nky; j++)
		{
			double k2 = kx[i]*kx[i] + ky[j]*ky[j];

			double Kx = eta*dt*kx[i];
			double Ky = eta*dt*ky[j];
			double K2 = Kx*Kx + Ky*Ky;

			double a = 1 + eta*dt*(T+D*k2);
			double b = 1 + eta*dt*K*k2;
			double det = a*(a*b + K2*E*v);
			double factor = dt / det;
			
			int k_index = i*Nky+j;

			mu_n00.value[k_index] = factor*(a*b + Ky*Ky*E*v);
			mu_n01.value[k_index] = -factor*(Kx*Ky*E*v);
			mu_n02.value[k_index] = -factor*(I*a*Kx*E);

			mu_n10.value[k_index] = -factor*(Kx*Ky*E*v);
			mu_n11.value[k_index] = factor*(a*b + Kx*Kx*E*v);
			mu_n12.value[k_index] = -factor*(I*a*Ky*E);

			mu_n20.value[k_index] = -factor*(I*a*Kx*v);
			mu_n21.value[k_index] = -factor*(I*a*Ky*v);
			mu_n22.value[k_index] = factor*a*a;

			factor = 1 / (det);

			mu_l00.value[k_index] = 1 + (-a*b*dt*(T+D*k2)+E*v*eta*dt*dt*(ky[j]*ky[j] - a*k2)) / det;
			mu_l01.value[k_index] = (-eta*dt*dt*kx[i]*ky[j]*E*v) / det;
			mu_l02.value[k_index] = (-I*a*dt*kx[i]*E) / det;

			mu_l10.value[k_index] = (-eta*dt*dt*kx[i]*ky[j]*E*v) / det;
			mu_l11.value[k_index] = 1 + (-a*b*dt*(T+D*k2)+E*v*eta*dt*dt*(kx[i]*kx[i] - a*k2)) / det;
			mu_l12.value[k_index] = (-I*a*dt*ky[j]*E) / det;

			mu_l20.value[k_index] = (-I*a*dt*kx[i]*v) / det;
			mu_l21.value[k_index] = (-I*a*dt*ky[j]*v) / det;
			mu_l22.value[k_index] = 1 + (-a*a*(dt*K*k2)-E*v*eta*dt*dt*a*k2) / det;

		}
	}
}


#endif
