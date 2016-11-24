#ifndef _CFIELD_
#define _CFIELD_

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <iostream>
#include <iomanip>

template <class a_type, int input_Nx, int input_Ny>
struct Field{
	static const int Nx = input_Nx;
	static const int Ny = input_Ny;
	a_type* value;

	Field();
	Field(const a_type input);
	Field(const Field&);
	~Field();
	a_type Get_Average();
	a_type Get_Max();
	a_type Get_Min();
	int Get_Min_Id();
	void Add_Noise(a_type amplitude, gsl_rng* gsl_r);
	void Write();
	void Sqrt();
	bool Check_Nan();
	Field& operator=(const Field& f);
	Field operator+(const Field& f);
	Field operator-(const Field& f);
	Field operator-(const a_type& shift);
	Field operator*(const Field& f);
	Field operator/(const Field& f);
	Field& operator+=(const Field& f);
	Field& operator-=(const Field& f);
	Field& operator*=(const Field& f);
	Field& operator/=(const Field& f);
	Field& operator=(const a_type& shift);
	Field& operator+=(const a_type& factor);
	Field& operator-=(const a_type& factor);
	Field& operator*=(const a_type& factor);
	Field& operator/=(const a_type& factor);
	friend std::ostream& operator<<(std::ostream& os, Field* f)
	{
		for (int j = 0; j < Ny; j++)
		{
			for (int i = 0; i < Nx; i++)
			{
//				os << creal(f->value[i*Ny+Ny-j-1]) << "+" << cimag(f->value[i*Ny+Ny-j-1]) << "i" << "\t";
//				os << sqrt(creal(f->value[i*Ny+Ny-j-1])*creal(f->value[i*Ny+Ny-j-1]) + cimag(f->value[i*Ny+Ny-j-1])*cimag(f->value[i*Ny+Ny-j-1])) << "\t";
				os << cimag(f->value[i*Ny+Ny-j-1]) << "i" << "\t";
			}
			os << "\n";
		}
		os << "\n";
		return (os);
	}
};

template <class a_type, int input_Nx, int input_Ny>
Field<a_type, input_Nx, input_Ny>::Field()
{
	value = new a_type[Nx*Ny];
	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
			value[i*Ny+j] = 0;
}

template <class a_type, int input_Nx, int input_Ny>
Field<a_type, input_Nx, input_Ny>::Field(const a_type input)
{
}

template <class a_type, int input_Nx, int input_Ny>
Field<a_type, input_Nx, input_Ny>::Field(const Field &f)
{
}

template <class a_type, int input_Nx, int input_Ny>
Field<a_type, input_Nx, input_Ny>::~Field()
{
	delete [] value;
}

template <class a_type, int input_Nx, int input_Ny>
a_type Field<a_type, input_Nx, input_Ny>::Get_Average()
{
	a_type result = 0;
	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
			result += value[i*Ny+j];
	result /= (Nx*Ny);
	return (result);
}

template <class a_type, int input_Nx, int input_Ny>
a_type Field<a_type, input_Nx, input_Ny>::Get_Max()
{
	a_type result = value[0];
	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
			if (value[i*Ny+j] > result)
				result = value[i*Ny+j];
	return (result);
}

template <class a_type, int input_Nx, int input_Ny>
a_type Field<a_type, input_Nx, input_Ny>::Get_Min()
{
	a_type result = value[0];
	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
			if (value[i*Ny+j] < result)
				result = value[i*Ny+j];
	return (result);
}

template <class a_type, int input_Nx, int input_Ny>
int Field<a_type, input_Nx, input_Ny>::Get_Min_Id()
{
	int id = 0;
	for (int i = 0; i < Nx*Ny; i++)
			if (value[i] < value[id])
				id = i;
	return (id);
}

template <class a_type, int input_Nx, int input_Ny>
void Field<a_type, input_Nx, input_Ny>::Add_Noise(const a_type amplitude, gsl_rng* gsl_r)
{
	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
		{
			a_type r = amplitude*gsl_ran_ugaussian(gsl_r);
			value[i*Ny+j] += r;
		}
}

template <class a_type, int input_Nx, int input_Ny>
void Field<a_type, input_Nx, input_Ny>::Write()
{
	std::cout << std::endl;
	for (int j = 0; j < Ny; j++)
	{
		for (int i = 0; i < Nx; i++)
		{
			std::cout << std::setprecision(4) << value[i*Ny+j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

template <class a_type, int input_Nx, int input_Ny>
void Field<a_type, input_Nx, input_Ny>::Sqrt()
{
	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
			value[i*Ny+j] = sqrt(value[i*Ny+j]);
}

template <class a_type, int input_Nx, int input_Ny>
bool Field<a_type, input_Nx, input_Ny>::Check_Nan()
{
	bool flag = false;
	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
			flag = flag || ((value[i*Ny+j] != value[i*Ny+j]) );
	return flag;
}

template <class a_type, int input_Nx, int input_Ny>
Field<a_type, input_Nx, input_Ny>& Field<a_type, input_Nx, input_Ny>::operator=(const Field& f)
{
	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
			value[i*Ny+j] = f.value[i*Ny+j];
	return (*this);
}

template <class a_type, int input_Nx, int input_Ny>
Field<a_type, input_Nx, input_Ny> Field<a_type, input_Nx, input_Ny>::operator+(const Field& f)
{
	Field result;
	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
			result.value[i*Ny+j] = value[i*Ny+j] + f.value[i*Ny+j];
	return (result);
}

template <class a_type, int input_Nx, int input_Ny>
Field<a_type, input_Nx, input_Ny> Field<a_type, input_Nx, input_Ny>::operator-(const Field& f)
{
	Field result;
	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
			result.value[i*Ny+j] = value[i*Ny+j] - f.value[i*Ny+j];
	return (result);
}

template <class a_type, int input_Nx, int input_Ny>
Field<a_type, input_Nx, input_Ny> Field<a_type, input_Nx, input_Ny>::operator-(const a_type& shift)
{
	Field result;
	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
			result.value[i*Ny+j] = value[i*Ny+j] - shift;
	return (result);
}

template <class a_type, int input_Nx, int input_Ny>
Field<a_type, input_Nx, input_Ny> Field<a_type, input_Nx, input_Ny>::operator*(const Field& f)
{
	Field result;
	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
			result.value[i*Ny+j] = value[i*Ny+j] * f.value[i*Ny+j];
	return (result);
}

template <class a_type, int input_Nx, int input_Ny>
Field<a_type, input_Nx, input_Ny> Field<a_type, input_Nx, input_Ny>::operator/(const Field& f)
{
	Field result;
	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
			result.value[i*Ny+j] = value[i*Ny+j] / f.value[i*Ny+j];
	return (result);
}

template <class a_type, int input_Nx, int input_Ny>
Field<a_type, input_Nx, input_Ny>& Field<a_type, input_Nx, input_Ny>::operator+=(const Field& f)
{
	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
			value[i*Ny+j] += f.value[i*Ny+j];
	return (*this);
}

template <class a_type, int input_Nx, int input_Ny>
Field<a_type, input_Nx, input_Ny>& Field<a_type, input_Nx, input_Ny>::operator-=(const Field& f)
{
	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
			value[i*Ny+j] -= f.value[i*Ny+j];
	return (*this);
}

template <class a_type, int input_Nx, int input_Ny>
Field<a_type, input_Nx, input_Ny>& Field<a_type, input_Nx, input_Ny>::operator*=(const Field& f)
{
	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
			value[i*Ny+j] *= f.value[i*Ny+j];
	return (*this);
}

template <class a_type, int input_Nx, int input_Ny>
Field<a_type, input_Nx, input_Ny>& Field<a_type, input_Nx, input_Ny>::operator/=(const Field& f)
{
	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
			if (f.value[i*Ny+j] != 0)
				value[i*Ny+j] /= f.value[i*Ny+j];
			else
				value[i*Ny+j] = 0;
	return (*this);
}

template <class a_type, int input_Nx, int input_Ny>
Field<a_type, input_Nx, input_Ny>& Field<a_type, input_Nx, input_Ny>::operator=(const a_type& shift)
{
	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
			value[i*Ny+j] = shift;
	return (*this);
}

template <class a_type, int input_Nx, int input_Ny>
Field<a_type, input_Nx, input_Ny>& Field<a_type, input_Nx, input_Ny>::operator+=(const a_type& shift)
{
	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
			value[i*Ny+j] += shift;
	return (*this);
}

template <class a_type, int input_Nx, int input_Ny>
Field<a_type, input_Nx, input_Ny>& Field<a_type, input_Nx, input_Ny>::operator-=(const a_type& shift)
{
	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
			value[i*Ny+j] -= shift;
	return (*this);
}

template <class a_type, int input_Nx, int input_Ny>
Field<a_type, input_Nx, input_Ny>& Field<a_type, input_Nx, input_Ny>::operator*=(const a_type& factor)
{
	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
			value[i*Ny+j] *= factor;
	return (*this);
}

template <class a_type, int input_Nx, int input_Ny>
Field<a_type, input_Nx, input_Ny>& Field<a_type, input_Nx, input_Ny>::operator/=(const a_type& factor)
{
	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
			value[i*Ny+j] /= factor;
	return (*this);
}


// Here I define operators between a field and a number with type a_type

template <class a_type, int input_Nx, int input_Ny>
Field<a_type, input_Nx, input_Ny> operator+(const Field<a_type, input_Nx, input_Ny>& f, const a_type shift)
{
	Field<a_type, input_Nx, input_Ny> result;
	for (int i = 0; i < f.Nx; i++)
		for (int j = 0; j < f.Ny; j++)
			result.value[i*f.Ny+j] = shift + f.value[i*f.Ny+j];
	return (result);
}

template <class a_type, int input_Nx, int input_Ny>
Field<a_type, input_Nx, input_Ny> operator+(const a_type shift, const Field<a_type, input_Nx, input_Ny>& f)
{
	Field<a_type, input_Nx, input_Ny> result;
	for (int i = 0; i < f.Nx; i++)
		for (int j = 0; j < f.Ny; j++)
			result.value[i*f.Ny+j] = shift + f.value[i*f.Ny+j];
	return (result);
}

template <class a_type, int input_Nx, int input_Ny>
Field<a_type, input_Nx, input_Ny> operator*(const Field<a_type, input_Nx, input_Ny>& f, const a_type factor)
{
	Field<a_type, input_Nx, input_Ny> result;
	for (int i = 0; i < f.Nx; i++)
		for (int j = 0; j < f.Ny; j++)
			result.value[i*f.Ny+j] = factor * f.value[i*f.Ny+j];
	return (result);
}

template <class a_type, int input_Nx, int input_Ny>
Field<a_type, input_Nx, input_Ny> operator*(const a_type factor, const Field<a_type, input_Nx, input_Ny>& f)
{
	Field<a_type, input_Nx, input_Ny> result;
	for (int i = 0; i < f.Nx; i++)
		for (int j = 0; j < f.Ny; j++)
			result.value[i*f.Ny+j] = factor * f.value[i*f.Ny+j];
	return (result);
}

template <class a_type, int input_Nx, int input_Ny>
Field<a_type, input_Nx, input_Ny> operator/(const Field<a_type, input_Nx, input_Ny>& f, const a_type factor)
{
	Field<a_type, input_Nx, input_Ny> result;
	for (int i = 0; i < f.Nx; i++)
		for (int j = 0; j < f.Ny; j++)
			result.value[i*f.Ny+j] = f.value[i*f.Ny+j] / factor;
	return (result);
}

#endif
