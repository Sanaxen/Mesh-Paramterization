#ifndef _MY_TYPES_H
#define _MY_TYPES_H

#include "qd/dd_real.h"

#define USE_DD

#ifdef USE_DD
#define MYdouble	dd_real
#else
#define MYdouble	double
#endif

#define GLfloat		double
#define GLdouble	MYdouble
#define GLuint		unsigned int

#ifndef GL_VERSION_1_3
#define GL_TEXTURE0                       0x84C0
#endif

#ifndef GL_ARB_multitexture
#define GL_TEXTURE0_ARB                   0x84C0
#endif


#if 0
class double_MYdouble
{
	MYdouble* b;
	double *A;
	int num;

public:

	inline MYdouble* ref()
	{
		return b;
	}
	inline MYdouble val()
	{
		return b[0];
	}

	double_MYdouble()
	{}
	double_MYdouble(int n, double* a)
	{
		num = n;
		A = a;
		b = new MYdouble[n];

		for ( int i = 0; i < n; i++ )
		{
			b[i] = a[i];
		}
	}
	double_MYdouble(int n, const double* a)
	{
		num = n;
		A = const_cast<double*>(a);
		b = new MYdouble[n];

		for ( int i = 0; i < n; i++ )
		{
			b[i] = a[i];
		}
	}
	double_MYdouble(const double a)
	{
		num = 1;
		b = new MYdouble[1];

		b[0] = a;
	}

	~double_MYdouble()
	{
		To_double(A);
		delete [] b;
	}

	inline void To_double(double* a)
	{
		for ( int i = 0; i < num; i++ )
		{
			a[i] = to_double(b[i]);
		}
	}
};
#endif

#endif