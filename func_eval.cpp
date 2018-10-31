#include <math.h>

#include <stdio.h>
#include "create_matrix.h"
#include "func_eval.h"

static double CalcNorm(int n, double* a)
{
	int i;
	int j;
	double tmp;
	double rezult;

	rezult = 0.0;
	for (i = 0; i < n; ++i)
	{
		tmp = 0.0;
		for (j = 0; j < n; ++j)
			tmp += fabs(a[i * n + j]);

		if (rezult < tmp)
			rezult = tmp;
	}

	return rezult;
}

static void Ref(int n, double* a)
{
	int i;
	int j;
	int k;
	double tmp1;
	double tmp2;

	for (i = 0; i < n - 2; ++i)
	{
		tmp1 = 0.0;
		for (j = i + 2; j < n; ++j)
			tmp1 += a[j * n + i] * a[j * n + i];

		tmp2 = sqrt(a[(i + 1) * n + i] * a[(i + 1) * n + i] + tmp1);

		if (tmp2 < 1e-100) // If all elements except diagonal are zero, go to next iter
		{
			a[(i + 1) * n + i] = 0.0;
			a[(i + 2) * n + i] = 0.0;

			continue;
		}

		if (tmp1 < 1e-100) // If all elements except diagonal and neighboor are zero, go to next iter
		{
			a[(i + 2) * n + i] = 0.0;

			continue;
		}
        
        // Evaluating vector x from lemma and reflection matrix
        
		a[(i + 1) * n + i] -= tmp2;

		tmp1 = 1.0 / sqrt(a[(i + 1) * n + i] * a[(i + 1) * n + i] + tmp1);
		for (j = i + 1; j < n; ++j)
			a[j * n + i] *= tmp1;
        
        // Evaluating unitar matrix A(i) = UA(i-1)U*
        
		for (j = i + 1; j < n; ++j)
		{
			tmp1 = 0.0;
			for (k = i + 1; k < n; ++k)
				tmp1 += a[k * n + i] * a[k * n + j];

			tmp1 *= 2.0;
			for (k = i + 1; k < n; ++k)
				a[k * n + j] -= tmp1 * a[k * n + i];
		}

		for (j = 0; j < n; ++j)
		{
			tmp1 = 0.0;
			for (k = i + 1; k < n; ++k)
				tmp1 += a[k * n + i] * a[j * n + k];

			tmp1 *= 2.0;
			for (k = i + 1; k < n; ++k)
				a[j * n + k] -= tmp1 * a[k * n + i];
		}

		a[(i + 1) * n + i] = tmp2;
		for (j = i + 2; j < n; ++j)
			a[j * n + i] = 0.0;
	}
}

static double GetShift(int n, double* a, int k)
{
	return a[(k - 1) * n + k - 1];
}

static void MakeShift(int n, double* a, int k, double s)
{
	int i;

	for (i = 0; i < k; ++i)
		a[i * n + i] -= s;
}

static void MakeDecomposition(int n, double* a, int k, double* x1, double* x2)
{
	int i;
	int j;
	double tmp1;
	double tmp2;

	for (i = 0; i < k - 1; ++i)
	{
        // Creating reflection vectors
        
		tmp1 = a[(i + 1) * n + i] * a[(i + 1) * n + i];

        if (tmp1 < 1e-100) // Only diagonal part of matrix
		{
			tmp2 = fabs(a[i * n + i]);
			x1[i] = (a[i * n + i] > 0.0 ? 1.0 : -1.0);
			x2[i] = 0.0;
		}
		else
		{
            tmp2 = sqrt(a[i * n + i] * a[i * n + i] + tmp1); // Norm

			a[i * n + i] -= tmp2;

			tmp1 = sqrt(a[i * n + i] * a[i * n + i] + tmp1);
			x1[i] = a[i * n + i] / tmp1;
			x2[i] = a[(i + 1) * n + i] / tmp1;
		}

        // Reflection (note, that a[(i + 1) * n + i] add up to zero)
		for (j = i + 1; j < k; ++j)
		{
			tmp1 = x1[i] * a[i * n + j];
			tmp1 += x2[i] * a[(i + 1) * n + j];

			tmp1 *= 2.0;

			a[i * n + j] -= tmp1 * x1[i];
			a[(i + 1) * n + j] -= tmp1 * x2[i];
		}

		a[i * n + i] = tmp2;
		a[(i + 1) * n + i] = 0.0;
	}
}

static void CalcProduction(int n, double* a, int k, double* x1, double* x2)
{
	int i;
	int j;
	double tmp;
    
    // Find R, and apply to A

	for (i = 0; i < k - 1; ++i)
		for (j = 0; j < i + 2; ++j)
		{
			tmp = a[j * n + i] * x1[i];
			tmp += a[j * n + i + 1] * x2[i];

			tmp *= 2.0;

			a[j * n + i] -= tmp * x1[i];
			a[j * n + i + 1] -= tmp * x2[i];
		}
}

void FindValues(int n, double* a, double* values, double eps, double* x1, double* x2, int* iterOut, int max_out, int v)
{
	int i;
	int k;
	int iter;
	double t;
	double s;

	iter = 0;

	t = CalcNorm(n, a) * eps;
    if (v == 1){
        printf("\nNorm * eps = %10.3e , where eps = %10.3e\n", t, eps);
    }

	Ref(n, a);

    if (v == 1){
        printf("\nMatrix A after reflection:\n");
        PrintMatrix(n, a, max_out);
    }

	for (k = n; k > 2; --k)
		while (fabs(a[(k - 1) * n + k - 2]) > t)
		{
			s = GetShift(n, a, k);
			MakeShift(n, a, k, s);


			MakeDecomposition(n, a, k, x1, x2);
			CalcProduction(n, a, k, x1, x2);

			MakeShift(n, a, k, -s);
            if (v == 1){
                printf("\nIteration = %d\n", iter);
                printf("\nMatrix x1:\n");
                PrintMatrix(n, x1, max_out);
                printf("\nMatrix x2:\n");
                PrintMatrix(n, x2, max_out);
                printf("\nMatrix R (=A(i)):\n");
                PrintMatrix(n, a, max_out);
            }
			++iter;
		}

	if (n > 1)
	{
		t = a[0 * n + 0] + a[1 * n + 1];
		s = a[0 * n + 0] * a[1 * n + 1] - a[0 * n + 1] * a[1 * n + 0];
		s = sqrt(t * t - 4.0 * s);

		a[0 * n + 0] = 0.5 * (t + s);
		a[1 * n + 1] = 0.5 * (t - s);
	}

	for (i = 0; i < n; ++i)
		values[i] = a[i * n + i];

	*iterOut = iter;
}
