#include <stdio.h>
#include <ctype.h>
#include <unistd.h>
#include <stdlib.h>
#include <fenv.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <time.h>
#include "args.h"
#include "create_matrix.h"
#include "func_eval.h"

#define EPS 1e-16

/* factorial supports the following command-line arguments:
 * 
 * -i input_file_name.txt - name of the input file (default = NULL)
 * -n number - number of elements (default = 0)
 * -v - option for debugging
 * -f formula - define formula
 */

struct timespec diff(struct timespec start, struct timespec end);

struct timespec diff(struct timespec start, struct timespec end)
{
    struct timespec temp;
    if ((end.tv_nsec-start.tv_nsec)<0) {
        temp.tv_sec = end.tv_sec-start.tv_sec-1;
        temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
    } else {
        temp.tv_sec = end.tv_sec-start.tv_sec;
        temp.tv_nsec = end.tv_nsec-start.tv_nsec;
    }
    return temp;
}

int main(int argc, char **argv){
    
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
    
    struct timespec time_start, time_end;
    
    int n = 0;
    int max_out = 0;
    char *inFileName = NULL;
    int verbose = 0;
    char *formula = NULL;
    double inv1;
    double inv2;
    double *x1 = NULL;
    double *x2 = NULL;
    double *A = NULL;
    double *values = NULL;
    FILE *fin;
    double eps = 1e-10;
    int iter;
    

    
    
    if ((get_args(&n, &inFileName, &verbose, &formula, &max_out, argc, argv)) != 0){
        fprintf (stderr, "Error: Can't get arguments!\n");
        return -1;
    }
    
    if (inFileName == NULL){
        A = (double *)malloc (n * n * sizeof(double));
	x1 = (double*)malloc(n * sizeof(double));
        x2 = (double*)malloc(n * sizeof(double));
        values = (double*)malloc(n * sizeof(double));
	
        if ((A && x1 && x2 && values) != 1) {
            printf ("Error: Not enough memory for matrice A!\n");
	    if (A == NULL)
	      free(A);
	    if (x1 == NULL)
	      free(x1);
            if (x2 == NULL)
	      free(x2);
            if (values == NULL)
	      free(values);
            return -1;
        }
        
        if ((create_matrix(A, n, formula)) != 0){
	    printf ("Error: Can't create matrix!\n");
	    free(A);
	    free(x1);
            free(x2);
            free(values);
	    return -1;
	}
    } else {
	fin = fopen(inFileName, "r");
	
	if (fin == NULL){
	    printf("Error: Can't open file!\n");
	    return -1;
	}
	
	if (fscanf(fin, "%d", &n) != 1){
	    printf("Error: Can't read dimension from file!\n");
	    fclose(fin);
	    return -1;
	}
	
	if (n < 1){
	    printf ("Error: Invalid matrix dimension!\n");
	    fclose(fin);
	    return -1;
	}
	
	A = (double *)malloc (n * n * sizeof(double));
	x1 = (double*)malloc(n * sizeof(double));
        x2 = (double*)malloc(n * sizeof(double));
        values = (double*)malloc(n * sizeof(double));
	
        if ((A && x1 && x2 && values) != 1) {
            printf ("Error: Not enough memory for matrice A!\n");
	    if (A == NULL)
	      free(A);
	    if (x1 == NULL)
	      free(x1);
            if (x2 == NULL)
	      free(x2);
            if (values == NULL)
	      free(values);
            return -1;
        }
        
        if (InputMatrix(n, A, fin) != 0){
	    printf("Error: Can't read matrix from file!\n");
	    fclose(fin);
	    free(A);
	    free(x1);
            free(x2);
            free(values);
	    return -1;
	}
    }
    
    printf("\nMatrix A:\n");
    
    PrintMatrix(n, A, max_out);
    
    printf("\n");
    
    printf("Calculating...\n");
    
    inv1 = 0.0;
    inv2 = 0.0;
    for (int i = 0; i < n; ++i){
        inv1 -= A[i * n + i];
        for (int j = 0; j < n; ++j)
            inv2 -= A[i * n + j] * A[j * n + i];
    }
    
    if( clock_gettime( CLOCK_MONOTONIC, &time_start) == -1 ) {
        perror( "clock gettime" );
        exit( EXIT_FAILURE );
    }
    
    FindValues(n, A, values, eps, x1, x2, &iter);
    
    if( clock_gettime( CLOCK_MONOTONIC, &time_end) == -1 ) {
        perror( "clock gettime" );
        exit( EXIT_FAILURE );
    }
    time_end = diff(time_start, time_end);
    
    for (int i = 0; i < n; ++i){
        inv1 += values[i];
        inv2 += values[i] * values[i];
    }
    
    printf("\nValues:\n");
    PrintVector(n, values, max_out);
    printf("\n");
    
    printf("Finding time\t\t= %f sec.\n\n", 
           (double)time_end.tv_sec + (double)time_end.tv_nsec/(double)1000000000);
    printf("Iterations\t= %d\n\n", iter);
    
    printf("Sum(x_i) - Sum(a_i)\t\t= %g\n", inv1);
    printf("Sum(x_i ^ 2) - Sum(a_ij * a_ji)\t= %g\n", inv2);

    
    

    

    free(A);
    free(x1);
    free(x2);
    free(values);
    
    return 0;
}
