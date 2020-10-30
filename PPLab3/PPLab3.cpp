#define CHUNK 100 
#define NMAX 100 
#include<omp.h>
#include <iostream>
#define DYNAMIC
int main() {
	int i;
	double* a = (double*)malloc(sizeof(double) * NMAX);
	double* sum = (double*)malloc(sizeof(double) * NMAX);
	double* b = (double*)malloc(sizeof(double) * NMAX);
	int chunk, n;
	chunk = CHUNK;
	n = NMAX;
	omp_set_num_threads(16);
	for (i = 0; i < NMAX; i++) {
		a[i] = i;
		b[i] = i+1;
	}
	double st_time, end_time;
	st_time = omp_get_wtime();
#ifdef GUIDED
#pragma omp for schedule(guided) 
#endif
#ifdef STATIC
#pragma omp for schedule(static) 
#endif
#ifdef DYNAMIC
#pragma omp for schedule(dynamic,chunk) 
#endif
	for (i = 0; i < NMAX; i++) {
		sum[i] = b[i] + a[i];
	}
	end_time = omp_get_wtime();
	end_time = end_time - st_time;
	for (i = 0; i < NMAX; i++) {
		printf("\nSum[%d] IS %lf ", i, sum[i]);
	}
	printf("\nTIME OF WORK IS %lf ", end_time);
	free(a);
	free(sum);
	free(b);
	return 0;
}
