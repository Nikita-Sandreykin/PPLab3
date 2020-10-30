#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <stdlib.h> // pulls in declaration of malloc, free
#include <string.h> // pulls in declaration for strlen.
#define ODD
int main(int argc, char* argv[]) {
	int *sendcounts;
	int *displs;
	#ifdef ODD
		const int NMAX = 7;
	#endif
	#ifdef EVEN
		const int NMAX = 8;
	#endif
	int ProcRank, ProcNum;
	int i;
	double* a = (double*)malloc(NMAX * sizeof(double));
	double* b = (double*)malloc(NMAX * sizeof(double));
	double* c = (double*)malloc(NMAX * sizeof(double));
	double st_time, end_time;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	if (ProcRank == 0) {
		for (i = 0; i < NMAX; i++) {
			a[i] =i;
			b[i] =i+1;
		}

	}
#ifdef ODD
	sendcounts = (int*)malloc(sizeof(int)*ProcNum);
	displs = (int*)malloc(sizeof(int)*ProcNum);
	int equalPortionSize = (NMAX) / ProcNum + ((NMAX%ProcNum == 0) ? 0 : 1);
	for (i = 0; i < ProcNum - 1; i++) {
		sendcounts[i] = equalPortionSize;
	}
	int lastPortionSize = NMAX - equalPortionSize * (ProcNum - 1);
	int portionSize = (ProcRank == ProcNum - 1) ? lastPortionSize : equalPortionSize;
	sendcounts[ProcNum - 1] = lastPortionSize;
	displs[0] = 0;
	for (i = 1; i < ProcNum; i++) {
		displs[i] = displs[i - 1] + sendcounts[i - 1];
	}
	double *abuf = (double*)malloc(portionSize * sizeof(double));
	double *bbuf = (double*)malloc(portionSize * sizeof(double));
	double *cbuf = (double*)malloc(portionSize * sizeof(double));
	MPI_Scatterv(a, sendcounts, displs, MPI_DOUBLE, abuf, sendcounts[ProcRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatterv(b, sendcounts, displs, MPI_DOUBLE, bbuf, sendcounts[ProcRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
	st_time = MPI_Wtime();
	for (i = 0; i < portionSize; i++) {
		cbuf[i] = abuf[i] + bbuf[i];
	}
	MPI_Gatherv(cbuf, sendcounts[ProcRank], MPI_DOUBLE, c, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
#ifdef EVEN 
	int count = NMAX / ProcNum;
	double *abuf = (double*)malloc(count * sizeof(double));
	double *bbuf = (double*)malloc(count * sizeof(double));
	double *cbuf = (double*)malloc(count * sizeof(double));
	MPI_Scatter(a, count, MPI_DOUBLE, abuf, count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatter(b, count, MPI_DOUBLE, bbuf, count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	st_time = MPI_Wtime();
	for (i = 0; i < count; i++) {
		cbuf[i] = abuf[i] + bbuf[i];
	}
	MPI_Gather(cbuf, count, MPI_DOUBLE, c, count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
	end_time = MPI_Wtime();
	end_time = end_time - st_time;
	if (ProcRank == 0) {
		for (i = 0; i < NMAX; i++) {
			printf("\nelement %d  = %10.2lf", i, c[i]);
		}
		printf("\nTIME OF WORK IS %lf ", end_time);
	}
	free(a);
	free(b);
	free(c);
	free(abuf);
	free(bbuf);
	free(cbuf);
#ifdef ODD
	free(displs);
	free(sendcounts);
#endif
	MPI_Finalize();
	return 0;
}
