#define _USE_MATH_DEFINES

#include <mpi.h>
#include <omp.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <iomanip>
using namespace std;

const int N = 1e9;
const double a = 0;
const double b = M_PI_4;
const double STEP = (b - a) / N;
inline double f(double x)
{
	return cos(x);
}

int main(int argc, char **argv)
{
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
	double partSum = 0;
	int pSize = rank != size - 1 ? 
			N / size : 
			N / size + N % size;
	int pOffset = rank * (N / size);
	double d0 = a + pOffset * STEP;

	double t = MPI_Wtime();
	#pragma omp parallel reduction(+ : partSum)
	{
		int nThreads = omp_get_num_threads();
		int tRank = omp_get_thread_num();
		int tSize = tRank != nThreads - 1 ? 
				pSize / nThreads : 
				pSize / nThreads + pSize % nThreads;
		int tOffset = tRank * (pSize / nThreads);
		double d = d0 + tOffset * STEP;
		for (int i = 0; i < nThreads; i++)
			if (i == tRank)
				for (int j = 0; j < tSize; j++, d += STEP)
					partSum += f(d);
	}
	t = MPI_Wtime() - t;
	for (int i = 0; i < size; i++)
	{
		if (i == rank)
			cout<<"rank"<<i<<": "<<t * 1000<<"ms"<<endl;
		MPI_Barrier(MPI_COMM_WORLD);
	}

	double sum;
	MPI_Allreduce(&partSum, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	if (rank == 0)
		cout<<"Result: "<<(sum + (f(b) - f(a)) / 2) * STEP<<endl;

    MPI_Finalize();
    return 0;
}