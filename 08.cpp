#include <omp.h>
#include <algorithm>
#include <functional>
#include <iostream>
#include <vector>
#include <iomanip>
using namespace std;

const int M = 5000;
const int N = M;
vector<vector<int> > A(M);
vector<int> b(M);
vector<int> x(N, 1);

void fill_A(vector<int> &row, int nRow)
{
	for (int i = 0; i < row.size(); i++)
		row[i] = nRow == i ? 2 : (nRow - 1 == i || nRow + 1 == i ? 1 : 0);
}
void mul1()
{
	#pragma omp parallel
	{
		#pragma omp for
		for (int i = 0; i < M; i++)
			for (int j = 0; j < N; j++)
				b[i] += A[i][j] * x[j];
	}
}
void mul2()
{
	#pragma omp parallel for
	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
			b[i] += A[i][j] * x[j];
}
void mul3()
{
	#pragma omp parallel
	{
		int size = omp_get_num_threads();
		int rank = omp_get_thread_num();
		int blockSize = M / size + (rank < M % size ? 1 : 0);
		int blockOffset = rank < M % size ? (M / size + 1) * rank : 
										 (M - (M / size) * (size - rank));

		for (int i = blockOffset; i < blockOffset + blockSize; i++)
			for (int j = 0; j < N; j++)
				b[i] += A[i][j] * x[j];
	}
}
void mul4()
{
	const int size = 3;
	#pragma omp parallel sections 
	{
		#pragma omp section
		{
			int rank = 0;
			int blockSize = M / size + (rank < M % size ? 1 : 0);
			int blockOffset = rank < M % size ? 
				(M / size + 1) * rank : 
				(M - (M / size) * (size - rank));

			for (int i = blockOffset; i < blockOffset + blockSize; i++)
				for (int j = 0; j < N; j++)
					b[i] += A[i][j] * x[j];
		}
		#pragma omp section
		{
			int rank = 1;
			int blockSize = M / size + (rank < M % size ? 1 : 0);
			int blockOffset = rank < M % size ? 
				(M / size + 1) * rank : 
				(M - (M / size) * (size - rank));

			for (int i = blockOffset; i < blockOffset + blockSize; i++)
				for (int j = 0; j < N; j++)
					b[i] += A[i][j] * x[j];
		}
		#pragma omp section
		{
			int rank = 2;
			int blockSize = M / size + (rank < M % size ? 1 : 0);
			int blockOffset = rank < M % size ? 
				(M / size + 1) * rank : 
				(M - (M / size) * (size - rank));

			for (int i = blockOffset; i < blockOffset + blockSize; i++)
				for (int j = 0; j < N; j++)
					b[i] += A[i][j] * x[j];
		}

	}
}
int main(int argc, char **argv)
{
	for (int i = 0; i < M; i++)
	{
		A[i].resize(N);
		fill_A(A[i], i);
	}
	double start = omp_get_wtime();
	mul4();
	cout<<"Time: "<<omp_get_wtime() - start<<" s"<<endl;
	for (int i = 0; i < M; i++)
		//cout<<setw(2)<<b[i]<<" ";
    return 0;
}