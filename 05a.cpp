#include <mpi.h>
#include <algorithm>
#include <functional>
#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
using namespace std;
int N = 100;

int getSize(int dataSize, int rank, int size)
{
	return dataSize / size + (rank < dataSize % size ? 1 : 0);
}
int getOffset(int dataSize, int rank, int size)
{
	return rank < dataSize % size ? (dataSize / size + 1) * rank : 
									(dataSize - (dataSize / size) * (size - rank));
}
void fill_A(vector<double> &row, int nRow)
{
	for (int i = 0; i < N; i++)
		row[i] = nRow == i ? 2 : (nRow - 1 == i || nRow + 1 == i ? 1 : 0);
}
void fill_b(double &xElem, int n)
{
	xElem = n == N - 1 ? 3 * (n + 1) - 1 : 4 * (n + 1);
}
struct Solver
{
	int rank;
	int size;
	vector<vector<double> > _A;
	vector<double> _x;
	int nCol;
	vector<int> sizes;
	vector<int> offsets;

	Solver()
	{
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &size);

		sizes.resize(size);
		offsets.resize(size);

		for (int i = 0; i < size; i++)
		{
			sizes[i] = getSize(N, i, size);
			offsets[i] = getOffset(N, i, size);
		}
		nCol = 0;

		_A.resize(sizes[rank]);
		_x.resize(sizes[rank]);

		for (int i = 0; i < sizes[rank]; i++)
		{
			_A[i].resize(N + 1);
			fill_A(_A[i], offsets[rank] + i);
			fill_b(_A[i][N], offsets[rank] + i);
		}
	}
	void solve()
	{
		for (int i = 0; i < size; i++)
		{
			if (i == rank)
				for (int j = 0; j < sizes[i]; j++, nCol++)
				{
					vector<double> &v = _A[j];
					transform(v.begin(), v.end(), v.begin(), bind1st(multiplies<double>(), 1 / v[nCol]));
					MPI_Bcast(&v[nCol], N - nCol + 1, MPI_DOUBLE, i, MPI_COMM_WORLD);
					for (int k = 0; k < j; k++)
					{
						vector<double> &u = _A[k];
						double d = u[nCol];
						for (int l = nCol; l <= N; l++)
							u[l] -= v[l] * d;
					}
					for (int k = j + 1; k < sizes[rank]; k++)
					{
						vector<double> &u = _A[k];
						double d = u[nCol];
						for (int l = nCol; l <= N; l++)
							u[l] -= v[l] * d;
					}
				}
			else
			{
				vector<double> v(N - nCol + 1);
				for (int j = 0; j < sizes[i]; j++, nCol++)
				{
					MPI_Bcast(&v[0], N - nCol + 1, MPI_DOUBLE, i, MPI_COMM_WORLD);
					for (int k = 0; k < sizes[rank]; k++)
					{
						vector<double> &u = _A[k];
						double d = u[nCol];
						for (int l = nCol; l <= N; l++)
							u[l] -= v[l - nCol] * d;
					}
				}
			}	
		}
	}
	void gather()
	{
		for (int i = 0; i < sizes[rank]; i++)
			_x[i] = _A[i][N];
		if (rank != 0)
			MPI_Gatherv(&_x[0], sizes[rank], MPI_DOUBLE, 
						NULL, NULL, NULL, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
		else
		{
			vector<double> x(N);
	
			MPI_Gatherv(&_x[0], sizes[rank], MPI_DOUBLE,  
				&x[0], &sizes[0], &offsets[0], MPI_DOUBLE, 0, MPI_COMM_WORLD); 
		
			cout<<"x:"<<endl;
			for (int i = 0; i < N; i++)
				cout<<setw(2)<<x[i]<<" ";
			cout<<endl;
		}
	}
};
void show()
{
	cout<<"A:\n";
	for (int i = 0; i < N; i++)
	{
		vector<double> row(N);
		fill_A(row, i);
		for (int j = 0; j < N; j++)
			cout<<setw(2)<<row[j]<<" ";
		cout<<endl;
	}
	cout<<"b:\n";
	vector<double> b(N);
	for (int i = 0; i < N; i++)
	{
		fill_b(b[i], i);
		cout<<setw(2)<<b[i]<<" ";
	}
	cout<<endl;
}
int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	if (argc > 1)
		N = atoi(argv[1]);
	Solver solver;
	//if (solver.rank == 0)
		//show();

	double t = MPI_Wtime();
	solver.solve();
	t = MPI_Wtime() - t;
	
	for (int i = 0; i < solver.size; i++)
	{
		if (i == solver.rank)
			cout<<"rank"<<i<<": "<<t * 1000<<"ms"<<endl;
		MPI_Barrier(MPI_COMM_WORLD);
	}
	//solver.gather();
	MPI_Finalize();
	return 0;
}