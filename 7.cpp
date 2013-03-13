#include <mpi.h>
#include <algorithm>
#include <functional>
#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
using namespace std;

const int Nx = 100;
const int Ny = 50;
const int N = Nx * Ny;
 
int getSize(int dataSize, int rank, int size)
{
	return dataSize / size + (rank < dataSize % size ? 1 : 0);
}
int getOffset(int dataSize, int rank, int size)
{
	return rank < dataSize % size ? (dataSize / size + 1) * rank : 
									(dataSize - (dataSize / size) * (size - rank));
}
inline int tIndex(int i, int j)
{
	return i * Nx + j;
}

struct Solver
{
	enum Border
	{
		Left, Right, Top, Bottom
	};

	vector<vector<double> > _A;
	vector<double> _x;
	vector<int> sizes;
	vector<int> offsets;
	int rank;
	int size;

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

		_A.resize(sizes[rank]);
		_x.resize(sizes[rank]);
		for (int i = 0; i < sizes[rank]; i++)
			_A[i].resize(N + 1);
		
		fillA();
	}
	int getLocalIndex(int index)
	{
		if (index % size != rank)
			return -1;
		return index / size;
	}
	void fillA()
	{
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				setAValue(i , j, j == i ? -4 : (j == i - 1 || j == i + 1 || j == i + Nx || j == i - Nx ? 1 : 0));
	}
	void setAValue(int i, int j, double value)
	{
		int index = getLocalIndex(i);
		if (index == -1)
			return;
		_A[index][j] = value;
	}
	void setFValue(int i, double value)
	{
		int index = getLocalIndex(i);
		if (index == -1)
			return;
		_A[index][N] = value;
	}

	void setTValue(int ti, double value)
	{ 
		for (int i = 0; i < N; i++)
			setAValue(ti, i, 0);
		setAValue(ti, ti, 1);
		setFValue(ti, value);
	}
	void setBoundaryCondition(Border border, double value)
	{
		switch (border)
		{
		case Left:
			for (int i = 0; i < Ny; i++)
				setTValue(tIndex(i, 0), value);
			break;
		case Right:
			for (int i = 0; i < Ny; i++)
				setTValue(tIndex(i, Nx - 1), value);
			break;
		case Top:
			for (int i = 0; i < Nx; i++)
				setTValue(tIndex(0, i), value);
			break;
		case Bottom:
			for (int i = 0; i < Nx; i++)
				setTValue(tIndex(Ny - 1, i), value);
			break;
		}
	}
	void setHeatSource(int i, int j, double value)
	{
		setFValue(tIndex(i, j), -value);
	}
	void setHeatFlow(Border border, double value)
	{
		switch (border)
		{
		case Left:
			for (int i = 0; i < Ny; i++)
			{
				int ti = tIndex(i, 0);
				for (int j = 0; j < N; j++)
					setAValue(ti, j, 0);
				setAValue(ti, ti, 1);
				setAValue(ti, ti + 1, -1);
				setFValue(ti, value);
			}
			break;
		case Right:
			for (int i = 0; i < Ny; i++)
			{
				int ti = tIndex(i, Nx - 1);
				for (int j = 0; j < N; j++)
					setAValue(ti, j, 0);
				setAValue(ti, ti, 1);
				setAValue(ti, ti - 1, -1);
				setFValue(ti, value);
			}
			break;
		case Top:
			for (int i = 0; i < Nx; i++)
			{
				int ti = tIndex(0, i);
				for (int j = 0; j < N; j++)
					setAValue(ti, j, 0);
				setAValue(ti, ti, 1);
				setAValue(ti, ti + Nx, -1);
				setFValue(ti, value);
			}
			break;
		case Bottom:
			for (int i = 0; i < Nx; i++)
			{
				int ti = tIndex(Ny - 1, i);
				for (int j = 0; j < N; j++)
					setAValue(ti, j, 0);
				setAValue(ti, ti, 1);
				setAValue(ti, ti - Nx, -1);
				setFValue(ti, value);
			}
			break;
		}
	}
	void solve()
	{
		for (int i = 0; i < N; i++)
		{
			int mod = i % size;
			int div = i / size;
			if (mod == rank)
			{
				vector<double> &v = _A[div];
				transform(v.begin() + i, v.end(), v.begin() + i, bind1st(multiplies<double>(), 1 / v[i]));
				MPI_Bcast(&v[i], N - i + 1, MPI_DOUBLE, mod, MPI_COMM_WORLD);
				for (int j = div + 1; j < sizes[rank]; j++)
				{
					double d = _A[j][i];
					vector<double> &u = _A[j];
					for (int k = i; k <= N; k++)
						u[k] -= v[k] * d;
				}
			}
			else
			{
				vector<double> v(N - i + 1);
				MPI_Bcast(&v[0], N - i + 1, MPI_DOUBLE, mod, MPI_COMM_WORLD);
				for (int j = mod < rank ? div : div + 1; j < sizes[rank]; j++)
				{
					double d = _A[j][i];
					vector<double> &u = _A[j];
					for (int k = i; k <= N; k++)
						u[k] -= v[k - i] * d;
				}
			}
		}
		for (int i = N - 1; i >= 0; i--)
		{
			int mod = i % size;
			int div = i / size;
			if (mod == rank)
			{
				_x[div] = _A[div][N];
				MPI_Bcast(& _x[div], 1, MPI_DOUBLE, mod, MPI_COMM_WORLD);
				for (int j = 0; j < div;  j++)
					_A[j][N] -= _A[j][i] * _A[div][N];
			}
			else
			{
				double xi;
				MPI_Bcast(&xi, 1, MPI_DOUBLE, mod, MPI_COMM_WORLD);
				for (int j = 0, to = mod < rank ? div : div + 1; j < to; j++)
					_A[j][N] -= _A[j][i] * xi;
			}
		}
	}
	void gather()
	{
		if (rank != 0)
			MPI_Gatherv(&_x[0], sizes[rank], MPI_DOUBLE, 
						NULL, NULL, NULL, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
		else
		{
			vector<double> buf(N);
			vector<double> x(N);
			MPI_Gatherv(&_x[0], sizes[rank], MPI_DOUBLE,  
				&buf[0], &sizes[0], &offsets[0], MPI_DOUBLE, 0, MPI_COMM_WORLD); 
			
			for (int i = 0; i < N; i++)
				x[i] = buf[offsets[i % size] + i / size];

			/*cout<<"x:"<<endl;
			for (int i = 0; i < N; i++)
				cout<<setw(2)<<x[i]<<" ";
			cout<<endl;*/
			ofstream out("plot", ios::binary);
			out.write((char*)&Nx, sizeof(Nx));
			out.write((char*)&Ny, sizeof(Ny));
			out.write((char*)&x[0], N * sizeof(double));
		}
	}
};
int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

	Solver solver;
	{
		solver.setHeatSource(Ny / 2, Nx / 4, 2.5);
		solver.setHeatSource(Ny / 2, 3 * Nx / 4, -2.5);
		solver.setBoundaryCondition(Solver::Right, 1);
		solver.setHeatFlow(Solver::Top, -0.5);
		solver.setHeatFlow(Solver::Bottom, 0);
		solver.setHeatFlow(Solver::Left, 0);

		/*solver.setBoundaryCondition(Solver::Top, -.5);
		solver.setBoundaryCondition(Solver::Bottom, .5);
		solver.setBoundaryCondition(Solver::Left, 0);
		solver.setBoundaryCondition(Solver::Right, 0);
		
		solver.setHeatFlow(Solver::Top, 0);
		solver.setHeatFlow(Solver::Bottom, 0);
		solver.setHeatFlow(Solver::Left, 0.5);
		solver.setHeatFlow(Solver::Right, -0.5);*/
	}

	//solver.outSLE();
	double t = MPI_Wtime();
	solver.solve();
	t = MPI_Wtime() - t;
	for (int i = 0; i < solver.size; i++)
	{
		if (i == solver.rank)
			cout<<"rank"<<i<<": "<<t * 1000<<"ms"<<endl;
		MPI_Barrier(MPI_COMM_WORLD);
	}
	solver.gather();

    MPI_Finalize();
    return 0;
}