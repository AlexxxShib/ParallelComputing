#include <mpi.h>
#include <iostream>
#include <vector>
#include <cstdlib>  
#include <iomanip>
#include <fstream>
using namespace std;

int Nx = 40;
int Ny = 40;
int N = Nx * Ny;
int MAX_ITERATIONS = 1e3;
const double epsilon = 1e-8;
double tau;
 
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
void mulPart(int rank, int size, vector<vector<double> > &A, vector<double> _x, vector<double> &res)
{
	int initSize = _x.size();
	int blockSize = getSize(N, rank, size);
	int blockOffset = getOffset(N, rank, size);
	const int rowCount = blockSize;
	_x.resize(N / size + 1);
	
	int nPrev = (rank - 1 + size) % size;
	int nNext = (rank + 1) % size;
	if (size == 1)
		for (int i = 0; i < size; i++)
			for (int j = 0; j < rowCount; j++)
			{
				vector<double> &v = A[j];
				for (int k = 0; k < blockSize; k++)
					res[j] += v[blockOffset + k] * _x[k];
			}
	else
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < rowCount; j++)
			{
				vector<double> &v = A[j];
				for (int k = 0; k < blockSize; k++)
					res[j] += v[blockOffset + k] * _x[k];
			}
			int curSize = blockSize;
			blockSize = getSize(N, (nPrev - i + size) % size, size);
			blockOffset = getOffset(N, (nPrev - i + size) % size, size);
			vector<double> __x(blockSize);
			MPI_Sendrecv(&_x[0], curSize, MPI_DOUBLE, nNext, 0, 
						 &__x[0], blockSize, MPI_DOUBLE, nPrev, 0,
						 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			_x = __x;
			_x.resize(N / size + 1);
		}
	_x.resize(initSize);
}
struct Solver
{
	enum Border
	{
		Left, Right, Top, Bottom
	};

	vector<vector<double> > A;
	vector<double> _t;
	vector<double> _f;
	int rank;
	int size;
	int blockSize;
	int blockOffset;
	double _normFSquared;
	double normFSquared;

	Solver()		
	{
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &size);

		blockSize = getSize(N, rank, size);
		blockOffset = getOffset(N, rank, size);

		A.resize(blockSize);
		//_t.resize(blockSize);
		_f.resize(blockSize);
		
		for (int i = 0; i < blockSize; i++)
			A[i].resize(N);
		
		fillA();
	}
	void outSLE()
	{
		if (rank == 0)
		{
			ofstream out("sle", ios::binary);
			out.write((char*)&N, sizeof(N));
		}
		for (int i = 0; i < size; i++)
		{
			if (i == rank)
			{
				ofstream out("sle", ios::binary | ios::app);
				for (int j = 0; j < blockSize; j++)
					out.write((char*)&A[j][0], N * sizeof(double));
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		if (rank != 0)
			gather(_f);
		else
		{
			vector<double> f(N);
			gather(_f, &f);
			ofstream out("sle", ios::binary | ios::app);
			out.write((char*)&f[0], N * sizeof(double));
		}
		//exit(0);
	}
	void fillA()
	{
		for (int i = blockOffset; i < blockOffset + blockSize; i++)
			for (int j = 0; j < N; j++)
				A[i - blockOffset][j] = j == i ? -4 : (j == i - 1 || j == i + 1 || j == i + Nx || j == i - Nx ? 1 : 0);
	}
	void setAValue(int i, int j, double value)
	{
		if (i >= blockOffset && i < blockOffset + blockSize)
			A[i - blockOffset][j] = value;
	}
	void setFValue(int i, double value)
	{
		if (i >= blockOffset && i < blockOffset + blockSize)
			_f[i - blockOffset] = value;
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
		int ti = tIndex(i, j);
		if (ti >= blockOffset && ti < blockOffset + blockSize)
			_f[ti - blockOffset] = -value;
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
	void solveSLE()
	{
		_normFSquared = .0;
		for (int i = 0; i < blockSize; i++)
			_normFSquared += _f[i] * _f[i];
		MPI_Allreduce(&_normFSquared, &normFSquared, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		_t = _f;
		//_t.resize(blockSize);
		//_t[0] = 1;

		for(int n = 0; n < MAX_ITERATIONS; n++)
		{
			vector<double> _r(blockSize);
			mulPart(rank, size, A, _t, _r);
			for (int i = 0; i < blockSize; i++)
				_r[i] -= _f[i];

			vector<double> _Ar(blockSize);
			mulPart(rank, size, A, _r, _Ar);

			double _normAr = .0, normAr, _Arr = .0, Arr;
			double _normRSquared = .0, normRSquared;
			for (int i = 0; i < blockSize; i++)
			{
				_normAr += _Ar[i] * _Ar[i];
				_Arr += _Ar[i] * _r[i];
				_normRSquared += _r[i] * _r[i];
			}
			MPI_Allreduce(&_Arr, &Arr, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&_normAr, &normAr, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&_normRSquared, &normRSquared, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

			if (Arr == 0)
				break;
			tau =  Arr / normAr;

			for (int i = 0; i < blockSize; i++)
				_t[i] -= _r[i] * tau;	
			
			//if (rank == 0 && n % 100 == 0)
			//	cout<<"\r"<<n<<" "<<normRSquared / normFSquared;
			//if (normRSquared / normFSquared < epsilon * epsilon)
			//	break;		
		}
		//if (rank == 0)
		//	cout<<endl;
	}
	void gather(vector<double> &_v, vector<double> *v = NULL)
	{
		if (rank != 0)
			MPI_Gatherv(&_v[0], blockSize, MPI_DOUBLE, 
						NULL, NULL, NULL, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
		else
		{
			vector<int> recvcounts(size);
			vector<int> displs(size);
			recvcounts[0] = getSize(N, 0, size);
			for (int i = 1; i < size; i++)
			{
				recvcounts[i] = getSize(N, i, size);
				displs[i] = displs[i - 1] + recvcounts[i - 1];
			}
			MPI_Gatherv(&_v[0], blockSize, MPI_DOUBLE,  
						&(*v)[0], &recvcounts[0], &displs[0], MPI_DOUBLE, 
						0, MPI_COMM_WORLD);
		}
	}
};
int main(int argc, char **argv)
{
	if (argc > 2)
	{
		Nx = atoi(argv[1]);
		Ny = Nx;
		MAX_ITERATIONS = atoi(argv[2]);
	}
	
	N = Nx * Ny;
    MPI_Init(&argc, &argv);

	Solver solver;
	{
		/*solver.setHeatSource(Ny / 2, Nx / 4, 250);
		solver.setHeatSource(Ny / 2, 3 * Nx / 4, -250);
		solver.setBoundaryCondition(Solver::Right, 1);
		solver.setHeatFlow(Solver::Top, -0.5);
		solver.setHeatFlow(Solver::Bottom, 0);
		solver.setHeatFlow(Solver::Left, 0);*/

		/*solver.setBoundaryCondition(Solver::Top, 0);
		solver.setBoundaryCondition(Solver::Bottom, .5);
		solver.setBoundaryCondition(Solver::Left, 0);
		solver.setBoundaryCondition(Solver::Right, 0);*/
		
		solver.setHeatFlow(Solver::Top, 10);
		solver.setHeatFlow(Solver::Bottom, 0);
		solver.setHeatFlow(Solver::Left, 0.5);
		solver.setHeatFlow(Solver::Right, -0.5);
	}

	solver.outSLE();
	//return 0;
	double t = MPI_Wtime();
	solver.solveSLE();
	t = MPI_Wtime() - t;
	for (int i = 0; i < solver.size; i++)
	{
		if (i == solver.rank)
			cout<<"rank"<<solver.rank<<": "<<t / MAX_ITERATIONS<<"s"<<endl;
		MPI_Barrier(MPI_COMM_WORLD);
	}

	if (solver.rank != 0)
		solver.gather(solver._t);
	else
	{
		vector<double> v(N);
		solver.gather(solver._t, &v);
		ofstream out("plot", ios::binary);
		out.write((char*)&Nx, sizeof(Nx));
		out.write((char*)&Ny, sizeof(Ny));
		out.write((char*)&v[0], N * sizeof(double));
	}
	
    MPI_Finalize();
    return 0;
}