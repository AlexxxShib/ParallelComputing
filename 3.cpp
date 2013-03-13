#include <mpi.h>
#include <iostream>
#include <vector>
#include <iomanip>
using namespace std;

const int N = 10;
const double tau = .5;
const double epsilon = 1e-4;
 
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
	for (int i = 0; i < row.size(); i++)
		row[i] = nRow == i ? 2 : (nRow - 1 == i || nRow + 1 == i ? 1 : 0);
}
void fill_b(double &xElem, int n)
{
	xElem = n == N - 1 ? 3 * (n + 1) - 1 : 4 * (n + 1);
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
				for (int k = 0; k < blockSize; k++)
					res[j] += A[j][blockOffset + k] * _x[k];
	else
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < rowCount; j++)
				for (int k = 0; k < blockSize; k++)
					res[j] += A[j][blockOffset + k] * _x[k];
			
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
void gather(int rank, int size, int blockSize, vector<double> &_x)
{
	if (rank != 0)
		MPI_Gatherv(&_x[0], blockSize, MPI_DOUBLE, 
					NULL, NULL, NULL, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
	else
	{
		vector<double> x(N);
		vector<int> recvcounts(size);
		vector<int> displs(size);
		recvcounts[0] = getSize(N, 0, size);
		for (int i = 1; i < size; i++)
		{
			recvcounts[i] = getSize(N, i, size);
			displs[i] = displs[i - 1] + recvcounts[i - 1];
		}
		MPI_Gatherv(&_x[0], blockSize, MPI_DOUBLE,  
					&x[0], &recvcounts[0], &displs[0], MPI_DOUBLE, 
					0, MPI_COMM_WORLD); 
		cout<<"_:\n";
		for (int i = 0; i < N; i++)
			cout<<setw(2)<<x[i]<<" ";
		cout<<endl;
	}
}
int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
	if (rank == 0)
		show();
	int blockSize = getSize(N, rank, size);
	int blockOffset = getOffset(N, rank, size);

	vector<vector<double> > A(blockSize);
	vector<double> _b(blockSize);
	vector<double> _x;
	double partNormBSquared = .0;
	double normBSquared;
	for (int i = 0; i < blockSize; i++)
	{
		A[i].resize(N);
		fill_A(A[i], i + blockOffset);
		fill_b(_b[i], i + blockOffset);
		partNormBSquared += _b[i] * _b[i];
	}
	_x = _b;
	MPI_Allreduce(&partNormBSquared, &normBSquared, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	for(;;)
	{
		vector<double> buf(blockSize);
		mulPart(rank, size, A, _x, buf);
		double partNormXSquared = .0;
		double normXSquared;
		for (int i = 0; i < blockSize; i++)
		{
			buf[i] -= _b[i];
			partNormXSquared += buf[i] * buf[i];
			_x[i] -= buf[i] * tau;
		}
		//gather(rank, size, blockSize, buf);
		MPI_Allreduce(&partNormXSquared, &normXSquared, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		//gather(rank, size, blockSize, _x);
		if (normXSquared / normBSquared < epsilon * epsilon)
			break;		
	}
	gather(rank, size, blockSize, _x);
    MPI_Finalize();
    return 0;
}