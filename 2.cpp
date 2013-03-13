#include <mpi.h>
#include <iostream>
#include <vector>
#include <iomanip>
using namespace std;
 
const int N = 1000, M = N;
int getSize(int dataSize, int rank, int size)
{
    return dataSize / size + (rank < dataSize % size ? 1 : 0);
}
int getOffset(int dataSize, int rank, int size)
{
    return rank < dataSize % size ? (dataSize / size + 1) * rank : 
									(dataSize - (dataSize / size) * (size - rank));
}
void fill_A(vector<int> &row, int nRow)
{
	for (int i = 0; i < row.size(); i++)
		row[i] = nRow == i ? 2 : (nRow - 1 == i || nRow + 1 == i ? 1 : 0);
}
void fill_x(int &xElem, int n)
{
	xElem = 1;
}
void show()
{
	cout<<"A:\n";
	for (int i = 0; i < M; i++)
	{
		vector<int> row(N);
		fill_A(row, i);
		for (int j = 0; j < N; j++)
			cout<<setw(2)<<row[j]<<" ";
		cout<<endl;
	}
	cout<<"x:\n";
	vector<int> x(N);
	for (int i = 0; i < N; i++)
	{
		fill_x(x[i], i);
		cout<<setw(2)<<x[i]<<" ";
	}
	cout<<endl;
}
int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
	//if (rank == 0)
	//	show();
	int rowCount = getSize(M, rank, size);
	int rowOffset = getOffset(M, rank, size);
	vector<vector<int> > A(rowCount);
	vector<int> _b(rowCount);
	for (int i = 0; i < rowCount; i++)
	{
			A[i].resize(N);
			fill_A(A[i], i + rowOffset);
	}
	
	int xSize = getSize(N, rank, size);
	int xOffset = getOffset(N, rank, size);
	vector<int> x(N);
	for (int i = 0; i < N; i++)
		fill_x(x[i], i);
	
	//int nPrev = (rank - 1 + size) % size;
	//int nNext = (rank + 1) % size;
	double t = MPI_Wtime();
	for (int i = 0; i < rowCount; i++)
		for (int j = 0; j < N; j++)
			_b[i] += A[i][j] * x[j];
	t = MPI_Wtime() - t;
	for (int i = 0; i < size; i++)
	{
		if (i == rank)
			cout<<"rank"<<i<<": "<<t * 1000<<"ms"<<endl;
		MPI_Barrier(MPI_COMM_WORLD);
	}
	vector<int> b(M);
	vector<int> recvcounts(size);
	vector<int> displs(size);
	recvcounts[0] = getSize(M, 0, size);
	for (int i = 1; i < size; i++)
	{
		recvcounts[i] = getSize(M, i, size);
		displs[i] = displs[i - 1] + recvcounts[i - 1];
	}
	MPI_Allgatherv(&_b[0], rowCount, MPI_INT, &b[0], &recvcounts[0], &displs[0], MPI_INT, MPI_COMM_WORLD); 
	if (rank == 0)
	{
		cout<<"Ax:"<<endl;
		for (int i = 0; i < M; i++)
			cout<<setw(2)<<b[i]<<" ";
		cout<<endl;
	}
	MPI_Finalize();
    return 0;
}
