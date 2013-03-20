#include <mpi.h>
#include <iostream>
#include <vector>
#include <iomanip>
using namespace std;
 
const int N = 2, M = 2;
int getSize(int totalSize, int rank, int size)
{
    return totalSize / size + (rank < totalSize % size ? 1 : 0);
}
int getOffset(int totalSize, int rank, int size)
{
    return rank < (totalSize % size) ? (totalSize / size + 1) * rank : (totalSize - (totalSize / size) * (size - rank));
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
	
	if (rank == 0)
		show();
	int rowCount = getSize(M, rank, size);
	int rowOffset = getOffset(M, rank, size);
	vector<vector<int> > A(rowCount);
	vector<int> _b(rowCount);
	
	for (int i = 0; i < rowCount; i++)
	{
			A[i].resize(N);
			fill_A(A[i], i + rowOffset);
	}
	vector<int> blockSizes(size);
	vector<int> blockOffsets(size);

	for (int i = 0; i < size; i++)
	{
		blockSizes[i] = getSize(N, i, size);
		blockOffsets[i] = getOffset(N, i, size);
	}
	int index = rank;
	vector<int> _x(N / size + 1);
	for (int i = 0; i < blockSizes[index]; i++)
			fill_x(_x[i], i + blockOffsets[index]);
			
	vector<int> __x(N / size + 1);
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < rowCount; j++)
			for (int k = 0; k < blockSizes[index]; k++)
				_b[j] += A[j][blockOffsets[index] + k] * _x[k];
				
		
		MPI_Sendrecv(&_x[0], blockSizes[index], MPI_INT, (rank + 1) % size, 0, 
					 &__x[0], blockSizes[(index + size - 1) % size], MPI_INT, (rank + size - 1) % size, 0,
					 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		index = (index + size - 1) % size;		
		_x = __x;
	}

	if (rank != 0)
		MPI_Gatherv(&_b[0], rowCount, MPI_INT,  NULL, NULL, NULL, MPI_INT, 0, MPI_COMM_WORLD); 
	else
	{
		vector<int> b(M);
		vector<int> recvcounts(size);
		vector<int> displs(size);
		recvcounts[0] = getSize(M, 0, size);
		for (int i = 1; i < size; i++)
		{
			recvcounts[i] = getSize(M, i, size);
			displs[i] = displs[i - 1] + recvcounts[i - 1];
		}
		MPI_Gatherv(&_b[0], rowCount, MPI_INT,  
					&b[0], &recvcounts[0], &displs[0], MPI_INT, 
					0, MPI_COMM_WORLD); 
		cout<<"Ax:"<<endl;
		for (int i = 0; i < M; i++)
			cout<<setw(2)<<b[i]<<" ";
		cout<<endl;
	}
    MPI_Finalize();
    return 0;
}