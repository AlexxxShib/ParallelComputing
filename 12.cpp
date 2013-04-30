#include <omp.h>
#include <mpi.h>
#include <algorithm>
#include <functional>
#include <iostream>
#include <vector>
#include <iomanip>
#include <cstring>
#include <fstream>
using namespace std;

const int Nx = 200;
const int Ny = 100;
const int MAXI = 1e4;
 
int getSize(int dataSize, int rank, int size)
{
	return dataSize / size + (rank < dataSize % size ? 1 : 0);
}

int getOffset(int dataSize, int rank, int size)
{
	return rank < dataSize % size ? (dataSize / size + 1) * rank : 
									(dataSize - (dataSize / size) * (size - rank));
}

enum Border
{
	Left,
	Right,
	Down,
	Up	
};

enum HandlePolicy
{
	Const,
	Flow,
};

struct BorderInfo
{
	HandlePolicy borderPolicy;
	double flowValue;
	double constValue;
	bool isContained;

	BorderInfo()
	{
		isContained = false;
	}
};

struct HeatSource
{
	int x;
	int y;
	double value;
};

struct ProblemData
{
	BorderInfo borderInfos[4];
	vector<HeatSource> heatSources;

	void setBoundaryCondition(Border border, double value)
	{
		borderInfos[border].borderPolicy = Const;
		borderInfos[border].constValue = value;
	}

	void setHeatSource(int i, int j, double value)
	{
		HeatSource heatSource = {i, j, value};
		heatSources.push_back(heatSource);
	}

	void setHeatFlow(Border border, double value)
	{
		borderInfos[border].borderPolicy = Flow;
		borderInfos[border].flowValue = value;
	}

};

struct BorderHandler
{
	function<void()> handlers[4];

	BorderHandler()
	{
	}

	BorderHandler(BorderInfo* infos, vector<vector<double> >& _A)
	{
		for (int i = 0; i < 4; i++)
		{
			double flowValue = infos[i].flowValue;
			double constValue = infos[i].constValue;

			if (infos[i].isContained)
			{
				if (infos[i].borderPolicy == Flow)
					switch((Border)i)
					{
					case Left:
						handlers[i] = [=, &_A]()
							{
								for (int i = 1; i < _A.size() - 1; i++)
									_A[i][1] = _A[i][2] + flowValue;
							};
						break;
					case Right:
						handlers[i] = [=, &_A]()
							{
								for (int i = 1; i < _A.size() - 1; i++)
									_A[i][_A.front().size() - 2] = _A[i][_A.front().size() - 3] + flowValue;
							};
						break;
					case Up:
						handlers[i] = [=, &_A]()
							{
								for (int i = 1; i < _A.front().size() - 1; i++)
									_A[1][i] = _A[2][i] + flowValue;
							};
						break;
					case Down:
						handlers[i] = [=, &_A]()
							{
								for (int i = 1; i < _A.front().size() - 1; i++)
									_A[_A.size() - 2][i] = _A[_A.size() - 3][i] + flowValue;
							};
						break;
					}
				else
					switch((Border)i)
					{
					case Left:
						handlers[i] = [=, &_A]()
							{
								for (int i = 1; i < _A.size() - 1; i++)
									_A[i][1] = constValue;
							};
						break;
					case Right:
						handlers[i] = [=, &_A]()
							{
								for (int i = 1; i < _A.size() - 1; i++)
									_A[i][_A.front().size() - 2] = constValue;
							};
						break;
					case Up:
						handlers[i] = [=, &_A]()
							{
								for (int i = 1; i < _A.front().size() - 1; i++)
									_A[1][i] = constValue;
							};
						break;
					case Down:
						handlers[i] = [=, &_A]()
							{
								for (int i = 1; i < _A.front().size() - 1; i++)
									_A[_A.size() - 2][i] = constValue;
							};
						break;
					}
			}
			else
				handlers[i] = []{};
		}
	}

	void handle()
	{
		for (int i = 0; i < 4; i++)
			handlers[i]();
	}
};

struct CommHelper
{
	int dest;
	vector<double> recvBuffer;
	vector<double> sendBuffer;

	MPI_Status recvStatus;
	MPI_Status sendStatus;
  	MPI_Request recvRequest;
  	MPI_Request sendRequest;

	function<void(vector<double>&)> readBorder;
	function<void(vector<double>&)> writeBorder;
	int rank;
	CommHelper()
	{
	}

	CommHelper(int dest, int borderSize,
		function<void(vector<double>&)> readBorder,
		function<void(vector<double>&)> writeBorder) :
		dest(dest),
		readBorder(readBorder),
		writeBorder(writeBorder)
	{
	
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		recvBuffer.resize(borderSize);
		sendBuffer.resize(borderSize);
	}

	void recv()
	{
		if (dest == -1)
			return;
		MPI_Irecv(&recvBuffer[0], recvBuffer.size(), MPI_DOUBLE, dest, 0, MPI_COMM_WORLD, &recvRequest);
	}

	void send()
	{
		if (dest == -1)
			return;
		readBorder(sendBuffer);
		MPI_Isend(&sendBuffer[0], sendBuffer.size(), MPI_DOUBLE, dest, 0, MPI_COMM_WORLD, &sendRequest);
	}


	void waitRecv()
	{
		if (dest == -1)
			return;
		MPI_Wait(&recvRequest, &recvStatus);
		writeBorder(recvBuffer);
	}

	void waitSend()
	{
		if (dest == -1)
			return;
		MPI_Wait(&sendRequest, &sendStatus);
	}
};

struct Solver
{
	vector<vector<double> > _A;
	vector<vector<double> > _f;

	MPI_Comm cartcomm;
	
	int size;
	int rank;

	int xLow;
	int xHi;
	int yLow;
	int yHi;
	int coords[2];

	int xSize;
	int ySize;
	int xOff;
	int yOff;

	int xFactor;
	int yFactor;

	BorderHandler borderHandler;

	CommHelper leftCommHelper;
	CommHelper rightCommHelper;
	CommHelper downCommHelper;
	CommHelper upCommHelper;

	Solver(ProblemData& problemData)
	{
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &size);
	   
		createCart(problemData);

		_A.resize(ySize + 2);
		_f.resize(ySize);

		for (int i = 0; i <= ySize + 1; i++)
			_A[i].resize(xSize + 2);

		for (int i = 0; i < ySize; i++)
			_f[i].resize(xSize);	

		for (int  i = 0; i < problemData.heatSources.size(); i++)
		{
			const HeatSource& src = problemData.heatSources[i];
			if (src.x >= xOff && src.x < xOff + xSize && src.y >= yOff && src.y < yOff + ySize) 
				_f[src.y - yOff + 1][src.x - xOff + 1] = src.value;
		}

		borderHandler = BorderHandler(problemData.borderInfos, _A);
		borderHandler.handle();

		initCommHelpers();
	}

	void createCart(ProblemData& problemData)
	{
		for (xFactor = sqrt(size); xFactor >= 1; xFactor--)
			if (size % xFactor == 0)
			{
				yFactor = size / xFactor;
				break;
			}
		int dims[2] = {yFactor, xFactor};
		int periods[2] = {0, 0};

	    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cartcomm);
		MPI_Cart_coords(cartcomm, rank, 2, coords);

		xSize = getSize(Nx, coords[1], xFactor);
		ySize = getSize(Ny, coords[0], yFactor);
		xOff = getOffset(Nx, coords[1], xFactor);
		yOff = getOffset(Ny, coords[0], yFactor);

		xLow = 1;
		yLow = 1;
		xHi = xSize;
		yHi = ySize;

		if (coords[1] == 0)
		{
			xLow++;
			problemData.borderInfos[Left].isContained = true;
		}
		if (coords[1] == xFactor - 1)
		{
			problemData.borderInfos[Right].isContained = true;
			xHi--;
		}
		if (coords[0] == 0)
		{
			problemData.borderInfos[Up].isContained = true;
			yLow++;
		}
		if (coords[0] == yFactor - 1)
		{
			problemData.borderInfos[Down].isContained = true;
			yHi--;
		}
		//cout<<rank<<" "<<xLow<<" "<<xHi<<" "<<yLow<<" "<<yHi<<" "<<coords[0]<<" "<<coords[1]<<endl;
	}

	void initCommHelpers()
	{
		int leftRank;
		int rightRank;
		int downRank;
		int upRank;

		MPI_Cart_shift(cartcomm, 1, -1, &rightRank, &leftRank); 
		MPI_Cart_shift(cartcomm, 0, -1, &downRank, &upRank); 

		//cout<<rank<<" "<<leftRank<<" "<<rightRank<<" "<<upRank<<" "<<downRank<<endl;

		leftCommHelper = CommHelper(leftRank, ySize, 
			[&](vector<double>& v)
			{
				for (int i = 0; i < ySize; i++)
					v[i] = _A[i + 1][1];
			},
			[&](vector<double>& v)
			{
				for (int i = 0; i < ySize; i++)
					_A[i + 1][0] = v[i];
			});

		rightCommHelper = CommHelper(rightRank, ySize, 
			[&](vector<double>& v)
			{
				for (int i = 0; i < ySize; i++)
					v[i] = _A[i + 1][xSize];
			},
			[&](vector<double>& v)
			{
				for (int i = 0; i < ySize; i++)
					_A[i + 1][xSize + 1] = v[i];
			});

		downCommHelper = CommHelper(downRank, xSize, 
			[&](vector<double>& v)
			{
				for (int i = 0; i < xSize; i++)
					v[i] = _A[ySize][i + 1];
			},
			[&](vector<double>& v)
			{
				for (int i = 0; i < xSize; i++)
					_A[ySize + 1][i + 1] = v[i];
			});

		upCommHelper = CommHelper(upRank, xSize, 
			[&](vector<double>& v)
			{
				for (int i = 0; i < xSize; i++)
					v[i] = _A[1][i + 1];
			},
			[&](vector<double>& v)
			{
				for (int i = 0; i < xSize; i++)
					_A[0][i + 1] = v[i];
			});
	}
	
	void solve()
	{
		rightCommHelper.send();
		downCommHelper.send();
		
		leftCommHelper.recv();
		upCommHelper.recv();
		for(int i = 0; i < MAXI; i++)
		{
			leftCommHelper.send();
			upCommHelper.send();

			leftCommHelper.waitRecv();
			upCommHelper.waitRecv();

			leftCommHelper.recv();
			rightCommHelper.recv();
			downCommHelper.recv();
			upCommHelper.recv();

			leftCommHelper.waitRecv();
			rightCommHelper.waitRecv();
			downCommHelper.waitRecv();
			upCommHelper.waitRecv();

			leftCommHelper.waitSend();
			downCommHelper.waitSend();
			rightCommHelper.waitSend();
			upCommHelper.waitSend();

			
			for (int d = xLow + yLow; d <= xHi + yHi; d++)
			{
				int x0 = (d > xHi + yLow) ? xHi : d - yLow;
				int y0 = d - x0;
				int offMax = min(x0 - xLow, yHi - y0);
				//#pragma omp parallel for num_threads(2)
				for (int off = 0; off <= offMax; off++)
				{
					int x = x0 - off;
					int y = y0 + off;
					_A[y][x] = 1. / 4 * (_A[y - 1][x] + _A[y + 1][x] + _A[y][x - 1] + _A[y][x + 1] + _f[y - 1][x - 1]);
				}
			}

			/*for (int y = yLow; y <= yHi; y++)
				for (int x = xLow; x <= xHi; x++)
					_A[y][x] = 1. / 4 * (_A[y - 1][x] + _A[y + 1][x] + _A[y][x - 1] + _A[y][x + 1] + _f[y - 1][x - 1]);
			*/

			rightCommHelper.send();
			downCommHelper.send();

			borderHandler.handle();

			if (rank == 0 && i % 10 == 0)
				cout<<"\r"<<i;
		}
	}

	void outData()
	{
		MPI_File fh;
		char fileName[256];
		strcpy(fileName,  "plot");
		MPI_File_open(MPI_COMM_WORLD, fileName, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

		MPI_Datatype blockType;
		MPI_Type_vector(ySize, xSize, Nx, MPI_DOUBLE, &blockType);
		MPI_Type_commit(&blockType);

		MPI_Offset disp = (MPI_Offset)(sizeof(double) * (xOff + Nx * yOff));
		char buffer[256];
		strcpy(buffer,  "native");
		MPI_File_set_view(fh, disp, MPI_DOUBLE, blockType, buffer, MPI_INFO_NULL);

		MPI_Datatype lineType;
		MPI_Type_contiguous(xSize, MPI_DOUBLE, &lineType);
		MPI_Type_commit(&lineType);

		vector<double> dummy(xSize);
		for (int i = 1; i <= ySize; i++)
			MPI_File_write(fh, &_A[i][1], 1, lineType, MPI_STATUS_IGNORE);
	}
};
int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

	ProblemData data;
	/*data.setHeatSource(Nx / 4, Ny / 2, 2.5);
	data.setHeatSource(3 * Nx / 4, Ny / 2, -0.5);
	data.setBoundaryCondition(Right, 1);
	data.setBoundaryCondition(Up, -0.5);
	data.setBoundaryCondition(Down, 0);
	data.setBoundaryCondition(Left, 0);*/

	/*data.setBoundaryCondition(Up, -.5);
	data.setBoundaryCondition(Down, .5);
	data.setBoundaryCondition(Left, 0);
	data.setBoundaryCondition(Right, 0);*/
	
	data.setHeatFlow(Up, 0);
	data.setHeatFlow(Down, 0);
	data.setHeatFlow(Left, 0.5);
	data.setHeatFlow(Right, -0.5);

	Solver solver(data);
	
	double t = MPI_Wtime();
	solver.solve();
	solver.outData();

	t = MPI_Wtime() - t;
	if (solver.rank == 0)
		cout<<endl;
	for (int i = 0; i < solver.size; i++)
	{
		if (i == solver.rank)
			cout<<"rank"<<i<<": "<<t * 1000<<"ms"<<endl;
		MPI_Barrier(MPI_COMM_WORLD);
	}

    MPI_Finalize();
    return 0;
}