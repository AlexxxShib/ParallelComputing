#include <omp.h>
#include <algorithm>
#include <functional>
#include <iostream>
#include <vector>
#include <iomanip>
using namespace std;

const int N = 10;
const double tau = .5;
const double epsilon = 1e-6;

vector<vector<double> > A(N);
vector<double> b(N);
vector<double> x(N);

double normBSquared;

void fill_A(vector<double> &row, int nRow)
{
	for (int i = 0; i < row.size(); i++)
		row[i] = nRow == i ? 2 : (nRow - 1 == i || nRow + 1 == i ? 1 : 0);
}
void fill_b(double &xElem, int n)
{
	xElem = n == N - 1 ? 3 * (n + 1) - 1 : 4 * (n + 1);
}
void solve1()
{
	vector<double> buf(N);
	double normXSquared = 0.;
	#pragma omp parallel
	{
		for(x = b;;)
		{
			#pragma omp for
			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
					buf[i] += A[i][j] * x[j];

			#pragma omp for reduction(+ : normXSquared)
			for (int i = 0; i < N; i++)
			{
				buf[i] -= b[i];
				normXSquared += buf[i] * buf[i];
				x[i] -= buf[i] * tau;
			}
			if (normXSquared / normBSquared < epsilon * epsilon)
				break;	
			#pragma omp barrier
			#pragma omp single
			{
				buf.assign(buf.size(), 0.);	
				normXSquared = 0;
			}
		}
	}
}
void solve2()
{
	vector<double> buf(N);
	for(x = b;;)
	{
		#pragma omp parallel for
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				buf[i] += A[i][j] * x[j];

		double normXSquared = 0.;
		#pragma omp parallel for reduction(+ : normXSquared)
		for (int i = 0; i < N; i++)
		{
			buf[i] -= b[i];
			normXSquared += buf[i] * buf[i];
			x[i] -= buf[i] * tau;
			
		}
		if (normXSquared / normBSquared < epsilon * epsilon)
			break;	
			
		buf.assign(buf.size(), 0.);	
		normXSquared = 0;
	}
}
int main(int argc, char **argv)
{
	for (int i = 0; i < N; i++)
	{
		A[i].resize(N);
		fill_A(A[i], i);
		fill_b(b[i], i);
		normBSquared += b[i] * b[i];
	}
	double start = omp_get_wtime();
	solve2();
	cout<<"Time: "<<omp_get_wtime() - start<<" s"<<endl;
	for (int i = 0; i < N; i++)
		cout<<setw(2)<<x[i]<<" ";
    return 0;
}