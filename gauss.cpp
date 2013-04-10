#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
using namespace std;
const double INF = numeric_limits<double>::infinity();

//mxn mx1
int gauss(vector<vector<double> > &a, vector<double> &f, vector<double> & ans)
{
	int n = (int)a.size();
	int m = (int)a[0].size();
 
	vector<int> where (m, -1);
	for (int col = 0, row = 0; col < m && row < n; ++col)
	{
		int sel = row;
		for (int i = row; i < n; ++i)
			if (abs(a[i][col]) > abs (a[sel][col]))
				sel = i;
		if (abs(a[sel][col]) < DBL_EPSILON)
			continue;
		for (int i = col; i < m; ++i)
			swap(a[sel][i], a[row][i]);
		swap(f[sel], f[row]);
		where[col] = row;
 
		for (int i = 0; i < n; ++i)
			if (i != row)
			{
				double c = a[i][col] / a[row][col];
				for (int j = col; j < m; ++j)
					a[i][j] -= a[row][j] * c;
				f[i] -= f[row] * c;
			}
		++row;
	}
 
	ans.assign(m, 0);
	for (int i = 0; i < m; ++i)
		if (where[i] != -1)
			ans[i] = f[where[i]] / a[where[i]][i];
	for (int i = 0; i < n; ++i)
	{
		double sum = 0;
		for (int j = 0; j < m; ++j)
			sum += ans[j] * a[i][j];
		if (abs(sum - f[i]) > DBL_EPSILON)
			return 0;
	}
 
	for (int i = 0; i < m; ++i)
		if (where[i] == -1)
			return INF;
	return 1;
}
int main(int argc, char **argv)
{
	ifstream in("sle", ios::binary);
	int N;
	in.read((char*)&N, sizeof(N));
	vector <vector<double> > A(N);
	for (int i = 0; i < N; ++i)
	{
		A[i].resize(N);
		in.read((char*)&A[i][0], N * sizeof(double));
	}
	vector<double> f(N);
	in.read((char*)&f[0], N * sizeof(double));

	vector<double> x(N);
	gauss(A, f, x);

	int Nx = 50;
	int Ny = N / Nx + .5;

	ofstream out("plot", ios::binary);
	out.write((char*)&Nx, sizeof(Nx));
	out.write((char*)&Ny, sizeof(Ny));
	out.write((char*)&x[0], N * sizeof(double));
    return 0;
}