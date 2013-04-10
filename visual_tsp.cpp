#pragma comment(linker, "/manifestdependency:\"type='win32' name='Microsoft.Windows.Common-Controls' version='6.0.0.0' processorArchitecture='x86' publicKeyToken='6595b64144ccf1df' language='*'\"")

#include <windows.h>
#include <algorithm>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
using namespace std;
LRESULT CALLBACK WndProc(HWND hwnd, UINT iMsg, WPARAM wParam, LPARAM lParam);

struct TSP
{
	int N;
	vector<vector<double> > c;
	vector<int> bestVisited;
	double best;

	vector<int> visited1;
	vector<int> rangeN;
	vector<pair<int, int> > points;
				
	int _count;

	TSP(int N = 10)
	{
		best = DBL_MAX;
		_count = 0;
		init(N);
	}

	void init(int N)
	{
		this->N = N;
		bestVisited.resize(N);
		visited1.resize(N);
		rangeN.resize(N);
		c.resize(N);
		points.resize(N);
		for (int i = 0; i < N; i++)
		{
			rangeN[i] = i;
			c[i].resize(N);
		}

	}
	short randx()
	{
		static unsigned state = 1;
		state = state * 1664525 + 1013904223;
		return (state >> 16) & 0x7FFF;
	}
	void generate(string outFile, string pointsFileName = string())
	{
		const int W = 500;
		const int H = 500;
		for (int i = 0; i < N; i++)
		{
			points[i].first = randx() % W;
			points[i].second = randx() % H;
		}
		if (!pointsFileName.empty())
		{
			ofstream out(pointsFileName);
			out<<N<<endl;
			for (int i = 0; i < N; i++)
				out<<points[i].first<<" "<<points[i].second<<endl;
			out.close();
		}
		c.resize(N);
		for (int i = 0; i < N; i++)
			c[i].resize(N);
		for (int i = 0; i < N; i++)
		{
			
			for (int j = i + 1; j < N; j++)
				c[i][j] = c[j][i] = sqrt(1. * (points[i].first - points[j].first) * (points[i].first - points[j].first) + 
					(points[i].second - points[j].second) * (points[i].second - points[j].second));
		}

		ofstream out(outFile);
		out<<N;
		out<<endl;
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
				out<<setw(15)<<c[i][j];
			out<<endl;
		}
		out.close();
	}

	void load(string inFile)
	{
		ifstream in(inFile);
		in>>N;
		c.resize(N);
		for (int i = 0; i < N; i++)
		{
			c[i].resize(N);
			for (int j = 0; j < N; j++)
				in>>c[i][j];
		}
		in.close();
		init(N);
	}
};

TSP tsp(15);
//int bestVisited[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
int bestVisited[] = {0, 2, 5, 7, 10, 8, 3, 1, 14, 9, 12, 13, 11, 4, 6};

int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, PSTR szCmdLine, int iCmdShow)
{
	tsp.generate("1.txt");

	static wchar_t szAppName[] = L"Window";
	HWND hwnd;
	MSG msg;
	WNDCLASSEX wndclass;

	wndclass.cbSize = sizeof(wndclass);
	wndclass.style = CS_HREDRAW | CS_VREDRAW;
	wndclass.lpfnWndProc = WndProc;
	wndclass.cbClsExtra = 0;
	wndclass.cbWndExtra = 0;
	wndclass.hInstance = hInstance;
	wndclass.hIcon = LoadIcon(NULL, IDI_APPLICATION);
	wndclass.hCursor = LoadCursor(NULL, IDC_ARROW);
	wndclass.hbrBackground = (HBRUSH)GetStockObject(WHITE_BRUSH);
	wndclass.lpszMenuName = NULL;
	wndclass.lpszClassName = szAppName;
	wndclass.hIconSm = LoadIcon(NULL, IDI_APPLICATION);

	RegisterClassEx(&wndclass);
	hwnd=CreateWindow(szAppName, L"Window", WS_OVERLAPPEDWINDOW, CW_USEDEFAULT,
		CW_USEDEFAULT, CW_USEDEFAULT, CW_USEDEFAULT, NULL, NULL, hInstance, NULL);
	ShowWindow(hwnd, iCmdShow);
	UpdateWindow(hwnd);
	while (GetMessage(&msg, NULL, 0, 0))
	{
		TranslateMessage(&msg);
		DispatchMessage(&msg);
	}
	return msg.wParam;
}

double dd(int x0, int x1, int y0, int y1)
{
	return sqrt(1. * (x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
}

LRESULT CALLBACK WndProc(HWND hwnd, UINT iMsg, WPARAM wParam, LPARAM lParam)
{
	HDC hdc;
	PAINTSTRUCT ps;
	switch(iMsg)
	{
	case WM_CREATE:
		return 0;
	case WM_PAINT:
		{
			hdc = BeginPaint(hwnd, &ps);
			for (int i = 0; i < tsp.N; i++)
				Ellipse(hdc, tsp.points[i].first - 2, tsp.points[i].second - 2, tsp.points[i].first + 2, tsp.points[i].second + 2);

			int first = bestVisited[0];
			double dist = 0;
			MoveToEx(hdc, tsp.points[first].first, tsp.points[first].second, NULL);
			for (int i = 1; i < tsp.N; i++)
			{
				LineTo(hdc, tsp.points[bestVisited[i]].first, tsp.points[bestVisited[i]].second);
				dist += dd(tsp.points[bestVisited[i - 1]].first, 
							tsp.points[bestVisited[i]].first,
							tsp.points[bestVisited[i - 1]].second,
							tsp.points[bestVisited[i]].second);

			}

			LineTo(hdc, tsp.points[first].first, tsp.points[first].second);
			dist += dd(tsp.points[first].first, 
							tsp.points[bestVisited[tsp.N - 1]].first,
							tsp.points[first].second,
							tsp.points[bestVisited[tsp.N - 1]].second);
			stringstream ss;
			ss<<dist;
			SetWindowTextA(hwnd, ss.str().c_str());
			EndPaint(hwnd, &ps);
			return 0;
		}
	case WM_DESTROY:
		PostQuitMessage(0);
		return 0;
	}
	return DefWindowProc(hwnd, iMsg, wParam, lParam);
}