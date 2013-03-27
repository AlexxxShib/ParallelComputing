#pragma comment(linker, "/manifestdependency:\"type='win32' name='Microsoft.Windows.Common-Controls' version='6.0.0.0' processorArchitecture='x86' publicKeyToken='6595b64144ccf1df' language='*'\"")

#include <windows.h>
#include <vector>
#include <fstream>
#include <algorithm>
using namespace std;
LRESULT CALLBACK WndProc(HWND hwnd, UINT iMsg, WPARAM wParam, LPARAM lParam);

vector<double> v;
double _max, _min;
int Nx, Ny;
const int cellLength = 8;
int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, PSTR szCmdLine, int iCmdShow)
{
	static wchar_t szAppName[] = L"Window";
	HWND hwnd;
	MSG msg;
	WNDCLASSEX wndclass;

	ifstream in("H:\\plot", ios::binary);
	in.read((char*)&Nx, sizeof(Nx));
	in.read((char*)&Ny, sizeof(Ny));
	v.resize(Nx * Ny);
	in.read((char*)&v[0], Nx * Ny * sizeof(double));
	_max = *max_element(v.begin(), v.end());
	_min = *min_element(v.begin(), v.end());
	/*vector<vector<double>> vv(Ny);
	for (int i = 0; i < Ny; i++)
		vv[i] = vector<double>(v.begin() + i * Nx, v.begin() + i * Nx + Nx);*/

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
		CW_USEDEFAULT, Nx * cellLength + 20,  Ny * cellLength + 40, NULL, NULL, hInstance, NULL);

	ShowWindow(hwnd, iCmdShow);
	UpdateWindow(hwnd);

	while (GetMessage(&msg, NULL, 0, 0))
	{
		TranslateMessage(&msg);
		DispatchMessage(&msg);
	}
	return msg.wParam;
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
		hdc = BeginPaint(hwnd, &ps);
		for (int i = 0; i < Ny; i++)
			for (int j = 0; j < Nx; j++)
			{
				int c = 255 - 255 * pow((v[i * Nx + j] - _min) / (_max - _min), 1);
				SelectObject(hdc, GetStockObject(NULL_PEN));
				HBRUSH br = CreateSolidBrush(RGB(255 - c, 0, c));
				SelectObject(hdc, br);
				Rectangle(hdc, cellLength * j, cellLength * i, cellLength * (j + 1) + 1, cellLength * (i + 1) + 1);
				DeleteObject(br);
			}
		EndPaint(hwnd, &ps);
		return 0;
	case WM_DESTROY:
		PostQuitMessage(0);
		return 0;
	}
	return DefWindowProc(hwnd, iMsg, wParam, lParam);
}