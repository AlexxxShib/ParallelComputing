#include <algorithm>
#include <fstream>
#include <vector>
#include <stack>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
using namespace std;

struct TSP
{
	int N;
	vector<vector<double> > c;
	vector<int> bestVisited;
	double best;

	vector<int> visited1;
	vector<int> rangeN;
				
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
		for (int i = 0; i < N; i++)
		{
			rangeN[i] = i;
			c[i].resize(N);
		}

	}
	void generate(string outFile, string pointsFileName = string())
	{
		const int W = 500;
		const int H = 500;
		srand(time(NULL));
		vector<pair<int, int> > points(N);
		for (int i = 0; i < N; i++)
		{
			points[i].first = rand() % W;
			points[i].second = rand() % H;
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

	void solve(void(TSP::*aux)(vector<int>&, double, int))
	{
		double initDist = visitedDistance();
		vector<int> visitedCopy = visited1;
		sort(visitedCopy.begin(), visitedCopy.end());
		vector<int> unvisited(N);
		vector<int>::iterator it = set_difference(rangeN.begin(), rangeN.end(),
			visitedCopy.begin(), visitedCopy.end(), unvisited.begin());
		unvisited.resize(it - unvisited.begin()); 
		(this->*aux)(unvisited, initDist, visited1.back());
		::_count = _count;
	}

	void execSolve(void(TSP::*aux)(vector<int>&, double, int))
	{
		greedy();
		for (int i = 1; i < N; i++)
		{
			visited1.clear();
			visited1.push_back(0);
			visited1.push_back(i);
			double dist = c[0][i];

			vector<int> visitedCopy = visited1;
			sort(visitedCopy.begin(), visitedCopy.end());
			vector<int> unvisited(N);
			vector<int>::iterator it = set_difference(rangeN.begin(), rangeN.end(),
				visitedCopy.begin(), visitedCopy.end(), unvisited.begin());
			unvisited.resize(it - unvisited.begin()); 
			(this->*aux)(unvisited, dist, i);
		}
	}

	void solveRec()
	{
		execSolve(&TSP::solveRecAux);
		cout<<_count<<endl;
		cout<<"best: "<<best<<endl;
	}

	void solveRecAux(vector<int>& unvisited, double dist, int lastVisited)
	{
		if (unvisited.empty())
		{
			/*cout<<endl;
			cout<<visited1[0]<<" "<<visited1[1]<<" "<<visited1[2]<<" "<<visited1[3]<<" "<<visited1[4]<<endl;*/
			_count++;
			best = min(best, dist + c[lastVisited][0]);
			bestVisited = visited1;
			return;
		}
		if (dist > best)
			return;
		for (int i = 0; i < unvisited.size(); i++)
		{
			int visiting = unvisited[i];
			dist += c[lastVisited][visiting];
			swap(unvisited[i], unvisited.back());
			unvisited.pop_back();

			visited1.push_back(visiting);
			solveRecAux(unvisited, dist, visiting);
			visited1.pop_back();
			dist -= c[lastVisited][visiting];
			unvisited.push_back(visiting);
			swap(unvisited[i], unvisited.back());
		}
	}

	void solveNoRec()
	{
		execSolve(&TSP::solveNoRecAux);
		cout<<_count<<endl;
		cout<<"best: "<<best<<endl;
	}

	struct State
	{
		int id;
		int index;
		int lastVisited;
		double dist;

		int visiting;
	};
	
	void solveNoRecAux(vector<int> &unvisited, double dist, int lastVisited)
	{
		int initVisited = visited1.size();
		stack<State> s;
		State state = {0, 0, lastVisited, dist};
		for(;;)
		{
			switch (state.id)
			{
			case 0:
				if (unvisited.empty())
				{
					_count++;
					if (state.dist + c[state.lastVisited][0] < best)
					{
						best = state.dist + c[state.lastVisited][0];
						bestVisited = visited1;
					}
					if (s.empty())
						return;
					state = s.top();
					s.pop();
					continue;
				}
				if (state.dist > best)
				{
					if (s.empty())
						return;
					state = s.top();
					s.pop();
					continue;
				}
			case 1:
				{
					int i = state.index;
					int lastVisited = state.lastVisited;
					int visiting = unvisited[i];
					double dist = state.dist;

					state.id = 2;
					state.visiting = visiting;

					s.push(state);

					swap(unvisited[i], unvisited.back());
					unvisited.pop_back();

					state.id = 0;
					state.index = 0;
					state.lastVisited = visiting;
					state.dist += c[lastVisited][visiting];

					visited1.push_back(visiting);
					continue;
				}
			case 2:
				{
					visited1.pop_back();

					int i = state.index;
					int lastVisited = state.lastVisited;
					int visiting = state.visiting;

					unvisited.push_back(visiting);
					swap(unvisited[i], unvisited.back());

					state.id = 1;
					state.index++;
					if (state.index != unvisited.size())
						continue;

					if (s.empty())
						return;
					state = s.top();
					s.pop();
					continue;
				}
			}
		}
	}

	double visitedDistance()
	{
		double dist = 0;
		for (int i = 1; i < visited1.size(); i++)
			dist += c[visited1[i - 1]][visited1[i]];
		return dist;
	}

	double greedy()
	{
		double best = 0;
		vector<bool> walked(N);
		walked[0] = true;
		bestVisited.push_back(0);
		int prev = 0;
		for (int i = 1; i < N; i++)
		{
			double minValue = DBL_MAX;
			int minIndex = -1;
			for (int j = 0; j < N; j++)
				if (!walked[j] && c[prev][j] < minValue)
				{
					minValue = c[prev][j];
					minIndex = j;
				}
			best += c[prev][minIndex];
			prev = minIndex;
			walked[minIndex] = true;
			bestVisited.push_back(minIndex);
			best += c[prev][0];
			return best;
		}
	}

	vector<int>::iterator skipPermut(vector<int>& v, int n)
	{
		int j;
		int index = -1;
		for (j = n - 1; j >= 0 && index == -1; j--)
		{
			int value = INT_MAX;
			for (int k = j + 1; k < N; k++)
				if (v[k] > v[j] && v[k] < value)
				{
					index = k;
					value = v[k];
				}
		}
		j++;
		if (index == -1)
			return v.end();
		swap(v[j], v[index]);
		sort(v.begin() + j + 1, v.end());
		return v.begin();
	}

	void solvePermut()
	{
		for (int i = 0; i < N; i++)
			visited1[i] = i;
		swap(visited1.front(), visited1.back());

		do
		{
begin:
			double dist = 0;
			for (int i = 1; i < N; i++)
			{
				if (dist > best)
				{
					//
					if (skipPermut(visited1, i) == visited1.end())
						goto end;
					goto begin;
				}
				dist += c[visited1[i - 1]][visited1[i]];
			}
			dist += c[visited1[N - 1]][visited1[0]];
			//log0<<visited[0]<<" "<<visited[1]<<" "<<visited[2]<<" "<<visited[3]<<" "<<visited[4]<<" "<<visited[5]<<endl;
			_count++;
			if (dist < best)
			{
				best = dist;
				bestVisited = visited1;
			}
		}
		while (next_permutation(visited1.begin(), visited1.end()));
end:
		cout<<_count<<endl;
		cout<<"best: "<<best<<endl;
	}

	void solveBruteForce()
	{
		for (int i = 0; i < N; i++)
			visited1[i] = i;
		swap(visited1.front(), visited1.back());

		do
		{
			_count++;
			double dist = c[visited1[N - 1]][visited1[0]];
			for (int i = 1; i < N; i++)
				dist += c[visited1[i - 1]][visited1[i]];
			if (dist < best)
			{
				best = dist;
				bestVisited = visited1;
			}
		}
		while (next_permutation(visited1.begin(), visited1.end()));
		cout<<_count<<endl;
		cout<<"best: "<<best<<endl;
	}
};

int main()
{
	TSP tsp(10);
	tsp.generate("tsp.txt");
	tsp.load("tsp.txt");
	time_t t = clock();
	tsp.greedy();
	tsp.solve(&TSP::solveRecAux);
	cout<<((clock() - t) / 1000.)<<endl;
	return 0;
}