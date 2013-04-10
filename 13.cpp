#include <sstream>
#include <unistd.h>
#include <mpi.h>
#include <map>
#include <stack>
#include <algorithm>
#include <functional>
#include <iostream>
#include <vector>
#include <queue>
#include <iomanip>
#include <fstream>
#include <cfloat>
#include <climits>
using namespace std;

int _count;

class Conditional;
class Mutex
{
	friend class Conditional;
	pthread_mutex_t mutex;
public:
	Mutex()
	{
		if (pthread_mutex_init(&mutex, NULL) != 0)
			cerr<<"error of call pthread_mutex_init"<<endl;
	}
	void lock()
	{
		if (pthread_mutex_lock(&mutex) != 0)
			cerr<<"error of call pthread_mutex_lock "<<endl;
		cerr.flush();
	}
	void unlock()
	{
		if (pthread_mutex_unlock(&mutex) != 0)
			cerr<<"error of call pthread_mutex_unlock "<<endl;
	}
	~Mutex()
	{
		if (pthread_mutex_destroy(&mutex) != 0)
			cerr<<"error of call pthread_mutex_destroy "<<endl;
	}
};

class Conditional
{
	pthread_cond_t cond;
public:
	Conditional()
	{
		if (pthread_cond_init(&cond, NULL) != 0)
			cerr<<"error of call pthread_cond_init"<<endl;
	}
	void wait(Mutex &mutex)
	{
		if (pthread_cond_wait(&cond, &mutex.mutex) != 0)
			cerr<<"error of call pthread_cond_wait"<<endl;
	}
	void signal()
	{
		if (pthread_cond_signal(&cond) != 0)
			cerr<<"error of call pthread_cond_signal"<<endl;
	}
	void broadcast()
	{
		if (pthread_cond_broadcast(&cond) != 0)
			cerr<<"error of call pthread_cond_broadcast"<<endl;
	}
	~Conditional()
	{
		if (pthread_cond_destroy(&cond) != 0)
			cerr<<"error of call pthread_cond_destroy"<<endl;
	}
};

class Runnable
{
protected:
	bool _abort;
public:
	Runnable() : _abort(false)
	{
	}
	virtual void abort()
	{
		_abort = true;
	}
	virtual void run() = 0;
	virtual ~Runnable()
	{
	}
};

struct Waiter
{
private:
	Mutex mutex;
	Conditional cond;
	bool waited;

public:
	Waiter()
	{
		waited = true;
	}
	void wait()
	{
		mutex.lock();
		while(waited)
			cond.wait(mutex);
		mutex.unlock();
	}
	void signal()
	{
		mutex.lock();
		waited = false;
		cond.signal();
		mutex.unlock();
	}
};

class Thread
{
	pthread_t tid;
	Runnable *runnable;
	Waiter* waiter;
public:
	Thread(Runnable *runnable, Waiter* waiter = NULL)
	{
		this->runnable = runnable;
		this->waiter = waiter;
	}
	void start()
	{
		int res = pthread_create(&tid, NULL, _start, this);
		while (res == EAGAIN)
		{
			usleep((useconds_t)1e5);
			res = pthread_create(&tid, NULL, _start, this);
		}

		if (res != 0)
			cerr<<"error of call pthread_create "<<res<<endl;
	}
private:
	static void* _start(void* param)
	{
		Thread* thread = reinterpret_cast<Thread*>(param);
		pthread_detach(thread->tid);
		thread->runnable->run();
		Waiter* waiter = thread->waiter;
		delete thread->runnable;
		delete thread;
		if (waiter != NULL)
			waiter->signal();
		return NULL;
	}
};

struct TSP
{
	static unsigned state;

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

	short randx()
	{
		state = state * 1664525 + 1013904223;
		return (state >> 16) & 0x7FFF;
	}

	void generate(string outFile, string pointsFileName = string())
	{
		const int W = 500;
		const int H = 500;
		//srand(1);
		vector<pair<int, int> > points(N);
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
			cout<<visited1[0]<<" "<<visited1[1]<<" "<<visited1[2]<<" "<<visited1[3]<<" "<<visited1[4]<<" - "<<best<<endl;
			cout<<"*"<<c[visited1[0]][visited1[1]]<<endl;
			cout<<"*"<<c[visited1[1]][visited1[2]]<<endl;
			cout<<"*"<<c[visited1[2]][visited1[3]]<<endl;
			cout<<"*"<<c[visited1[3]][visited1[4]]<<endl;
			cout<<"*"<<c[visited1[4]][visited1[0]]<<endl;
			cout<<"*"<<lastVisited<<endl;*/
			_count++;
			if (dist + c[lastVisited][0] < best)
			{
				best = dist + c[lastVisited][0];
				bestVisited = visited1;
			}
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

unsigned TSP::state = 1;

template <class T>
struct MessageKeeper
{
	int size;
	T* body;
	
	MessageKeeper(int N) : 
		size (sizeof(T) + N * sizeof(int)),
		body((T *)new char[size])
	{
	}

	~MessageKeeper()
	{
		delete[] (char*)body;
	}
};

struct BestMessage
{
	double updatedBest;
	int visited[0];
};

struct TaskMessage
{
	double currentBest;
	int nVisited;
	int visited[0];
};

struct Task
{
	vector<int> visited;
	double currentBest;

	Task()
	{
		currentBest = 0;
	}

	Task& operator=(const TaskMessage& msg)
	{
		currentBest = msg.currentBest;
		visited.assign(msg.visited, msg.visited + msg.nVisited);
	}
};

struct MessageManager
{
	enum MessageType
	{
		NEW_BEST,
		TASK,
		NO_TASKS
	};
	const int N;
	MessageKeeper<TaskMessage> taskMsg;
	MessageKeeper<BestMessage> bestMsg;

	MessageManager(int N) : N(N), taskMsg(N), bestMsg(N)
	{
	}

	MessageType receive(int src, MPI_Status* pst = NULL)
	{
		MPI_Status st;
		if (pst == NULL)
			pst = &st;
		MPI_Probe(src, MPI_ANY_TAG, MPI_COMM_WORLD, pst);
		switch (pst->MPI_TAG)
		{
		case TASK:
			if (pst->count != 0)
				MPI_Recv(taskMsg.body, taskMsg.size, MPI_BYTE, src, TASK, MPI_COMM_WORLD, pst);
			else
				MPI_Recv(NULL, 0, MPI_BYTE, src, TASK, MPI_COMM_WORLD, pst);
			return TASK;
		case NO_TASKS:
			MPI_Recv(NULL, 0, MPI_BYTE, src, NO_TASKS, MPI_COMM_WORLD,pst);
			return NO_TASKS;
		case NEW_BEST:
			MPI_Recv(bestMsg.body, bestMsg.size, MPI_BYTE, src, NEW_BEST, MPI_COMM_WORLD, pst);
			return NEW_BEST;
		}
	}
};

struct MasterSender : public Runnable
{
	double best;
	queue<Task> tasks;
	queue<int> requests;
	int nWorkers;
	MessageKeeper<TaskMessage> taskMsg;

	Mutex mutex;
	Conditional cond;

	MasterSender(int N, int nWorkers) : nWorkers(nWorkers), taskMsg(N)
	{
	}

	bool addRequest(int rank, double currentBest)
	{
		mutex.lock();
		best = currentBest;
		requests.push(rank);
		cond.signal();
		bool ret = !(tasks.empty() && nWorkers == requests.size());
		mutex.unlock();
		return ret;
	}

	void run()
	{
		int nTasks = tasks.size();
		for(; nWorkers > 0;)
		{
			mutex.lock();
			while (requests.empty())
				cond.wait(mutex);
			int iWorker = requests.front();
			requests.pop();
			
			if (tasks.empty())
			{
				nWorkers--;
				mutex.unlock();
				MPI_Send(NULL, 0, MPI_BYTE, iWorker, MessageManager::NO_TASKS, MPI_COMM_WORLD);
			}
			else
			{
				Task &task = tasks.front();
				taskMsg.body->nVisited = task.visited.size();
				taskMsg.body->currentBest = best;
				copy(task.visited.begin(), task.visited.end(), taskMsg.body->visited);
				tasks.pop();
				mutex.unlock();
				{
					stringstream ss;
					ss<<"master: "<<tasks.size()<<"/"<<nTasks<<" tasks\n";
					cout<<ss.str();
				}
				MPI_Send(taskMsg.body, taskMsg.size, MPI_BYTE, iWorker, MessageManager::TASK, MPI_COMM_WORLD);
			}
		}
	}
};

struct Master
{
	static const int DEPTH = 3;
	
	int size;
	TSP tsp;

	Master(string dataFile, int size)
	{
		this->size = size;
		tsp.load(dataFile);
		MPI_Bcast(&tsp.N, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);

		for (int i = 0; i < tsp.N; i++)
			MPI_Bcast(&tsp.c[i][0], tsp.N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}

	void createTasks(const Task& root, int depth, queue<Task>& tasks)
	{
		if (--depth == 0)
		{
			tasks.push(root);
			return;
		}
		for (int i = 0; i < tsp.N; i++)
		{
			if (find(root.visited.begin(), root.visited.end(), i) == root.visited.end())
			{
				Task task = root;
				task.visited.push_back(i);
				createTasks(task, depth, tasks);
			}
		}
	}

	void start()
	{
		if (DEPTH - 1 >= tsp.N)
		{
			tsp.solveBruteForce();
			return;
		}
		MasterSender* sender = new MasterSender(tsp.N, size - 1);
		queue<Task>& tasks = sender->tasks;
		Task root;
		root.visited.push_back(0);
		
		createTasks(root, DEPTH, tasks);
		tsp.greedy();

		Waiter waiter;
		Thread* thread = new Thread(sender, &waiter);
		thread->start();

		MessageManager msgManager(tsp.N);
		for(bool stop = false; !stop;)
		{
			MPI_Status st;
			switch(msgManager.receive(MPI_ANY_SOURCE, &st))
			{
			case MessageManager::NEW_BEST:
				{
					const BestMessage& msg = *msgManager.bestMsg.body;
					if (msg.updatedBest < tsp.best)
					{
						tsp.best = msg.updatedBest;
						tsp.bestVisited.assign(msg.visited, msg.visited + tsp.N);
					}
				}
				break;
			case MessageManager::TASK:
				stop = !sender->addRequest(st.MPI_SOURCE, tsp.best);
				break;
			}
		}
		waiter.wait();
		cout<<"ANSWER "<<tsp.best<<"\n";
		for (int i = 0; i < tsp.bestVisited.size(); i++)
			cout<<" - "<<tsp.bestVisited[i];
		cout<<endl;
	}
};

struct WorkerSender : public Runnable
{
	const int N;
	Mutex mutex;
	Conditional cond;
	bool requestTask;
	bool updateBest;
	MessageKeeper<BestMessage> bestMsg;

	WorkerSender(int N) : N(N), bestMsg(N)
	{
		requestTask = true;
		updateBest = false;
	}

	void abort()
	{
		mutex.lock();
		_abort = true;
		cond.signal();
		mutex.unlock();
	}

	void sendBest(double newBest, const vector<int>& visited)
	{
		mutex.lock();
		bestMsg.body->updatedBest = newBest;
		copy(visited.begin(), visited.end(), bestMsg.body->visited);
		updateBest = true;
		cond.signal();
		mutex.unlock();
	}

	void sendRequest()
	{
		mutex.lock();
		requestTask = true;
		cond.signal();
		mutex.unlock();
	}

	void run()
	{
		MessageKeeper<BestMessage> bestMsgLocal(N);
		for(;;)
		{
			mutex.lock();
			while (!requestTask && !updateBest && !_abort)
				cond.wait(mutex);

			if (_abort)
			{
				mutex.unlock();
				return;
			}
			bool requestTaskLocal = requestTask;
			bool updateBestLocal = updateBest;
			requestTask = false;
			updateBest = false;
			if (updateBestLocal)
			{
				bestMsgLocal.body->updatedBest = bestMsg.body->updatedBest;
				copy(bestMsg.body->visited, bestMsg.body->visited + N, bestMsgLocal.body->visited);
			}
			mutex.unlock();

			if (requestTaskLocal)
				MPI_Send(NULL, 0, MPI_BYTE, 0, MessageManager::TASK, MPI_COMM_WORLD);
			if (updateBestLocal)
				MPI_Send(bestMsgLocal.body, bestMsgLocal.size, MPI_BYTE, 0, MessageManager::NEW_BEST, MPI_COMM_WORLD);
		}
	}
};

struct Worker
{
	TSP tsp;
	Task task;
	WorkerSender* sender;
	double t;
	double tIdle;

	Worker()
	{
		int N;
		MPI_Bcast(&N, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
		tsp.init(N);

		for (int i = 0; i < N; i++)
			MPI_Bcast(&tsp.c[i][0], N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		sender = new WorkerSender(N);
	}

	void start()
	{
		tIdle = 0;
		t = MPI_Wtime();
		Waiter waiter;
		Thread* thread = new Thread(sender, &waiter);
		thread->start();
		
		MessageManager msgManager(tsp.N);
		for(;;)
		{
			tIdle -= MPI_Wtime();
			MessageManager::MessageType msgType = msgManager.receive(0);
			tIdle += MPI_Wtime();
			switch (msgType)
			{
			case MessageManager::TASK:
				task = *msgManager.taskMsg.body;
				break;
			case MessageManager::NO_TASKS:
				goto end;
			}
			tsp.best = task.currentBest;
			tsp.visited1 = task.visited;
			tsp.solve(&TSP::solveRecAux);
			if (tsp.best < task.currentBest)
				sender->sendBest(tsp.best, tsp.bestVisited);

			sender->sendRequest();
		}
end:
		sender->abort();
		waiter.wait();
		t = MPI_Wtime() - t;
	}
};

int main(int argc, char **argv)
{
	if (argc <= 3)
		return 0;
	int provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE,  &provided);

	int rank;
	int size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	if (rank == 0)
	{
		if (argc == 4)
			TSP::state = atoi(argv[3]);
		TSP tsp(atoi(argv[2]));
		tsp.generate(argv[1], "pf.txt");
	}
	if (size < 2)
	{
		MPI_Finalize();
		return 0;
	}
	

	if (rank == 0)
	{
		Master(argv[1], size).start();
		MPI_Barrier(MPI_COMM_WORLD);
	}
	else
	{
		Worker worker;
		worker.start();
		MPI_Barrier(MPI_COMM_WORLD);
		stringstream ss;
		ss<<"worker"<<rank<<" time: "<<worker.t<<"s\n\tidle time: "<<worker.tIdle<<"s\t"<<
			worker.tIdle * 100 / worker.t<<"%"<<endl;
		cout<<ss.str();
	}
    MPI_Finalize();
    return 0;
}