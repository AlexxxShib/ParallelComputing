#include <unistd.h>
#include <mpi.h>
#include <map>
#include <sstream>
#include <algorithm>
#include <functional>
#include <iostream>
#include <vector>
#include <queue>
#include <iomanip>
#include <fstream>
#include <cfloat>
#include <climits>
#include <ctime>
using namespace std;

const int MAX_DURATION = 1e5;

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
		int res = pthread_mutex_destroy(&mutex);
		if (res != 0)
			cerr<<"error of call pthread_mutex_destroy "<<res<<endl;
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
	Runnable* runnable;
	Waiter* waiter;
public:
	Thread(Runnable* runnable, Waiter* waiter = NULL)
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

struct Task
{
private:
	static int ID;
public:
	static Task EMPTY_TASK;
	int duration;
	int id;

	Task() : duration(-1)
	{
	}
	
	Task(int rank, int duration) : duration(duration), id(rank * 1000 + ID++)
	{
	}

	void execute(int rank)
	{
		usleep(duration);
		{
			stringstream ss;
			ss<<"Task "<<id<<" was executed by "<<rank<<" process"<<endl;
			cout<<ss.str();
		}
	}
	bool isEmpty()
	{
		return duration == -1;
	}
};
int Task::ID = 0;
Task Task::EMPTY_TASK;

struct TaskPool
{
	Mutex mutex;
	Conditional cond;
	bool closed;

	TaskPool()
	{
		closed = false;
	}
	void enqueue(const Task& task)
	{
		mutex.lock();
		tasks.push(task);
		cond.signal();
		mutex.unlock();
	}
	Task waitTask()
	{
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		mutex.lock();
		while(tasks.empty() && !closed )
			cond.wait(mutex);

		if (tasks.empty())
		{
			mutex.unlock();
			return Task::EMPTY_TASK;
		}
		Task task = tasks.front();
		tasks.pop();
		mutex.unlock();
		return task;
	}

	Task dequeue()
	{
		mutex.lock();
		if (tasks.empty())
		{
			mutex.unlock();
			return Task::EMPTY_TASK;
		}
		Task task = tasks.front();
		tasks.pop();
		mutex.unlock();
		return task;
	}

	bool isEnds()
	{
		return tasks.size() < 2;
	}

	void closePool()
	{
		mutex.lock();
		closed = true;
		cond.signal();
		mutex.unlock();
	}

private:
	queue<Task> tasks;
};

struct MessageManager
{
	enum MessageType
	{
		REQUEST_TASK,
		RESPONSE_TASK,
		NO_TASKS,
		FINISH
	};

	Task task;
	MessageType waitResponseMessage(int src, MPI_Status* pst = NULL)
	{
		MPI_Status st;
		if (pst == NULL)
			pst = &st;
		MPI_Probe(src, RESPONSE_TASK, MPI_COMM_WORLD, pst);
		if (pst->count != 0)
		{
			MPI_Recv(&task, sizeof(Task), MPI_BYTE, src, RESPONSE_TASK, MPI_COMM_WORLD, pst);
			return RESPONSE_TASK;
		}
		MPI_Recv(NULL, 0, MPI_BYTE, src, RESPONSE_TASK, MPI_COMM_WORLD, pst);
		return NO_TASKS;
	}

	MessageType waitRequestMessage(int src, MPI_Status* pst = NULL)
	{
		MPI_Status st;
		if (pst == NULL)
			pst = &st;
		MPI_Probe(src, REQUEST_TASK, MPI_COMM_WORLD, pst);
		if (pst->count != 0)
		{
			int dummy;
			MPI_Recv(&dummy, sizeof(dummy), MPI_BYTE, src, REQUEST_TASK, MPI_COMM_WORLD, pst);
			return REQUEST_TASK;
		}
		MPI_Recv(NULL, 0, MPI_BYTE, src, REQUEST_TASK, MPI_COMM_WORLD, pst);
		return FINISH;
	}

	void sendMessage(int dest, MessageType msgType, const Task& task = Task())
	{
		switch (msgType)
		{
		case REQUEST_TASK:
			MPI_Send(&dest, sizeof(dest), MPI_BYTE, dest, REQUEST_TASK, MPI_COMM_WORLD);
			break;
		case RESPONSE_TASK:
			MPI_Send((void*)&task, sizeof(task), MPI_BYTE, dest, RESPONSE_TASK, MPI_COMM_WORLD);
			break;
		case NO_TASKS:
			MPI_Send(NULL, 0, MPI_BYTE, dest, RESPONSE_TASK, MPI_COMM_WORLD);
			break;
		case FINISH:
			MPI_Send(NULL, 0, MPI_BYTE, dest, REQUEST_TASK, MPI_COMM_WORLD);
			break;

		}
	}
};

struct Sender : public Runnable
{
	int size;
	int rank;
	TaskPool& taskPool;
	bool requested;
	Mutex mutex;
	Conditional cond;

	Sender(TaskPool& taskPool, int rank, int size) : taskPool(taskPool), rank(rank), size(size)
	{
		requested = false;
	}

	void abort()
	{
		mutex.lock();
		_abort = true;
		cond.signal();
		mutex.unlock();
	}

	void sendRequest()
	{
		//if (taskPool.closed)
		//	return;
		mutex.lock();
		requested = true;
		cond.signal();
		mutex.unlock();
	}

	void sendFinish()
	{
		MessageManager msgManager;
		for (int i = 0; i < size; i++)
			if (i != rank)
				msgManager.sendMessage(i, MessageManager::FINISH);
	}

	void run()
	{
		for(;;)
		{
			MessageManager msgManager;
			mutex.lock();
			while (!requested && !_abort)
				cond.wait(mutex);
			if (_abort)
			{
				mutex.unlock();
				sendFinish();
				return;
			}
			requested = false;
			mutex.unlock();
			bool taskReceived = false;
			for (int i = 0; i < size && !taskReceived; i++)
				if (i != rank)
				{
					msgManager.sendMessage(i, MessageManager::REQUEST_TASK);
					switch (msgManager.waitResponseMessage(i))
					{
					case MessageManager::RESPONSE_TASK:
						taskPool.enqueue(msgManager.task);
						taskReceived = true;
						break;
					case MessageManager::NO_TASKS:
						break;
					}
					
				}
			if (!taskReceived)
			{
				taskPool.closePool();
				sendFinish();
				mutex.lock();
				while (!_abort)
					cond.wait(mutex);
				mutex.unlock();
				return;
			}
		}
	}
};

struct Receiver : public Runnable
{
	int rank;
	int size;
	TaskPool& taskPool;

	Receiver(TaskPool& taskPool, int rank, int size) : rank(rank), size(size), taskPool(taskPool)
	{
	}

	void run()
	{
		MessageManager msgManager;
		for(int nLived = size - 1; nLived > 0;)
		{
			MPI_Status st;
			switch (msgManager.waitRequestMessage(MPI_ANY_SOURCE, &st))
			{
			case MessageManager::REQUEST_TASK:
				{
					Task task = taskPool.dequeue();
					if (task.isEmpty())
						msgManager.sendMessage(st.MPI_SOURCE, MessageManager::NO_TASKS);
					else
						msgManager.sendMessage(st.MPI_SOURCE, MessageManager::RESPONSE_TASK, task);
					break;
				}
			case MessageManager::FINISH:
				nLived--;
				break;
			}
		}
		taskPool.closePool();
	}
};

struct Worker
{
	const int N;
	TaskPool taskPool;

	int nOwnTasks;
	int nTasks;

	int rank;
	int size;

	Worker() : N(rand() % 1000 + 10)
	{
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &size);

		for (int i = 0; i < N; i++)
			taskPool.enqueue(Task(rank, rand() % MAX_DURATION));

		nOwnTasks = N;
		nTasks = 0;
	}

	void start()
	{
		Waiter recvWaiter;
		Waiter sendWaiter;
		Sender* sender = new Sender(taskPool, rank, size);
		Receiver* receiver = new Receiver(taskPool, rank, size);

		Thread* senderThread = new Thread(sender, &sendWaiter);
		Thread* receiverThread = new Thread(receiver, &recvWaiter);
		senderThread->start(); 
		receiverThread->start();

		for(;;)
		{
			Task task = taskPool.waitTask();
			if (task.isEmpty())
				break;

			nTasks++;
			task.execute(rank);
			if (taskPool.isEnds())
				sender->sendRequest();
		}
		sender->abort();
		sendWaiter.wait();
		recvWaiter.wait();
	}
};

int main(int argc, char **argv)
{
	int provided;
	int rank;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE,  &provided);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	srand((rank + 1) * time(NULL));

	Worker worker;
	worker.start();
	MPI_Barrier(MPI_COMM_WORLD);

	{
		stringstream ss;
		ss<<"worker"<<rank<<"\n\texecute "<<worker.nTasks<<" tasks\n\town tasks: "<<worker.nOwnTasks<<endl;
		cout<<ss.str();
	}

    MPI_Finalize();
    return 0;
}