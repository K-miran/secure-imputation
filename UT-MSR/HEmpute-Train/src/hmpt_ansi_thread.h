#ifndef _ANSI_THREAD_
#define _ANSI_THREAD_

/*
t_socket_thread class: Represents a thread for any type of connection oriented socket, this can be on either side of connection, i.e., it can be
used in client or server. 
*/

enum{THREAD_INITING, THREAD_RUNNING, THREAD_TERMINAL}; 

#ifdef __unix__
// For linux multithreading
#include <pthread.h>
#elif defined _WIN32
#include <windows.h>
#else 
#error "The environment could not be determined, must specify manually."
#endif

class t_ansi_thread
{
public:
	static int thread_count;
	t_ansi_thread(void* (*_thread_function)(void*), void* _thread_params);
	~t_ansi_thread();

	void* (*thread_function)(void*);

	void* thread_params;

	// State of the thread, necessary to keep track of detached pthread's. 
	volatile int thread_state; 

	int id;

	bool run_thread(); // Start the socket thread. (Start protocol_implementing_function as a thread)
	void wait_thread(); // Wait for the socket thread to finish.
	void exit_thread(); // exit this thread.

#ifdef _WIN32
	// Real (Win32) handle of accepted client.
	// Have to find a corresponding element to refer to the created thread in linux platform.
	// Use WaitForSingleObject to wait for the threads to finish.
	HANDLE thread_handle;
#endif

#ifdef __unix__
	// pthread handle of accepted client.
	// Have to find a corresponding element to refer to the created thread in linux platform.
	// Use pthread_join (http://www.yolinux.com/TUTORIALS/LinuxTutorialPosixThreads.html#BASICS)
	// to wait for the threads to finish.
	pthread_t thread_handle;
#endif
};

#endif // _ANSI_THREAD_
