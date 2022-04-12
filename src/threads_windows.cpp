
struct Thread_Windows {
  HANDLE handle;
  int (*thread_proc)(void*);
  void* thread_data;
};

bool start_thread_windows(Thread* thread, int (*proc)(void*), void* data) {
  Thread_Windows* t = (Thread_Windows*) thread;

  t->thread_proc = proc;
  t->thread_data = data;
  t->handle = CreateThread(NULL,
                           0,
                           (LPTHREAD_START_ROUTINE) t->thread_proc,
                           t->thread_data,
                           0,
                           NULL);
  return t->handle;
}

bool is_thread_running_windows(Thread* thread) {
  Thread_Windows* t = (Thread_Windows*) thread;

  bool thread_is_initialized = t->handle != NULL;
  bool thread_is_running     = WaitForSingleObject(t->handle, 0) != WAIT_OBJECT_0;

  return thread_is_initialized && thread_is_running;
}

bool suspend_thread_windows(Thread* thread) {
  Thread_Windows* t = (Thread_Windows*) thread;
  return SuspendThread(t->handle) != -1;
}

bool resume_thread_windows(Thread* thread) {
  Thread_Windows* t = (Thread_Windows*) thread;
  return ResumeThread(t->handle) != -1;
}

bool kill_thread_windows(Thread* thread) {
  Thread_Windows* t = (Thread_Windows*) thread;
  return TerminateThread(t->handle, 0);
}

void check_threads_api() {
  assert(is_thread_running);
  assert(start_thread);
  assert(suspend_thread);
  assert(resume_thread);
  assert(kill_thread);
}

void init_threads_api() {
  is_thread_running = is_thread_running_windows;
  start_thread   = start_thread_windows;
  suspend_thread = suspend_thread_windows;
  resume_thread  = resume_thread_windows;
  kill_thread    = kill_thread_windows;
}
