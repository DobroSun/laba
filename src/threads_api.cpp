
struct Thread {
  unsigned char blob[32];
};

bool (*start_thread)(Thread*, int (*proc)(void*), void* data);
bool (*is_thread_running)(Thread*);
bool (*suspend_thread)(Thread*);
bool (*resume_thread)(Thread*);
bool (*kill_thread)(Thread*);
void check_threads_api();
void init_threads_api();
