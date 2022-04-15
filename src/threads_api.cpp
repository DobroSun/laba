
struct Thread {
  unsigned char blob[32];
};

struct Mutex {
  unsigned char blob[16];
};

bool (*start_thread)(Thread*, int (*proc)(void*), void* data);
bool (*is_thread_running)(Thread*);
bool (*suspend_thread)(Thread*);
bool (*resume_thread)(Thread*);
bool (*kill_thread)(Thread*);

void check_threads_api();
void init_threads_api();

Mutex (*create_mutex)();
bool  (*lock_mutex)(Mutex*);
bool  (*release_mutex)(Mutex*);

struct Scoped_Lock {
  Mutex* mutex;
  Scoped_Lock(Mutex* m): mutex(m) { lock_mutex(mutex); }
  ~Scoped_Lock()                  { release_mutex(mutex); }
};

