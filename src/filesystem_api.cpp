
struct File {
  unsigned char blob[8];
};

bool   (*file_open)(File*, const char*);
bool   (*file_close)(File*);
bool   (*file_read)(File*, void*, size_t, size_t*);
bool   (*file_write)(File*, const char*, size_t);
size_t (*file_get_size)(File*);
void check_filesystem_api();
void init_filesystem_api();
