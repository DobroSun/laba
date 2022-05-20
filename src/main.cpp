// Dear ImGui: standalone example application for DirectX 11
// If you are new to Dear ImGui, read documentation from the docs/ folder + read the top of imgui.cpp.
// Read online: https://github.com/ocornut/imgui/tree/master/docs

#include "../imgui/imgui.cpp"
#include "../imgui/imgui_demo.cpp"
#include "../imgui/imgui_draw.cpp"
#include "../imgui/imgui_tables.cpp"
#include "../imgui/imgui_widgets.cpp"
#include "../imgui/backends/imgui_impl_win32.cpp"
#include "../imgui/backends/imgui_impl_dx11.cpp"
#include "../implot/implot.cpp"
#include "../implot/implot_demo.cpp"
#include "../implot/implot_items.cpp"

#include "../std/pch.h"

#include "threads_api.cpp"
#include "threads_windows.cpp"

#pragma comment(lib, "d3d11.lib")
#pragma comment(lib, "d3dcompiler.lib")

#define clamp(w, mi, ma) ( min(max((w), (mi)), (ma)) )

typedef uint32_t uint;
typedef uint32_t uint32;
typedef uint64_t uint64;
typedef  int64_t  int64;

// Forward declarations.
bool  ImGui_ImplWin32_Init(void* hwnd);
void  ImGui_ImplWin32_Shutdown();
void  ImGui_ImplWin32_NewFrame();
void  ImGui_ImplWin32_EnableDpiAwareness();
float ImGui_ImplWin32_GetDpiScaleForHwnd(void* hwnd);       // HWND hwnd
float ImGui_ImplWin32_GetDpiScaleForMonitor(void* monitor); // HMONITOR monitor
void  ImGui_ImplWin32_EnableAlphaCompositing(void* hwnd);   // HWND hwnd
bool  ImGui_ImplDX11_Init(ID3D11Device* device, ID3D11DeviceContext* device_context);
void  ImGui_ImplDX11_Shutdown();
void  ImGui_ImplDX11_NewFrame();
void  ImGui_ImplDX11_RenderDrawData(ImDrawData* draw_data);
void  ImGui_ImplDX11_InvalidateDeviceObjects();
bool  ImGui_ImplDX11_CreateDeviceObjects();

bool CreateDeviceD3D(HWND hWnd);
void CleanupDeviceD3D();
void CreateRenderTarget();
void CleanupRenderTarget();
LRESULT WINAPI WndProc(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam);
// end.



static ID3D11Device*            g_pd3dDevice = NULL;
static ID3D11DeviceContext*     g_pd3dDeviceContext = NULL;
static IDXGISwapChain*          g_pSwapChain = NULL;
static ID3D11RenderTargetView*  g_mainRenderTargetView = NULL;

#define min(x, y) ( ((x) < (y)) ? (x) : (y) )
#define max(x, y) ( ((x) > (y)) ? (x) : (y) )
#define square(x) ( (x)*(x) )
#define abs(x)    ( (x>0) ? (x) : -(x) )
#define array_size(x) ( sizeof(x)/sizeof(*(x)) )


// 
// p (x, t=0) := ...;
// u (x, t=0) := 0;
// ro(x, t=0) := ro0;
// 
double p_t0(double x, double p0, double x0, double r0) {
  return p0 * exp(-square(x - x0) / square(r0));
}

double u_t0(double x) {
  return 0;
}

double ro_t0(double x, double ro0) {
  return ro0;
}

// 
// u(x=0, t) := 0;
// u(x=1, t) := 0;
//
double u_x0(double t) {
  return 0;
}

struct Parameters_Laba1 {
  float xmin;
  float xmax;

  float p0;
  float x0;
  float r0;
  float ro0;
  float gamma;

  float dt;
  float dx;
};

struct Parameters_Laba2 {
  float xmin;
  float xmax;

  float p0;
  float x0;
  float r0;
  float ro0;
  float gamma;
  float alpha;

  float a;

  float dt;
  float dx;
};

struct Parameters_Laba3 {
  float xmin;
  float xmax;
  float ymin;
  float ymax;
  float center_x;
  float center_y;
  float a;
  float b;
  float Re;
  float U0;
  float dx;
  float dy;
  float dt;
};

struct Thread_Data1 {
  Parameters_Laba1* parameters;
  double* global_t;
  double* dt;
};

struct Thread_Data2 {
  Parameters_Laba2* parameters;
  double* global_t;
  double* dt;
};

struct Thread_Data3 {
  Parameters_Laba3* parameters;
  double* global_t;
  double* dt;
};



struct Vertex {
  double r;
  double u;
  union {
    double p;
    double E;
  };
};

struct Voxel {
  double X;
  double u;
  double r;
  double p;
  double E;
};

struct Velocity_And_Euler {
  double X;
  double u;
};


struct Voxel2d {
  double eps;
  double psi;
  double u;
  double v;
};

struct Result {
  Mutex data_mutex = {};

  array<double> time; // time domain;
  array<double> grid; // for now we only have 1 dimension.
  array<Vertex> data; // this is going to evolve in time. size := NUMBER_OF_LAYERS * grid.size()
};

struct Result_Lagrange {
  Mutex data_mutex = {};

  array<double> time;
  array<double>   points_grid; 
  array<double> internal_grid;

  array<Velocity_And_Euler> points_data;
  array<Voxel>            internal_data;
};

struct Vec2 {
  double x;
  double y;
};

struct Result2d {
  Mutex data_mutex = {};

  array<Voxel2d> data;
  //array<double>  time; // @CleanUp: actually we never use time, because our dt is always const. So we can just remove time from Result, Result_Lagrange, Result2d.

  size_t size_x, size_y;
};

Vertex get_data(Result* r, size_t i, size_t j) {
  return r->data[i * r->grid.size + j];
}

double get_grid(Result_Lagrange* l, size_t j) {
  return l->points_grid[j];
}

Velocity_And_Euler get_point(Result_Lagrange* l, size_t i, size_t j) {
  return l->points_data[i * l->points_grid.size + j];
}

Voxel get_internal(Result_Lagrange* l, float i, float j) {
  size_t ii = i;
  size_t jj = j - 1/2.0f;
  return l->internal_data[ii * l->internal_grid.size + jj];
}

Voxel2d get_data(Result2d* r, size_t i, size_t j, size_t t) {
  return r->data[t * (r->size_x * r->size_y) + i * r->size_y + j];
}

Voxel2d* get_pointer(Result2d* r, array<Voxel2d>* v, size_t i, size_t j) {
  return &v->data[i * r->size_y + j];
}

Voxel2d get_data(Result2d* r, array<Voxel2d>* v, size_t i, size_t j) {
  return *get_pointer(r, v, i, j);
}

static Result result;
static Result_Lagrange lagrange;
static Result2d result2d;


static int euler_upwind_method(void* param) {
  Thread_Data1* data = (Thread_Data1*) param;
  Parameters_Laba1 params = *data->parameters; // copy

  double xmin = params.xmin; 
  double xmax = params.xmax;
  double p0   = params.p0;
  double x0   = params.x0;
  double r0   = params.r0;
  double ro0  = params.ro0;
  double gamma = params.gamma;
  double dt   = params.dt;
  double dx   = params.dx;

  size_t NUMBER_OF_POINTS_IN_X = (xmax - xmin) / dx;

  double t = 0.0; // real time.

  InterlockedExchange64((int64*) data->dt,       *(int64*) &dt);
  InterlockedExchange64((int64*) data->global_t, *(int64*) &t);

  {
    // 
    // create a grid.
    //
    array<double> grid = temp_array(double, NUMBER_OF_POINTS_IN_X);
    for (size_t i = 0; i < grid.size; i++) {
      grid[i] = xmin + i*dx;
    }

    Scoped_Lock lock(&result.data_mutex);
    result.grid = array_copy(&grid);
  }

  array<Vertex> solution = temp_array(Vertex, NUMBER_OF_POINTS_IN_X);

  {
    // 
    // apply initial conditions on t=0
    //
    for (size_t j = 0; j < solution.size; j++)  {
      double x = result.grid[j];

      Vertex* v = &solution[j];
      v->r = ro_t0(x, ro0);
      v->u =  u_t0(x);
      v->p =  p_t0(x, p0, x0, r0);
    }

    Scoped_Lock lock(&result.data_mutex);
    array_add(&result.data, &solution);
    array_add(&result.time, t);
  }

  size_t i = 0;
  while (1) {
    i++;
    t += dt;

    for (size_t j = 1; j < solution.size-1; j++) {
      // 
      // internal points.
      //
      Vertex curr = get_data(&result, i-1, j);
      Vertex prev = get_data(&result, i-1, j-1);
      Vertex next = get_data(&result, i-1, j+1);

      double x = result.grid[j];
      double a = (curr.u >= 0) ? 0 : 1;

      Vertex* v = &solution[j];
      v->r = curr.r - (1 - a) * curr.u * dt / dx * (curr.r - prev.r) - a * curr.u * dt / dx * (next.r - curr.r) - dt / (2*dx) * curr.r * (next.u - prev.u);
      v->u = curr.u - (1 - a) * curr.u * dt / dx * (curr.u - prev.u) - a * curr.u * dt / dx * (next.u - curr.u) - dt / (2*dx) / curr.r * (next.p - prev.p);;
      v->p = curr.p - (1 - a) * curr.u * dt / dx * (curr.p - prev.p) - a * curr.u * dt / dx * (next.p - curr.p) - dt / (2*dx) * gamma * curr.p * (next.u - prev.u);
    }

    size_t j = 0;
    solution[j]   = solution[j+1];
    solution[j].u = u_x0(0);

    j = solution.size-1;
    solution[j]   = solution[j-1];
    solution[j].u = u_x0(0);

    {
      Scoped_Lock lock(&result.data_mutex);
      array_add(&result.data, &solution);
      array_add(&result.time, t);
    }
    InterlockedExchange64((int64*) data->global_t, *(int64*) &t);
  }
  return 0;
}

static int non_conservative_lax_method(void* param) {
  Thread_Data1* data = (Thread_Data1*) param;
  Parameters_Laba1 params = *data->parameters; // copy

  double xmin = params.xmin; 
  double xmax = params.xmax;
  double p0   = params.p0;
  double x0   = params.x0;
  double r0   = params.r0;
  double ro0  = params.ro0;
  double gamma = params.gamma;
  double dt   = params.dt;
  double dx   = params.dx;

  size_t NUMBER_OF_POINTS_IN_X = (xmax - xmin) / dx;

  double t = 0.0;

  InterlockedExchange64((int64*) data->dt,       *(int64*) &dt);
  InterlockedExchange64((int64*) data->global_t, *(int64*) &t);

  {
    // 
    // create a grid.
    //
    array<double> grid = temp_array(double, NUMBER_OF_POINTS_IN_X);
    for (size_t i = 0; i < grid.size; i++) {
      grid[i] = xmin + i*dx;
    }

    Scoped_Lock lock(&result.data_mutex);
    result.grid = array_copy(&grid);
  }

  array<Vertex> solution = temp_array(Vertex, NUMBER_OF_POINTS_IN_X);

  {
    // 
    // apply initial conditions on t=0
    //
    for (size_t j = 0; j < solution.size; j++)  {
      double x = result.grid[j];

      Vertex* v = &solution[j];
      v->r = ro_t0(x, ro0);
      v->u =  u_t0(x);
      v->p =  p_t0(x, p0, x0, r0);
    }

    {
      Scoped_Lock lock(&result.data_mutex);
      array_add(&result.data, &solution);
      array_add(&result.time, t);
    }
  }

  size_t i = 0;
  while (1) {
    i++;
    t += dt;

    for (size_t j = 1; j < solution.size-1; j++) {
      // 
      // internal points.
      //
      Vertex curr = get_data(&result, i-1, j);
      Vertex prev = get_data(&result, i-1, j-1);
      Vertex next = get_data(&result, i-1, j+1);

      double x = result.grid[j];

      Vertex* v = &solution[j];
      v->r = 1/2.0f * (next.r + prev.r) - curr.u * dt / (2.0f * dx) * (next.r - prev.r) - dt / (2.0f * dx) * curr.r * (next.u - prev.u);
      v->u = 1/2.0f * (next.u + prev.u) - curr.u * dt / (2.0f * dx) * (next.u - prev.u) - dt / (2.0f * dx) / curr.r * (next.p - prev.p);
      v->p = 1/2.0f * (next.p + prev.p) - curr.u * dt / (2.0f * dx) * (next.p - prev.p) - dt / (2.0f * dx) * gamma * curr.p * (next.u - prev.u);

      // Courant–Friedrichs–Lewy condition:
      double frac = dx / (abs(curr.u) + sqrt(gamma * curr.p / curr.r));
      if (dt > frac) {
        printf("Exiting the loop on %d layer!\n", i);
        goto end;
      }
    }

    size_t j = 0;
    solution[j]   = solution[j+1];
    solution[j].u = u_x0(0);

    j = solution.size-1;
    solution[j]   = solution[j-1];
    solution[j].u = u_x0(0);

    {
      Scoped_Lock lock(&result.data_mutex);
      array_add(&result.data, &solution);
      array_add(&result.time, t);
    }
    InterlockedExchange64((int64*) data->global_t, *(int64*) &t);
  }
end:
  return 0;
}

static double compute_p(Vertex v, double gamma) {
  return v.E * v.r * (gamma - 1);
}

static int conservative_lax_method(void* param) {
  Thread_Data1* data = (Thread_Data1*) param;
  Parameters_Laba1 params = *data->parameters; // copy

  double xmin = params.xmin; 
  double xmax = params.xmax;
  double p0   = params.p0;
  double x0   = params.x0;
  double r0   = params.r0;
  double ro0  = params.ro0;
  double gamma = params.gamma;
  double dt   = params.dt;
  double dx   = params.dx;

  size_t NUMBER_OF_POINTS_IN_X = (xmax - xmin) / dx;

  double t = 0.0; // real time.

  InterlockedExchange64((int64*) data->dt,       *(int64*) &dt);
  InterlockedExchange64((int64*) data->global_t, *(int64*) &t);

  {
    // 
    // create a grid.
    //
    array<double> grid = temp_array(double, NUMBER_OF_POINTS_IN_X);
    for (size_t i = 0; i < grid.size; i++) {
      grid[i] = xmin + i*dx;
    }

    Scoped_Lock lock(&result.data_mutex);
    result.grid = array_copy(&grid);
  }

  array<Vertex> solution = temp_array(Vertex, NUMBER_OF_POINTS_IN_X);

  {
    // 
    // apply initial conditions on t=0
    //
    for (size_t j = 0; j < solution.size; j++)  {
      double x = result.grid[j];

      Vertex* v = &solution[j];
      v->r = ro_t0(x, ro0);
      v->u =  u_t0(x);
      v->p =  p_t0(x, p0, x0, r0);
      v->E =  v->p / (v->r * (gamma-1));
    }

    Scoped_Lock lock(&result.data_mutex);
    array_add(&result.data, &solution);
    array_add(&result.time, t);
  }

  size_t i = 0;
  while (1) {
    i++;
    t += dt;

    for (size_t j = 1; j < solution.size-1; j++) {
      // 
      // internal points.
      //
      Vertex curr = get_data(&result, i-1, j);
      Vertex prev = get_data(&result, i-1, j-1);
      Vertex next = get_data(&result, i-1, j+1);

      double x = result.grid[j];

      Vertex* v = &solution[j];
      v->r = 1/2.0f * (next.r + prev.r) - dt / (2*dx) * (next.r*next.u - prev.r*prev.u);
      v->u = 1/v->r * 1/2.0f * (next.r*next.u + prev.r*prev.u) - dt / (2*dx) * (next.r*square(next.u) + next.p - prev.r*square(prev.u) - prev.p);
      v->E = -1/2.0f * square(v->u) + 1/v->r * (1/2.0f * (next.r*(next.E + 1/2.0f*square(next.u)) + prev.r*(prev.E + 1/2.0f*square(prev.u))) - dt / (2*dx) * ((next.r*next.E + compute_p(next, gamma) + 1/2.0f*next.r*square(next.u))*next.u - (prev.r*prev.E + compute_p(prev, gamma) + 1/2.0f*prev.r*square(prev.u))*prev.u));
    }

    size_t j = 0;
    solution[j]   = solution[j+1];
    solution[j].u = u_x0(0);

    j = solution.size-1;
    solution[j]   = solution[j-1];
    solution[j].u = u_x0(0);

    {
      Scoped_Lock lock(&result.data_mutex);
      array_add(&result.data, &solution);
      array_add(&result.time, t);
    }
    InterlockedExchange64((int64*) data->global_t, *(int64*) &t);
  }

  return 0;
}

static int lagrange_method(void* param) {
  Thread_Data2* data = (Thread_Data2*) param;
  Parameters_Laba2 params = *data->parameters; // copy

  double xmin = params.xmin; 
  double xmax = params.xmax;
  double p0   = params.p0;
  double x0   = params.x0;
  double r0   = params.r0;
  double ro0  = params.ro0;
  double gamma = params.gamma;
  double alpha = params.alpha;
  double a    = params.a;
  double dt   = params.dt;
  double dx   = params.dx;

  size_t NUMBER_OF_POINTS_IN_X = (xmax - xmin) / dx;
  size_t MAX_ALLOCATED_MEMORY = NUMBER_OF_POINTS_IN_X* (sizeof(Velocity_And_Euler) + sizeof(Voxel) + sizeof(double));

  Memory_Arena arena;
  begin_memory_arena(&arena, MAX_ALLOCATED_MEMORY);
  defer { end_memory_arena(&arena); };


  double t = 0.0;
  InterlockedExchange64((int64*) data->dt,       *(int64*) &dt);
  InterlockedExchange64((int64*) data->global_t, *(int64*) &t);


  {
    // 
    // create a grid.
    //
    array<double> points_grid;
    array<double> internal_grid;

    points_grid.allocator   = arena.allocator;
    internal_grid.allocator = arena.allocator;

    array_resize(&points_grid,   NUMBER_OF_POINTS_IN_X);
    array_resize(&internal_grid, points_grid.size-1);

    for (size_t i = 0; i < points_grid.size; i++) {
      points_grid[i] = i;
    }

    for (size_t i = 0; i < points_grid.size-1; i++) {
      float curr = points_grid[i];
      float next = points_grid[i+1];
      internal_grid[i] = (curr + next) / 2.0f;
    }

    {
      Scoped_Lock lock(&lagrange.data_mutex);
      lagrange.points_grid   = array_copy(&points_grid);
      lagrange.internal_grid = array_copy(&internal_grid);
    }
  }

  reset_memory_arena(&arena);

  array<Velocity_And_Euler> points_solution;
  array<Voxel>              internal_solution;
  array<double>             ques;
  points_solution.allocator   = arena.allocator;
  internal_solution.allocator = arena.allocator;
  ques.allocator = arena.allocator;

  array_resize(&points_solution,   NUMBER_OF_POINTS_IN_X);
  array_resize(&internal_solution, points_solution.size-1);
  array_resize(&ques, internal_solution.size);

  {
    //
    // apply initial conditions on t=0
    //
    for (size_t j = 0; j < points_solution.size; j++) {
      size_t i = lagrange.points_grid[j];

      Velocity_And_Euler* v = &points_solution[j];
      v->X = xmin + i*dx;
      v->u = u_t0(v->X);
    }

    for (size_t j = 0; j < internal_solution.size; j++) {
      size_t i = lagrange.internal_grid[j];

      Voxel* v = &internal_solution[j];
      v->X = xmin + i*dx;
      v->u = u_t0(v->X);
      v->r = ro_t0(v->X, ro0);
      v->u =  u_t0(v->X);
      v->p =  p_t0(v->X, p0, x0, r0);
      v->E =  v->p / (v->r * (gamma-1));
    }

    {
      Scoped_Lock lock(&lagrange.data_mutex);
      array_add(&lagrange.points_data,   &points_solution);
      array_add(&lagrange.internal_data, &internal_solution);
      array_add(&lagrange.time, t);
    }
  }

  size_t i = 0;
  while (1) {
    i++;
    t += dt;

    for (size_t j = 0; j < points_solution.size; j++) {
      double curr_X = get_point(&lagrange, i-1, j).X;

      if (j == 0 || j == points_solution.size-1) {
        // 
        // boundary points
        //  

      } else {
        // 
        // internal points
        //

        double next_x = get_grid(&lagrange, j+1);
        double curr_x = get_grid(&lagrange, j);

        double next_u = get_point(&lagrange, i-1, j+1).u;
        double curr_u = get_point(&lagrange, i-1, j).u;
        double prev_u = get_point(&lagrange, i-1, j-1).u;

        auto prev = get_internal(&lagrange, i-1, j-1/2.0);
        auto next = get_internal(&lagrange, i-1, j+1/2.0);

        double V_nc;
        double V_np;
        double V_pc;
        double V_pp;

        if (i > 1) {
          auto prev_prev = get_internal(&lagrange, i-2, j-1/2.0);
          auto next_next = get_internal(&lagrange, i-2, j+1/2.0);
          V_nc = 1 / next.r;
          V_np = 1 / next_next.r;
          V_pc = 1 / prev.r;
          V_pp = 1 / prev_prev.r;

        } else {
          V_nc = 1 / next.r;
          V_np = 1 / ro0;
          V_pc = 1 / prev.r;
          V_pp = 1 / ro0;
        }

        double du1 = next_u - curr_u;
        double du2 = curr_u - prev_u;

        double dq_next = (du1 < 0) ? (2*square(a) / (V_nc + V_np) * square(du1)) : 0;
        double dq_prev = (du2 < 0) ? (2*square(a) / (V_pc + V_pp) * square(du2)) : 0;
        double dq = dq_next - dq_prev;

        double dp = next.p - prev.p;
        ques[j] = dq_next;
        
        Velocity_And_Euler* v = &points_solution[j];
        v->u = curr_u + dt * (-1/ro0 * (dp + dq) / (next_x - curr_x) * pow((curr_X / curr_x), alpha-1.0));
        v->X = curr_X + dt * v->u; 
      }
    }

    for (size_t j = 0; j < internal_solution.size; j++) {
      auto next_x = get_grid(&lagrange, j+1);
      auto curr_x = get_grid(&lagrange, j);

      auto next_u = points_solution[j+1].u;
      auto curr_u = points_solution[j].u;

      auto next_X = points_solution[j+1].X;
      auto curr_X = points_solution[j].X;

      auto next = get_internal(&lagrange, i-1, j+1/2.0);

      if (j == 0 || j == internal_solution.size-1) {
        // 
        // boundary points
        //

      } else {
        // 
        // internal points.
        //

        Voxel* v = &internal_solution[j];
        v->X = (next_X + curr_X) / 2.0f; // @Incomplete: lerp();
        v->u = (next_u + curr_u) / 2.0f; // @Incomplete: lerp();
        v->r = ro0 * (pow(next_x, alpha) - pow(curr_x, alpha)) / (pow(next_X, alpha) - pow(curr_X, alpha));
        v->E = next.E - (next.p + ques[j]) * (1 / v->r - 1 / next.r);
        v->p = v->E * v->r * (gamma-1);
      }
    }

    size_t j = 0;
    points_solution[j]   = points_solution[j+1];
    points_solution[j].u = u_x0(0);

    j = points_solution.size-1;
    points_solution[j] = points_solution[j-1];
    points_solution[j].u = u_x0(0);

    j = 0;
    internal_solution[j] = internal_solution[j+1];
    internal_solution[j].u = u_x0(0);
    internal_solution[j].p = 0;

    j = internal_solution.size-1;
    internal_solution[j] = internal_solution[j-1];
    internal_solution[j].u = u_x0(0);
    internal_solution[j].p = 0;

    {
      Scoped_Lock lock(&lagrange.data_mutex);
      array_add(&lagrange.points_data,   &points_solution);
      array_add(&lagrange.internal_data, &internal_solution);
      array_add(&lagrange.time, t);
    }
    InterlockedExchange64((int64*) data->global_t, *(int64*) &t);
  }

  return 0;
}

static int method_2d(void* param) {
  Thread_Data3* data = (Thread_Data3*) param;
  Parameters_Laba3 params = *data->parameters; // copy

  double xmin = params.xmin;
  double xmax = params.xmax;
  double ymin = params.ymin;
  double ymax = params.ymax;
  double center_x = params.center_x;
  double center_y = params.center_y;
  double b_width  = params.a;
  double b_height = params.b;
  double Re = params.Re;
  double U0 = params.U0;
  double dx = params.dx;
  double dy = params.dy;
  double dt = params.dt;
  double t = 0.0;

  size_t NUMBER_OF_POINTS_IN_X      = (xmax - xmin) / dx;
  size_t NUMBER_OF_POINTS_IN_Y      = (ymax - ymin) / dy;
  size_t NUMBER_OF_POINTS_PER_LAYER = NUMBER_OF_POINTS_IN_X * NUMBER_OF_POINTS_IN_Y;
  size_t MAX_ALLOCATED_MEMORY = NUMBER_OF_POINTS_PER_LAYER * sizeof(Voxel2d);

  size_t center_x_index = (center_x - xmin) / (xmax - xmin) * NUMBER_OF_POINTS_IN_X;
  size_t center_y_index = (center_y - ymin) / (ymax - ymin) * NUMBER_OF_POINTS_IN_Y;

  size_t width_index  = b_width  / dx;
  size_t height_index = b_height / dy;

  size_t left   = center_x_index - width_index;
  size_t right  = center_x_index + width_index;
  size_t top    = center_y_index - height_index;;
  size_t bottom = center_y_index + height_index;

  // @Incomplete: x,y and i,j are inverted. How do I deal with that?

  InterlockedExchange64((int64*)  data->dt,        *(int64*) &dt);
  InterlockedExchange64((int64*)  data->global_t,  *(int64*) &t); 
  InterlockedExchange64((int64*) &result2d.size_x, *(int64*) &NUMBER_OF_POINTS_IN_X);
  InterlockedExchange64((int64*) &result2d.size_y, *(int64*) &NUMBER_OF_POINTS_IN_Y);

  Memory_Arena arena;
  begin_memory_arena(&arena, MAX_ALLOCATED_MEMORY);
  defer { end_memory_arena(&arena); }; // @MemoryLeak: @Incomplete: what happens to this memory if we kill a thread? 

  reset_memory_arena(&arena);

  array<Voxel2d> temp_data;
  temp_data.allocator = arena.allocator;

  array_resize(&temp_data, NUMBER_OF_POINTS_PER_LAYER);

  { // initial conditions (t = 0).
    for (size_t i = 0; i < NUMBER_OF_POINTS_PER_LAYER; i++) {
      temp_data[i].eps = 0;
      temp_data[i].psi = 0;
      temp_data[i].u   = 0;
      temp_data[i].v   = 0;
    }
    {
      Scoped_Lock lock(&result2d.data_mutex);
      array_add(&result2d.data, &temp_data);
    }

    t += dt;
    InterlockedExchange64((int64*) data->global_t, *(int64*) &t);
  }

  size_t k = 0;
  while (1) {

    memcpy(temp_data.data, result2d.data.data + k*NUMBER_OF_POINTS_PER_LAYER, NUMBER_OF_POINTS_PER_LAYER * sizeof(Voxel2d));

    for (size_t i = 1; i < NUMBER_OF_POINTS_IN_X-1; i++) {
      for (size_t j = 1; j < NUMBER_OF_POINTS_IN_Y-1; j++) {
        size_t idx = i * NUMBER_OF_POINTS_IN_Y + j;

        Voxel2d v_cur = get_data(&result2d, i,   j, k);
        Voxel2d v_ip1 = get_data(&result2d, i+1, j, k);
        Voxel2d v_im1 = get_data(&result2d, i-1, j, k);
        Voxel2d v_jp1 = get_data(&result2d, i, j+1, k);
        Voxel2d v_jm1 = get_data(&result2d, i, j-1, k);

        double u = v_cur.u;
        double v = v_cur.v;

        double a = u > 0 ? 0 : 1;
        double b = v > 0 ? 0 : 1;

        double eps = v_cur.eps - dt / dx * u * (1 - a) * (v_cur.eps - v_im1.eps) - dt / dx * a * u * (v_ip1.eps - v_cur.eps) - dt / dy * v * (1 - b) * (v_cur.eps - v_jm1.eps) - dt / dy * b * v * (v_jp1.eps - v_cur.eps) + dt / Re * ((v_ip1.eps - 2*v_cur.eps + v_im1.eps) / (dx*dx) + (v_jp1.eps - 2*v_cur.eps + v_jm1.eps) / (dy*dy)); 

        if (i >= left && i <= right && j >= top && j <= bottom) {
          continue;
        } else {
          temp_data[idx].eps = eps;
        }
      }
    }

    for (size_t k = 0; k < 5; k++) { // @Incomplete: do something about this.
      // 
      // One iteration by Jacobi.
      // 
      for (size_t i = 1; i < NUMBER_OF_POINTS_IN_X-1; i++) {
        for (size_t j = 1; j < NUMBER_OF_POINTS_IN_Y-1; j++) {
          size_t idx = i * NUMBER_OF_POINTS_IN_Y + j;

          Voxel2d v_cur = get_data(&result2d, &temp_data, i,   j);
          Voxel2d v_ip1 = get_data(&result2d, &temp_data, i+1, j);
          Voxel2d v_im1 = get_data(&result2d, &temp_data, i-1, j);
          Voxel2d v_jp1 = get_data(&result2d, &temp_data, i, j+1);
          Voxel2d v_jm1 = get_data(&result2d, &temp_data, i, j-1);

          double q   = (dx*dx) / (dy*dy);
          double psi = (v_ip1.psi + v_im1.psi + q*(v_jp1.psi + v_jm1.psi)) / (2.0 + 2.0*q) - dx*dx / (2.0 + 2.0*q) * v_cur.eps;

          if (i >= left && i <= right && j >= top && j <= bottom) {
            continue;
          } else {
            temp_data[idx].psi = psi;
          }
        }
      }
    }

    for (size_t i = 0; i < NUMBER_OF_POINTS_IN_X; i++) {
      for (size_t j = 0; j < NUMBER_OF_POINTS_IN_Y-1; j++) {
        size_t idx = i * NUMBER_OF_POINTS_IN_Y + j;

        Voxel2d v_cur = get_data(&result2d, &temp_data, i,   j);
        Voxel2d v_jp1 = get_data(&result2d, &temp_data, i, j+1);

        double u = (v_jp1.psi - v_cur.psi) / dy;

        if (i >= left && i <= right && j >= top && j <= bottom) {
          continue;
        } else {
          temp_data[idx].u = u;
        }
      }
    }

    for (size_t i = 0; i < NUMBER_OF_POINTS_IN_X-1; i++) {
      for (size_t j = 0; j < NUMBER_OF_POINTS_IN_Y; j++) {
        size_t idx = i * NUMBER_OF_POINTS_IN_Y + j;

        Voxel2d v_cur = get_data(&result2d, &temp_data, i,   j);
        Voxel2d v_ip1 = get_data(&result2d, &temp_data, i+1, j);

        double v = -(v_ip1.psi - v_cur.psi) / dx;

        if (i >= left && i <= right && j >= top && j <= bottom) {
          continue;
        } else {
          temp_data[idx].v = v;
        }
      }
    }

    for (size_t i = 0; i < NUMBER_OF_POINTS_IN_X-1; i++) {
      size_t idx;
      size_t jd1;

      idx = i * NUMBER_OF_POINTS_IN_Y + 0; 
      jd1 = i * NUMBER_OF_POINTS_IN_Y + 1;

      auto r = temp_data.data;
      temp_data[idx].eps = 2 * (r[jd1].psi - r[idx].psi - U0*dy) / (dy*dy);
      temp_data[idx].psi = r[jd1].psi - U0*dy;
      temp_data[idx].u   = U0;
      temp_data[idx].v   = 0;

      idx = i * NUMBER_OF_POINTS_IN_Y + (NUMBER_OF_POINTS_IN_Y-1); 
      jd1 = i * NUMBER_OF_POINTS_IN_Y + (NUMBER_OF_POINTS_IN_Y-2);
      temp_data[idx].eps = 2.0 * (r[jd1].psi - r[idx].psi + U0*dy) / (dy*dy);
      temp_data[idx].psi = r[jd1].psi + U0*dy;
      temp_data[idx].u   = U0;
      temp_data[idx].v   = 0;
    }

    for (size_t j = 0; j < NUMBER_OF_POINTS_IN_Y; j++) {
      size_t idx;
      size_t id1;

      idx = 0 * NUMBER_OF_POINTS_IN_Y + j; 
      id1 = 1 * NUMBER_OF_POINTS_IN_Y + j;
      temp_data[idx].eps = temp_data[id1].eps;
      temp_data[idx].psi = temp_data[id1].psi;

      idx = (NUMBER_OF_POINTS_IN_X-1) * NUMBER_OF_POINTS_IN_Y + j; 
      id1 = (NUMBER_OF_POINTS_IN_X-2) * NUMBER_OF_POINTS_IN_Y + j;
      temp_data[idx].eps = temp_data[id1].eps;
      temp_data[idx].psi = temp_data[id1].psi;
    }

    for (size_t j = left; j <= right; j++) {
      Voxel2d* jt = get_pointer(&result2d, &temp_data, j, top);
      Voxel2d* jb = get_pointer(&result2d, &temp_data, j, bottom);

      jt->u   = 0;
      jt->v   = 0;
      jt->psi = 0;
      jt->eps = 2.0 * get_data(&result2d, &temp_data, j, top-1).psi / (dy*dy);

      jb->u   = 0;
      jb->v   = 0;
      jb->psi = 0;
      jb->eps = 2.0 * get_data(&result2d, &temp_data, j, bottom+1).psi / (dy*dy);
    }

    for (size_t i = top; i <= bottom; i++) {
      Voxel2d* il = get_pointer(&result2d, &temp_data, left,  i);
      Voxel2d* ir = get_pointer(&result2d, &temp_data, right, i);

      il->u   = 0;
      il->v   = 0;
      il->psi = 0;
      il->eps = 2.0 * get_data(&result2d, &temp_data, left-1, i).psi / (dx*dx);;

      ir->u   = 0;
      ir->v   = 0;
      ir->psi = 0;
      ir->eps = 2.0 * get_data(&result2d, &temp_data, right+1, i).psi / (dx*dx);
    }

    {
      get_pointer(&result2d, &temp_data, left, bottom)->eps  = 2.0 * get_pointer(&result2d, &temp_data, left, bottom+1)->psi / (dy*dy)
                                                             + 2.0 * get_pointer(&result2d, &temp_data, left-1, bottom)->psi / (dx*dx);
      get_pointer(&result2d, &temp_data, right, bottom)->eps = 2.0 * get_pointer(&result2d, &temp_data, right, bottom+1)->psi / (dy*dy)
                                                             + 2.0 * get_pointer(&result2d, &temp_data, right+1, bottom)->psi / (dx*dx);

      get_pointer(&result2d, &temp_data, left, top)->eps  = 2.0 * get_pointer(&result2d, &temp_data, left, top-1)->psi / (dy*dy)
                                                          + 2.0 * get_pointer(&result2d, &temp_data, left-1, top)->psi / (dx*dx);
      get_pointer(&result2d, &temp_data, right, top)->eps = 2.0 * get_pointer(&result2d, &temp_data, right, top-1)->psi / (dy*dy)
                                                          + 2.0 * get_pointer(&result2d, &temp_data, right+1, top)->psi / (dx*dx);

    }

    {
      Scoped_Lock lock(&result2d.data_mutex);
      array_add(&result2d.data, &temp_data);
    }

    k += 1;
    t += dt;

    InterlockedExchange64((int64*) data->global_t, *(int64*) &t);
  }

  return 0;
}

typedef int (*Thread_Proc)(void*);

struct Input_Float_Settings {
  const char* name;
  float* pointer;
};

struct Method_Spec {
  Thread_Proc proc;
  const char* name;
};

struct The_Thing {
  Thread      thread = {};
  Thread_Proc thread_proc = {};
  void*       thread_parameter = {};

  const char* name = {};
  const char* method_name = {}; 

  Input_Float_Settings* inputs       = NULL;
  size_t                inputs_count = 0;

  Method_Spec* methods       = NULL;
  size_t       methods_count = 0;

  bool thread_is_paused = false; // @Incomplete: move this shit out of here.
  bool clear_the_thing = false;

  bool auto_fit = true;
  bool replay = false;
  int  replay_multiplier = 0;

  float time = 0;

  double global_t = 0;
  double dt = 0;
};

void render_laba1(Memory_Arena* arena, The_Thing* thing) {
  if (thing->clear_the_thing) {
    thing->clear_the_thing = false;

    result.grid.size = 0;
    result.data.size = 0;
    result.time.size = 0;

    thing->global_t = 0;
    thing->dt       = 0;
  }


  array<double> grid_to_be_drawn_this_frame;
  array<Vertex> data_to_be_drawn_this_frame;
  array<double> p_array;
  array<double> u_array;
  array<double> r_array;


  grid_to_be_drawn_this_frame.allocator = arena->allocator;
  data_to_be_drawn_this_frame.allocator = arena->allocator;
  p_array.allocator = arena->allocator;
  u_array.allocator = arena->allocator;
  r_array.allocator = arena->allocator;


  {
    Scoped_Lock mutex(&result.data_mutex);

    size_t n = thing->time / thing->dt;
    size_t needed_data_start =  n    * result.grid.size;
    size_t needed_data_end   = (n+1) * result.grid.size;

    if (result.data.size && result.grid.size) {
      array_copy_range(&data_to_be_drawn_this_frame, &result.data, needed_data_start, needed_data_end);
      array_copy      (&grid_to_be_drawn_this_frame, &result.grid);
    }
  }

  size_t data_size = data_to_be_drawn_this_frame.size;

  array_resize(&p_array, data_size);
  array_resize(&u_array, data_size);
  array_resize(&r_array, data_size);

  for (size_t i = 0; i < data_size; i++) {
    p_array[i] = data_to_be_drawn_this_frame[i].p;
    u_array[i] = data_to_be_drawn_this_frame[i].u;
    r_array[i] = data_to_be_drawn_this_frame[i].r;
  }

         const bool auto_fit = thing->auto_fit;
  static const ImVec2 plot_rect_size = ImVec2(300, 300);
  static const auto   plot_flags = thing->auto_fit ? ImPlotAxisFlags_AutoFit : ImPlotAxisFlags_None;
  if (ImPlot::BeginPlot("Pressure", plot_rect_size)) {
    ImPlot::SetupAxes("", "", plot_flags, plot_flags);
    ImPlot::PlotLine("p(x, t)", grid_to_be_drawn_this_frame.data, p_array.data, p_array.size);
    ImPlot::EndPlot();
  }

  ImGui::SameLine();

  if (ImPlot::BeginPlot("Velocity", plot_rect_size)) {
    ImPlot::SetupAxes("", "", plot_flags, plot_flags);
    ImPlot::PlotLine("u(x, t)", grid_to_be_drawn_this_frame.data, u_array.data, u_array.size);
    ImPlot::EndPlot();
  }

  ImGui::SameLine();

  if (ImPlot::BeginPlot("Density", plot_rect_size)) {
    ImPlot::SetupAxes("", "", plot_flags, plot_flags);
    ImPlot::PlotLine("r(x, t)", grid_to_be_drawn_this_frame.data, r_array.data, r_array.size);
    ImPlot::EndPlot();
  }
}

void render_laba2(Memory_Arena* arena, The_Thing* thing) {
  if (thing->clear_the_thing) {
    thing->clear_the_thing = false;

    lagrange.points_grid.size = 0;
    lagrange.internal_grid.size = 0;
    lagrange.time.size = 0;
    lagrange.points_data.size = 0;
    lagrange.internal_data.size = 0;

    thing->global_t = 0;
    thing->dt       = 0;
  }


  array<double> grid_to_be_drawn_this_frame;
  array<Voxel> data_to_be_drawn_this_frame;
  array<double> X_array;
  array<double> E_array;
  array<double> p_array;
  array<double> u_array;
  array<double> r_array;

  grid_to_be_drawn_this_frame.allocator = arena->allocator;
  data_to_be_drawn_this_frame.allocator = arena->allocator;
  X_array.allocator = arena->allocator;
  E_array.allocator = arena->allocator;
  p_array.allocator = arena->allocator;
  u_array.allocator = arena->allocator;
  r_array.allocator = arena->allocator;

  {
    Scoped_Lock mutex(&lagrange.data_mutex);

    size_t n = thing->time / thing->dt;
    size_t needed_data_start =  n    * lagrange.internal_grid.size;
    size_t needed_data_end   = (n+1) * lagrange.internal_grid.size;

    if (lagrange.internal_data.size && lagrange.internal_grid.size) {
      array_copy_range(&data_to_be_drawn_this_frame, &lagrange.internal_data, needed_data_start, needed_data_end);
      array_copy      (&grid_to_be_drawn_this_frame, &lagrange.internal_grid);
    }
  }

  size_t data_size = data_to_be_drawn_this_frame.size;

  array_resize(&X_array, data_size);
  array_resize(&E_array, data_size);
  array_resize(&p_array, data_size);
  array_resize(&u_array, data_size);
  array_resize(&r_array, data_size);

  for (size_t i = 0; i < data_size; i++) {
    X_array[i] = data_to_be_drawn_this_frame[i].X;
    E_array[i] = data_to_be_drawn_this_frame[i].E;
    p_array[i] = data_to_be_drawn_this_frame[i].p;
    u_array[i] = data_to_be_drawn_this_frame[i].u;
    r_array[i] = data_to_be_drawn_this_frame[i].r;
  }

         const bool auto_fit = thing->auto_fit;
  static const ImVec2 plot_rect_size = ImVec2(300, 300);
  static const auto   plot_flags = auto_fit ? ImPlotAxisFlags_AutoFit : ImPlotAxisFlags_None;
  if (ImPlot::BeginPlot("Euler coordinate R(r, t)", plot_rect_size)) {
    ImPlot::SetupAxes("", "", plot_flags, plot_flags);
    ImPlot::PlotLine("R(r, t)", grid_to_be_drawn_this_frame.data, X_array.data, X_array.size);
    ImPlot::EndPlot();
  }

  if (ImPlot::BeginPlot("Velocity : u(r, t)", plot_rect_size)) {
    ImPlot::SetupAxes("", "", plot_flags, plot_flags);
    ImPlot::PlotLine("u(r, t)", grid_to_be_drawn_this_frame.data, u_array.data, u_array.size);
    ImPlot::EndPlot();
  }

  ImGui::SameLine();

  if (ImPlot::BeginPlot("Pressure : p(r, t)", plot_rect_size)) {
    ImPlot::SetupAxes("", "", plot_flags, plot_flags);
    ImPlot::PlotLine("p(r, t)", grid_to_be_drawn_this_frame.data, p_array.data, p_array.size);
    ImPlot::EndPlot();
  }

  ImGui::SameLine();

  if (ImPlot::BeginPlot("Energy : E(r, t)", plot_rect_size)) {
    ImPlot::SetupAxes("", "", plot_flags, plot_flags);
    ImPlot::PlotLine("E(r, t)", grid_to_be_drawn_this_frame.data, E_array.data, E_array.size);
    ImPlot::EndPlot();
  }

  ImGui::SameLine();

  if (ImPlot::BeginPlot("Density : ro(r, t)", plot_rect_size)) {
    ImPlot::SetupAxes("", "", plot_flags, plot_flags);
    ImPlot::PlotLine("ro(r, t)", grid_to_be_drawn_this_frame.data, r_array.data, r_array.size);
    ImPlot::EndPlot();
  }

  if (ImPlot::BeginPlot("Velocity : u(R, t)", plot_rect_size)) {
    ImPlot::SetupAxes("", "", plot_flags, plot_flags);
    ImPlot::PlotLine("u(R, t)", X_array.data, u_array.data, u_array.size);
    ImPlot::EndPlot();
  }

  ImGui::SameLine();

  if (ImPlot::BeginPlot("Pressure : p(R, t)", plot_rect_size)) {
    ImPlot::SetupAxes("", "", plot_flags, plot_flags);
    ImPlot::PlotLine("p(R, t)", X_array.data, p_array.data, p_array.size);
    ImPlot::EndPlot();
  }

  ImGui::SameLine();

  if (ImPlot::BeginPlot("Energy : E(R, t)", plot_rect_size)) {
    ImPlot::SetupAxes("", "", plot_flags, plot_flags);
    ImPlot::PlotLine("E(R, t)", X_array.data, E_array.data, E_array.size);
    ImPlot::EndPlot();
  }

  ImGui::SameLine();

  if (ImPlot::BeginPlot("Density : ro(R, t)", plot_rect_size)) {
    ImPlot::SetupAxes("", "", plot_flags, plot_flags);
    ImPlot::PlotLine("ro(R, t)", X_array.data, r_array.data, r_array.size);
    ImPlot::EndPlot();
  }
}

void render_laba3(Memory_Arena* arena, The_Thing* thing) {
  if (thing->clear_the_thing) {
    thing->clear_the_thing = false;

    result2d.data.size = 0;
    result2d.size_x = 0;
    result2d.size_y = 0;

    thing->global_t = 0;
    thing->dt       = 0;
  }


  ImPlotColormap  map = ImPlotColormap_Plasma;
  ImPlotAxisFlags axes_flags = ImPlotAxisFlags_Lock | ImPlotAxisFlags_NoGridLines | ImPlotAxisFlags_NoTickMarks;

  ImPlot::PushColormap(map);

  ImPlot::ColormapScale("##HeatScale", 0.0, 1.0, ImVec2(60,600));

  ImGui::SameLine();

  array<Voxel2d> heatmap_to_draw_this_frame;
  array<double> velocity;

  size_t size_x;
  size_t size_y;

  heatmap_to_draw_this_frame.allocator = arena->allocator;
  velocity.allocator = arena->allocator;


  {
    Scoped_Lock mutex(&result2d.data_mutex);

    size_x = result2d.size_x;
    size_y = result2d.size_y;

    size_t n = thing->time / thing->dt;
    size_t start =  n    * size_x * size_y;
    size_t end   = (n+1) * size_x * size_y;

    if (result2d.data.size && start <= result2d.data.size && end <= result2d.data.size) { // @Incomplete: wtf is that.
      array_copy_range(&heatmap_to_draw_this_frame, &result2d.data, start, end);
    }
  }

  array_resize(&velocity, heatmap_to_draw_this_frame.size);

  for (size_t i = 0; i < velocity.size; i++) {
    double u = heatmap_to_draw_this_frame[i].u;
    double v = heatmap_to_draw_this_frame[i].v;

    velocity[i] = sqrt(u*u + v*v);
  }

  if (ImPlot::BeginPlot("##Heatmap2", ImVec2(600,600))) {
    ImPlot::SetupAxes(NULL, NULL, ImPlotAxisFlags_NoDecorations, ImPlotAxisFlags_NoDecorations);
    ImPlot::PlotHeatmap("heat1", velocity.data, size_x, size_y, 0.0f, 1.0f, NULL);
    ImPlot::EndPlot();
  }

  ImPlot::PopColormap();

}

int main(int, char**) {
  WNDCLASSEX wc = { sizeof(WNDCLASSEX), CS_CLASSDC, WndProc, 0L, 0L, GetModuleHandle(NULL), NULL, NULL, NULL, NULL, _T("ImGui Example"), NULL };
  ::RegisterClassEx(&wc);
  HWND hwnd = ::CreateWindow(wc.lpszClassName, _T("Dear ImGui DirectX11 Example"), WS_OVERLAPPEDWINDOW ^ WS_THICKFRAME, 100, 100, 1280, 800, NULL, NULL, wc.hInstance, NULL);

  if (!CreateDeviceD3D(hwnd))
  {
    CleanupDeviceD3D();
    ::UnregisterClass(wc.lpszClassName, wc.hInstance);
    return 1;
  }

  ::ShowWindow(hwnd, SW_SHOWDEFAULT);
  ::UpdateWindow(hwnd);

  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImPlot::CreateContext();
  ImGuiIO& io = ImGui::GetIO();

  // Setup Dear ImGui style
  ImGui::StyleColorsLight();

  // Setup Platform/Renderer backends
  ImGui_ImplWin32_Init(hwnd);
  ImGui_ImplDX11_Init(g_pd3dDevice, g_pd3dDeviceContext);


  bool show_imgui_demo  = false;
  bool show_implot_demo = false;
  bool done             = false;

  init_threads_api();
  check_threads_api();

  Parameters_Laba1 parameters_laba1 = {};
  Parameters_Laba2 parameters_laba2 = {};
  Parameters_Laba3 parameters_laba3 = {};

  Input_Float_Settings input_parameters_laba1[] = {
    { "xmin", &parameters_laba1.xmin },
    { "xmax", &parameters_laba1.xmax },
    { "p0",   &parameters_laba1.p0 },
    { "x0",   &parameters_laba1.x0 },
    { "r0",   &parameters_laba1.r0 },
    { "ro0",  &parameters_laba1.ro0 },
    { "gamma", &parameters_laba1.gamma },
    { "dx",   &parameters_laba1.dx },
    { "dt",   &parameters_laba1.dt },
  };

  Input_Float_Settings input_parameters_laba2[] = {
    { "xmin", &parameters_laba2.xmin },
    { "xmax", &parameters_laba2.xmax },
    { "p0",   &parameters_laba2.p0 },
    { "x0",   &parameters_laba2.x0 },
    { "r0",   &parameters_laba2.r0 },
    { "ro0",  &parameters_laba2.ro0 },
    { "gamma", &parameters_laba2.gamma },
    { "alpha", &parameters_laba2.alpha },
    { "a",     &parameters_laba2.a },
    { "dx",   &parameters_laba2.dx },
    { "dt",   &parameters_laba2.dt },
  };

  Input_Float_Settings input_parameters_laba3[] = {
    { "xmin",   &parameters_laba3.xmin },
    { "xmax",   &parameters_laba3.xmax },
    { "ymin",   &parameters_laba3.ymin },
    { "ymax",   &parameters_laba3.ymax },
    { "barrier center x", &parameters_laba3.center_x },
    { "barrier center y", &parameters_laba3.center_y },
    { "barrier width",  &parameters_laba3.a },
    { "barrier height", &parameters_laba3.b },
    { "Re",     &parameters_laba3.Re },
    { "U0",     &parameters_laba3.U0 },
    { "dx",     &parameters_laba3.dx },
    { "dy",     &parameters_laba3.dy },
    { "dt",     &parameters_laba3.dt },
  };


  Method_Spec methods_laba1[] = {
    { non_conservative_lax_method, "Non Conservative Lax Method" },
    { euler_upwind_method,        "Euler Upwind Method"         },
    { conservative_lax_method,    "Conservative Lax Method"     },
  };

  Method_Spec methods_laba2[] = {
    { lagrange_method, "Lagrange Method" },
  };

  Method_Spec methods_laba3[] = {
    { method_2d, "2D Method" },
  };

  The_Thing things[3];

  Thread_Data1 thread_data1 = { &parameters_laba1, &things[0].global_t, &things[0].dt };
  Thread_Data2 thread_data2 = { &parameters_laba2, &things[1].global_t, &things[1].dt };
  Thread_Data3 thread_data3 = { &parameters_laba3, &things[2].global_t, &things[2].dt };


  { // init methods

    things[0].thread_proc = methods_laba1[0].proc;
    things[0].thread_parameter = &thread_data1;
    things[0].name = "Laba 1";
    things[0].method_name = methods_laba1[0].name;
    things[0].inputs        = input_parameters_laba1;
    things[0].inputs_count  = array_size(input_parameters_laba1);
    things[0].methods       = methods_laba1;
    things[0].methods_count = array_size(methods_laba1);


    things[1].thread_proc = methods_laba2[0].proc;
    things[1].thread_parameter = &thread_data2;
    things[1].name        = "Laba 2";
    things[1].method_name = methods_laba2[0].name;
    things[1].inputs        = input_parameters_laba2;
    things[1].inputs_count  = array_size(input_parameters_laba2);
    things[1].methods       = methods_laba2;
    things[1].methods_count = array_size(methods_laba2);

    things[2].thread_proc = methods_laba3[0].proc;
    things[2].thread_parameter = &thread_data3;
    things[2].name        = "Laba 3";
    things[2].method_name = methods_laba3[0].name;
    things[2].inputs        = input_parameters_laba3;
    things[2].inputs_count  = array_size(input_parameters_laba3);
    things[2].methods       = methods_laba3;
    things[2].methods_count = array_size(methods_laba3);

    // init parameters
    parameters_laba1.xmin = -3.0f;
    parameters_laba1.xmax = 3.0f;
    parameters_laba1.p0   = 1.0f;
    parameters_laba1.x0   = 0.0f;
    parameters_laba1.r0   = 0.5f;
    parameters_laba1.ro0  = 1.0f;
    parameters_laba1.gamma = 1.67f;
    parameters_laba1.dt   = 0.0001f;
    parameters_laba1.dx   = 0.01f;

    parameters_laba2.xmin = -3.0f;
    parameters_laba2.xmax = 3.0f;
    parameters_laba2.p0   = 1.0f;
    parameters_laba2.x0   = 0.0f;
    parameters_laba2.r0   = 0.5f;
    parameters_laba2.ro0  = 1.0f;
    parameters_laba2.gamma = 1.67f;
    parameters_laba2.alpha = 1;
    parameters_laba2.a     = 2.3;
    parameters_laba2.dt   = 0.0001f;
    parameters_laba2.dx   = 0.01f;

    parameters_laba3.xmin = -3;
    parameters_laba3.xmax =  3;
    parameters_laba3.ymin = -3;
    parameters_laba3.ymax =  3;
    parameters_laba3.center_x = -1.5;
    parameters_laba3.center_y = 0;
    parameters_laba3.a = 0.2;
    parameters_laba3.b = 0.2;
    parameters_laba3.Re = 500;
    parameters_laba3.U0 = 1;
    parameters_laba3.dx = 0.05;
    parameters_laba3.dy = 0.05;
    parameters_laba3.dt = 0.01;

    result.data_mutex = create_mutex();
    lagrange.data_mutex = create_mutex();
  }

  Memory_Arena temporary_storage;
  begin_memory_arena(&temporary_storage, MB(200));

  ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);
  while (!done) {
    MSG msg;
    while (::PeekMessage(&msg, NULL, 0U, 0U, PM_REMOVE)) {
      ::TranslateMessage(&msg);
      ::DispatchMessage(&msg);
      if (msg.message == WM_QUIT) done = true;
    }
    if (done) { break; }
        
    auto framerate = io.Framerate;

    reset_memory_arena(&temporary_storage);

    // Start the Dear ImGui frame
    ImGui_ImplDX11_NewFrame();
    ImGui_ImplWin32_NewFrame();
    ImGui::NewFrame();

    {
      ImGui::Begin("Control window");

      ImGui::Checkbox("ImGui Demo",  &show_imgui_demo);
      ImGui::Checkbox("ImPlot Demo", &show_implot_demo);

      for (size_t i = 0; i < array_size(things); i++) {
        The_Thing* thing = &things[i];

        if (ImGui::CollapsingHeader(thing->name)) {
          for (size_t i = 0; i < thing->inputs_count; i++) {
            Input_Float_Settings* s = &thing->inputs[i];

            static const float step      = 0.0f;
            static const float step_fast = 0.0f;
            static const char* format    = "%g";
            static const ImGuiInputTextFlags flags = ImGuiInputTextFlags_CharsScientific;
            ImGui::InputFloat(s->name, s->pointer, step, step_fast, format, flags);
          }

          ImGui::SliderFloat("Time", &thing->time, 0, thing->global_t);
          ImGui::Checkbox("Auto Fit", &thing->auto_fit);

          const char* start_or_continue = thing->thread_is_paused ? "Continue" : "Start";
          if (ImGui::Button(start_or_continue)) {
            if (thing->thread_is_paused) {
              resume_thread(&thing->thread);
              thing->thread_is_paused = false;
            } else {
              if (!is_thread_running(&thing->thread)) {
                start_thread(&thing->thread, thing->thread_proc, thing->thread_parameter);
              }
            }
          }

          ImGui::SameLine();

          if (ImGui::Button("Pause")) {
            if (!thing->thread_is_paused && is_thread_running(&thing->thread)) {
              suspend_thread(&thing->thread);
            }
            thing->thread_is_paused = true;
          }

          ImGui::SameLine();

          if (ImGui::Button("Stop")) {
            if (is_thread_running(&thing->thread)) {
              kill_thread(&thing->thread);
            }
            thing->thread_is_paused = false;
          }

          ImGui::SameLine();

          if (ImGui::Button("Clear")) {
            if (is_thread_running(&thing->thread)) { // just kill a thread.
              kill_thread(&thing->thread);
              thing->thread_is_paused = false;
            }
            thing->clear_the_thing = true;
          }

          if (ImGui::Button("Start Playing")) {
            thing->replay = true;
          }

          ImGui::SameLine();

          if (ImGui::Button("Stop Playing")) {
            thing->replay = false;
          }

          ImGui::SameLine();

          if (ImGui::Button("x2")) {
            thing->replay_multiplier += 1;
          }

          ImGui::SameLine();

          if (ImGui::Button("/2")) {
            thing->replay_multiplier -= 1;
          }

          ImGui::SameLine();

          ImGui::Text("%s%g", (thing->replay_multiplier >= 0) ? "x" : "/", pow(2.0, abs(thing->replay_multiplier)));

          for (size_t i = 0; i < thing->methods_count; i++) {
            Method_Spec* method = &thing->methods[i];
            if (i != 0) { ImGui::SameLine(); }
            if (ImGui::Button(method->name)) {
              thing->thread_proc = method->proc;
              thing->method_name = method->name;
            }
          }

          ImGui::Text("Using %s", thing->method_name);

          if (thing->replay) {
            double m = pow(2.0, thing->replay_multiplier);
            thing->time += m * 1/framerate;
          }

          thing->time = (thing->time <= thing->global_t) ? thing->time : 0;

          switch(i) {
          case 0: render_laba1(&temporary_storage, thing); break;
          case 1: render_laba2(&temporary_storage, thing); break;
          case 2: render_laba3(&temporary_storage, thing); break;
          }
        }
      }

      ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f/framerate, framerate);
      ImGui::End();
    }

    if (show_imgui_demo) {
      ImGui::ShowDemoWindow(&show_imgui_demo);
    }

    if (show_implot_demo) {
      ImPlot::ShowDemoWindow();
    }


    // Rendering
    ImGui::Render();
    const float clear_color_with_alpha[4] = { clear_color.x * clear_color.w, clear_color.y * clear_color.w, clear_color.z * clear_color.w, clear_color.w };
    g_pd3dDeviceContext->OMSetRenderTargets(1, &g_mainRenderTargetView, NULL);
    g_pd3dDeviceContext->ClearRenderTargetView(g_mainRenderTargetView, clear_color_with_alpha);
    ImGui_ImplDX11_RenderDrawData(ImGui::GetDrawData());

    g_pSwapChain->Present(1, 0); // Present with vsync
    //g_pSwapChain->Present(0, 0); // Present without vsync
  }

  end_memory_arena(&temporary_storage);

  // Cleanup
  ImGui_ImplDX11_Shutdown();
  ImGui_ImplWin32_Shutdown();
  ImPlot::DestroyContext();
  ImGui::DestroyContext();

  CleanupDeviceD3D();
  ::DestroyWindow(hwnd);
  ::UnregisterClass(wc.lpszClassName, wc.hInstance);
  return 0;
}

// Helper functions

bool CreateDeviceD3D(HWND hWnd)
{
    // Setup swap chain
    DXGI_SWAP_CHAIN_DESC sd;
    ZeroMemory(&sd, sizeof(sd));
    sd.BufferCount = 2;
    sd.BufferDesc.Width = 0;
    sd.BufferDesc.Height = 0;
    sd.BufferDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
    sd.BufferDesc.RefreshRate.Numerator = 60;
    sd.BufferDesc.RefreshRate.Denominator = 1;
    sd.Flags = DXGI_SWAP_CHAIN_FLAG_ALLOW_MODE_SWITCH;
    sd.BufferUsage = DXGI_USAGE_RENDER_TARGET_OUTPUT;
    sd.OutputWindow = hWnd;
    sd.SampleDesc.Count = 1;
    sd.SampleDesc.Quality = 0;
    sd.Windowed = TRUE;
    sd.SwapEffect = DXGI_SWAP_EFFECT_DISCARD;

    UINT createDeviceFlags = 0;
    //createDeviceFlags |= D3D11_CREATE_DEVICE_DEBUG;
    D3D_FEATURE_LEVEL featureLevel;
    const D3D_FEATURE_LEVEL featureLevelArray[2] = { D3D_FEATURE_LEVEL_11_0, D3D_FEATURE_LEVEL_10_0, };
    if (D3D11CreateDeviceAndSwapChain(NULL, D3D_DRIVER_TYPE_HARDWARE, NULL, createDeviceFlags, featureLevelArray, 2, D3D11_SDK_VERSION, &sd, &g_pSwapChain, &g_pd3dDevice, &featureLevel, &g_pd3dDeviceContext) != S_OK)
        return false;

    CreateRenderTarget();
    return true;
}

void CleanupDeviceD3D()
{
    CleanupRenderTarget();
    if (g_pSwapChain) { g_pSwapChain->Release(); g_pSwapChain = NULL; }
    if (g_pd3dDeviceContext) { g_pd3dDeviceContext->Release(); g_pd3dDeviceContext = NULL; }
    if (g_pd3dDevice) { g_pd3dDevice->Release(); g_pd3dDevice = NULL; }
}

void CreateRenderTarget()
{
    ID3D11Texture2D* pBackBuffer;
    g_pSwapChain->GetBuffer(0, IID_PPV_ARGS(&pBackBuffer));
    g_pd3dDevice->CreateRenderTargetView(pBackBuffer, NULL, &g_mainRenderTargetView);
    pBackBuffer->Release();
}

void CleanupRenderTarget()
{
    if (g_mainRenderTargetView) { g_mainRenderTargetView->Release(); g_mainRenderTargetView = NULL; }
}

// Forward declare message handler from imgui_impl_win32.cpp
extern IMGUI_IMPL_API LRESULT ImGui_ImplWin32_WndProcHandler(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam);

// Win32 message handler
// You can read the io.WantCaptureMouse, io.WantCaptureKeyboard flags to tell if dear imgui wants to use your inputs.
// - When io.WantCaptureMouse is true, do not dispatch mouse input data to your main application, or clear/overwrite your copy of the mouse data.
// - When io.WantCaptureKeyboard is true, do not dispatch keyboard input data to your main application, or clear/overwrite your copy of the keyboard data.
// Generally you may always pass all inputs to dear imgui, and hide them from your application based on those two flags.
LRESULT WINAPI WndProc(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam)
{
    if (ImGui_ImplWin32_WndProcHandler(hWnd, msg, wParam, lParam))
        return true;

    switch (msg)
    {
    case WM_SIZE:
        if (g_pd3dDevice != NULL && wParam != SIZE_MINIMIZED)
        {
            CleanupRenderTarget();
            g_pSwapChain->ResizeBuffers(0, (UINT)LOWORD(lParam), (UINT)HIWORD(lParam), DXGI_FORMAT_UNKNOWN, 0);
            CreateRenderTarget();
        }
        return 0;
    case WM_SYSCOMMAND:
        if ((wParam & 0xfff0) == SC_KEYMENU) // Disable ALT application menu
            return 0;
        break;
    case WM_DESTROY:
        ::PostQuitMessage(0);
        return 0;
    }
    return ::DefWindowProc(hWnd, msg, wParam, lParam);
}
