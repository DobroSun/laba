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
  float xmin = 0.0f;
  float xmax = 0.0f;

  float p0  = 0.0f;
  float x0  = 0.0f;
  float r0  = 0.0f;
  float ro0 = 0.0f;
  float gamma = 0.0f;

  float dt = 0.0f;
  float dx = 0.0f;
};

struct Parameters_Laba2 {
  float xmin = 0.0f;
  float xmax = 0.0f;

  float p0  = 0.0f;
  float x0  = 0.0f;
  float r0  = 0.0f;
  float ro0 = 0.0f;
  float gamma = 0.0f;
  float alpha = 0;

  float dt = 0.0f;
  float dx = 0.0f;
};

struct Vertex {
  float r;
  float u;
  union {
    float p;
    float E;
  };
};

struct Input_Float_Settings {
  const char* name;
  float* pointer;
};

struct Input_Parameter_Settings {
  Input_Float_Settings* inputs;
  size_t count;
};

struct Result {
  Mutex data_mutex = {};
  float dt;
  float  t;
  array<float> grid; // for now we only have 1 dimension.
  array<float> time; // time domain;
  array<Vertex> data; // this is going to evolve in time. size := NUMBER_OF_LAYERS * grid.size()
};

Vertex get_data(Result* r, size_t t, size_t j) {
  return r->data[t * r->grid.size + j];
}

static Result result;

static int euler_upwind_method(void* param) {
  Parameters_Laba1 data = *(Parameters_Laba1*) param; // copy

  float xmin = data.xmin; 
  float xmax = data.xmax;
  float p0   = data.p0;
  float x0   = data.x0;
  float r0   = data.r0;
  float ro0  = data.ro0;
  float gamma = data.gamma;
  float dt   = data.dt;
  float dx   = data.dx;

  size_t NUMBER_OF_POINTS_IN_X = (xmax - xmin) / dx;

  float t = 0.0; // real time.

  {
    // 
    // create a grid.
    //
    array<float> grid = temp_array(float, NUMBER_OF_POINTS_IN_X);
    for (size_t i = 0; i < grid.size; i++) {
      grid[i] = xmin + i*dx;
    }

    {
      Scoped_Lock lock(&result.data_mutex);
      result.grid = array_copy(&grid);
    }
    InterlockedExchange((uint*) &result.dt, *(uint*) &dt);
  }

  array<Vertex> solution = temp_array(Vertex, NUMBER_OF_POINTS_IN_X);

  {
    // 
    // apply initial conditions on t=0
    //
    for (size_t j = 0; j < solution.size; j++)  {
      float x = result.grid[j];

      Vertex* v = &solution[j];
      v->r = ro_t0(x, ro0);
      v->u =  u_t0(x);
      v->p =  p_t0(x, p0, x0, r0);
    }

    {
      Scoped_Lock lock(&result.data_mutex);
      array_add(&result.data, &solution);
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

      float x = result.grid[j];
      float a = (curr.u >= 0) ? 0 : 1;

      Vertex* v = &solution[j];
      v->r = curr.r - (1 - a) * curr.u * dt / dx * (curr.r - prev.r) - a * curr.u * dt / dx * (next.r - curr.r) - dt / (2*dx) * curr.r * (next.u - prev.u);
      v->u = curr.u - (1 - a) * curr.u * dt / dx * (curr.u - prev.u) - a * curr.u * dt / dx * (next.u - curr.u) - dt / (2*dx) / curr.r * (next.p - prev.p);;
      v->p = curr.p - (1 - a) * curr.u * dt / dx * (curr.p - prev.p) - a * curr.u * dt / dx * (next.p - curr.p) - dt / (2*dx) * gamma * curr.p * (next.u - prev.u);
    }


    size_t j;
    {
      j = 0;
      Vertex* v = &solution[j];
      v->r = solution[j+1].r;
      v->u = u_x0(t);
      v->p = solution[j+1].p;
    }
    {
      j = solution.size-1;
      Vertex* v = &solution[j];
      v->r = solution[j-1].r;
      v->u = u_x0(t);
      v->p = solution[j-1].p;
    }

    {
      Scoped_Lock lock(&result.data_mutex);
      array_add(&result.data, &solution);
      array_add(&result.time, t);
    }
    InterlockedExchange((uint*) &result.t, *(uint*) &t);
  }
  return 0;
}

static int nonconservative_lax_method(void* param) {
  Parameters_Laba1 data = *(Parameters_Laba1*) param; // copy

  float xmin = data.xmin; 
  float xmax = data.xmax;
  float p0   = data.p0;
  float x0   = data.x0;
  float r0   = data.r0;
  float ro0  = data.ro0;
  float gamma = data.gamma;
  float dt   = data.dt;
  float dx   = data.dx;

  size_t NUMBER_OF_POINTS_IN_X = (xmax - xmin) / dx;

  float t = 0.0; // real time.

  {
    // 
    // create a grid.
    //
    array<float> grid = temp_array(float, NUMBER_OF_POINTS_IN_X);
    for (size_t i = 0; i < grid.size; i++) {
      grid[i] = xmin + i*dx;
    }

    {
      Scoped_Lock lock(&result.data_mutex);
      result.grid = array_copy(&grid);
    }
    InterlockedExchange((uint*) &result.dt, *(uint*) &dt);
  }

  array<Vertex> solution = temp_array(Vertex, NUMBER_OF_POINTS_IN_X);

  {
    // 
    // apply initial conditions on t=0
    //
    for (size_t j = 0; j < solution.size; j++)  {
      float x = result.grid[j];

      Vertex* v = &solution[j];
      v->r = ro_t0(x, ro0);
      v->u =  u_t0(x);
      v->p =  p_t0(x, p0, x0, r0);
    }

    {
      Scoped_Lock lock(&result.data_mutex);
      array_add(&result.data, &solution);
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

      float x = result.grid[j];

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

    size_t j;
    {
      j = 0;
      Vertex* v = &solution[j];
      v->r = solution[j+1].r;
      v->u = u_x0(t);
      v->p = solution[j+1].p;
    }
    {
      j = solution.size-1;
      Vertex* v = &solution[j];
      v->r = solution[j-1].r;
      v->u = u_x0(t);
      v->p = solution[j-1].p;
    }

    {
      Scoped_Lock lock(&result.data_mutex);
      array_add(&result.data, &solution);
      array_add(&result.time, t);
    }
    InterlockedExchange((uint*) &result.t, *(uint*) &t);
  }
end:
  return 0;
}

static double compute_p(Vertex v, double gamma) {
  return v.E * v.r * (gamma - 1);
}

static int conservative_lax_method(void* param) {
  Parameters_Laba1 data = *(Parameters_Laba1*) param; // copy

  float xmin = data.xmin; 
  float xmax = data.xmax;
  float p0   = data.p0;
  float x0   = data.x0;
  float r0   = data.r0;
  float ro0  = data.ro0;
  float gamma = data.gamma;
  float dt   = data.dt;
  float dx   = data.dx;

  size_t NUMBER_OF_POINTS_IN_X = (xmax - xmin) / dx;

  float t = 0.0; // real time.

  {
    // 
    // create a grid.
    //
    array<float> grid = temp_array(float, NUMBER_OF_POINTS_IN_X);
    for (size_t i = 0; i < grid.size; i++) {
      grid[i] = xmin + i*dx;
    }

    {
      Scoped_Lock lock(&result.data_mutex);
      result.grid = array_copy(&grid);
    }
    InterlockedExchange((uint*) &result.dt, *(uint*) &dt);
  }

  array<Vertex> solution = temp_array(Vertex, NUMBER_OF_POINTS_IN_X);

  {
    // 
    // apply initial conditions on t=0
    //
    for (size_t j = 0; j < solution.size; j++)  {
      float x = result.grid[j];

      Vertex* v = &solution[j];
      v->r = ro_t0(x, ro0);
      v->u =  u_t0(x);
      v->p =  p_t0(x, p0, x0, r0);
      v->E =  v->p / (v->r * (gamma-1));
    }

    {
      Scoped_Lock lock(&result.data_mutex);
      array_add(&result.data, &solution);
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

      float x = result.grid[j];

      Vertex* v = &solution[j];
      v->r = 1/2.0f * (next.r + prev.r) - dt / (2*dx) * (next.r*next.u - prev.r*prev.u);
      v->u = 1/v->r * 1/2.0f * (next.r*next.u + prev.r*prev.u) - dt / (2*dx) * (next.r*square(next.u) + next.p - prev.r*square(prev.u) - prev.p);
      v->E = -1/2.0f * square(v->u) + 1/v->r * (1/2.0f * (next.r*(next.E + 1/2.0f*square(next.u)) + prev.r*(prev.E + 1/2.0f*square(prev.u))) - dt / (2*dx) * ((next.r*next.E + compute_p(next, gamma) + 1/2.0f*next.r*square(next.u))*next.u - (prev.r*prev.E + compute_p(prev, gamma) + 1/2.0f*prev.r*square(prev.u))*prev.u));
    }

    size_t j;
    {
      j = 0;
      Vertex* v = &solution[j];
      v->r = solution[j+1].r;
      v->u = u_x0(t);
      v->E = solution[j+1].E;
    }
    {
      j = solution.size-1;
      Vertex* v = &solution[j];
      v->r = solution[j-1].r;
      v->u = u_x0(t);
      v->E = solution[j-1].E;
    }

    {
      Scoped_Lock lock(&result.data_mutex);
      array_add(&result.data, &solution);
      array_add(&result.time, t);
    }
    InterlockedExchange((uint*) &result.t, *(uint*) &t);
  }

  return 0;
}

struct Voxel {
  float X;
  float u;
  float r;
  float E;
  float p;
};

struct Result_Lagrange {
  Mutex data_mutex = {};

  float t;  // @Incomplete: move this shit out of here, keep it in The_Thing;
  float dt; // @Incomplete: same

  array<float> grid; 
  array<float> internal_points_grid;
  array<Voxel> data;
};

static Result_Lagrange lagrange;

static int lagrange_method(void* param) {
  Parameters_Laba2 data = *(Parameters_Laba2*) param; // copy

  float xmin = data.xmin; 
  float xmax = data.xmax;
  float p0   = data.p0;
  float x0   = data.x0;
  float r0   = data.r0;
  float ro0  = data.ro0;
  float gamma = data.gamma;
  float alpha = data.alpha;
  float dt   = data.dt;
  float dx   = data.dx;

  size_t NUMBER_OF_POINTS_IN_X = (xmax - xmin) / dx;

  float t = 0;

  {
    // 
    // create a grid.
    //
    array<float> grid = temp_array(float, NUMBER_OF_POINTS_IN_X);
    for (size_t i = 0; i < grid.size; i++) {
      grid[i] = xmin + i*dx;
    }

    {
      Scoped_Lock lock(&lagrange.data_mutex);
      lagrange.grid = array_copy(&grid);
    }
    InterlockedExchange((uint*) &lagrange.dt, *(uint*) &dt);
  }

  array<Voxel> solution = temp_array(Voxel, NUMBER_OF_POINTS_IN_X);

#if 0
  {
    // 
    // apply initial conditions on t=0
    //
    for (size_t j = 0; j < solution.size; j++)  {
      double x = result.grid[j];

      Vertex* v = &data[j];
      v->r = ro_t0(x, ro0);
      v->u =  u_t0(x);
      v->p =  p_t0(x, p0, x0, r0);
    }

    {
      Scoped_Lock lock(&lagrange.data_mutex);
      array_add(&lagrange.data, &solution);
    }
  }
#endif


  size_t i = 0;
  while (1) {
    i++;
    t += dt;

    InterlockedExchange((uint*) &lagrange.t, *(uint*) &t);
  }

  return 0;
}

typedef int (*Thread_Proc)(void*);

struct Method_Spec {
  Thread_Proc proc;
  const char* name;
};

struct The_Thing {
  Thread      thread = {};
  Thread_Proc thread_proc = {};
  void*       thread_parameter = {};

  const char* name = {};
  const char* method_name = "";

  Input_Parameter_Settings settings = {};

  bool thread_is_paused = false; // @Incomplete: move this shit out of here.
  bool clear_the_thing = false;

  bool auto_fit = false;
  bool replay = false;
  int  replay_multiplier = 0;
  float time = 0;

  Method_Spec* methods       = NULL;
  size_t       methods_count = 0;
};

void render_laba1(Memory_Arena* arena, The_Thing* thing) {
  if (thing->clear_the_thing) {
    thing->clear_the_thing = false;

    array_free(&result.grid);
    array_free(&result.time);
    array_free(&result.data);
    result.grid = {};
    result.data = {};
    result.time = {};
  }


  array<float> grid_to_be_drawn_this_frame;
  array<Vertex> data_to_be_drawn_this_frame;
  array<float> p_array;
  array<float> u_array;
  array<float> r_array;


  grid_to_be_drawn_this_frame.allocator = arena->allocator;
  data_to_be_drawn_this_frame.allocator = arena->allocator;
  p_array.allocator = arena->allocator;
  u_array.allocator = arena->allocator;
  r_array.allocator = arena->allocator;


  {
    Scoped_Lock mutex(&result.data_mutex);

    size_t n = thing->time / result.dt;
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

  static const ImVec2 plot_rect_size = ImVec2(300, 300);
  static const auto   plot_flags = thing->auto_fit ? ImPlotAxisFlags_AutoFit : 0;
  if (ImPlot::BeginPlot("Davlenie", plot_rect_size)) {
    ImPlot::SetupAxes("", "", plot_flags, plot_flags);
    ImPlot::PlotLine("p(x)", grid_to_be_drawn_this_frame.data, p_array.data, p_array.size);
    ImPlot::EndPlot();
  }

  ImGui::SameLine();

  if (ImPlot::BeginPlot("Velocity", plot_rect_size)) {
    ImPlot::SetupAxes("", "", plot_flags, plot_flags);
    ImPlot::PlotLine("u(x)", grid_to_be_drawn_this_frame.data, u_array.data, u_array.size);
    ImPlot::EndPlot();
  }

  ImGui::SameLine();

  if (ImPlot::BeginPlot("Ro", plot_rect_size)) {
    ImPlot::SetupAxes("", "", plot_flags, plot_flags);
    ImPlot::PlotLine("r(x)", grid_to_be_drawn_this_frame.data, r_array.data, r_array.size);
    ImPlot::EndPlot();
  }
}

void render_laba2(Memory_Arena* arena, The_Thing* thing) {
  if (thing->clear_the_thing) {
    thing->clear_the_thing = false;
  }
}

// Main code
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

  Method_Spec methods_laba1[] = {
    { nonconservative_lax_method, "Non Conservative Lax Method" },
    { euler_upwind_method,        "Euler Upwind Method"         },
    { conservative_lax_method,    "Conservative Lax Method"     },
  };

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
    { "dx",   &parameters_laba2.dx },
    { "dt",   &parameters_laba2.dt },
  };

  The_Thing things[2];
  things[0].thread_proc = methods_laba1[0].proc;
  things[0].thread_parameter = &parameters_laba1;
  things[0].settings = { input_parameters_laba1, array_size(input_parameters_laba1) }; 

  things[0].name = "Laba 1";
  things[0].method_name = methods_laba1[0].name;

  things[0].methods       = methods_laba1;
  things[0].methods_count = array_size(methods_laba1);

#if 0
#endif

  things[1].name = "Laba 2";
  things[1].thread_proc = lagrange_method;
  things[1].thread_parameter = &parameters_laba2;
  things[1].settings = { input_parameters_laba2, array_size(input_parameters_laba2) };

  { // init_program();
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
    parameters_laba2.dt   = 0.0001f;
    parameters_laba2.dx   = 0.01f;


    result.data_mutex = create_mutex();
    lagrange.data_mutex = create_mutex();
  }

  Memory_Arena temporary_storage;
  begin_memory_arena(&temporary_storage, KB(200));

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
          for (size_t i = 0; i < thing->settings.count; i++) {
            Input_Float_Settings* s = &thing->settings.inputs[i];

            static const float step      = 0.0f;
            static const float step_fast = 0.0f;
            static const char* format    = "%.4f";
            static const ImGuiInputTextFlags flags = ImGuiInputTextFlags_CharsScientific;
            ImGui::InputFloat(s->name, s->pointer, step, step_fast, format, flags);
          }

          ImGui::SliderFloat("Time", &thing->time, 0, result.t);
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
            thing->time = (thing->time <= result.t) ? thing->time : 0; // @Incomplete: result.t;
          }

          switch(i) {
          case 0: render_laba1(&temporary_storage, thing); break;
          case 1: render_laba2(&temporary_storage, thing); break;
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
