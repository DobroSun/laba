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

#include <cmath>
#include <cstdint>

#include "array.cpp"
#include "threads_api.cpp"
#include "threads_windows.cpp"

#pragma comment(lib, "d3d11.lib")
#pragma comment(lib, "d3dcompiler.lib")

#define clamp(w, mi, ma) ( min(max((w), (mi)), (ma)) )

typedef uint32_t uint;
typedef uint32_t uint32;
typedef uint64_t uint64;

typedef  int32_t  int32;
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

  float dt = 0.0f;
  float dx = 0.0f;
};

struct Thread_Data1 {
  Parameters_Laba1* parameters;
             float* global_t;
};

struct Thread_Data2 {
  Parameters_Laba2* parameters;
             float* global_t;
};

struct Input_Float_Settings {
  const char* name;
  float* pointer;
};

struct Input_Parameters_Settings {
  Input_Float_Settings* inputs;
  size_t count;
};

typedef int (*Thread_Proc)(void*);

struct The_Thing {
  const char* name;

  Thread      thread;
  Thread_Proc thread_proc;
  void*       thread_param = NULL;

  Input_Parameters_Settings input_parameters_settings = {};

  float global_t = 0.0f; // this is going to be passed in thread_param.
  float target_t = 0.0f;

  bool use_auto_fit = false;
  bool thread_is_paused = false; // @Incomplete: refactor that into threads_api
  bool should_clear_solution = false;
  bool in_replay    = false;

  int32 replay_modifier = 0;
};

struct Vertex {
  double t;
  double r;
  double u;
  double p;
};

struct Voxel {
  double x_euler; // euler coordinate, will recompute them.
  double p;
  double r;
  double u;
};

struct Result_Laba2 {
  Mutex mutex;
  array<double> grid_lagrange; // lagrange coordinate, those are const.
  array<Voxel> data;
};

#if 0
struct Laba_Results {
  Mutex mutex;

  double t;
  double dt;

  array<double> grid;
  union {
    array<Vertex> data1;
    array<Voxel>  data2;
  };
};
#endif

struct Result {
  Mutex mutex;
  double dt;
  double  t;
  array<double> grid; // for now we only have 1 dimension.
  array<Vertex> data; // this is going to evolve in time. size := NUMBER_OF_LAYERS * grid.size()
};

Vertex get_data(Result* r, size_t t, size_t j) {
  return r->data[t * r->grid.size + j];
}

static Result      result = {};
static Result_Laba2 laba2 = {};

static int computation_thread_proc(void* param) {
  //Thread_Data1* data = (Thread_Data1*) param;

  Parameters_Laba1 inputs = *(Parameters_Laba1*) param; // copy
  //float* global_t = data->global_t;


  float xmin = inputs.xmin; 
  float xmax = inputs.xmax;
  float p0   = inputs.p0;
  float x0   = inputs.x0;
  float r0   = inputs.r0;
  float ro0  = inputs.ro0;
  float gamma = inputs.gamma;
  float dt_f = inputs.dt;
  float dx   = inputs.dx;

  size_t NUMBER_OF_POINTS_IN_X = (xmax - xmin) / dx;

  double dt = dt_f;
  double t = 0.0; // real time.

  {
    // 
    // create a grid.
    //
    double* grid = (double*) alloca(sizeof(double) * NUMBER_OF_POINTS_IN_X);
    for (size_t i = 0; i < NUMBER_OF_POINTS_IN_X; i++) {
      grid[i] = xmin + i*dx;
    }

    {
      Scoped_Lock lock(&result.mutex);
      result.grid = array_copy(grid, NUMBER_OF_POINTS_IN_X);
    }
    InterlockedExchange64((int64*) &result.dt, *(int64*) &dt);
    //InterlockedExchange((uint*) global_t, *(uint*) &dt_f);
  }
  {
    // 
    // apply initial conditions on t=0
    //
    Vertex* data = (Vertex*) alloca(sizeof(Vertex) * NUMBER_OF_POINTS_IN_X);
    for (size_t j = 0; j < NUMBER_OF_POINTS_IN_X; j++)  {
      double x = result.grid[j];

      Vertex* v = &data[j];
      v->t = t;
      v->r = ro_t0(x, ro0);
      v->u =  u_t0(x);
      v->p =  p_t0(x, p0, x0, r0);
    }

    {
      Scoped_Lock lock(&result.mutex);
      array_add(&result.data, data, NUMBER_OF_POINTS_IN_X);
    }
  }

  Vertex* solution = (Vertex*) alloca(sizeof(Vertex) * NUMBER_OF_POINTS_IN_X);

  size_t i = 0;
  while (1) {
    i++;
    t += dt;

    for (size_t j = 1; j < NUMBER_OF_POINTS_IN_X-1; j++) {
      // 
      // internal points.
      //
      Vertex curr = get_data(&result, i-1, j);
      Vertex prev = get_data(&result, i-1, j-1);
      Vertex next = get_data(&result, i-1, j+1);

      double x = result.grid[j];

      Vertex* v = &solution[j];
      v->t = t;
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
      v->t = t;
      v->r = solution[j+1].r;
      v->u = u_x0(t);
      v->p = solution[j+1].p;
    }
    {
      j = NUMBER_OF_POINTS_IN_X-1;
      Vertex* v = &solution[j];
      v->t = t;
      v->r = solution[j-1].r;
      v->u = u_x0(t);
      v->p = solution[j-1].p;
    }

    {
      Scoped_Lock lock(&result.mutex);
      array_add(&result.data, solution, NUMBER_OF_POINTS_IN_X);
    }
    InterlockedExchange64((int64*) &result.t, *(int64*) &t);

    float t_f = t;
    //InterlockedExchange((uint*) global_t, *(uint*) &t_f);
  }
end:
  return 0;
}


static int computation_thread_proc_laba2(void* param) {
  Parameters_Laba2 data = *(Parameters_Laba2*) param; // copy

  float xmin = data.xmin; 
  float xmax = data.xmax;
  float dt_f = data.dt;
  float dx   = data.dx;

  size_t NUMBER_OF_POINTS_IN_X = (xmax - xmin) / dx;

  double dt = dt_f;
  puts("Yey!");


end:
  return 0;
}

void render_laba2(ImGuiIO& io, The_Thing* thing) {}

void render_laba1(ImGuiIO& io, The_Thing* thing) {
  auto framerate = io.Framerate;

  array<double> grid_to_be_drawn_this_frame;
  array<Vertex> data_to_be_drawn_this_frame;

  {
    Scoped_Lock mutex(&result.mutex);
    size_t n = thing->target_t / thing->global_t;
    size_t needed_data_start =  n    * result.grid.size;
    size_t needed_data_end   = (n+1) * result.grid.size;

    array_copy_range(&data_to_be_drawn_this_frame, &result.data, needed_data_start, needed_data_end);
    array_copy      (&grid_to_be_drawn_this_frame, &result.grid);
  }

  if (thing->in_replay) {
    double m = pow(2.0, thing->replay_modifier);
    thing->target_t += m * 1/framerate;
    thing->target_t  = (thing->target_t <= thing->global_t) ? thing->target_t : 0;
  }

  size_t data_size = data_to_be_drawn_this_frame.size;

  array<double> p_array;
  array<double> u_array;
  array<double> r_array;

  array_resize(&p_array, data_size);
  array_resize(&u_array, data_size);
  array_resize(&r_array, data_size);

  for (size_t i = 0; i < data_size; i++) {
    p_array[i] = data_to_be_drawn_this_frame[i].p;
    u_array[i] = data_to_be_drawn_this_frame[i].u;
    r_array[i] = data_to_be_drawn_this_frame[i].r;
  }

  static const ImVec2 plot_rect_size = ImVec2(300, 300);
  static const auto   plot_flags = thing->use_auto_fit ? ImPlotAxisFlags_AutoFit : 0;
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



  init_threads_api();
  check_threads_api();

  bool show_imgui_demo  = false;
  bool show_implot_demo = false;
  bool show_laba1       = true;
  bool show_laba2       = false;
  bool show_laba3       = false;
  bool done             = false;
  bool thread_is_paused = false;

  bool auto_fit = false;
  bool replay = false;
  int  replay_multiplier = 0;
  float time = 0;


  Parameters_Laba1 parameters_laba1 = {};
  Parameters_Laba2 parameters_laba2 = {};

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
    { "dx",   &parameters_laba2.dx   },
    { "dt",   &parameters_laba2.dt   },
  };

  The_Thing the_things[2];
#if 0
  Thread_Data1 thread_data1 = {};
  thread_data1.parameters = &parameters_laba1;
  thread_data1.global_t   = &the_things[0].global_t;

  Thread_Data2 thread_data2 = {};
  thread_data2.parameters = &parameters_laba2;
  thread_data2.global_t   = &the_things[1].global_t;

  {
    the_things[0].name = "Laba 1";
    the_things[0].thread_proc  = computation_thread_proc;
    the_things[0].thread_param = &thread_data1;
    the_things[0].input_parameters_settings = { input_parameters_laba1, array_size(input_parameters_laba1) };

    the_things[1].name = "Laba 2";
    the_things[1].thread_proc  = computation_thread_proc_laba2;
    the_things[1].thread_param = &thread_data2;
    the_things[1].input_parameters_settings = { input_parameters_laba2, array_size(input_parameters_laba2) };
  }
#endif

  { // default parameters.
    parameters_laba1.xmin = -3.0f;
    parameters_laba1.xmax = 3.0f;
    parameters_laba1.p0   = 1.0f;
    parameters_laba1.x0   = 0.0f;
    parameters_laba1.r0   = 0.5f;
    parameters_laba1.ro0  = 1.0f;
    parameters_laba1.gamma = 1.67f;
    parameters_laba1.dt   = 0.0001f;
    parameters_laba1.dx   = 0.01f;


    parameters_laba2.dt = 0.0001f;
    parameters_laba2.dx = 0.01f;

    result.mutex = create_mutex();
    laba2.mutex  = create_mutex();
  }

  const ImVec4 clear_color               = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);
  const float  clear_color_with_alpha[4] = { clear_color.x * clear_color.w, clear_color.y * clear_color.w, clear_color.z * clear_color.w, clear_color.w };
  while (!done) {
    MSG msg;
    while (::PeekMessage(&msg, NULL, 0U, 0U, PM_REMOVE)) {
      ::TranslateMessage(&msg);
      ::DispatchMessage(&msg);
      if (msg.message == WM_QUIT) done = true;
    }
    if (done) { break; }
        
    auto framerate = io.Framerate;


    // Start the Dear ImGui frame
    ImGui_ImplDX11_NewFrame();
    ImGui_ImplWin32_NewFrame();
    ImGui::NewFrame();

    {
      ImGui::Begin("Control window");
      ImGui::Checkbox("ImGui Demo",  &show_imgui_demo);
      ImGui::Checkbox("ImPlot Demo", &show_implot_demo);

#if 0
      for (size_t i = 0; i < array_size(the_things); i++) {
        The_Thing* thing = &the_things[i];
        const char* name = thing->name;

        Thread*     thread       = &thing->thread;
        Thread_Proc thread_proc  = thing->thread_proc;
        void*       thread_param = thing->thread_param;

        Input_Parameters_Settings* settings = &thing->input_parameters_settings;

        const float*  global_t = &thing->global_t;
              float*  target_t = &thing->target_t;
        bool*     use_auto_fit = &thing->use_auto_fit;
        bool* thread_is_paused = &thing->thread_is_paused;
        bool* should_clear_solution = &thing->should_clear_solution;
        bool* in_replay        = &thing->in_replay;
        int32* replay_modifier = &thing->replay_modifier;

        bool thread_is_running = is_thread_running(thread);

        if (ImGui::CollapsingHeader(name)) {
          for (size_t k = 0; k < settings->count; k++) {
            Input_Float_Settings* input = &settings->inputs[k];

            static const float step      = 0.0f;
            static const float step_fast = 0.0f;
            static const char* format    = "%.4f";
            static const ImGuiInputTextFlags flags = ImGuiInputTextFlags_CharsScientific;
            ImGui::InputFloat(input->name, input->pointer, step, step_fast, format, flags);
          }

          ImGui::SliderFloat("Time", target_t, 0, thing->global_t);
          ImGui::Checkbox("Auto Fit", use_auto_fit);

          const char* start_or_continue = *thread_is_paused ? "Continue" : "Start";
          if (ImGui::Button(start_or_continue)) {
            if (*thread_is_paused) {
              resume_thread(thread);
              *thread_is_paused = false;
            } else {
              if (!thread_is_running) {
                start_thread(thread, thread_proc, thread_param);
              }
            }
          }

          ImGui::SameLine();

          if (ImGui::Button("Pause")) {
            if (!*thread_is_paused && thread_is_running) {
              suspend_thread(thread);
              *thread_is_paused = true;
            }
          }

          ImGui::SameLine();

          if (ImGui::Button("Stop")) {
            if (thread_is_running) {
              kill_thread(thread);
              *thread_is_paused = false;
            }
          }

          ImGui::SameLine();

          if (ImGui::Button("Clear")) {
            if (thread_is_running) {
              kill_thread(thread);
              *thread_is_paused = false;
            }
            *should_clear_solution = true;
          }

          if (ImGui::Button("Start Playing")) {
            *in_replay = true;
          }

          ImGui::SameLine();

          if (ImGui::Button("x2")) {
            *replay_modifier += 1;
          }

          ImGui::SameLine();

          if (ImGui::Button("/2")) {
            *replay_modifier -= 1;
          }

          ImGui::SameLine();

          ImGui::Text("%s%g", (*replay_modifier >= 0) ? "x" : "/", pow(2.0, abs(*replay_modifier)));
          if (ImGui::Button("Stop Playing")) {
            *in_replay = false;
          }

          switch (i) {
          case 0: render_laba1(io, thing); break;
          case 1: render_laba2(io, thing); break;
          }
        }
      }
#else
      if (ImGui::CollapsingHeader("Laba 1")) {
        ImGui::Checkbox("Show", &show_laba1);

        for (size_t i = 0; i < array_size(input_parameters_laba1); i++) {
          Input_Float_Settings s = input_parameters_laba1[i];

          static const float step      = 0.0f;
          static const float step_fast = 0.0f;
          static const char* format    = "%.4f";
          static const ImGuiInputTextFlags flags = ImGuiInputTextFlags_CharsScientific;
          ImGui::InputFloat(s.name, s.pointer, step, step_fast, format, flags);
        }

        ImGui::SliderFloat("Time", &time, 0, result.t);
        ImGui::Checkbox("Auto Fit", &auto_fit);

        const char* start_or_continue = thread_is_paused ? "Continue" : "Start";
        if (ImGui::Button(start_or_continue)) {
          if (thread_is_paused) {
            resume_thread(&the_things[0].thread);
            thread_is_paused = false;
          } else {
            if (!is_thread_running(&the_things[0].thread)) {
              start_thread(&the_things[0].thread, computation_thread_proc, &parameters_laba1);
            }
          }
        }

        ImGui::SameLine();

        if (ImGui::Button("Pause")) {
          if (!thread_is_paused && is_thread_running(&the_things[0].thread)) {
            suspend_thread(&the_things[0].thread);
            thread_is_paused = true;
          }
        }

        ImGui::SameLine();

        if (ImGui::Button("Stop")) {
          if (is_thread_running(&the_things[0].thread)) {
            kill_thread(&the_things[0].thread);
            thread_is_paused = false;
          }
        }

        ImGui::SameLine();

        if (ImGui::Button("Clear")) {
          if (is_thread_running(&the_things[0].thread)) { // just kill a thread.
            kill_thread(&the_things[0].thread);
            thread_is_paused = false;
          }
          result = {}; // @MemoryLeak: 
        }

        if (ImGui::Button("Start Playing")) {
          replay = true;
        }

        ImGui::SameLine();

        if (ImGui::Button("x2")) {
          replay_multiplier += 1;
        }

        ImGui::SameLine();

        if (ImGui::Button("/2")) {
          replay_multiplier -= 1;
        }

        ImGui::SameLine();

        ImGui::Text("%s%g", (replay_multiplier >= 0) ? "x" : "/", pow(2.0, abs(replay_multiplier)));
        if (ImGui::Button("Stop Playing")) {
          replay = false;
        }

        if (show_laba1) {
          array<double> grid_to_be_drawn_this_frame;
          array<Vertex> data_to_be_drawn_this_frame;
          {
            Scoped_Lock mutex(&result.mutex);

            size_t n = time / result.dt;
            size_t needed_data_start =  n    * result.grid.size;
            size_t needed_data_end   = (n+1) * result.grid.size;

            array_copy_range(&data_to_be_drawn_this_frame, &result.data, needed_data_start, needed_data_end);
            array_copy      (&grid_to_be_drawn_this_frame, &result.grid);
          }

          if (replay) {
            double m = pow(2.0, replay_multiplier);
            time += m * 1/framerate;
            time = (time <= result.t) ? time : 0;
          }

          size_t data_size = data_to_be_drawn_this_frame.size;

          array<double> p_array;
          array<double> u_array;
          array<double> r_array;

          array_resize(&p_array, data_size);
          array_resize(&u_array, data_size);
          array_resize(&r_array, data_size);

          for (size_t i = 0; i < data_size; i++) {
            p_array[i] = data_to_be_drawn_this_frame[i].p;
            u_array[i] = data_to_be_drawn_this_frame[i].u;
            r_array[i] = data_to_be_drawn_this_frame[i].r;
          }

          static const ImVec2 plot_rect_size = ImVec2(300, 300);
          static const auto   plot_flags = auto_fit ? ImPlotAxisFlags_AutoFit : 0;
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
      }

      if (ImGui::CollapsingHeader("Laba 2")) {

      }

      if (ImGui::CollapsingHeader("Laba 3")) {}
      if (ImGui::CollapsingHeader("Laba 4")) {}
#endif

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
    g_pd3dDeviceContext->OMSetRenderTargets(1, &g_mainRenderTargetView, NULL);
    g_pd3dDeviceContext->ClearRenderTargetView(g_mainRenderTargetView, clear_color_with_alpha);
    ImGui_ImplDX11_RenderDrawData(ImGui::GetDrawData());
    g_pSwapChain->Present(1, 0); // Present with vsync
    // g_pSwapChain->Present(0, 0); // Present without vsync
  }

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
