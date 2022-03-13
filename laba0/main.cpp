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

#define max(x, y) ( ((x) > (y)) ? (x) : (y) )

struct Thread_Data {
  float* ts = NULL;
  float* xs = NULL;
  float* ys = NULL; 

  float w1 = 1;
  float w2 = 1;
  float d1 = 0;
  float d2 = 0;
  unsigned count = 0;

  float dt = 0;
};


static bool show_imgui_demo  = false;
static bool show_implot_demo = false;

static bool show_laba1 = true;
static bool show_laba2 = false;
static bool show_laba3 = false;
static bool show_laba4 = false;

static bool auto_cursor     = true;
static bool processing_data = false;
static bool want_user_confirmation_for_reset = false;
static bool done            = false;


static const int NUMBER_OF_POINTS = 100000;
static float ts[NUMBER_OF_POINTS];
static float xs[NUMBER_OF_POINTS];
static float ys[NUMBER_OF_POINTS];

static float o1 = 1;
static float o2 = 1;
static float b1 = 0;
static float b2 = 0;

static Thread_Data thread_data;

static int cursor = 0;

static HANDLE computation_thread = NULL;

static bool already_computed_some_data() { return thread_data.count != 0; }

static DWORD computation_thread_proc(LPVOID parameter) {
  Thread_Data* data = (Thread_Data*) parameter;

  float w1 = data->w1;
  float w2 = data->w2;
  float d1 = data->d1;
  float d2 = data->d2;
  float dt = data->dt;

  float t = 0;
  while (data->count < NUMBER_OF_POINTS) {
    float x = sin(w1 * t + d1);
    float y = sin(w2 * t + d2);
    unsigned count = data->count;

    {
      InterlockedExchange((unsigned*) &data->ts[count], *(unsigned*) &t);
      InterlockedExchange((unsigned*) &data->xs[count], *(unsigned*) &x);
      InterlockedExchange((unsigned*) &data->ys[count], *(unsigned*) &y);
      InterlockedIncrement(&data->count);
    }

    t += dt;

    Sleep(20); // @RemoveMe: 
  }

  return 0;
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


  thread_data.ts = ts;
  thread_data.xs = xs;
  thread_data.ys = ys;
  thread_data.dt = io.DeltaTime;


  ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);
  while (!done) {
    MSG msg;
    while (::PeekMessage(&msg, NULL, 0U, 0U, PM_REMOVE)) {
      ::TranslateMessage(&msg);
      ::DispatchMessage(&msg);
      if (msg.message == WM_QUIT) done = true;
    }
    if (done) { break; }
        

    // Start the Dear ImGui frame
    ImGui_ImplDX11_NewFrame();
    ImGui_ImplWin32_NewFrame();
    ImGui::NewFrame();


  
    if (auto_cursor) cursor = max(cursor, thread_data.count);

    {
      ImGui::Begin("Control window");

      ImGui::Checkbox("ImGui Demo", &show_imgui_demo);
      ImGui::Checkbox("ImPlot Demo", &show_implot_demo);

      if (ImGui::CollapsingHeader("Laba 1")) {
        ImGui::Checkbox("Show", &show_laba1);

        static const float step      = 0.0f;
        static const float step_fast = 0.0f;
        static const char* format    = "%.4f";
        static const ImGuiInputTextFlags flags = ImGuiInputTextFlags_CharsScientific;

        ImGui::InputFloat("omega 1", &o1, step, step_fast, format, flags);
        ImGui::InputFloat("delta 1", &b1, step, step_fast, format, flags);
        ImGui::InputFloat("omega 2", &o2, step, step_fast, format, flags);
        ImGui::InputFloat("delta 2", &b2, step, step_fast, format, flags);

        ImGui::SliderInt("Time Domain: ", &cursor, 0, thread_data.count);
        ImGui::Checkbox("Auto Cursor", &auto_cursor);
      }

      if (ImGui::CollapsingHeader("Laba 2")) {}
      if (ImGui::CollapsingHeader("Laba 3")) {}
      if (ImGui::CollapsingHeader("Laba 4")) {}

      if (ImGui::Button("Start")) {
        // If there is no data computed, we want to reset everything automatically without pressing 'Reset' button.
        if (!already_computed_some_data()) {
          thread_data.w1 = o1;
          thread_data.w2 = o2;
          thread_data.d1 = b1;
          thread_data.d2 = b2;
        }

        if (!processing_data) {
          if (already_computed_some_data()) {
            // resume thread execution.
            ResumeThread(computation_thread);

          } else {
            // destroy current thread if there is one.
            // and start a new thread.
            if (computation_thread) TerminateThread(computation_thread, 0);
            computation_thread = CreateThread(NULL, 0, computation_thread_proc, &thread_data, 0, NULL);
          }
        }
        processing_data = true;
      }

      ImGui::SameLine();

      if (ImGui::Button("Reset")) {
        want_user_confirmation_for_reset = true;
      }

      ImGui::SameLine();

      if (ImGui::Button("Stop")) {
        // suspend thread execution.
        SuspendThread(computation_thread);
        processing_data = false;
      }

      auto framerate = io.Framerate;
      ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f/framerate, framerate);
      ImGui::End();
    }

    if (show_imgui_demo) {
      ImGui::ShowDemoWindow(&show_imgui_demo);
    }

    if (show_implot_demo) {
      ImPlot::ShowDemoWindow();
    }

    if (show_laba1) {
      if (ImPlot::BeginPlot("Sin")) {
        ImPlot::SetupAxisLimits(ImAxis_X1,  0.0, cursor * thread_data.dt, ImGuiCond_Always);
        ImPlot::SetupAxisLimits(ImAxis_Y1, -1.0, 1.0,                     ImGuiCond_Always);

        ImPlot::PlotLine("Data 1", thread_data.ts, thread_data.xs, cursor);
        ImPlot::PlotLine("Data 2", thread_data.ts, thread_data.ys, cursor);
        ImPlot::EndPlot();
      }

      if (ImPlot::BeginPlot("Phase Portrait")) {
        ImPlot::SetupAxisLimits(ImAxis_X1, -1.0, 1.0, ImGuiCond_Always);
        ImPlot::SetupAxisLimits(ImAxis_Y1, -1.0, 1.0, ImGuiCond_Always);
        
        ImPlot::PlotLine("Data", thread_data.xs, thread_data.ys, cursor);
        ImPlot::EndPlot();
      }
    }

    if (want_user_confirmation_for_reset) {
      if (already_computed_some_data()) {
        ImGui::Begin("Are you sure want to reset everything?");
        if (ImGui::Button("Reset")) {


          if (processing_data) {
            if (computation_thread) TerminateThread(computation_thread, 0);
          }

          thread_data.w1 = o1;
          thread_data.w2 = o2;
          thread_data.d1 = b1;
          thread_data.d2 = b2;
          thread_data.count = 0;

          if (processing_data) {
            computation_thread = CreateThread(NULL, 0, computation_thread_proc, &thread_data, 0, NULL);
          }

          cursor = 0;
          want_user_confirmation_for_reset = false;
        }

        ImGui::SameLine();

        if (ImGui::Button("Cancel")) {
          want_user_confirmation_for_reset = false;
        }
        ImGui::End();
      } else {

        if (processing_data) {
          if (computation_thread) TerminateThread(computation_thread, 0);
        }

        thread_data.w1 = o1;
        thread_data.w2 = o2;
        thread_data.d1 = b1;
        thread_data.d2 = b2;
        thread_data.count = 0;

        if (processing_data) {
          computation_thread = CreateThread(NULL, 0, computation_thread_proc, &thread_data, 0, NULL);
        }

        cursor = 0;
        want_user_confirmation_for_reset = false;
      }
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
