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

  bool show_laba1 = false;
  bool show_laba2 = false;
  bool show_laba3 = false;
  bool show_laba4 = false;

  bool start_processing_data = false;
  bool want_user_confirmation_for_reset = false;
  bool done = false;


  static float dt = 0;
  static float ts[100000];
  static float xs[100000];
  static float ys[100000];
  static int count = 0;

  static float w1 = 1;
  static float w2 = 1;
  static float d1 = 0;
  static float d2 = 0;

  float o1 = w1;
  float o2 = w2;
  float b1 = d1;
  float b2 = d2;




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
      }

      if (ImGui::CollapsingHeader("Laba 2")) {}
      if (ImGui::CollapsingHeader("Laba 3")) {}
      if (ImGui::CollapsingHeader("Laba 4")) {}

      if (ImGui::Button("Start")) {
        // If there is no data computed, we want to reset everything automatically without pressing 'Reset' button.
        if (!count) {
          w1 = o1;
          w2 = o2;
          d1 = b1;
          d2 = b2;
        }
        start_processing_data = true;
      }

      ImGui::SameLine();

      if (ImGui::Button("Reset")) {
        want_user_confirmation_for_reset = true;
      }

      ImGui::SameLine();

      if (ImGui::Button("Stop")) {
        start_processing_data = false;
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

    if (start_processing_data) {
      dt += io.DeltaTime;
      ts[count] = dt;
      xs[count] = sin(w1 * dt + d1);
      ys[count] = sin(w2 * dt + d2);
      count++;
    }

    if (show_laba1) {
      if (ImPlot::BeginPlot("Sin")) {
        ImPlot::SetupAxisLimits(ImAxis_X1,  0.0,  dt, ImGuiCond_Always);
        ImPlot::SetupAxisLimits(ImAxis_Y1, -1.0, 1.0, ImGuiCond_Always);

        ImPlot::PlotLine("Data 1", ts, xs, count);
        ImPlot::PlotLine("Data 2", ts, ys, count);
        ImPlot::EndPlot();
      }

      if (ImPlot::BeginPlot("Phase Portrait")) {
        ImPlot::SetupAxisLimits(ImAxis_X1, -1.0, 1.0, ImGuiCond_Always);
        ImPlot::SetupAxisLimits(ImAxis_Y1, -1.0, 1.0, ImGuiCond_Always);
        
        ImPlot::PlotLine("Data", xs, ys, count);
        ImPlot::EndPlot();
      }
    }

    if (want_user_confirmation_for_reset) {

      // If there is some data computed, we want to ask for a permission to reset everything, otherwise just reset everything.
      if (count) {
        ImGui::Begin("Are you sure want to reset everything?");
        if (ImGui::Button("Reset")) {
          w1 = o1;
          w2 = o2;
          d1 = b1;
          d2 = b2;

          dt    = 0;
          count = 0;

          want_user_confirmation_for_reset = false;
        }

        ImGui::SameLine();

        if (ImGui::Button("Cancel")) {
          want_user_confirmation_for_reset = false;
        }
        ImGui::End();
      } else {
        w1 = o1;
        w2 = o2;
        d1 = b1;
        d2 = b2;

        dt    = 0;
        count = 0;
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
