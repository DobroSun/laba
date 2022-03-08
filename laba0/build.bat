@REM Build for Visual Studio compiler. Run your copy of vcvars32.bat or vcvarsall.bat to setup command-line compiler.
@set OUT_DIR=Binaries
@set OUT_EXE=main.exe
@set INCLUDES=/I..\imgui /I..\imgui\backends /I..\implot /I "%WindowsSdkDir%Include\um" /I "%WindowsSdkDir%Include\shared" /I "%DXSDK_DIR%Include"
@set SOURCES=main.cpp 
@set LIBS=/LIBPATH:"%DXSDK_DIR%/Lib/x86" d3d11.lib d3dcompiler.lib
mkdir %OUT_DIR%
cl /nologo /Zi /MD %INCLUDES% /D UNICODE /D _UNICODE %SOURCES% /Fe%OUT_DIR%/%OUT_EXE%.exe /Fo%OUT_DIR%/ /link %LIBS%


@REM // @set SOURCES = ..\imgui\backends\imgui_impl_dx11.cpp ..\imgui\backends\imgui_impl_win32.cpp ..\imgui\imgui*.cpp ..\implot\implot*.cpp
