mkdir build-release

cd build-release

IF "%ARCHITECTURE%" == "x64" (
  set ARCH=x64
) else (
  set ARCH=Win32
)

IF "%BUILD_PLATFORM%" == "VS2019" (
    set LIBPATH=E:\libs\VS2019
    set GENERATOR=Visual Studio 16
    set VS_PATH="C:\Program Files (x86)\Microsoft Visual Studio\2019\Professional\Common7\IDE\devenv.com"
)

ECHO "==============================================================="
ECHO "==============================================================="
ECHO "Building with :"
ECHO "ARCHITECTURE        : %ARCHITECTURE%"
ECHO "BUILD_PLATFORM      : %BUILD_PLATFORM%"
ECHO "GENERATOR           : %GENERATOR%"
ECHO "CMAKE_CONFIGURATION : %CMAKE_CONFIGURATION%"
ECHO "VS_PATH             : %VS_PATH%"
ECHO "LIBPATH             : %LIBPATH%"
ECHO "==============================================================="
ECHO "==============================================================="

"C:\Program Files\CMake\bin\cmake.exe" -G "%GENERATOR%"  -A %ARCH% ..

%VS_PATH% /Build "Release" OpenVolumeMesh.sln /Project "ALL_BUILD"

IF %errorlevel% NEQ 0 exit /b %errorlevel%

"C:\Program Files\CMake\bin\ctest.exe" -C Release .

IF %errorlevel% NEQ 0 exit /b %errorlevel%

cd ..

mkdir build-debug

cd build-debug

"C:\Program Files\CMake\bin\cmake.exe" -G "%GENERATOR%" -A %ARCH% ..

%VS_PATH% /Build "Debug" OpenVolumeMesh.sln /Project "ALL_BUILD"

IF %errorlevel% NEQ 0 exit /b %errorlevel%

"C:\Program Files\CMake\bin\ctest.exe" -C Debug .

IF %errorlevel% NEQ 0 exit /b %errorlevel%
