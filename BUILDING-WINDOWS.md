# Building CAFloodPro for Windows

The GitHub Actions workflow builds the Pro model in this repository as a
64-bit, double-precision OpenCL executable. It packages the executable, the
OpenCL loader when supplied by vcpkg, GPU configuration examples, source code,
dependency licences, and the CMake configuration used for the build.

## Run the online build

1. Open **Actions** and select **Build CAFloodPro for Windows**.
2. Select the branch containing these build fixes and click **Run workflow**.
3. Open a successful run and download **cafloodpro-windows-gpu** under
   **Artifacts**. Extract the archive before running the executable.

The workflow also runs automatically on pushes to `codex/windows-gpu-build`.

## Framework source

The older Exeter component archives omit the implementation of the raster
reader expected by the Pro source. The build therefore downloads the complete
CADDIES framework mirrored at [WaterDesk/caflood](https://github.com/WaterDesk/caflood),
whose README identifies the original University of Exeter framework sources.
The snapshot is pinned to commit
`ef408ee3b6126fadd85cb7c3a2bc73adfc3d8011` and checked against SHA256
`10936983b1440122217440a42513c86bfdad703c06090ca50ce06d82c3bbdad3`.

Only its framework is used. Its bundled flood and example applications are
removed from the temporary build tree and replaced with the Pro application
from this repository. The Pro hydraulic equations are unchanged.

The application explicitly depends on the framework's generated common CA
headers so that selecting the executable target also builds those headers.
MSVC definitions are set on the application target, preserving the compiler's
default flags, including C++ exception handling.

## Scope of validation

The workflow checks that the Windows executable starts successfully with
`/version`, `/help`, and `/test`. Those checks do not validate an OpenCL kernel
or a flood simulation. Test an example simulation on the intended computer
before relying on model outputs.

The initial build disables optional GDAL and NetCDF support and uses ASCII
grids. Running the model requires a compatible OpenCL device and its vendor
driver. Double precision requires device support for FP64. The Microsoft
Visual C++ x64 Redistributable may also be needed.

The workflow currently disables OpenCL events following the older NVIDIA
guidance. Review that setting for other OpenCL implementations.

## Local compilation

The same CMake configuration can be used on a Windows desktop with Visual
Studio 2022's **Desktop development with C++** workload, CMake 3, the complete
framework above, and the OpenCL headers and library. Place this repository's
contents under `Source/apps/cafloodpro`, remove the other application folders,
and configure from the framework's `Source` directory. The workflow records
all required configuration arguments; its build command is:

```powershell
cmake --build BUILD_DIRECTORY --config Release --target cafloodpro_GPU_64
```

The `_64` suffix describes the executable architecture; numerical precision is
selected separately using `CAAPI_REAL_PRECISION=double`.
