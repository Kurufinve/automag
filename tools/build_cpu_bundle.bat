@echo off
rem Windows counterpart of tools/build_cpu_bundle.sh.
rem
rem The 5.7 GB bundle is big because the CUDA build of torch brings ~3 GB of GPU
rem libraries.  A CPU wheel of the SAME torch version gives identical numbers and
rem takes a fraction of that, so the bundle stays small once the checkpoints are
rem kept outside.  Run this from the "Anaconda Prompt" of the automag environment
rem (Python 3.12): the prophet packages are copied from that environment, because
rem prophet-mlip is distributed as a wheel inside the Prophet repository and not
rem from PyPI.
rem
rem Usage:
rem   tools\build_cpu_bundle.bat                 both bundles (checkpoints inside + outside)
rem   tools\build_cpu_bundle.bat --no-models     only the small one (checkpoints outside)
rem   tools\build_cpu_bundle.bat --setup-only    just prepare the CPU environment
rem   tools\build_cpu_bundle.bat --help
rem
rem Environment:
rem   AUTOMAG_CPU_ENV      venv to build from (default %LOCALAPPDATA%\automag\cpuenv)
rem   AUTOMAG_CPU_PYTHON   interpreter used to create that venv (default python)
rem   AUTOMAG_TORCH_CPU    torch version to install (default 2.14.1)
rem   AUTOMAG_SOURCE_ENV   environment to copy prophet/prophet_spin from
rem                        (default: the active conda environment)
setlocal EnableDelayedExpansion

set "ROOT=%~dp0.."
pushd "%ROOT%" || (echo ERROR: cannot enter "%ROOT%" & exit /b 1)

set "BUILD_EMBEDDED=1"
set "SETUP_ONLY=0"

:parse
if "%~1"=="" goto parsed
if /i "%~1"=="--no-models"  set "BUILD_EMBEDDED=0" & shift & goto parse
if /i "%~1"=="--setup-only" set "SETUP_ONLY=1"     & shift & goto parse
if /i "%~1"=="--help"       goto usage
if /i "%~1"=="-h"           goto usage
echo unknown option: %~1
goto usage_error

:parsed
if "%AUTOMAG_CPU_ENV%"==""    set "AUTOMAG_CPU_ENV=%LOCALAPPDATA%\automag\cpuenv"
if "%AUTOMAG_TORCH_CPU%"==""  set "AUTOMAG_TORCH_CPU=2.14.1"
if "%AUTOMAG_CPU_PYTHON%"=="" set "AUTOMAG_CPU_PYTHON=python"
set "ENV_DIR=%AUTOMAG_CPU_ENV%"
set "VENV_PY=%ENV_DIR%\Scripts\python.exe"

rem --- the CPU build must match the interpreter of the reference environment ---
"%AUTOMAG_CPU_PYTHON%" -c "import sys;raise SystemExit(0 if sys.version_info[:2]==(3,12) else 1)"
if errorlevel 1 (
    echo ERROR: this script expects Python 3.12.
    echo        Activate the automag environment first ^(conda activate automag^),
    echo        or set AUTOMAG_CPU_PYTHON to a 3.12 interpreter.
    popd & exit /b 1
)

if not exist "%VENV_PY%" (
    echo === creating the CPU build environment: %ENV_DIR%
    "%AUTOMAG_CPU_PYTHON%" -m venv "%ENV_DIR%"
    if errorlevel 1 (echo ERROR: venv creation failed & popd & exit /b 1)
)

rem The wheel name depends on the ABI tag, so read it from the venv interpreter
rem through a temp file (a for /f loop would have to survive cmd's % and ' rules).
"%VENV_PY%" -c "import sys;sys.stdout.write('cp'+str(sys.version_info[0])+str(sys.version_info[1]))" > "%TEMP%\automag_tag.txt"
set /p TAG=<"%TEMP%\automag_tag.txt"
del "%TEMP%\automag_tag.txt" 2>nul
echo === build environment: %ENV_DIR% (!TAG!)

"%VENV_PY%" -m pip install -q --upgrade pip

echo === installing torch %AUTOMAG_TORCH_CPU%+cpu
"%VENV_PY%" -m pip install -q --force-reinstall "torch==%AUTOMAG_TORCH_CPU%+cpu" --index-url https://download.pytorch.org/whl/cpu
if errorlevel 1 (
    echo     the CPU index has no %AUTOMAG_TORCH_CPU%+cpu for !TAG!, downloading the wheel
    set "WHL_DIR=%ENV_DIR%\wheels"
    if not exist "!WHL_DIR!" mkdir "!WHL_DIR!"
    set "WHL=!WHL_DIR!\torch-%AUTOMAG_TORCH_CPU%+cpu-!TAG!-!TAG!-win_amd64.whl"
    if not exist "!WHL!" (
        curl -L --fail --max-time 900 -o "!WHL!" "https://download-r2.pytorch.org/whl/cpu/torch-%AUTOMAG_TORCH_CPU%%2Bcpu-!TAG!-!TAG!-win_amd64.whl"
        if errorlevel 1 (
            echo ERROR: could not download the CPU wheel.
            echo        Open https://download.pytorch.org/whl/cpu/torch/ , fetch
            echo        torch-%AUTOMAG_TORCH_CPU%+cpu-!TAG!-!TAG!-win_amd64.whl into
            echo        !WHL_DIR! and run this script again.
            popd & exit /b 1
        )
    )
    "%VENV_PY%" -m pip install -q --force-reinstall --no-deps "!WHL!"
    if errorlevel 1 (echo ERROR: installing the CPU wheel failed & popd & exit /b 1)
)

echo === installing the project requirements and the build tools
set "REQ=requirements.txt"
if exist "requirements-windows.txt" set "REQ=requirements-windows.txt"
"%VENV_PY%" -m pip install -q -r "!REQ!" pyinstaller ninja e3nn matscipy torch_geometric
if errorlevel 1 (echo ERROR: pip install failed & popd & exit /b 1)

rem --- prophet-mlip is not on PyPI: reuse the copy of the reference environment ---
if not "%AUTOMAG_SOURCE_ENV%"=="" (
    set "REF_SITE=%AUTOMAG_SOURCE_ENV%\Lib\site-packages"
) else if not "%CONDA_PREFIX%"=="" (
    set "REF_SITE=%CONDA_PREFIX%\Lib\site-packages"
) else if not "%VIRTUAL_ENV%"=="" (
    set "REF_SITE=%VIRTUAL_ENV%\Lib\site-packages"
) else (
    set "REF_SITE="
)
set "VENV_SITE=%ENV_DIR%\Lib\site-packages"

if not exist "!REF_SITE!\prophet_spin" (
    echo WARNING: !REF_SITE! has no prophet_spin - the bundle will not bundle it.
    echo          Activate the automag environment or set AUTOMAG_SOURCE_ENV.
) else (
    for %%p in (prophet prophet_spin) do (
        if exist "!REF_SITE!\%%p" (
            if not exist "!VENV_SITE!\%%p" (
                echo     copying %%p from !REF_SITE!
                xcopy /E /I /Y /Q "!REF_SITE!\%%p" "!VENV_SITE!\%%p" >nul
            )
        )
    )
    for /d %%d in ("!REF_SITE!\prophet_mlip-*.dist-info") do (
        if not exist "!VENV_SITE!\%%~nxd" xcopy /E /I /Y /Q "%%d" "!VENV_SITE!\%%~nxd" >nul
    )
)

"%VENV_PY%" -c "import torch, prophet_spin; print('    torch', torch.__version__, '- cuda available:', torch.cuda.is_available())"
if errorlevel 1 (
    echo ERROR: torch/prophet_spin do not import in %ENV_DIR%.
    echo        Install prophet-mlip there, or set AUTOMAG_SOURCE_ENV to an
    echo        environment that has prophet_spin.
    popd & exit /b 1
)

if "%SETUP_ONLY%"=="1" (
    echo === environment ready: %ENV_DIR%
    popd & exit /b 0
)

if "%BUILD_EMBEDDED%"=="1" (
    echo === building dist\run_automag-cpu ^(checkpoints inside^)
    if exist "dist\run_automag-cpu" rmdir /S /Q "dist\run_automag-cpu"
    "%VENV_PY%" tools\build_executable.py --onedir --name run_automag-cpu
    if errorlevel 1 (echo ERROR: the build failed & popd & exit /b 1)
)

echo === building dist\run_automag-cpu-extmodels ^(checkpoints outside^)
if exist "dist\run_automag-cpu-extmodels" rmdir /S /Q "dist\run_automag-cpu-extmodels"
"%VENV_PY%" tools\build_executable.py --onedir --name run_automag-cpu-extmodels --no-models
if errorlevel 1 (echo ERROR: the build failed & popd & exit /b 1)

echo.
echo === result
for %%n in (run_automag-cpu run_automag-cpu-extmodels) do (
    if exist "dist\%%n" (
        for /f "delims=" %%m in ('powershell -NoProfile -Command "if (Test-Path 'dist\%%n') { [math]::Round((Get-ChildItem -Recurse -File 'dist\%%n' ^| Measure-Object -Property Length -Sum).Sum / 1MB) }"') do echo     dist\%%n: %%m MB
    )
)

echo.
echo Use the small bundle from the project directory you want to run in:
echo     set PROPHET_MODEL=C:\path\to\models\prophet-v1-mag.pt
echo     set PROPHET_HEAD=C:\path\to\models\prophet-v1-mag-head.pt
echo     "%ROOT%\dist\run_automag-cpu-extmodels\run_automag-cpu-extmodels.exe" --wait
echo.
echo Without bundling torch at all, keep the 359 MB bundle and point
echo AUTOMAG_HELPER_PYTHON at any environment that has torch + prophet_spin.
popd
exit /b 0

:usage
echo Usage: tools\build_cpu_bundle.bat [--no-models] [--setup-only]
echo.
echo   --no-models    build only the small bundle ^(checkpoints outside^)
echo   --setup-only   prepare the CPU build environment and stop
echo.
echo Environment: AUTOMAG_CPU_ENV, AUTOMAG_CPU_PYTHON, AUTOMAG_TORCH_CPU,
echo              AUTOMAG_SOURCE_ENV
popd
exit /b 0

:usage_error
popd
exit /b 2
