@echo off

set BUILD_DIR=%~dp0..\build\exe

if not exist %BUILD_DIR% mkdir %BUILD_DIR%

cd %BUILD_DIR%

pyinstaller --onefile ../../eoriver/cli.py

cd dist
move cli.exe eoriver.exe
cd ..