@echo off
REM Compile CLR managed code. Dump .dll, .obj and .pdb to out directory.
call "C:\Program Files (x86)\Microsoft Visual Studio\2019\BuildTools\VC\Auxiliary\Build\vcvarsall.bat" x64
cd ../out/debug
set compilerflags=/clr /Od /Zi
set linkerflags=/DLL /OUT:QPFloat.dll
set includePaths=/I INCLUDE /I..\..\..\..\src\
cl.exe %includePaths% %compilerflags% ../../src/*.cpp /link %linkerflags% 