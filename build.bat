@echo off
REM Compile CLR managed code. Dump .dll, .obj and .pdb to out directory.
cd ./out
call "C:\Program Files (x86)\Microsoft Visual Studio\2019\BuildTools\VC\Auxiliary\Build\vcvarsall.bat" x64
set compilerflags=/Clr /Od /Zi
set linkerflags=/DLL /OUT:./QPfloat.dll
cl.exe %compilerflags% ../src/QPfloat/*.cpp /link %linkerflags% 