
::  run this file under vim by  ': !%'

:: set path=%path%;C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\bin\amd64\


:: 要不然找不到fortran compiler
set path=%path%;C:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2016.0.110\windows\bin\intel64\
set LIB=%LIB%;C:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2016.0.110\windows\compiler\lib\intel64_win\

: : set path=%path%;C:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2016.0.110\windows\bin\ia32\
: : set LIB=%LIB%;C:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2016.0.110\windows\compiler\lib\ia32_win\
        
: :  下面一行很重要，f2py 默认需要vs2008，下面指向了 vs 2015 即 14.0 
:: set VS90COMNTOOLS=%VS140COMNTOOLS%
:: set VS90COMNTOOLS=C:\Program Files (x86)\Microsoft Visual Studio 14.0\Common7\Tools\
:: set VS90COMNTOOLS=C:\Program Files\Microsoft Visual Studio 12.0\Common7\Tools\
set VS90COMNTOOLS=C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\bin\

:: f2py -c -m random_64_ifort ../random4py.f90 --fcompiler=intelvem
:: f2py -c -m random_64_ifort ../random4py.f90 
:: f2py -c -m TTT test.c --fcompiler=intelvem  --build-dir build 
::  f2py -c -m aaa test.c --build-dir build 
:: f2py -c -m ttt test.f90 --build-dir build 

::  f2py -c -m random_64_ifort ../random4py.f90 --fcompiler=msvc 



