
f2py -c -m random_gfort ../random4py.f90  --fcompiler=gfortran --compiler=mingw32
f2py -c -m array_permutation_gfort ../array_permutation4py.f90 --fcompiler=gfortran --compiler=mingw32  --f90flags="-openmp" -lgomp  
::set path=%path%;C:\Users\zhihua\Dropbox\My-documents\My-code\quantum-many-body\mera-algorithms\python\merapy\lib\win
::set LIB=%LIB%;C:\Users\zhihua\Dropbox\My-documents\My-code\quantum-many-body\mera-algorithms\python\merapy\lib\win
::f2py -c -m common_gfort ../common4py.f90  -llapack -lblas  --f90flags="-openmp" --fcompiler=gfortran  --compiler=mingw32  --f90flags="-openmp" -lgomp   -LC:\Users\zhihua\Dropbox\My-documents\My-code\quantum-many-body\mera-algorithms\python\merapy\lib\win

:: note I put liblapack.dll etc under currect dir so -L. and -llapack 
f2py -c -m common_gfort ../common4py.f90  -llapack -lblas  -L.  --f90flags="-openmp" --fcompiler=gfortran  --compiler=mingw32  --f90flags="-openmp" -lgomp   

:: instruction 
    :: without -lgomp then, error:  set_num_thre... not defined 
    :: 



