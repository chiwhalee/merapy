!fortran
!f90

!#define DATATYPE real*8
!q:
    
    ! why header.f90 don't work? 
    !1. why define not work?
    !1. how to use params
    ! XGEMM 等宏定义怎么处理？


subroutine test_add_inplace(data1, data2, totDim)
    integer rec(6), totDim 
    real*8:: data1(totDim)
    real*8:: data2(totDim)
!f2py intent(inout):: data3
    do i = 1, totDim
        data2(i)= data1(i) + data2(i)
    end do

end subroutine test_add_inplace


!RESHAPE is a function -- usually used in an assignment statement so there will be a copy operation. Perhaps it can be done without copying using pointers. Unless the array is huge, it is probably better to use RESHAPE.


subroutine contract_core_player_fort(data1, data2, data3, totDim1,totDim2, totDim3, rec, num_rec)
    !USE data_buffer
    integer  totDim1, totDim2, totDim3, num_rec
    integer rec(num_rec, 6)
    integer p1, p2, p3, Dim1, Dimc, Dim2
    real*8, target:: data1(totDim1), data2(totDim2), data3(totDim3)
!f2py intent(inout):: data3
    real*8:: alpha,beta
    real*8, pointer:: matrix1(:), matrix2(:), matrix3(:)

    !rec(:, 1:3)=rec(:, 1:3) + 1
        
    data3 = 0.d0
    matrix1 => data1
    matrix2 => data2
    matrix3 => data3
    !print "('T%data=', 200(1x,',', E10.4))", data3(:)
    alpha = 1.d0; beta=1.d0
    do i = 1, num_rec 
 
        !call Matrix_Multiply_inplace(Dim1, Dimc, Dim2, matrix1(p1), matrix2(p2), matrix3(p3), alpha, beta)
        !change new line so that older comiler like gfortran can pass 
        call Matrix_Multiply_inplace(rec(i,4), rec(i,6), rec(i,5), matrix1(rec(i,1) + 1), & 
          matrix2(rec(i,2) + 1), matrix3(rec(i,3) + 1), alpha, beta)
        !call Matrix_Multiply_inplace(rec(i,4), rec(i,6), rec(i,5), matrix1(rec(i,1)), matrix2(rec(i,2)), matrix3(rec(i,3)), alpha, beta)
     end do
end subroutine contract_core_player_fort

!MODULE data_buffer
!    real*8, pointer:: matrix_buf(:)
!    integer:: allocateddd
!    !allocate(matrix_buf(1000))
!
!END MODULE data_buffer

subroutine set_num_of_threads(n)
    integer n
    call OMP_set_num_threads(n)
end subroutine set_num_of_threads

subroutine get_max_threads(n) 
    integer n
!f2py intent(out)::n
    n=OMP_get_max_threads()
end subroutine get_max_threads

subroutine get_num_of_threads(n) 
    integer n
!f2py intent(out)::n
    n=OMP_get_num_threads()
end subroutine get_num_of_threads

subroutine get_num_of_procs(n) 
    integer n
!f2py intent(out)::n
    n=OMP_get_num_procs()
end subroutine get_num_of_procs


subroutine contract_core_player_fort_paralell_bac(data1, data2, data3, totDim1,totDim2, totDim3, rec, num_rec)
        !USE data_buffer
        integer  totDim1, totDim2, totDim3, num_rec
        integer rec(num_rec, 6)
        integer p1, p2, p3, Dim1, Dimc, Dim2
        real*8, target:: data1(totDim1), data2(totDim2) !, data3(totDim3)
        real*8:: data3(totDim3)
        real*8:: alpha,beta
        real*8, pointer:: matrix1(:), matrix2(:), matrix3(:)
        
        !real*8, pointer:: matrix11(:), matrix22(:), matrix33(:)
        real*8, pointer:: matrix11, matrix22, matrix33
        real*8, pointer:: matrix_temp(:)
            
        data3 = 0.d0
        matrix1 => data1;   matrix2 => data2;   
        !print "('T%data=', 200(1x,',', E10.4))", data3(:)
        alpha = 1.d0; beta=1.d0
    
    !!reduction(+:data3)
    !$OMP PARALLEL shared(data3) 
    !$OMP DO private(p1, p2, p3, dim1, dim2, dimc, matrix_temp) 
    
        do i = 1, num_rec 
            p1 = rec(i, 1) + 1
            p2 = rec(i, 2) + 1
            p3 = rec(i, 3) + 1
            Dim1 = rec(i, 4)
            Dim2 = rec(i, 5)
            Dimc = rec(i, 6)
            
            allocate(matrix_temp(Dim1*Dim2))
            matrix_temp = 0.d0
            !matrix_temp(:) = data3(p3:Dim1*Dim2)
            call Matrix_Multiply_inplace(Dim1, Dimc, Dim2, matrix1(p1), matrix2(p2), matrix_temp, alpha, beta)
            !$FLUSH(data3)
            data3(p3:p3 + Dim1*Dim2-1) = data3(p3:p3 + Dim1*Dim2-1)  +  Matrix_temp(:)
            !data3(p3:Dim1*Dim2) = data3(p3:Dim1*Dim2)  +  Matrix_temp(:)
            !$OMP critical
            !print *, p3, data3(:)
            !print *, Matrix_temp(:)
            !$OMP end critical
            !data3(p3:Dim1*Dim2) = data3(p3:Dim1*Dim2)  +  Matrix_temp(:)
            deallocate(Matrix_temp)
            
         end do
    !$OMP END DO 
    !$OMP END PARALLEL 
end subroutine contract_core_player_fort_paralell_bac

!!!!
    !OMP atomic is useless, since it only support scalars
    !one can use any of OMP ordered,  OMP critical, OMP reduction to solve the race
    !condition problem

subroutine contract_core_player_fort_paralell_ordered(data1, data2, data3, totDim1,totDim2, totDim3, rec, num_rec)
        !USE data_buffer
        integer  totDim1, totDim2, totDim3, num_rec
        integer rec(num_rec, 6)
        integer p1, p2, p3, Dim1, Dimc, Dim2
        real*8, target:: data1(totDim1), data2(totDim2) !, data3(totDim3)
        real*8:: data3(totDim3)
        real*8:: alpha,beta
        real*8, pointer:: matrix1(:), matrix2(:), matrix3(:)
        
        !real*8, pointer:: matrix11(:), matrix22(:), matrix33(:)
        real*8, pointer:: matrix11, matrix22, matrix33
        real*8, pointer:: matrix_temp(:)
            
        data3 = 0.d0
        matrix1 => data1;   matrix2 => data2;   
        !print "('T%data=', 200(1x,',', E10.4))", data3(:)
        alpha = 1.d0; beta=1.d0
    
    !!reduction(+:data3)
    !$OMP PARALLEL shared(data3) 
    !!!$OMP PARALLEL 
    !$OMP DO ordered schedule(runtime) private(p1, p2, p3, dim1, dim2, dimc, matrix_temp) 
    
        do i = 1, num_rec 
            p1 = rec(i, 1) + 1
            p2 = rec(i, 2) + 1
            p3 = rec(i, 3) + 1
            Dim1 = rec(i, 4)
            Dim2 = rec(i, 5)
            Dimc = rec(i, 6)
            
            allocate(matrix_temp(Dim1*Dim2))
            matrix_temp = 0.d0
            !matrix_temp(:) = data3(p3:Dim1*Dim2)
            call Matrix_Multiply_inplace(Dim1, Dimc, Dim2, matrix1(p1), matrix2(p2), matrix_temp, alpha, beta)
            !$OMP ordered
            
            !print *, i, p3, Dim1*Dim2 !, data3(:8) 
            
            !print "('T%data=', 200(1x,',', E10.2))", data3(p3:p3 + Dim1*Dim2)
            !data3(p3:p3 + Dim1*Dim2) = data3(p3:p3 + Dim1*Dim2)  +  Matrix_temp(1:1+Dim1*Dim2)
            !data3(p3:p3 + Dim1*Dim2-1) = data3(p3:p3 + Dim1*Dim2-1)  +  Matrix_temp(1:1+Dim1*Dim2-1)
            data3(p3:p3 + Dim1*Dim2-1) = data3(p3:p3 + Dim1*Dim2-1)  +  Matrix_temp(:)
            
            !do j = 1, Dim1*Dim2
            !    data3(p3 + j-1)=data3(p3 + j-1) + Matrix_temp(j)
            !end do

            !print "('T%data=', 200(1x,',', E10.2))", data3(p3:p3 + Dim1*Dim2)
            !print "('T%data=', 200(1x,',', E10.2))", Matrix_temp(:)
            !print "('T%data=', 200(1x,',', E10.2))", data3(:8)
            
            !$OMP end ordered
            
            deallocate(Matrix_temp)
            
         end do
    !$OMP END DO 
    !$OMP END PARALLEL 
end subroutine contract_core_player_fort_paralell_ordered

subroutine contract_core_player_fort_paralell_critical(data1, data2, data3, totDim1,totDim2, totDim3, rec, num_rec)
        !USE data_buffer
        integer  totDim1, totDim2, totDim3, num_rec
        integer rec(num_rec, 6)
        integer p1, p2, p3, Dim1, Dimc, Dim2
        real*8, target:: data1(totDim1), data2(totDim2) !, data3(totDim3)
        real*8:: data3(totDim3)
        real*8:: alpha,beta
        real*8, pointer:: matrix1(:), matrix2(:), matrix3(:)
        
        !real*8, pointer:: matrix11(:), matrix22(:), matrix33(:)
        real*8, pointer:: matrix11, matrix22, matrix33
        real*8, pointer:: matrix_temp(:)
            
        data3 = 0.d0
        matrix1 => data1;   matrix2 => data2;   
        !print "('T%data=', 200(1x,',', E10.4))", data3(:)
        alpha = 1.d0; beta=1.d0
    
    !!reduction(+:data3)
    !$OMP PARALLEL shared(data3) 
    !!!$OMP PARALLEL 
    !$OMP DO schedule(runtime) private(p1, p2, p3, dim1, dim2, dimc, matrix_temp) 
    
        do i = 1, num_rec 
            p1 = rec(i, 1) + 1
            p2 = rec(i, 2) + 1
            p3 = rec(i, 3) + 1
            Dim1 = rec(i, 4)
            Dim2 = rec(i, 5)
            Dimc = rec(i, 6)
            
            allocate(matrix_temp(Dim1*Dim2))
            matrix_temp = 0.d0
            !matrix_temp(:) = data3(p3:Dim1*Dim2)
            call Matrix_Multiply_inplace(Dim1, Dimc, Dim2, matrix1(p1), matrix2(p2), matrix_temp, alpha, beta)
            !$OMP critical
            
            !print *, i, p3, Dim1*Dim2 !, data3(:8) 
            
            !print "('T%data=', 200(1x,',', E10.2))", data3(p3:p3 + Dim1*Dim2)
            !data3(p3:p3 + Dim1*Dim2) = data3(p3:p3 + Dim1*Dim2)  +  Matrix_temp(1:1+Dim1*Dim2)
            !data3(p3:p3 + Dim1*Dim2-1) = data3(p3:p3 + Dim1*Dim2-1)  +  Matrix_temp(1:1+Dim1*Dim2-1)
            data3(p3:p3 + Dim1*Dim2-1) = data3(p3:p3 + Dim1*Dim2-1)  +  Matrix_temp(:)
            
            !do j = 1, Dim1*Dim2
            !    data3(p3 + j-1)=data3(p3 + j-1) + Matrix_temp(j)
            !end do

            !print "('T%data=', 200(1x,',', E10.2))", data3(p3:p3 + Dim1*Dim2)
            !print "('T%data=', 200(1x,',', E10.2))", Matrix_temp(:)
            !print "('T%data=', 200(1x,',', E10.2))", data3(:8)
            
            !$OMP end critical
            
            deallocate(Matrix_temp)
            
         end do
    !$OMP END DO 
    !$OMP END PARALLEL 
end subroutine contract_core_player_fort_paralell_critical

subroutine contract_core_player_fort_paralell_reduction(data1, data2, data3, totDim1,totDim2, totDim3, rec, num_rec)
        !USE data_buffer
        integer  totDim1, totDim2, totDim3, num_rec
        integer rec(num_rec, 6)
        integer p1, p2, p3, Dim1, Dimc, Dim2
        real*8, target:: data1(totDim1), data2(totDim2) !, data3(totDim3)
        real*8:: data3(totDim3)
        real*8:: alpha,beta
        real*8, pointer:: matrix1(:), matrix2(:), matrix3(:)
        
        !real*8, pointer:: matrix11(:), matrix22(:), matrix33(:)
        real*8, pointer:: matrix11, matrix22, matrix33
        real*8, pointer:: matrix_temp(:)
        
        !this line may be not needed, since reduction default is 0        
        data3 = 0.d0
        matrix1 => data1;   matrix2 => data2;   
        !print "('T%data=', 200(1x,',', E10.4))", data3(:)
        alpha = 1.d0; beta=1.d0
    
    !$OMP PARALLEL shared(data3) 
    !$OMP DO schedule(runtime) private(p1, p2, p3, dim1, dim2, dimc, matrix_temp) &
    !$OMP reduction(+:data3)
    
        do i = 1, num_rec 
            p1 = rec(i, 1) + 1
            p2 = rec(i, 2) + 1
            p3 = rec(i, 3) + 1
            Dim1 = rec(i, 4)
            Dim2 = rec(i, 5)
            Dimc = rec(i, 6)
            
            allocate(matrix_temp(Dim1*Dim2))
            matrix_temp = 0.d0
            call Matrix_Multiply_inplace(Dim1, Dimc, Dim2, matrix1(p1), matrix2(p2), matrix_temp, alpha, beta)
            
            data3(p3:p3 + Dim1*Dim2-1) = data3(p3:p3 + Dim1*Dim2-1)  +  Matrix_temp(:)
            deallocate(Matrix_temp)
            
         end do
    !$OMP END DO 
    !$OMP END PARALLEL 
end subroutine contract_core_player_fort_paralell_reduction

subroutine contract_core_player_fort_paralell_reduction_1(data1, data2, data3, totDim1,totDim2, totDim3, rec, num_rec)
        !USE data_buffer
        integer  totDim1, totDim2, totDim3, num_rec
        integer rec(num_rec, 6)
        integer p1, p2, p3, Dim1, Dimc, Dim2, d12
        real*8, target:: data1(totDim1), data2(totDim2) !, data3(totDim3)
        real*8:: data3(totDim3)
        real*8:: alpha,beta
        real*8, pointer:: matrix1(:), matrix2(:), matrix3(:)
        
        !real*8, pointer:: matrix11(:), matrix22(:), matrix33(:)
        real*8, pointer:: matrix11, matrix22, matrix33
        real*8, pointer:: matrix_temp(:)
        logical::  alloced = .False.
        
        !this line may be not needed, since reduction default is 0        
        data3 = 0.d0
        matrix1 => data1;   matrix2 => data2;   
        !print "('T%data=', 200(1x,',', E10.4))", data3(:)
        alpha = 1.d0; beta=1.d0
    
        alloced = .False.
        d12 = 0
    !$OMP PARALLEL shared(data3) 
    !$OMP DO schedule(runtime) private(p1, p2, p3, dim1, dim2, dimc, matrix_temp, alloced, d12) &
    !$OMP reduction(+:data3)
    
        do i = 1, num_rec 
            p1 = rec(i, 1) + 1
            p2 = rec(i, 2) + 1
            p3 = rec(i, 3) + 1
            Dim1 = rec(i, 4)
            Dim2 = rec(i, 5)
            Dimc = rec(i, 6)
            if (d12<Dim1*Dim2) alloced=.False.
            !if(not(alloced)) then
            if(.not.alloced) then
                allocate(matrix_temp(Dim1*Dim2))
                alloced = .true.
            end if
            d12 = Dim1*Dim2
            matrix_temp = 0.d0
            call Matrix_Multiply_inplace(Dim1, Dimc, Dim2, matrix1(p1), matrix2(p2), matrix_temp, alpha, beta)
            
            data3(p3:p3 + d12-1) = data3(p3:p3 + d12-1)  +  Matrix_temp(:)
            !deallocate(Matrix_temp)
            
         end do
    !$OMP END DO 
    !$OMP END PARALLEL 
end subroutine contract_core_player_fort_paralell_reduction_1

subroutine contract_core_player_fort_paralell_test(data1, data2, data3, totDim1,totDim2, totDim3, rec, num_rec)
        !USE data_buffer
        integer  totDim1, totDim2, totDim3, num_rec
        integer rec(num_rec, 6)
        integer p1, p2, p3, Dim1, Dimc, Dim2
        real*8, target:: data1(totDim1), data2(totDim2), data3(totDim3)
        real*8:: alpha,beta
        real*8, pointer:: matrix1(:), matrix2(:), matrix3(:)
        
        real*8, pointer:: matrix11, matrix22, matrix33
            
        data3 = 0.d0
        matrix1 => data1;   matrix2 => data2;   matrix3 => data3

        !print "('T%data=', 200(1x,',', E10.4))", data3(:)
        alpha = 1.d0; beta=1.d0
    
    !$OMP PARALLEL shared(data3, matrix1, matrix2, matrix3) 
    !$OMP DO schedule(runtime) private(p1, p2, p3, dim1, dim2, dimc) 
    
        do i = 1, num_rec 
            p1 = rec(i, 1) + 1
            p2 = rec(i, 2) + 1
            p3 = rec(i, 3) + 1
            Dim1 = rec(i, 4)
            Dim2 = rec(i, 5)
            Dimc = rec(i, 6)
            !!!$OMP FLUSH (matrix3)
            !!!$OMP FLUSH 
            call Matrix_Multiply_inplace(Dim1, Dimc, Dim2, matrix1(p1), matrix2(p2), matrix3(p3), alpha, beta)
         end do
    !$OMP END DO 
    !$OMP END PARALLEL 
end subroutine contract_core_player_fort_paralell_test




function vector_distance(n,A,B)
   implicit none
   integer n
   real*8:: A(n), B(n)
   real*8 vector_distance
   integer i
   
   vector_distance = 0.d0
   do i = 1, n
      vector_distance = vector_distance+abs(A(i)-B(i))
   end do
end function vector_distance


subroutine transpose4py(A,nA,At)
!
   integer nA
   real*8::A(nA,nA)
   real*8::At(nA,nA)
!f2py intent(in,out)::At   
   At=transpose(A)

end subroutine transpose4py

subroutine Unit_Matrix(A, nA)
! this is to be replaced by  A=np.identity(nA)
   integer nA
   !DATATYPE:: A(nA, nA)
   !real*8,intent(inout):: A(nA, nA)   ! this is needless 
   real*8:: A(nA, nA)
!f2py intent(inout):: A
   integer i
   
   A = 0.d0
   do i = 1, nA
      A(i,i) = 1.d0
   end do
end subroutine Unit_Matrix

subroutine Matrix_Direct_Product(A,nA,mA, B,nB,mB, C)
! this maybe replaced by numpy' method
   implicit none
   integer nA,mA, nB,mB
   !DATATYPE A(nA,mA), B(nB,mB), C(nA*nB, mA*mB)
   real*8 A(nA,mA), B(nB,mB), C(nA*nB, mA*mB)
!f2py intent(out)::C
   !real*8, intent(in,out)::C(nA*nB, mA*mB)

   integer i,j, pi, pj
!$OMP PARALLEL DO private(i,j,pi,pj)
   do i = 1, nB
      pi = nA*(i-1)
      do j = 1, mB
         pj = mA*(j-1)
         !整个矩阵块A乘以B(i, j)
         !C^{i1, i2}_{j1, j2} = A^{i1}_{j1} * B^{i2}_{j2}
         C(pi+1:pi+nA, pj+1:pj+mA) = A(:,:)*B(i,j)
      end do
   end do
!$OMP END PARALLEL DO   
end subroutine Matrix_Direct_Product

subroutine Matrix_Get_Position(p, rank, pos, Dims)
   !: pos-->p
   !: array start from zero
   integer p, rank
!f2py intent(out):: p   
   integer pos(rank), Dims(rank)
   integer i
   
   if (rank.eq.0) then
      p = 1
      return
   end if
   
   p = 0
   do i = rank, 2, -1
      !p = (p-1+pos(i))*Dims(i-1)
      p = (p+pos(i))*Dims(i-1) 
      !this line is modified by lzh to be compatible  with numpy
   end do
   p = p+pos(1)
end subroutine Matrix_Get_Position

subroutine Matrix_Get_Position_rev(p, rank, pos, Dims)
   !: p-->pos
   integer p, rank
   integer pos(rank), Dims(rank)
!f2py intent(out)::pos
   integer i, res
   
   pos = 0
   !res = p-1
   !modified by lzh
   res=p
   do i = 1, rank
      !pos(i) = mod(res, Dims(i))+1
      !modified by lzh
      pos(i) = mod(res, Dims(i))
      res = res/Dims(i)
      if (res.eq.0) return
   end do
end subroutine Matrix_Get_Position_rev

subroutine Matrix_Multiply_inplace(n,l,m, A,B,C,alpha,beta)
    !use params
    !attention_omitted_something  seems params is not used at all
   implicit none
   integer n,l,m
   real*8 A(n,l), B(l,m)
   real*8,intent(inout)::C(n,m)
!f2py intent(out):: C
   real*8 :: alpha,beta

!!$   C=matmul(A,B)
!see http://www.math.utah.edu/software/lapack/lapack-blas/zgemm.html
   !call XGEMM('N', 'N', n, m, l, alpha, A, n, B, l, beta, C, n)
   call dgemm('N', 'N', n, m, l, alpha, A, n, B, l, beta, C, n)
   
end subroutine Matrix_Multiply_inplace


subroutine Matrix_Multiply(n,l,m, A,B,C,alpha,beta)
    !use params
    !attention_omitted_something  seems params is not used at all
    implicit none
    integer n,l,m
    real*8 A(n,l), B(l,m), C(n, m)
!f2py intent(out):: C
    real*8 :: alpha,beta

!!$   C=matmul(A,B)
!see http://www.math.utah.edu/software/lapack/lapack-blas/zgemm.html
    !call XGEMM('N', 'N', n, m, l, alpha, A, n, B, l, beta, C, n)
    call dgemm('N', 'N', n, m, l, alpha, A, n, B, l, beta, C, n)
end subroutine Matrix_Multiply

subroutine Matrix_Multiply_complex(n,l,m, A,B,C,alpha,beta)
    implicit none
    integer n,l,m
    complex*16 A(n,l), B(l,m), C(n, m)
!f2py intent(out):: C
    real*8 :: alpha,beta

!!$   C=matmul(A,B)
!see http://www.math.utah.edu/software/lapack/lapack-blas/zgemm.html
    !call XGEMM('N', 'N', n, m, l, alpha, A, n, B, l, beta, C, n)
    call zgemm('N', 'N', n, m, l, alpha, A, n, B, l, beta, C, n)
end subroutine Matrix_Multiply_complex 


subroutine Matrix_eigen_vector_bac(nV, A, E)
      integer nV
      real*8 A(nV,nV)
!f2py intent(inout):: A
      integer,parameter:: E_size=1024*64
      integer:: lwork, info
      !real*8:: E(E_size)
      real*8:: E(nV)
!f2py intent(inout):: E
      
      !real*8:: work(5*E_size)
      !lwork = 5*E_size
    
      !lzh:  issue: I dont understand how to set lwork, and work 
      real*8:: work(5*nV)
      !lwork = 5*E_size
      lwork = 5*nV
      !issue: should take use of work to speed up 
      call dsyev('V', 'U', nV, A, nV, E, work, lwork, info)   

end subroutine Matrix_eigen_vector_bac 


subroutine Matrix_eigen_vector(nV, A, E)
      integer nV
      real*8 A(nV,nV)
!f2py intent(inout):: A
      integer:: lwork, info
      real*8:: E(nV)
!f2py intent(inout):: E
      
      !integer,parameter:: E_size=1024*64
      !real*8:: work(5*E_size)
      !lwork = 5*E_size
      !lzh:  issue: I dont understand how to set lwork, and work 
      real*8:: work(6*nV)
      lwork = 6*nV

      call dsyev('V', 'U', nV, A, nV, E, work, lwork, info)   

    end subroutine Matrix_eigen_vector




subroutine Matrix_Multiply_inplace_del(n,l,m, A,B,C,alpha,beta)
    !use params
    !attention_omitted_something  seems params is not used at all
   implicit none
   integer n,l,m
   real*8 A(n,l), B(l,m)
   real*8 C(n,m)
!f2py intent(inout):: C
   real*8 :: alpha,beta

!!$   C=matmul(A,B)
!see http://www.math.utah.edu/software/lapack/lapack-blas/zgemm.html
   !call XGEMM('N', 'N', n, m, l, alpha, A, n, B, l, beta, C, n)
   call dgemm('N', 'N', n, m, l, alpha, A, n, B, l, beta, C, n)
   
end subroutine Matrix_Multiply_inplace_del



subroutine Matrix_SVD(m,n, min_mn, A, U, S, VT)
    !svd of A, a m by n matrix,   S is singular values
   implicit none
   integer m,n, min_mn
   real*8:: A(m,n), U(m,min_mn), VT(min_mn,n)
!f2py intent(out):: U, VT
   real*8:: S(min_mn), work(10*(m+n))
!f2py intent(out):: S 
   integer info, lwork
   
   lwork = 10*(m+n)
!!$   print *, 'm,n,min_mn=', m, n, min_mn
   !call XGESVD('S', 'S', m, n, A, m, S, U, m, VT, min_mn, work, lwork, info)
   call dgesvd('S', 'S', m, n, A, m, S, U, m, VT, min_mn, work, lwork, info)
!!$   print *, 'svd_info=', info
end subroutine Matrix_SVD

subroutine Matrix_SVD_Unitilize(m,n, A,A1, S)
! svd of A,  then let A1 = U*VT
   implicit none
   integer m,n, min_mn
   real*8:: A(m,n), U(m,m), VT(n,n), A1(m, n)
!f2py intent(out):: A1   
   real*8:: S(m+n), work(10*(m+n))
!f2py intent(out)::S 
   integer info, lwork
   
   lwork = 10*(m+n)
   min_mn = min(m,n)
!!$   print *, 'm,n,min_mn=', m, n, min_mn
!attention_this_may_be_wrong   !XGESVD
! XGESVD 在header.f90 中定义        
   call dgesvd('S', 'S', m, n, A, m, S, U, m, VT, n, &
        work, lwork, info)
    
   call dgemm('N', 'N', m, min_mn, n, 1.d0, U, m, Vt, n, &
        0.d0, A1, m)
!!$   print "('SVD=', 200(1x,E15.8))", S(1:min_mn)
end subroutine Matrix_SVD_Unitilize


function Matrix_Trace(n, A) RESULT(X)
   integer n
   real*8:: A(n,n), X
   integer i
   X = 0.d0
   do i = 1, n
      X = X+A(i,i)
   end do
end function Matrix_Trace

!!$ shift one dimension array 
!!$ A(1:n-i+1) = A(i:n)
!!$ A(n-i:) = b
!lzh:  EOSHIFT 在老版本的ifort上有bug，干脆把整个函数注释掉
!subroutine iShift(n,A,i,b)
!   implicit none
!   integer n,i, b, A(n)
!!f2py intent(inout):: A
!   integer X(n)
!   X=EOSHIFT(A, shift=i, boundary=b)
!   A=X
!end subroutine iShift
 
subroutine Set_Matrix(nA,mA,A, nB,mB,B, x,y, forward)
!this maybe replaced by some numpy methods
! 把B塞到A中, 起始坐标为x, y; 或反过来
   integer nA,mA, nB,mB, x,y
   logical forward
   real*8:: A(nA,mA), B(nB,mB)   
!f2py intent(inout):: A, B
   if (forward) then
      A(x+1:x+nB,y+1:y+mB) = B(1:nB,1:mB)
   else
      B(1:nB,1:mB) = A(x+1:x+nB,y+1:y+mB)
   end if
end subroutine Set_Matrix




