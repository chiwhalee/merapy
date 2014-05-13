
subroutine Array_Permutation_fort_parallel(rank, Dims, order_, totDim, A,B)
    ! rank是腿的个数，
    !Dims是每条腿的维度，
    !order是permutation的次序，
    !A,B是存放数据的array, 
    !totDim是总维度(各条腿维度之积，fortran子程序需要指定array的长度）

       implicit none
       integer:: rank, Dims(32), order_(32), totDim
       real*8,intent(IN):: A(totDim)
       real*8,intent(inout)::B(totDim)
!f2py intent(inout) :: B
       integer p, q, idx(32), order(32), rorder(32)
       integer nDims(32), nidx(32)
       integer strides(32), rstrides(32)
       integer ntotDim, inc
       integer i,j, np, p1,p2, pp, res
       integer nths, nth
       integer, external:: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
       character*128 fmt
       !these are added by lzh, to make it compatible with python index which
       !starts at 0
       do i = 1, rank
            order(i)=order_(i) + 1
       end do
       
       ntotDim = 1
       do i = 1, rank
          ntotDim = ntotDim*Dims(i)
          nDims(i) = Dims(order(i))
          rorder(order(i)) = i
       end do
       if (ntotDim.ne.totDim) then 
           print *, 'error in ArrayPermutation, size not match'
           print *, "ntotDim", ntotDim, "totDim", totDim
           stop
        end if
       
       strides = 1
       do i = 1, rank
          do j = 1, rorder(i)-1
             strides(i) = strides(i)*nDims(j)
          end do
          rstrides(i) = (Dims(i)-1)*strides(i)
       end do
!!$   write(fmt,*) "(2(I5,1x,'| '),", rank, "(1x,I5), ' | ',", &
!!$        rank, "(1x,I5), ' | ', 200(1x,I5))"
!$OMP PARALLEL private(p,q,i,np,inc,idx,p1,p2,pp,res)   
       nths = OMP_GET_NUM_THREADS()
       nth = OMP_GET_THREAD_NUM()
       
       np = totDim/nths
       p1 = nth*np+1
       p2 = (nth+1)*np
       if (nth.eq.nths-1) p2 = totDim
       
       pp = p1-1
       do i = 1, rank-1
          res = pp/Dims(i)
          idx(i) = pp-res*Dims(i)+1      
          pp = res
       end do
       idx(rank) = pp+1   
       
       q = idx(order(rank))-1
       do i = rank-1, 1, -1
          q = nDims(i)*q+idx(order(i))-1
       end do
       q = q+1
       
       do p = p1, p2
          B(q) = A(p)
          
          inc = 1
          i = 1
          do while ((inc.eq.1).and.(i.le.rank))
             idx(i) = idx(i)+1
             if (idx(i).le.Dims(i)) then
                q = q+strides(i)
                inc = 0            
             else
                q = q-rstrides(i)
                idx(i) = 1
                i = i+1            
             end if
          end do
       end do
!$OMP END PARALLEL   
end subroutine Array_Permutation_fort_parallel


subroutine Array_Permutation_fort(rank, Dims, order_, totDim, A,B)
       implicit none
       integer,intent(IN):: rank, Dims(32), order_(32), totDim
       real*8,intent(IN):: A(totDim)
       real*8,intent(inout)::B(totDim)
!f2py intent(inout) :: B
       integer p, q, idx(32), rorder(32), order(32)
       integer nDims(32), nidx(32)
       integer strides(32), rstrides(32)
       integer ntotDim, inc
       integer i,j, np, p1, p2, pp, res
       character*128 fmt
        
        
       ntotDim = 1
       do i = 1, rank
            order(i)=order_(i) + 1
       end do
       
        do i = 1, rank
          ntotDim = ntotDim*Dims(i)
          nDims(i) = Dims(order(i))
          !rorder(order(i)+1) = i-1
          rorder(order(i)) = i
       end do
       if (ntotDim.ne.totDim) print *, 'error in ArrayPermutation, size not match'
       
       !print *, "order", order(0:rank)
       !print *, "rorder", rorder(0:rank)
       
       strides = 1
       do i = 1, rank
          !do j = 1, i-1
          do j = 1, rorder(i)-1
             strides(i) = strides(i)*nDims(j)
          end do
          rstrides(i) = (Dims(i)-1)*strides(i)
       end do

       !print *, "strides", strides(0:rank)
       !print *, "rstrides", rstrides(0:rank)

       
       p1 = 1; p2 = totDim   
       q=1; idx(1:rank) = 1
       do p = p1, p2
          B(q) = A(p)
          !print  * , "q, ", q,"B(q)", B(q)
          !print  * , "A", "(E5.5)", A(p), B(p)
          
          inc = 1
          i = 1
          do while ((inc.eq.1).and.(i.le.rank))
             !print  * , "q" 
             idx(i) = idx(i)+1
             if (idx(i).le.Dims(i)) then
                q = q+strides(i)
                inc = 0            
             else
                q = q-rstrides(i)
                idx(i) = 1
                i = i+1            
             end if
          end do
       end do
end subroutine Array_Permutation_fort


subroutine Array_Permutation_fort_test(rank, Dims, order_, totDim, A,B)
       implicit none
       integer rank, totDim
       real*8,intent(IN):: A(totDim)
       real*8,intent(inout)::B(totDim)
!f2py intent(inout) :: B
       integer, intent(IN):: Dims(rank), order_(rank)
       integer p, q, idx(32), rorder(32), order(32)
       integer nDims(32), nidx(32)
       integer strides(32), rstrides(32)
       integer ntotDim, inc
       integer i,j, np, p1, p2, pp, res
       character*128 fmt
        
        
       ntotDim = 1
       do i = 1, rank
            order(i)=order_(i) + 1
       end do
       
        do i = 1, rank
          ntotDim = ntotDim*Dims(i)
          nDims(i) = Dims(order(i))
          !rorder(order(i)+1) = i-1
          rorder(order(i)) = i
       end do
       if (ntotDim.ne.totDim) print *, 'error in ArrayPermutation, size not match'
       
       !print *, "order", order(0:rank)
       !print *, "rorder", rorder(0:rank)
       
       strides = 1
       do i = 1, rank
          !do j = 1, i-1
          do j = 1, rorder(i)-1
             strides(i) = strides(i)*nDims(j)
          end do
          rstrides(i) = (Dims(i)-1)*strides(i)
       end do

       !print *, "strides", strides(0:rank)
       !print *, "rstrides", rstrides(0:rank)

       
       p1 = 1; p2 = totDim   
       q=1; idx(1:rank) = 1
       do p = p1, p2
          B(q) = A(p)
          !print  * , "q, ", q,"B(q)", B(q)
          !print  * , "A", "(E5.5)", A(p), B(p)
          
          inc = 1
          i = 1
          do while ((inc.eq.1).and.(i.le.rank))
             !print  * , "q" 
             idx(i) = idx(i)+1
             if (idx(i).le.Dims(i)) then
                q = q+strides(i)
                inc = 0            
             else
                q = q-rstrides(i)
                idx(i) = 1
                i = i+1            
             end if
          end do
       end do
end subroutine Array_Permutation_fort_test


subroutine permute_player_fort(rank, nidx, tape_ind, tape_dim, tape_ord, totDim, data1, data2)
       integer:: rank, nidx, totDim
       integer:: pidx, qidx, length
       integer:: tape_ind(nidx, 3), tape_dim(nidx, 32), tape_ord(nidx, 32)
       real*8,intent(IN):: data1(totDim)
       real*8,intent(inout)::data2(totDim)
!f2py intent(inout)::data2
    do i  = 1, nidx 
        pidx = tape_ind(i, 1) + 1
        qidx = tape_ind(i, 2) + 1
        length = tape_ind(i, 3)
        call Array_Permutation_fort_parallel(rank, tape_dim(i, 1:32),&
            tape_ord(i, 1:32), length, data1(pidx:pidx + length), data2(qidx:qidx + length))
    end do
end subroutine permute_player_fort



subroutine permute_player_fort_parallel(rank, nidx, tape_ind, tape_dim, tape_ord, totDim, data1, data2)
       integer:: rank, nidx, totDim
       integer:: pidx, qidx, length
       integer:: tape_ind(nidx, 3), tape_dim(nidx, 32), tape_ord(nidx, 32)
       real*8,intent(IN):: data1(totDim)
       real*8,intent(inout)::data2(totDim)
!f2py intent(inout)::data2

!$OMP PARALLEL private(pidx, qidx, length)
!$OMP do
    do i  = 1, nidx 
        pidx = tape_ind(i, 1) + 1
        qidx = tape_ind(i, 2) + 1
        length = tape_ind(i, 3)
        call Array_Permutation_fort(rank, tape_dim(i, 1:32),&
            tape_ord(i, 1:32), length, data1(pidx:pidx + length), data2(qidx:qidx + length))
    end do
!$OMP END do nowait
!$OMP END PARALLEL
end subroutine permute_player_fort_parallel



subroutine permute_player_fort_parallel_dynamic(rank, nidx, tape_ind, tape_dim, tape_ord, totDim, data1, data2)
       integer:: rank, nidx, totDim
       integer:: pidx, qidx, length
       integer:: tape_ind(nidx, 3), tape_dim(nidx, 32), tape_ord(nidx, 32)
       real*8,intent(IN):: data1(totDim)
       real*8,intent(inout)::data2(totDim)
!f2py intent(inout)::data2

!$OMP PARALLEL private(pidx, qidx, length)
!$OMP do schedule(dynamic)
    do i  = 1, nidx 
        pidx = tape_ind(i, 1) + 1
        qidx = tape_ind(i, 2) + 1
        length = tape_ind(i, 3)
        call Array_Permutation_fort(rank, tape_dim(i, 1:32),&
            tape_ord(i, 1:32), length, data1(pidx:pidx + length), data2(qidx:qidx + length))
    end do
!$OMP END do 
!$OMP END PARALLEL
end subroutine permute_player_fort_parallel_dynamic


subroutine permute_player_fort_parallel_guided(rank, nidx, tape_ind, tape_dim, tape_ord, totDim, data1, data2)
       integer:: rank, nidx, totDim
       integer:: pidx, qidx, length
       integer:: tape_ind(nidx, 3), tape_dim(nidx, 32), tape_ord(nidx, 32)
       real*8,intent(IN):: data1(totDim)
       real*8,intent(inout)::data2(totDim)
!f2py intent(inout)::data2

!$OMP PARALLEL private(pidx, qidx, length)
!$OMP do schedule(guided)
    do i  = 1, nidx 
        pidx = tape_ind(i, 1) + 1
        qidx = tape_ind(i, 2) + 1
        length = tape_ind(i, 3)
        call Array_Permutation_fort(rank, tape_dim(i, 1:32),&
            tape_ord(i, 1:32), length, data1(pidx:pidx + length), data2(qidx:qidx + length))
    end do
!$OMP END do nowait
!$OMP END PARALLEL
end subroutine permute_player_fort_parallel_guided


subroutine permute_player_fort_parallel_runtime(rank, nidx, tape_ind, tape_dim, tape_ord, totDim, data1, data2)
       integer:: rank, nidx, totDim
       integer:: pidx, qidx, length
       integer:: tape_ind(nidx, 3), tape_dim(nidx, 32), tape_ord(nidx, 32)
       real*8,intent(IN):: data1(totDim)
       real*8,intent(inout)::data2(totDim)
!f2py intent(inout)::data2

!$OMP PARALLEL private(pidx, qidx, length)
!$OMP do schedule(runtime)
    do i  = 1, nidx 
        pidx = tape_ind(i, 1) + 1
        qidx = tape_ind(i, 2) + 1
        length = tape_ind(i, 3)
        call Array_Permutation_fort(rank, tape_dim(i, 1:32),&
            tape_ord(i, 1:32), length, data1(pidx:pidx + length), data2(qidx:qidx + length))
    end do
!$OMP END do 
!$OMP END PARALLEL
end subroutine permute_player_fort_parallel_runtime

