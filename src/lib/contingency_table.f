
        subroutine contingency_table(o,f,ni,nj,nk             ! i
     1                              ,thresh_o,thresh_f        ! i
     1                              ,lun_out                  ! i
     1                              ,ilow,ihigh,jlow,jhigh    ! i
     1                              ,lmask_rqc_3d             ! i
     1                              ,contable)                ! o

        real o(ni,nj,nk)
        real f(ni,nj,nk)

        logical lmask_rqc_3d(ni,nj,nk)

!       first index is observed, second index is forecast
!       0 is yes, 1 is no
        integer,parameter :: k12 = selected_int_kind(12)
        integer (kind=k12) :: contable(0:1,0:1)

        contable = 0 ! initialize

        do k = 1,nk
        do i = ilow,ihigh
        do j = jlow,jhigh

          if(lmask_rqc_3d(i,j,k))then

            if(o(i,j,k) .ge. thresh_o)then
                index1 = 0
            else
                index1 = 1
            endif

            if(f(i,j,k) .ge. thresh_f)then
                index2 = 0
            else
                index2 = 1
            endif

            contable(index1,index2) = contable(index1,index2) + 1

          endif

        enddo ! j
        enddo ! i
        enddo ! k

!       write contingency table
        write(lun_out,*)
        write(lun_out,1)
        write(lun_out,2)
        write(lun_out,3)contable(0,0),contable(1,0) 
        write(lun_out,4)contable(0,1),contable(1,1) 

 1      format('                    obs')
 2      format('                  y         n     ')
 3      format('  fcst  y',  i11,         i11)
 4      format('        n',  i11,         i11)

        return
        end
