  subroutine smooth2 (nx,ny,dist,data)

    implicit none
    integer, intent(in)  :: nx,ny,dist
    real,    intent(inout) :: data(nx,ny)
    real                   :: tempdata(nx,ny)

    integer i,j,ii,jj,ipt,jpt
    real  :: np, datasum

    tempdata(:,:) = 0.
    do j = 1, ny
      do i = 1, nx
        
        np = 0.
        datasum = 0.
        do jj = -dist,dist,1
          do ii = -dist,dist,1
            ipt = i + ii
            jpt = j + jj
            if ((ipt .ge. 1) .and. (ipt .le. nx) .and. &
                (jpt .ge. 1) .and. (jpt .le. ny) ) then
               datasum = datasum + data(ipt,jpt)
               np = np + 1.
            endif
          enddo
        enddo

        tempdata(i,j) = datasum / np
      enddo
    enddo
    data = tempdata
    end subroutine smooth2
      

