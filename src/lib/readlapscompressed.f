
       subroutine read_laps_compressed(i4time,dir,ext,
     1                    iimax,jjmax,kdim,
     1                    var_req,lvl_req,lvl_coord_req,
     1                    units_req_o,comment_req_o,data_o,istatus)

c**********************************************************************
c
c      subroutine read_laps_compressed
c
c      author:    steve albers
c
c      reads data requested by arrays var_req and lvl_req for the
c      i4time, dir and ext specified.  returns lvl_coord-req,
c      units_req, comment_req, data and istatus
c
c      based on 'read_laps_data' by linda wharton
c
c**********************************************************************
c
c

      integer nf
      parameter (nf=3)

      integer       i4time,               !input i4time of data
     1                iimax,jjmax,          !input # cols, # rows
     1                kdim,                 !input k dimension of data array
     1                lvl_req(kdim*nf),     !input requested levels
     1                istatus               !output
      character*(*)   dir                   !input directory to read data from
      character*(*)   ext                   !input file name ext 
      character*(*)   var_req(kdim)         !input 3 letter id of requested fields
      character*(*)   lvl_coord_req(kdim)   !output vertical coordinate of fields

      character*(*)   units_req_o(kdim)     !output units of requested fields
      character*(10)  units_req_l(kdim,nf)  !local  units of requested fields

      character*(*)   comment_req_o(kdim)   !output comments for requested fields
      character*(125) comment_req_l(kdim,nf)!local  comments for requested fields

      real        data_o(iimax,jjmax,kdim)               !output data
!     real        data_l(iimax,jjmax,kdim,nf)            !local data

      real,    allocatable, dimension(:,:,:,:) :: data_l !local data

      integer, allocatable, dimension(:) :: array1       !local compressed array
      real,    allocatable, dimension(:) :: array2       !local compressed array
c
      integer fn_length,
     1          flag,                   !print flag (1 = off)
     1          error(3),
     1          comm_len,
     1          ext_len
  
c
      character*5       fcst_hh_mm
      character*9       gtime
      character*150     file_name
c
      common            /prt/flag
c
c-------------------------------------------------------------------------------
c
      error(1)=1
      error(2)=0
      error(3)=-2
c
c ****  create file name.
c
      call make_fnam_lp(i4time,gtime,istatus)
      if (istatus .ne. 1) then
        write (6,*) 
     1'error converting i4time to file name...read aborted.'
        istatus=error(2)
        return
      endif

      call s_len(ext, ext_len)

c **** get actual reftime from gtime...
c
      comm_len = len(comment_req_o(1))
c
c **** read in compressed file
c
      lun = 65

      itry = 0 ! retry if file is only partially written out
5     if(.false.)then ! old way
          call open_lapsprd_file_read(lun,i4time,ext,istatus)
          if(istatus .ne. 1)goto 950
      else            ! new way using 'dir'

c fcst_hh_mm: hard wired as a place holder - will be used in filename only
c   if read_laps_compressed is called on lga, lgb, fua, fsf, ram, rsf
c to fix this, call read_laps instead
          fcst_hh_mm = '0000'

          call cvt_fname_v3(dir,gtime,fcst_hh_mm,ext,ext_len,
     1                      file_name,fn_length,istatus)
          if (istatus .eq. error(2)) goto 930

          open(lun,file=file_name(1:fn_length),status='old',err=950)
      endif

10    read(lun,*,err=11,end=11)kkdim
      go to 12
11    write(6,*)' warning: could not read kkdim'
      call sleep(2)
      itry = itry + 1
      if(itry .le. 2)then
          close(lun) 
          write(6,*)' retrying the read'
          go to 5 
      else
          go to 960
      endif

12    if(kkdim .ne. kdim * nf)then
          write(6,*)kkdim,kdim*nf
          go to 980
      endif

      write(6,*)' reading comments, length = ',comm_len
      do l = 1,nf
          do k = 1,kdim
              read(lun,1)comment_req_l(k,l)
              ! temporarily comment this following write statement out for saving output time and file july 2013
              ! steve and yuanfu will add a parameter later controling the output.
              ! write(6,1)comment_req_l(k,l)
1             format(a)
          enddo ! k
      enddo ! l

      read(lun,*)n_cmprs

      allocate(array1(n_cmprs),stat=istat_alloc)
      if(istat_alloc .ne. 0)then
          write(6,*)' error: could not allocate array1'
          stop
      endif

      allocate(array2(n_cmprs),stat=istat_alloc)
      if(istat_alloc .ne. 0)then
          write(6,*)' error: could not allocate array2'
          stop
      endif

      icheck_sum = 0

      do i = 1,n_cmprs
20        read(lun,*,err=21,end=21)array1(i),array2(i)
          go to 22
21        write(6,*)' warning: could not read array1/array2'
          call sleep(5)
          itry = itry + 1
          if(itry .le. 5)then
              close(lun)
              write(6,*)' retrying the read'
              go to 5 
          else
              go to 960
          endif

22        icheck_sum = icheck_sum + array1(i)
      enddo ! i

      close(lun)

      ngrids = iimax*jjmax*kdim*nf

      if(icheck_sum .ne. ngrids)then
          write(6,*)icheck_sum, ngrids
          go to 980
      endif

!     decode the data 
      write(6,*)' decoding the data'

      n_cmprs_max = ngrids

      allocate(data_l(iimax,jjmax,kdim,nf),stat=istat_alloc)
      if(istat_alloc .ne. 0)then
          write(6,*)' error: could not allocate data_l'
          stop
      endif

      call runlength_decode(ngrids,n_cmprs,array1,array2        ! i
     1                     ,data_l                              ! o
     1                     ,istatus)                            ! o
      deallocate(array1)
      deallocate(array2)
      if(istatus .ne. 1)goto 970

!     position the data into the proper 3-d array
      if(var_req(1) .eq. 'ref')then
          ifield = 1
      elseif(var_req(1) .eq. 'vel')then
          ifield = 2
      elseif(var_req(1) .eq. 'nyq')then
          ifield = 3
      else
          write(6,*)' error: unknown variable requested ',var_req
          goto930
      endif

      do k = 1,kdim
!         assign comment
          comment_req_o(k) = comment_req_l(k,ifield)

!         assign level coord?

!         assign units?

!         do i = 1,iimax
!         do j = 1,jjmax
!             data_o(i,j,k) = data_l(i,j,k,ifield)
!         enddo ! j
!         enddo ! i

      enddo ! k      

      data_o(:,:,:) = data_l(:,:,:,ifield)
c
c ****  return normally.
c
        istatus=error(1)
999     if(allocated(data_l))deallocate(data_l)
        return
c
c ****  error trapping.
c
930     if (flag .ne. 1)
     1    write (6,*) 'file_name variable too short...read aborted.'
        istatus=error(2)
        goto 999
c
950     if (flag .ne. 1)
     1    write (6,*) 'error opening compressed file...read aborted.'       
        istatus=error(2)
        goto 999
c
960     if (flag .ne. 1)
     1    write (6,*) 'invalid read as compressed file...read aborted.'
        istatus=error(2)
        goto 999
c
970     if (flag .ne. 1)
     1    write (6,*) 'error decoding data...read aborted.'
        istatus=error(2)
        goto 999
c
980     if (flag .ne. 1)
     1    write (6,*) 'error in array dimensioning...read aborted.'
        istatus=error(2)
        goto 999
c
990     if (flag .ne. 1)
     1    write (6,*) 'error in version, not a valid laps file... '
     1                ,'read aborted.'
        istatus=error(2)
        goto 999
c
992     continue
        istatus=error(2)
        goto 999
c
        end


        subroutine runlength_decode(ngrids,n_cmprs,array1,array2        ! i
     1                     ,data                                        ! o
     1                     ,istatus)                                    ! o

        integer array1(n_cmprs)
        real array2(n_cmprs)
        real data(ngrids)

!       setup for first point
        i_end = 0

        do i = 1,n_cmprs ! loop through array list
            i_start = i_end + 1
            i_end = i_start + array1(i) - 1

            do ii = i_start,i_end
                data(ii) = array2(i)
            enddo ! ii

        enddo ! i

        if(i_end .ne. ngrids)then
            write(6,*)' error in runlength_decode',i_end,ngrids
            istatus = 0
            return
        endif

        istatus = 1
        return
        end
