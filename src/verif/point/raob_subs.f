      subroutine get_laps(nx,ny,nk,dir_in,i4time_syn,a9_time,
     1     ulapsgs,vlapsgs,tlapsgs,rhlapsgs,htlapsgs,
     1     ulapsgp,vlapsgp,tlapsgp,rhlapsgp,
     1     htlapsgp,htlgags,htlgagp,
     1     balance,laps_levels,status)

      implicit none
      include 'netcdf.inc'

      integer       nx, ny, nk
      character*(*)   dir_in
      integer       i4time_syn,status(2,6)
      real          ulapsgs(nx,ny,nk), vlapsgs(nx,ny,nk), 
     1                tlapsgs(nx,ny,nk), rhlapsgs(nx,ny,nk),
     1                htlapsgs(nx,ny,nk),ulapsgp(nx,ny,nk), 
     1                vlapsgp(nx,ny,nk),tlapsgp(nx,ny,nk), 
     1                rhlapsgp(nx,ny,nk), htlapsgp(nx,ny,nk),
     1                htlgags(nx,ny,nk), htlgagp(nx,ny,nk)
      real          laps_levels(nk) !laps pressure levels in mb
      real          make_rh

      character*4     lvl_coord(nk)
      character*3     var(nk), ext
      character*10    units(nk)
      character*20    cmdltype
      character*125   comment(nk)
      character*180   fdir_in
      character*9     fnamesyn,fnameprev,fnamepprev,a9_time
      integer       i, j, k, lvl(nk), fdir_len, i4time_pprev,
     1                i4time_prev, len_dir_in, istatus,
     1                balance
      integer         i4time_init,i4time_fcst

! begin
      i4time_prev = i4time_syn - 3600
      call make_fnam_lp(i4time_syn,fnamesyn,istatus)
      call make_fnam_lp(i4time_prev,fnameprev,istatus)

! status(1,x) = syn  (2,x) = prev
! status(x,1)=u (x,2)=v (x,3)=t (x,4)=ht (x,5)=rh (x,6)=lgaht

      status = 1

      call s_len(dir_in,len_dir_in)

      j = 6

      var='ht'
      lvl=laps_levels
c
c js: this needs to consider nest7grid.parms fdda_model_source
c     since that can be our background too. done: 10-11-02 with
c     mod to wind3d.pl
c     ext = 'lga'

      fdir_in = dir_in(1:len_dir_in)//'lga/'
      call s_len(fdir_in,fdir_len)

c     write(6,*) dir_in
c     write(6,*) fdir_in

c js: need to determine what background was actually used in the
c     analysis and configure i4time_syn appropriately (likely will
c     need an i4time_syn_valid too as second argument to this sub)
c     read lga ht valid i4time_syn

      call bkgd_wgi(a9_time,i4time_init,i4time_fcst
     +,ext,cmdltype,balance,istatus)
      if(istatus.ne.1)then
         print*,'failure in bkgd_wgi to get model bkgd time'
         status=0
         return
      endif

      call read_laps(i4time_init,i4time_fcst,fdir_in,ext,
     1                    nx,ny,nk,nk,var,lvl,lvl_coord,
     1                    units,comment,htlgags,istatus)

c     call get_modelfg_3d(i4time_syn,var(1),nx,ny,nk,htlgags
c    1,istatus)

      if (istatus .ne. 1) then

        write (6,*) ' error in readlapsdata for lga ht syn',
     1             fdir_in(1:fdir_len)//fnamesyn//'0000.'//ext(1:3)
      endif

c       try previous hr run, 1 hr fcst
c       call read_laps(i4time_prev,i4time_syn,fdir_in,ext,
c    1                      nx,ny,nk,nk,var,lvl,lvl_coord,
c    1                      units,comment,htlgags,istatus)
c       call get_modelfg_3d(i4time_prev,var(1),nx,ny,nk,htlgags
c    1,istatus)
c       if (istatus .ne. 1) then
c         write (6,*) ' error in readlapsdata for lga ht syn',
c    1            fdir_in(1:fdir_len)//fnameprev//'0100.'//ext(1:3)
c         status(1,j) = 0
c       else
c         write(6,*) 'success reading lga ht from ',
c    1            fdir_in(1:fdir_len)//fnameprev//'0100.'//ext(1:3)
c       endif
c     else
c       write(6,*) 'success reading lga ht from ',
c    1             fdir_in(1:fdir_len)//fnamesyn//'0000.'//ext(1:3)
c     endif
c     read lga ht valid i4time_prev
c js: probably don't need to get the previous hour lga.  if the
c     current one doesn't exist, then we cannot verify. 
c     call read_laps(i4time_prev,i4time_prev,fdir_in,ext,
c    1                    nx,ny,nk,nk,var,lvl,lvl_coord,
c    1                    units,comment,htlgagp,istatus)
c     call get_modelfg_3d(i4time_prev,var(1),nx,ny,nk,htlgagp
c    1,istatus)
c     if (istatus .ne. 1) then
c       write (6,*) ' error in readlapsdata for lga ht prev',
c    1             fdir_in(1:fdir_len)//fnameprev//'0000.'//ext(1:3)
c       try previous hr run, 1 hr fcst
c       i4time_pprev = i4time_prev-3600
c       call make_fnam_lp(i4time_pprev,fnamepprev,istatus)
c       call read_laps(i4time_pprev,i4time_prev,fdir_in,ext,
c    1                      nx,ny,nk,nk,var,lvl,lvl_coord,
c    1                      units,comment,htlgagp,istatus)
c       call get_modelfg_3d(i4time_pprev,var(1),nx,ny,nk,htlgagp
c    1,istatus)
c       if (istatus .ne. 1) then
c         write (6,*) ' error in readlapsdata for lga ht pprev',
c    1            fdir_in(1:fdir_len)//fnamepprev//'0100.'//ext(1:3)
c         status(2,j) = 0
c       else
c         write(6,*) 'success reading lga ht from ',
c    1            fdir_in(1:fdir_len)//fnamepprev//'0100.'//ext(1:3)
c       endif
c     else
c       write(6,*) 'success reading lga ht from ',
c    1             fdir_in(1:fdir_len)//fnameprev//'0000.'//ext(1:3)
c     endif

c js: new code for verifying background with sounding or profiler
c   we'll use balance code = 2 to indicate verification of background.

      if(balance.eq.2)then

         ext = 'lga'
         var = 't3'
         call read_laps(i4time_init,i4time_fcst,fdir_in,ext,
     1                    nx,ny,nk,nk,var,lvl,lvl_coord,
     1                    units,comment,tlapsgs,istatus)
         if(istatus.ne.1)then
            write (6,*) ' error in readlapsdata for lga ht pprev',
     1fdir_in(1:fdir_len)//fnamepprev//'0100.'//ext(1:3)
            status(1,3)=0
            return
         endif
         var='ht'
         call read_laps(i4time_init,i4time_fcst,fdir_in,ext,
     1                    nx,ny,nk,nk,var,lvl,lvl_coord,
     1                    units,comment,htlapsgs,istatus)
         if(istatus.ne.1)then
            write (6,*) ' error in readlapsdata for lga ht pprev',
     1fdir_in(1:fdir_len)//fnamepprev//'0100.'//ext(1:3)
            status(1,4)=0
            return
         endif
         var='u3'
         call read_laps(i4time_init,i4time_fcst,fdir_in,ext,
     1                    nx,ny,nk,nk,var,lvl,lvl_coord,
     1                    units,comment,ulapsgs,istatus)
         if(istatus.ne.1)then
            write (6,*) ' error in readlapsdata for lga ht pprev',
     1fdir_in(1:fdir_len)//fnamepprev//'0100.'//ext(1:3)
            status(1,1)=0
            return
         endif
         var='v3'
         call read_laps(i4time_init,i4time_fcst,fdir_in,ext,
     1                    nx,ny,nk,nk,var,lvl,lvl_coord,
     1                    units,comment,vlapsgs,istatus)
         if(istatus.ne.1)then
            write (6,*) ' error in readlapsdata for lga ht pprev',
     1fdir_in(1:fdir_len)//fnamepprev//'0100.'//ext(1:3)
            status(1,2)=0
            return
         endif
         var='sh'
         call read_laps(i4time_init,i4time_fcst,fdir_in,ext,
     1                    nx,ny,nk,nk,var,lvl,lvl_coord,
     1                    units,comment,rhlapsgs,istatus)
         if(istatus.ne.1)then
            write (6,*) ' error in readlapsdata for lga ht pprev',
     1fdir_in(1:fdir_len)//fnamepprev//'0100.'//ext(1:3)
            status(1,5)=0
            return
         endif

         do k=1,nk
          do j=1,ny
           do i=1,nx
              rhlapsgs(i,j,k)=make_rh(laps_levels(k)
     1,tlapsgs(i,j,k)-273.16,rhlapsgs(i,j,k)*1000.,-132.)*100.
           enddo
          enddo
         enddo

      else

      j = 1
      var = 'u3'
      ext = 'lw3'

      if (balance .eq. 1) then
        fdir_in = dir_in(1:len_dir_in)//'balance/lw3/'
      elseif(balance .eq. 0)then
        fdir_in = dir_in(1:len_dir_in)//'lw3/'
      endif

      call s_len(fdir_in,fdir_len)

      call read_laps_data(i4time_syn,fdir_in,ext,nx,ny,nk,
     1                    nk,var,lvl,lvl_coord,
     1                    units,comment,ulapsgs,istatus)

      if (istatus .ne. 1) then
        write (6,*) ' error in readlapsdata for lw3 u syn'
        status(1,j) = 0
      else
        write(6,*) 'success reading u3 from ',
     1             fdir_in(1:fdir_len)//fnamesyn//'.'//ext(1:3)
      endif

      call read_laps_data(i4time_prev,fdir_in,ext,nx,ny,nk,
     1                    nk,var,lvl,lvl_coord,
     1                    units,comment,ulapsgp,istatus)

      if (istatus .ne. 1) then
        write (6,*) ' error in readlapsdata for lw3 u prev'
        status(2,j) = 0
      else
        write(6,*) 'success reading u3 from ',
     1             fdir_in(1:fdir_len)//fnameprev//'.'//ext(1:3)
      endif

      j = 2
      do i = 1, nk
        var(i) = 'v3'
      enddo

      call read_laps_data(i4time_syn,fdir_in,ext,nx,ny,nk,
     1                    nk,var,lvl,lvl_coord,
     1                    units,comment,vlapsgs,istatus)

      if (istatus .ne. 1) then
        write (6,*) ' error in readlapsdata for lw3 v syn'
        status(1,j) = 0
      else
        write(6,*) 'success reading v3 from ',
     1             fdir_in(1:fdir_len)//fnamesyn//'.'//ext(1:3)
      endif

      call read_laps_data(i4time_prev,fdir_in,ext,nx,ny,nk,
     1                    nk,var,lvl,lvl_coord,
     1                    units,comment,vlapsgp,istatus)

      if (istatus .ne. 1) then
        write (6,*) ' error in readlapsdata for lw3 v prev'
        status(2,j) = 0
      else
        write(6,*) 'success reading v3 from ',
     1             fdir_in(1:fdir_len)//fnameprev//'.'//ext(1:3)
      endif

      ext = 'lt1'

      j = 3
      do i = 1, nk
        var(i) = 't3'
      enddo
      if (balance .eq. 1) then
        fdir_in = dir_in(1:len_dir_in)//'balance/lt1/'
      else
        fdir_in = dir_in(1:len_dir_in)//'lt1/'
      endif
      call s_len(fdir_in,fdir_len)

      call read_laps_data(i4time_syn,fdir_in,ext,nx,ny,nk,
     1                    nk,var,lvl,lvl_coord,
     1                    units,comment,tlapsgs,istatus)

      if (istatus .ne. 1) then
        write (6,*) ' error in readlapsdata for lt1 t3 syn '
        status(1,j) = 0
        status = 0
      else
        write(6,*) 'success reading t3 from ',
     1             fdir_in(1:fdir_len)//fnamesyn//'.'//ext(1:3)
      endif

      call read_laps_data(i4time_prev,fdir_in,ext,nx,ny,nk,
     1                    nk,var,lvl,lvl_coord,
     1                    units,comment,tlapsgp,istatus)

      if (istatus .ne. 1) then
        write (6,*) ' error in readlapsdata for lt1 t3 prev '
        status(2,j) = 0
      else
        write(6,*) 'success reading t3 from ',
     1             fdir_in(1:fdir_len)//fnameprev//'.'//ext(1:3)
      endif

      j = 4
      do i = 1, nk
        var(i) = 'ht'
      enddo

      call read_laps_data(i4time_syn,fdir_in,ext,nx,ny,nk,
     1                    nk,var,lvl,lvl_coord,
     1                    units,comment,htlapsgs,istatus)

      if (istatus .ne. 1) then
        write (6,*) ' error in readlapsdata for lt1 ht syn'
        status(1,j) = 0
      else
        write(6,*) 'success reading ht from ',
     1             fdir_in(1:fdir_len)//fnamesyn//'.'//ext(1:3)
      endif

      call read_laps_data(i4time_prev,fdir_in,ext,nx,ny,nk,
     1                    nk,var,lvl,lvl_coord,
     1                    units,comment,htlapsgp,istatus)

      if (istatus .ne. 1) then
        write (6,*) ' error in readlapsdata for lt1 ht prev'
        status(2,j) = 0
      else
        write(6,*) 'success reading ht from ',
     1             fdir_in(1:fdir_len)//fnameprev//'.'//ext(1:3)
      endif

      ext = 'lh3'

      j = 5
c     do i = 1, nk
c       var(i) = 'rh3'
c     enddo

      var = 'rhl'

      if (balance .eq. 1) then
        fdir_in = dir_in(1:len_dir_in)//'balance/lh3/'
      else
        fdir_in = dir_in(1:len_dir_in)//'lh3/'
      endif
      call s_len(fdir_in,fdir_len)

      call read_laps_data(i4time_syn,fdir_in,ext,nx,ny,nk,
     1                    nk,var,lvl,lvl_coord,
     1                    units,comment,rhlapsgs,istatus)

      if (istatus .ne. 1) then
        write (6,*) ' error in readlapsdata for lh3 rh3 syn'
        status(1,j) = 0
      else
        write(6,*) 'success reading rh from ',
     1             fdir_in(1:fdir_len)//fnamesyn//'.'//ext(1:3)
      endif

      call read_laps_data(i4time_prev,fdir_in,ext,nx,ny,nk,
     1                    nk,var,lvl,lvl_coord,
     1                    units,comment,rhlapsgp,istatus)

      if (istatus .ne. 1) then
        write (6,*) ' error in readlapsdata for lh3 rh3 prev'
        status(2,j) = 0
      else
        write(6,*) 'success reading rh from ',
     1             fdir_in(1:fdir_len)//fnameprev//'.'//ext(1:3)
      endif

      endif

      return
      end
!1.....................................................................................
      subroutine ascend_w(timesyn, timerel, numsigw, numw,fileavail,
     1                    wmostanum, stalat, stalon, staelev, maxw,
     1                    max_ht_m_proc,typew,
     1                    maxraob, nx, ny, nk, raobid, statusl,
     1                    htsigw, wdsigw, wssigw,  !remember (maxw,maxraob)
     1                    ulapsgs,vlapsgs,ulapsgp,vlapsgp,  !(nx,ny,nk)
     1                    htlapsgs, htlapsgp, raob_missing_data,
     1                    verif_missing_data,
     1                    lat, lon, laps_levels,  !(nx,ny)
!................................variables below returned...................
     1                    riw,rjw,rkw,latw,lonw,htw, uiw,viw,
     1                    upw,vpw,timelapsw,status)

      include 'trigd.inc'

      implicit none

      real          height_to_zcoord3
      integer       timesyn, timerel, numsigw, numw, fileavail,
     1                wmostanum
      real          stalat, stalon, staelev,max_ht_m_proc
      character*1     typew(maxw,maxraob),typewm(maxw,maxraob)
      integer       maxw, maxraob,nx,ny,nk, raobid,statusl(2,6)
      real          htsigw(maxw,maxraob),
     1                wdsigw(maxw,maxraob), wssigw(maxw,maxraob),
     1                ulapsgs(nx,ny,nk), vlapsgs(nx,ny,nk),
     1                ulapsgp(nx,ny,nk), vlapsgp(nx,ny,nk),
     1                htlapsgs(nx,ny,nk), htlapsgp(nx,ny,nk),
     1                raob_missing_data, verif_missing_data,
     1                lat(nx,ny), lon(nx,ny),
     1                laps_levels(nk), 
     1                riw(maxw,maxraob),rjw(maxw,maxraob),
     1                rkw(maxw,maxraob),latw(maxw,maxraob),
     1                lonw(maxw,maxraob),htw(maxw,maxraob),
     1		      uiw(maxw,maxraob),viw(maxw,maxraob),
     1		      upw(maxw,maxraob),vpw(maxw,maxraob)
      integer       timelapsw(maxw,maxraob), status

      integer       j,lvl,timeprev,index, time_tot, istatus,
     1                timeprevcomp, timesyncomp,int_ri,int_rj
      real          prsigw, usigw(maxw),vsigw(maxw),
     1                u_disp, v_disp, delta_h,delta_t, delta_u,
     1                delta_v, rise_rate, maxhtproc, rilaps, 
     1                rjlaps, rklaps, u_laps, v_laps, delta_s

! begin
      status = 1   !assume a good return
      timeprev = timesyn - 3600
      timesyncomp = timesyn
      timeprevcomp = timeprev
      if (fileavail .eq. 1) timeprevcomp = 0
      if (fileavail .eq. 2) timesyncomp = 0

      if ((timeprevcomp .eq. 0) .and. (timesyncomp .eq. 0)) then
        write(6,*) ' no wind data to process '
        istatus = 0
        return
      endif

      maxhtproc = max_ht_m_proc  !max htsigw value to process (in meters)

!     setup for ascention
      index = 0
      u_disp = 0.
      v_disp = 0.
      time_tot = timerel
      rise_rate = 300. / 60. ! m/s

!     loop thru levels
      numw = -1
      lvl = 1
      do while (lvl .le. numsigw)
        if ((htsigw(lvl,raobid) .gt. 0) .and.
     1      (htsigw(lvl,raobid) .le. maxhtproc) .and.
     1      (htsigw(lvl,raobid) .ne. raob_missing_data)) then
          index = index + 1

          if (index .ge. 1) then
            if(index .gt. 1)then ! lower level is above 1st (sfc) level
              delta_h = htsigw(lvl,raobid) - htsigw(lvl-1,raobid)
            else                 ! lower level is at the sfc
              delta_h = htsigw(lvl,raobid) - staelev
            endif

            htw(index,raobid) = htsigw(lvl,raobid)
            typewm(index,raobid) = typew(index,raobid)

            if ((wdsigw(lvl,raobid) .eq. raob_missing_data) .or.
     1          (wssigw(lvl,raobid) .eq. raob_missing_data)) then
              index = index - 1  ! don't use this lvl ob...go to next
              lvl = lvl + 1
              goto 777
            else
              call disp_to_uv(wdsigw(lvl,raobid),wssigw(lvl,raobid),
     1                        upw(index,raobid),vpw(index,raobid))
            endif

            delta_t = delta_h / rise_rate
            time_tot = time_tot + nint(delta_t + 0.5)
            if (abs(time_tot - timesyncomp) .gt. 
     1          abs(time_tot - timeprevcomp)) then
              timelapsw(index,raobid) = timeprev
            else
              timelapsw(index,raobid) = timesyn
            endif
            delta_s = wssigw(lvl,raobid) * delta_t
            delta_u = -sind(wdsigw(lvl,raobid)) * delta_s
            delta_v = -cosd(wdsigw(lvl,raobid)) * delta_s
            u_disp = u_disp + delta_u
            v_disp = v_disp + delta_v
            latw(index,raobid) = stalat + v_disp / 111000.
            lonw(index,raobid) = stalon + u_disp / 
     1          (111000.* cosd((stalat+latw(index,raobid))/2))

!           determine laps i,j,k, calc deltau and deltav
!           index has number of elements to correlate

            call latlon_to_rlapsgrid(latw(index,raobid),
     1                               lonw(index,raobid),
     1                               lat,lon,nx,ny,
     1                               riw(index,raobid),
     1                               rjw(index,raobid),istatus)

            if (istatus .ne. 1) then  !if this level is out, rest probably will be too
              numw = index -1
              lvl = numsigw + 1
              if (index .gt. 1) then  !keep what have
                goto 777
              else
                write(6,*) 'raob ',wmostanum,' not in laps domain  ',
     1                   latw(index,raobid),' ',lonw(index,raobid)
                goto 777
              endif
            endif
             
c           get integral ri, rj for height_to_zcoord3
            int_ri = nint(riw(index,raobid))
            int_rj = nint(rjw(index,raobid))

            if (timelapsw(index,raobid) .eq. timesyn) then   !use htlapsgs
              rkw(index,raobid) = height_to_zcoord3(htsigw(lvl,raobid),
     1                                       htlapsgs,laps_levels,
     1                                       nx,ny,nk,int_ri,
     1                                       int_rj,istatus)
            else  !use htlapsgp
              rkw(index,raobid) = height_to_zcoord3(htsigw(lvl,raobid),
     1                                       htlapsgp,laps_levels,
     1                                       nx,ny,nk,int_ri,
     1                                       int_rj,istatus)
            endif

            if ((istatus .eq. 1) .and. (rkw(index,raobid) .le. nk)) then

!             determine which laps grid to use
              if (timelapsw(index,raobid) .eq. timesyn) then
                call trilinear_laps(riw(index,raobid),rjw(index,raobid),
     1                              rkw(index,raobid),nx,ny,nk,
     1                              ulapsgs, uiw(index,raobid))

                call trilinear_laps(riw(index,raobid),rjw(index,raobid),
     1                              rkw(index,raobid),nx,ny,nk,
     1                              vlapsgs, viw(index,raobid))
              else
                call trilinear_laps(riw(index,raobid),rjw(index,raobid),
     1                              rkw(index,raobid),nx,ny,nk,
     1                              ulapsgp, uiw(index,raobid))

                call trilinear_laps(riw(index,raobid),rjw(index,raobid),
     1                              rkw(index,raobid),nx,ny,nk,
     1                              vlapsgp, viw(index,raobid))
              endif

            else  !height is above laps domain
              numw = index - 1
              lvl = numsigw + 1
              if (index .gt. 1) then  !save what have
                goto 777
              else
                write(6,*) 'raob ',wmostanum, ' above laps domain  ',
     1                     rkw(index,raobid)
                goto 777
              endif
            endif

          endif
        endif  !if data for lvl not missing

        lvl = lvl + 1
777     continue

      enddo

      if (numw .eq. -1) numw = index

c     pass back modified typew
      do j = 1, numw
        typew(j,raobid) = typewm(j,raobid)
      enddo
 
      return
      end
!2.....................................................................................
      subroutine ascend_t(timesyn, timerel, numsigt,
     1                    numt,fileavailuv,fileavail,wmostanum,
     1                    stalat, stalon, staelev, typet,
     1                    maxt, maxw, max_ht_m_proc,
     1                    maxraob, nx, ny, nk, raobid, statusl,
     1                    prsigt, tsigt, tdsigt,htsigt,
     1                    tlapsgs,rhlapsgs,htlapsgs,tlapsgp,  !(nx,ny,nk)
     1                    rhlapsgp, htlapsgp,
     1                    htlgags, htlgagp, raob_missing_data,
     1                    verif_missing_data, lat,lon, !(nx,ny,nk)
     1                    laps_levels, numsigw,
     1                    htsigw, wdsigw, wssigw, !remember (maxw,maxraob)
!................................variables below returned...................
     1                    rit,rjt,rkt,latt,lont,htt,prit,tit,tdit,
     1                    prpt,tpt,tdpt,timelapst,status)

      include 'trigd.inc'

      implicit none

      real          ztopsa, zcoord_of_logpressure, psatoz,
     1                make_ssh, k_to_c, make_td, c_to_k

      integer       timesyn, timerel, numsigt, numt,
     1                fileavailuv, fileavail,wmostanum
      real          stalat, stalon, staelev,max_ht_m_proc
      character*1     typet(maxt,maxraob),typetm(maxt,maxraob)
      integer       maxt, maxw, maxraob,nx,ny,nk
      integer       raobid, statusl(2,6)
      real          prsigt(maxt,maxraob),
     1                tsigt(maxt,maxraob), 
     1                tdsigt(maxt,maxraob), 
     1                htsigt(maxt,maxraob), 
     1                tlapsgs(nx,ny,nk), htlapsgs(nx,ny,nk),
     1                tlapsgp(nx,ny,nk), htlapsgp(nx,ny,nk),
     1                rhlapsgs(nx,ny,nk), rhlapsgp(nx,ny,nk),
     1                htlgags(nx,ny,nk), htlgagp(nx,ny,nk),
     1                raob_missing_data,verif_missing_data,
     1                lat(nx,ny), lon(nx,ny), laps_levels(nk)
      integer	      numsigw
      real          htsigw(maxw,maxraob), wdsigw(maxw,maxraob),
     1                wssigw(maxw,maxraob),
     1                rit(maxt,maxraob),rjt(maxt,maxraob),
     1                rkt(maxt,maxraob),latt(maxt,maxraob),
     1                lont(maxt,maxraob),htt(maxt,maxraob),
     1		      prit(maxt,maxraob),tit(maxt,maxraob),
     1		      tdit(maxt,maxraob),prpt(maxt,maxraob),
     1                tpt(maxt,maxraob),tdpt(maxt,maxraob)
      integer       timelapst(maxt,maxraob), status

      integer       timeprevcomp, timesyncomp,istatus,
     1                j,lvl,timeprev,index, time_tot,
     1                bkgprev, bkgsyn,int_ri,int_rj
      real          wd, ws, rh, q, t_ref, tc,
     1                u_disp, v_disp,delta_h,delta_t,delta_u,
     1                delta_v, rise_rate, maxhtproc,htm,htz, 
     1                t_laps, ht, htprev,delta_s,pres_pa

! begin
      status = 1   !assume a good return

      timeprev = timesyn - 3600
      timesyncomp = timesyn
      timeprevcomp = timeprev
      if ((fileavail .eq. 1) .or. (fileavailuv .eq. 1)) timeprevcomp = 0
      if ((fileavail .eq. 2) .or. (fileavailuv .eq. 2)) timesyncomp = 0
      if (statusl(1,6) .eq. 1) then
        bkgsyn = 1 
      else
        bkgsyn = 0 
      endif
      if (statusl(2,6) .eq. 1) then
        bkgprev = 1 
      else
        bkgprev = 0 
      endif

      if ((bkgsyn .eq. 0) .and. (bkgprev .eq. 0)) then
        write(6,*) 'no lga background available...',
     1             'using psatoz for pressure/height calcs for temp'
      endif
      if ((timeprevcomp .eq. 0) .and. (timesyncomp .eq. 0)) then
        write(6,*) ' no temp/wind data to process '
        status = 0
        return
      endif

      maxhtproc = max_ht_m_proc  !max htsigw value to process (in meters)

!     setup for ascention
      index = 0
      u_disp = 0.
      v_disp = 0.
      time_tot = timerel
      rise_rate = 300. / 60. ! m/s
      htprev = staelev

!     loop thru levels
      numt = -1
      lvl = 1

      do while (lvl .le. numsigt)

c       write(6,*) 'prsigt(lvl,raobid)=',prsigt(lvl,raobid), lvl, raobid

        if (htsigt(lvl,raobid) .ne. verif_missing_data) then 
          write(6,*) 'height from htman used'
          ht = htsigt(lvl,raobid)
        else
          if (prsigt(lvl,raobid) .ne. raob_missing_data) then
c           quick and dirty for now...need ri,rj,timelaps for pressure_to_height
            ht = psatoz(prsigt(lvl,raobid))
          else
            ht = raob_missing_data
          endif
        endif

c       write(6,*) 'height from ztopsa=',ht, prsigt(lvl,raobid)

c       stop below htsigw(max) because of ascend
        if ((ht .gt. 0) .and. (ht .le. maxhtproc) .and.
     1      (ht .ne. raob_missing_data) .and.
     1      (ht .le. htsigw(numsigw,raobid))) then

          call interpwind2t(maxw, maxraob, numsigw, wdsigw, 
     1                      wssigw, htsigw, raobid, maxhtproc, 
     1                      ht, ws, wd, istatus)
          if (istatus .eq.1) then
            index = index + 1

            if (index .ge. 1) then

c             fill in prpt, tpt, tdpt from sigt vars
              prpt(index,raobid) = prsigt(lvl,raobid)
              tpt(index,raobid) = tsigt(lvl,raobid)
              tdpt(index,raobid) = tdsigt(lvl,raobid)
              typetm(index,raobid) = typet(lvl,raobid)

              delta_h = ht - htprev

              delta_t = delta_h / rise_rate
              time_tot = time_tot + nint(delta_t + 0.5)
              if (abs(time_tot - timesyncomp) .gt. 
     1            abs(time_tot - timeprevcomp)) then
                timelapst(index,raobid) = timeprev
              else
                timelapst(index,raobid) = timesyn
              endif

              delta_s = ws * delta_t
              delta_u = -sind(wd) * delta_s
              delta_v = -cosd(wd) * delta_s
              u_disp = u_disp + delta_u
              v_disp = v_disp + delta_v
              latt(index,raobid) = stalat + v_disp / 111000.
              lont(index,raobid) = stalon + u_disp / 
     1                 (111000.* cosd((stalat+latt(index,raobid))/2))
              htprev = ht

!             determine laps i,j,k

              call latlon_to_rlapsgrid(latt(index,raobid),
     1                                 lont(index,raobid),
     1                                 lat,lon,nx,ny,
     1                                 rit(index,raobid),
     1                                 rjt(index,raobid),
     1                                 istatus)

              if (istatus .ne. 1) then  !out of laps domain...rest probably too
                 numt = index - 1
                 lvl = numsigt + 1
                if (index .gt. 1) then  !save what have
                  goto 777
                else
                  write(6,*) 'raob ',wmostanum,' not in laps domain  ',
     1                      latt(index,raobid),' ',lont(index,raobid)
                  goto 777
                endif
              endif

c             re-calc height more accurately now that have ri, rj
              if (htsigt(lvl,raobid) .eq. verif_missing_data) then  !not man height...recalc
              else
                htm = ht
              endif

                pres_pa = prsigt(lvl,raobid) * 100
                int_ri = nint(rit(index,raobid))
                int_rj = nint(rjt(index,raobid))

                if ((bkgsyn.eq.1).and.(bkgprev.eq.1)) then  !use bkg closest to timelaps(index,raobid)
                  if (timelapst(index,raobid) .eq. timesyn) then
                    call pressure_to_height(pres_pa,htlgags,nx,ny,
     1                                      nk,int_ri,int_rj,ht,istatus)
                  else
                    call pressure_to_height(pres_pa,htlgagp,nx,ny,
     1                                      nk,int_ri,int_rj,ht,istatus)
                  endif
                else
                  if (bkgsyn .eq. 1) then  ! use htlgags
                    call pressure_to_height(pres_pa,htlgags,nx,ny,
     1                                      nk,int_ri,int_rj,ht,istatus)
                  else
                    if (bkgprev .eq. 1) then  ! use htlgagp
                      call pressure_to_height(pres_pa,htlgagp,nx,ny,
     1                                      nk,int_ri,int_rj,ht,istatus)
                    else
                      ! leave calc from psatoz
                    endif
                  endif
                endif
c             endif

              if (htsigt(lvl,raobid) .eq. verif_missing_data) then  !not man height...recalc
              htt(index,raobid) = ht
              else
              htt(index,raobid) = htm
              htz = psatoz(prsigt(lvl,raobid))
              write(6,*) 'ht=',ht,'  htm=',htm,'  htz=',htz

              endif
              rkt(index,raobid) = 
     1          zcoord_of_logpressure(prsigt(lvl,raobid)*100.)

              if (rkt(index,raobid) .le. nk) then

!               determine which laps grid to use and interpolate t in the vertical 
                if (timelapst(index,raobid) .eq. timesyn) then

                  call trilinear_laps(rit(index,raobid),
     1                                rjt(index,raobid),
     1                                rkt(index,raobid),nx,ny,nk,
     1                                tlapsgs, tit(index,raobid))

                  call trilinear_laps(rit(index,raobid),
     1                                rjt(index,raobid),
     1                                rkt(index,raobid),nx,ny,nk,
     1                                htlapsgs, ht)

                  call trilinear_laps(rit(index,raobid),
     1                                rjt(index,raobid),
     1                                rkt(index,raobid),nx,ny,nk,
     1                                rhlapsgs, rh)

                  rh = rh/100.0
                  prit(index,raobid) = ztopsa(ht)
                  t_ref = -132.0
                  tc = k_to_c(tit(index,raobid))
                  q = make_ssh(prit(index,raobid),tc,rh,t_ref) 
                  tdit(index,raobid) = 
     1              make_td(prit(index,raobid),tc,q,t_ref)
                  tdit(index,raobid) = c_to_k(tdit(index,raobid))

                else

                  call trilinear_laps(rit(index,raobid),
     1                                rjt(index,raobid),
     1                                rkt(index,raobid),nx,ny,nk,
     1                                tlapsgp, tit(index,raobid))

                  call trilinear_laps(rit(index,raobid),
     1                                rjt(index,raobid),
     1                                rkt(index,raobid),nx,ny,nk,
     1                                htlapsgs, ht)

                  call trilinear_laps(rit(index,raobid),
     1                                rjt(index,raobid),
     1                                rkt(index,raobid),nx,ny,nk,
     1                                rhlapsgp, rh)

                  prit(index,raobid) = ztopsa(ht)
                  t_ref = -132.0
                  tc = k_to_c(tit(index,raobid))
                  q = make_ssh(prit(index,raobid),tc,rh,t_ref) 
                  tdit(index,raobid) = 
     1               make_td(prit(index,raobid),tc,q,t_ref)

                endif

!               check tit(index,raobid) for validity
                if (tit(index,raobid) .ge. 350) then
                  index = index - 1
                endif
              else
                numt = index - 1
                lvl = numsigt + 1
                if (index .gt. 1) then  !save what have
                  goto 777
                else
                 write(6,*) 'raob ',wmostanum,' is above laps domain  ',
     1                       rkt(index,raobid)
                  goto 777
                endif
              endif

            endif
          endif
        endif
 
        lvl = lvl + 1
777     continue

      enddo

      if (numt .eq. -1) numt = index

c     pass back modified typet
      do j = 1, numt
        typet(j,raobid) = typetm(j,raobid)
      enddo
 
      return
      end
!3.....................................................................................
      subroutine interpwind2t(maxw, maxraob, numsigw, wdsigw, 
     1                        wssigw, htsigw, raobid, maxhtproc, 
     1                        ht, ws, wd, istatus)

      implicit none
      integer    maxw, maxraob, numsigw, raobid, maxhtproc, 
     1             istatus, u, l, i, true, false, found, last
      real       htsigw(maxw,maxraob), wdsigw(maxw,maxraob),
     1             wssigw(maxw,maxraob), ht, ws, wd, deltaws

! begin
      istatus = 1
      true = 1
      false = 0

!     we know that ht < maxhtproc and ht <= htsigw(numsigw,raobid)
        l = 1
      if (htsigw(1,raobid) .eq. 0.0) then
        if (numsigw .ge. 3) then
          l = 2
          u = 3
        else
          istatus = 0
          return
        endif
      elseif (numsigw .ge. 2) then
        l = 1
        u = 2
      endif

      if (ht .lt. htsigw(l,raobid)) then
        istatus = 0
        return
      endif

!     we know that ht >= htsigw(l,raobid) and ht <= htsigw(numsigw,raobid)

      found = false
      last = false
      do while ((u .le. numsigw) .and. (found .eq. false) .and.
     1          (last .eq. false)) 
        if ((ht .ge. htsigw(l,raobid)) .and.
     1      (ht .le. htsigw(u,raobid))) then
          found = true
        elseif (u .eq. numsigw) then
          last = true
        elseif (ht .gt. htsigw(u,raobid)) then
          l = l + 1
          u = u + 1
        endif
      enddo

      if (found .eq. true) then
        ws = wssigw(l,raobid) + ((ht-htsigw(l,raobid)) *
     1       ((wssigw(u,raobid)-wssigw(l,raobid)) /
     1        (htsigw(u,raobid)-htsigw(l,raobid))))
        wd = wdsigw(l,raobid) + ((ht-htsigw(l,raobid)) *
     1       ((wdsigw(u,raobid)-wdsigw(l,raobid)) /
     1        (htsigw(u,raobid)-htsigw(l,raobid))))
        istatus = 1
      else
        istatus = 0
      endif

      return
      end
!4....................................................................................
      subroutine mergetw(max_hts,max_raobs,maxw,maxt,numw,numt,
     1                   nraobs,use_raob,n_raobs_use,
     1                   ri,rj,rk,lat,lon,hts,type,up,ui,vp,vi,
     1                   tp,ti,tdp,tdi,prp,pri,timelaps,nhts,
     1                   rit,rjt,rkt,latt,lont,htt,typet,
     1                   riw,rjw,rkw,latw,lonw,htw,typew,
     1                   uiw,viw,upw,vpw,timelapsw,
     1                   prit,tit,tdit,prpt,tpt,tdpt,timelapst,
     1                   verif_missing_data, raob_missing_data,
     1                   istatus)

      implicit none

      integer		max_hts,max_raobs,maxw,maxt,
     1                  nraobs,use_raob(max_raobs),
     1                  n_raobs_use

c     data written out
      real            ri(max_hts,max_raobs),
     1                  rj(max_hts,max_raobs),
     1                  rk(max_hts,max_raobs),
     1                  lat(max_hts,max_raobs),
     1                  lon(max_hts,max_raobs),
     1                  hts(max_hts,max_raobs)
      character*1       type(max_hts,max_raobs)
      real            up(max_hts,max_raobs),
     1                  ui(max_hts,max_raobs),
     1                  vp(max_hts,max_raobs),
     1                  vi(max_hts,max_raobs),
     1                  tp(max_hts,max_raobs),
     1                  ti(max_hts,max_raobs),
     1                  tdp(max_hts,max_raobs),
     1                  tdi(max_hts,max_raobs),
     1                  prp(max_hts,max_raobs),
     1                  pri(max_hts,max_raobs)
      integer         timelaps(max_hts,max_raobs),
     1                  nhts(max_raobs)

c     temp variables
      integer         numt(max_raobs),
     1			timelapst(maxt,max_raobs)
      character*1       typet(maxt,max_raobs)
      real            rit(maxt,max_raobs),
     1                  rjt(maxt,max_raobs),
     1                  rkt(maxt,max_raobs),
     1                  latt(maxt,max_raobs),
     1                  lont(maxt,max_raobs),
     1                  htt(maxt,max_raobs),
     1                  prit(maxt,max_raobs),
     1                  tit(maxt,max_raobs),
     1                  tdit(maxt,max_raobs),
     1                  prpt(maxt,max_raobs),
     1                  tpt(maxt,max_raobs),
     1                  tdpt(maxt,max_raobs)

c     wind variables
      integer         numw(max_raobs),
     1 			timelapsw(maxw,max_raobs)
      character*1       typew(maxw,max_raobs)
      real            riw(maxw,max_raobs),
     1                  rjw(maxw,max_raobs),
     1                  rkw(maxw,max_raobs),
     1                  latw(maxw,max_raobs),
     1                  lonw(maxw,max_raobs),
     1                  htw(maxw,max_raobs),
     1                  uiw(maxw,max_raobs),
     1                  viw(maxw,max_raobs),
     1                  upw(maxw,max_raobs),
     1                  vpw(maxw,max_raobs)

      real		verif_missing_data,
     1                  raob_missing_data
      integer		istatus

c     local variables

      integer		tptr,wptr,index,i,j,ndiff
      real		sumdiff,avgdiff,pctdiff,
     1                  sumpctdiff,avgpctdiff,
     1                  rkdiff

c
c     begin
c
      istatus = 1   !assume good return
      ndiff = 0
      sumdiff = 0

c     loop through raobs
      do i = 1, nraobs
        if (use_raob(i) .eq. 1) then

          write(6,*) 'raob ',i,' before merge'
          write(6,*) numw(i)
          write(6,*) 'htw, typew, upw, vpw'
          do j = 1, numw(i)
            write(6,*) j,htw(j,i),typew(j,i),upw(j,i),
     1                 vpw(j,i)
          enddo

          write(6,*) numt(i)
          write(6,*) 'htt, typet, prpt, tpt, tdpt'
          do j = 1, numt(i)
            write(6,*) j,htt(j,i),typet(j,i),prpt(j,i),
     1                 tpt(j,i),tdpt(j,i)
          enddo

c         setup pointers
          tptr = 1
          wptr = 1
          index = 1

c         both htt and htw should have no missing data...if it does, don't use it
c         drive by sigt because it has more levels

          do while (tptr .le. numt(i))
            if (wptr .le. numw(i)) then

c             verify ri,rj,rk,lat,lon,ht are valid

              if (htw(wptr,i) .lt. htt(tptr,i)) then  ! w ob
                hts(index,i) = htw(wptr,i)
                ri(index,i) = riw(wptr,i)
                rj(index,i) = rjw(wptr,i)
                rk(index,i) = rkw(wptr,i)
                lat(index,i) = latw(wptr,i)
                lon(index,i) = lonw(wptr,i)
                up(index,i) = upw(wptr,i)
                ui(index,i) = uiw(wptr,i)
                vp(index,i) = vpw(wptr,i)
                vi(index,i) = viw(wptr,i)
                tp(index,i) = verif_missing_data
                ti(index,i) = verif_missing_data
                tdp(index,i) = verif_missing_data
                tdi(index,i) = verif_missing_data
                prp(index,i) = verif_missing_data
                pri(index,i) = verif_missing_data
                type(index,i) = 'w'
                if (typew(wptr,i) .eq. 'm') then 
                  write(6,*)'raob ',i,' has mismatch man htw.lt.htt',
     1                      wptr,' ',tptr,' ',htw(wptr,i),' ',
     1                      htt(tptr,i)
                endif
                wptr = wptr + 1
                index = index + 1
              else
                if (htt(tptr,i) .lt. htw(wptr,i)) then
                  hts(index,i) = htt(tptr,i)
                  ri(index,i) = rit(tptr,i)
                  rj(index,i) = rjt(tptr,i)
                  rk(index,i) = rkt(tptr,i)
                  lat(index,i) = latt(tptr,i)
                  lon(index,i) = lont(tptr,i)
                  up(index,i) = verif_missing_data
                  ui(index,i) = verif_missing_data
                  vp(index,i) = verif_missing_data
                  vi(index,i) = verif_missing_data
                  tp(index,i) = tpt(tptr,i)
                  ti(index,i) = tit(tptr,i)
                  tdp(index,i) = tdpt(tptr,i)
                  tdi(index,i) = tdit(tptr,i)
                  prp(index,i) = prpt(tptr,i)
                  pri(index,i) = prit(tptr,i)
                  type(index,i) = 't'
                  if (typet(tptr,i) .eq. 'm') then 
                    write(6,*)'raob ',i,' has mismatch man htt.lt.htw',
     1                        wptr,' ',tptr,' ',htw(wptr,i),' ',
     1                        htt(tptr,i)
                  endif
                  tptr = tptr + 1
                  index = index + 1
                else
                  if (htw(wptr,i) .eq. htt(tptr,i)) then
                    if ((typew(wptr,i) .eq. 'm') .and.  !find corresponding entry in t data
     1                  (riw(wptr,i) .eq. rit(tptr,i)) .and.
     1                  (rjw(wptr,i) .eq. rjt(tptr,i)) .and.
     1                  (latw(wptr,i) .eq. latt(tptr,i)) .and.
     1                  (lonw(wptr,i) .eq. lont(tptr,i))) then
c                     if (abs(rkdiff) .lt. 0.02) then
                        hts(index,i) = htw(wptr,i)
                        ri(index,i) = riw(wptr,i)
                        rj(index,i) = rjw(wptr,i)
                        rk(index,i) = rkw(wptr,i)
                        lat(index,i) = latw(wptr,i)
                        lon(index,i) = lonw(wptr,i)
                        up(index,i) = upw(wptr,i)
                        ui(index,i) = uiw(wptr,i)
                        vp(index,i) = vpw(wptr,i)
                        vi(index,i) = viw(wptr,i)
                        tp(index,i) = tpt(tptr,i)
                        ti(index,i) = tit(tptr,i)
                        tdp(index,i) = tdpt(tptr,i)
                        tdi(index,i) = tdit(tptr,i)
                        prp(index,i) = prpt(tptr,i)
                        pri(index,i) = prit(tptr,i)
                        type(index,i) = 'm'
c                     else   !mismatch on nav...write out wind only
c                       hts(index,i) = htw(wptr,i)
c                       ri(index,i) = riw(wptr,i)
c                       rj(index,i) = rjw(wptr,i)
c                       rk(index,i) = rkw(wptr,i)
c                       lat(index,i) = latw(wptr,i)
c                       lon(index,i) = lonw(wptr,i)
c                       up(index,i) = upw(wptr,i)
c                       ui(index,i) = uiw(wptr,i)
c                       vp(index,i) = vpw(wptr,i)
c                       vi(index,i) = viw(wptr,i)
c                       tp(index,i) = verif_missing_data
c                       ti(index,i) = verif_missing_data
c                       tdp(index,i) = verif_missing_data
c                       tdi(index,i) = verif_missing_data
c                       prp(index,i) = verif_missing_data
c                       pri(index,i) = verif_missing_data
c                       type(index,i) = 'm'

c                     endif
                    else
                      hts(index,i) = htw(wptr,i)
                      ri(index,i) = riw(wptr,i)
                      rj(index,i) = rjw(wptr,i)
                      rk(index,i) = rkw(wptr,i)
                      lat(index,i) = latw(wptr,i)
                      lon(index,i) = lonw(wptr,i)
                      up(index,i) = upw(wptr,i)
                      ui(index,i) = uiw(wptr,i)
                      vp(index,i) = vpw(wptr,i)
                      vi(index,i) = viw(wptr,i)
                      tp(index,i) = verif_missing_data
                      ti(index,i) = verif_missing_data
                      tdp(index,i) = verif_missing_data
                      tdi(index,i) = verif_missing_data
                      prp(index,i) = verif_missing_data
                      pri(index,i) = verif_missing_data
                      type(index,i) = 'm'

                      rkdiff = rkw(wptr,i) - rkt(tptr,i)
                      write(6,*) 'rkdiff = ',rkdiff
  
                      pctdiff = ((rkw(wptr,i)-rkt(tptr,i))/
     1                         rkw(wptr,i))*100
                      sumpctdiff = sumpctdiff + pctdiff
                      sumdiff = sumdiff + (rkw(wptr,i)-rkt(tptr,i))
                      ndiff = ndiff + 1
                      write(6,*)'raob ',i, ' has man data'
     1                         , ' rkw rkt pres %diff(w-t)/w ht'
                      write(6,*)rkw(wptr,i),rkt(tptr,i),prpt(tptr,i),
     1                          pctdiff,hts(index,i)
                      write(6,*) 'latw/t lonw/t riw/t rjw/t'
                      write(6,*) latw(wptr,i), latt(tptr,i),
     1                           lonw(wptr,i), lont(tptr,i),
     1                           riw(wptr,i), rit(tptr,i),
     1                           rjw(wptr,i), rjt(tptr,i)
                      write(6,*)
                    endif

                    wptr = wptr + 1
                    tptr = tptr + 1
                    index = index + 1

                  endif   !htw(wptr,i) .eq. htt(tptr,i)
                endif   !htt(tptr,i) .lt. htw(wptr,i)
              endif   !htw(wptr,i) .lt. htt(tptr,i) 

            else  !wptr .gt numw(i)

              do while (tptr .le. numt(i))
                hts(index,i) = htt(tptr,i)
                ri(index,i) = rit(tptr,i)
                rj(index,i) = rjt(tptr,i)
                rk(index,i) = rkt(tptr,i)
                lat(index,i) = latt(tptr,i)
                lon(index,i) = lont(tptr,i)
                up(index,i) = verif_missing_data
                ui(index,i) = verif_missing_data
                vp(index,i) = verif_missing_data
                vi(index,i) = verif_missing_data
                tp(index,i) = tpt(tptr,i)
                ti(index,i) = tit(tptr,i)
                tdp(index,i) = tdpt(tptr,i)
                tdi(index,i) = tdit(tptr,i)
                prp(index,i) = prpt(tptr,i)
                pri(index,i) = prit(tptr,i)
                type(index,i) = 't'

                tptr = tptr + 1
                index = index + 1
              enddo
            endif  !wptr .le. numw(i)

          enddo    !tptr .le. numt(i)

c         save number of heights for raob(i)
          nhts(i) = index - 1

c         do final sweep on data to make sure all missing values equal verif_missing_data
          do j = 1, nhts(i)
            if (up(j,i).eq.raob_missing_data)  
     1          up(j,i) = verif_missing_data
            if (ui(j,i).eq.raob_missing_data)  
     1          ui(j,i) = verif_missing_data
            if (vp(j,i).eq.raob_missing_data)  
     1          vp(j,i) = verif_missing_data
            if (vi(j,i).eq.raob_missing_data)  
     1          vi(j,i) = verif_missing_data
            if (tp(j,i).eq.raob_missing_data)  
     1          tp(j,i) = verif_missing_data
            if (ti(j,i).eq.raob_missing_data)  
     1          ti(j,i) = verif_missing_data
            if (tdp(j,i).eq.raob_missing_data)  
     1          tdp(j,i) = verif_missing_data
            if (tdi(j,i).eq.raob_missing_data)  
     1          tdi(j,i) = verif_missing_data
            if (prp(j,i).eq.raob_missing_data)  
     1          prp(j,i) = verif_missing_data
            if (pri(j,i).eq.raob_missing_data)  
     1          pri(j,i) = verif_missing_data
          enddo 

          write(6,*) 'raob ',i,' after merge'
          write(6,*) nhts(i)
          write(6,*) 'hts, type, prp, tp, up, vp'
          do j = 1, nhts(i)
            write(6,*) j,hts(j,i),type(j,i), prp(j,i),
     1                 tp(j,i), up(j,i), vp(j,i)
          enddo

        endif    !use_raob(i) .eq. 1
      enddo    !loop through raobs

      if (ndiff .gt. 0) then
        avgpctdiff = sumpctdiff/ndiff
        avgdiff = sumdiff/ndiff
        write(6,*) 'avg diff on rk = ',avgdiff,'  ndiff = ',ndiff
        write(6,*) 'avg pct diff = ',avgpctdiff
        write(6,*)
      endif

      return
      end
!5....................................................................................
