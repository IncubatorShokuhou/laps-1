      subroutine write_lga(nx_laps,ny_laps,nz_laps,bgtime,bgvalid,
     .cmodel,missingflag,pr,ht,tp,sh,uw,vw,ww,istatus)

      implicit none

      integer nx_laps,ny_laps,nz_laps
      integer i,j,k,ic,lendir
      integer bgtime,bgvalid
      integer istatus
      integer warncnt
      integer ip(nz_laps)

      real    missingflag 
      real    pr(nz_laps)                  !pressure or sigma-p grid levels
      real    ht(nx_laps,ny_laps,nz_laps), !height (m)
     .        tp(nx_laps,ny_laps,nz_laps), !temperature (k)
     .        sh(nx_laps,ny_laps,nz_laps), !specific humidity (kg/kg)
     .        uw(nx_laps,ny_laps,nz_laps), !u-wind (m/s)
     .        vw(nx_laps,ny_laps,nz_laps), !v-wind (m/s)
     .        ww(nx_laps,ny_laps,nz_laps)  !w-wind (pa/s)

      character*256 outdir
      character*31  ext
      character*3   var(nz_laps)
      character*4   lvl_coord(nz_laps)
      character*10  units(nz_laps)
      character*125 comment(nz_laps)
      character*5 fcst_hh_mm
      character*9 a9time
      character*256 outfile
      character*(*) cmodel

      character*40   v_g

      logical l_exist

      call get_vertical_grid(v_g,istatus)
      call upcase(v_g, v_g)

      ext='lga'
      call get_directory('lga',outdir,lendir)

c     test for existence of output file
      call make_fnam_lp(bgtime,a9time,istatus)
      if (istatus .ne. 1) then
        write (6,*)
     1'error converting i4time to file name...write aborted.'
        istatus=0
        return
      endif

      call make_fcst_time(bgvalid,bgtime,fcst_hh_mm,istatus)

      outfile = trim(outdir)//'/'//a9time//trim(fcst_hh_mm)//'.lga'
      inquire(file=outfile,exist=l_exist)
      if(l_exist .eqv. .false.)then
          write(6,*)' lga file = ',trim(outfile)
      else
          write(6,*)' lga file (already exists) = ',trim(outfile)
          write(6,*)
     1         ' returning from write_lga without writing new lga file'     
          istatus = 1
          return
      endif
c

      call s_len(cmodel,ic)
      do k=1,nz_laps                 ! height
         do j=1,ny_laps
         do i=1,nx_laps
            if(ht(i,j,k) .ge. missingflag .and. warncnt.lt.100) 
     +              then
               print*,'missingflag at ',i,j,k,' in ht'
               warncnt=warncnt+1
            endif
         enddo
         enddo
c only need to do some of these once.
         if (v_g .eq. 'pressure') then     
            ip(k)=int(pr(k)) ! integer mb
            lvl_coord(k)='mb  '
            comment(k)=cmodel(1:ic)//' interpolated to laps isobaric.'
         elseif (v_g .eq. 'sigma_p') then
            ip(k)=nint(pr(k)*1000.) ! integer sigma * 1000
            lvl_coord(k)='    '
            comment(k)=cmodel(1:ic)//' interpolated to laps sigma p'
         endif

         var(k)='ht '
         units(k)='meters'
      enddo
      print*,'ht'
      call write_laps(bgtime,bgvalid,outdir,ext,
     .                nx_laps,ny_laps,nz_laps,nz_laps,var,
     .                ip,lvl_coord,units,comment,ht,istatus)
      if (istatus .ne. 1) then
          print *,'error writing interpolated data to laps database.'
          return
      endif

      do k=1,nz_laps                 ! temperature
         do j=1,ny_laps
         do i=1,nx_laps
            if(tp(i,j,k) .ge. missingflag .and. warncnt.lt.100) 
     +              then
               print*,'missingflag at ',i,j,k,' in tp'
               warncnt=warncnt+1
            endif
         enddo
         enddo
         var(k)='t3 '
         units(k)='kelvin'
      enddo
      print*,'t3'
      call write_laps(bgtime,bgvalid,outdir,ext,
     .                nx_laps,ny_laps,nz_laps,nz_laps,var,
     .                ip,lvl_coord,units,comment,tp,istatus)
      if (istatus .ne. 1) then
          print *,'error writing interpolated data to laps database.'
          return
      endif

      do k=1,nz_laps                 ! specific humidity
         do j=1,ny_laps
         do i=1,nx_laps
            if(sh(i,j,k) .ge. missingflag .and. warncnt.lt.100) 
     +           then
               print*,'missingflag at ',i,j,k,' in sh'
               warncnt=warncnt+1
            endif
         enddo
         enddo
         var(k)='sh '
         units(k)='kg/kg'
      enddo
      print*,'sh'
      call write_laps(bgtime,bgvalid,outdir,ext,
     .                nx_laps,ny_laps,nz_laps,nz_laps,var,
     .                ip,lvl_coord,units,comment,sh,istatus)
      if (istatus .ne. 1) then
          print *,'error writing interpolated data to laps database.'
          return
      endif

      do k=1,nz_laps                 ! u-component wind
         do j=1,ny_laps
         do i=1,nx_laps
            if(uw(i,j,k) .ge. missingflag .and. warncnt.lt.100) 
     +              then
               print*,'missingflag at ',i,j,k,' in uw'
               warncnt=warncnt+1
            endif
         enddo
         enddo
         var(k)='u3 '
         units(k)='m/s'
      enddo
      print*,'u3'
      call write_laps(bgtime,bgvalid,outdir,ext,
     .                nx_laps,ny_laps,nz_laps,nz_laps,var,
     .                ip,lvl_coord,units,comment,uw,istatus)
      if (istatus .ne. 1) then
          print *,'error writing interpolated data to laps database.'
          return
      endif

      do k=1,nz_laps                 ! v-component wind
         do j=1,ny_laps
         do i=1,nx_laps
            if(vw(i,j,k) .ge. missingflag .and. warncnt.lt.100) 
     +              then
               print*,'missingflag at ',i,j,k,' in vw'
               warncnt=warncnt+1
            endif
         enddo
         enddo
         var(k)='v3 '
      enddo
      print*,'v3'
      call write_laps(bgtime,bgvalid,outdir,ext,
     .                nx_laps,ny_laps,nz_laps,nz_laps,var,
     .                ip,lvl_coord,units,comment,vw,istatus)
      if (istatus .ne. 1) then
          print *,'error writing interpolated data to laps database.'
          return
      endif
c
c if laps is background then we'll allow writing missing ww data.
c --------------------------------------------------------------
      if(cmodel.eq.'laps')then
         do k=1,nz_laps                 ! w-component wind
            do j=1,ny_laps
            do i=1,nx_laps
             if(ww(i,j,k) .ge. missingflag .and. warncnt.lt.50)
     +              then
               print*,'missingflag at ',i,j,k,' in ww'
               warncnt=warncnt+1
             endif
            enddo
            enddo
            var(k)='om '
            units(k)='pa/s'
         enddo


      else

         do k=1,nz_laps                 ! w-component wind
            do j=1,ny_laps
            do i=1,nx_laps
             if(ww(i,j,k) .ge. missingflag .and. warncnt.lt.100)
     +              then
               print*,'missingflag at ',i,j,k,' in ww'
               warncnt=warncnt+1
             endif
            enddo
            enddo
            var(k)='om '
            units(k)='pa/s'
         enddo
      endif

      print*,'om'
      call write_laps(bgtime,bgvalid,outdir,ext,
     .                nx_laps,ny_laps,nz_laps,nz_laps,var,
     .                ip,lvl_coord,units,comment,ww,istatus)
      if (istatus .ne. 1) then
          print *,'error writing interpolated data to laps database.'
          return
      endif

      return
      end
c
c -------------------------------------------------------------------------------
c
      subroutine write_lgap(nx_laps,ny_laps,nz_laps,bgtime,bgvalid,
     .cmodel,missingflag,ht,pr,tp,sh,uw,vw,ww,istatus)

      implicit none

      integer nx_laps,ny_laps,nz_laps
      integer i,j,k,ic,lendir
      integer bgtime,bgvalid
      integer istatus
      integer warncnt
      integer ip(nz_laps)

      real    missingflag
      real    ht(nz_laps)                  !sigma_ht height levels
      real    pr(nx_laps,ny_laps,nz_laps), !pressure (pa)
     .        tp(nx_laps,ny_laps,nz_laps), !temperature (k)
     .        sh(nx_laps,ny_laps,nz_laps), !specific humidity (kg/kg)
     .        uw(nx_laps,ny_laps,nz_laps), !u-wind (m/s)
     .        vw(nx_laps,ny_laps,nz_laps), !v-wind (m/s)
     .        ww(nx_laps,ny_laps,nz_laps)  !w-wind (m/s)

      real scale_height

      character*256 outdir
      character*31  ext
      character*3   var(nz_laps)
      character*4   lvl_coord(nz_laps)
      character*10  units(nz_laps)
      character*125 comment(nz_laps)
      character*(*) cmodel

      ext='lga'
      call get_directory('lga',outdir,lendir)

      call s_len(cmodel,ic)
      do k=1,nz_laps                 ! pressure
         do j=1,ny_laps
         do i=1,nx_laps
            if(pr(i,j,k) .ge. missingflag .and. warncnt.lt.100) 
     +              then
               print*,'missingflag at ',i,j,k,' in pr'
               warncnt=warncnt+1
            endif
         enddo
         enddo
c only need to do some of these once.
         ip(k)=nint(ht(k))
         var(k)='p3 '
         lvl_coord(k)='    '
         units(k)='pascals'
         comment(k)=cmodel(1:ic)//' interpolated to laps sigma_ht.'
      enddo
      print*,'p3'
      call write_laps(bgtime,bgvalid,outdir,ext,
     .                nx_laps,ny_laps,nz_laps,nz_laps,var,
     .                ip,lvl_coord,units,comment,pr,istatus)
      if (istatus .ne. 1) then
          print *,'error writing interpolated data to laps database.'
          return
      endif

      do k=1,nz_laps                 ! temperature
         do j=1,ny_laps
         do i=1,nx_laps
            if(tp(i,j,k) .ge. missingflag .and. warncnt.lt.100) 
     +              then
               print*,'missingflag at ',i,j,k,' in tp'
               warncnt=warncnt+1
            endif
         enddo
         enddo
         var(k)='t3 '
         units(k)='kelvin'
      enddo
      print*,'t3'
      call write_laps(bgtime,bgvalid,outdir,ext,
     .                nx_laps,ny_laps,nz_laps,nz_laps,var,
     .                ip,lvl_coord,units,comment,tp,istatus)
      if (istatus .ne. 1) then
          print *,'error writing interpolated data to laps database.'
          return
      endif

      do k=1,nz_laps                 ! specific humidity
         do j=1,ny_laps
         do i=1,nx_laps
            if(sh(i,j,k) .ge. missingflag .and. warncnt.lt.100) 
     +           then
               print*,'missingflag at ',i,j,k,' in sh'
               warncnt=warncnt+1
            endif
         enddo
         enddo
         var(k)='sh '
         units(k)='kg/kg'
      enddo
      print*,'sh'
      call write_laps(bgtime,bgvalid,outdir,ext,
     .                nx_laps,ny_laps,nz_laps,nz_laps,var,
     .                ip,lvl_coord,units,comment,sh,istatus)
      if (istatus .ne. 1) then
          print *,'error writing interpolated data to laps database.'
          return
      endif

      do k=1,nz_laps                 ! u-component wind
         do j=1,ny_laps
         do i=1,nx_laps
            if(uw(i,j,k) .ge. missingflag .and. warncnt.lt.100) 
     +              then
               print*,'missingflag at ',i,j,k,' in uw'
               warncnt=warncnt+1
            endif
         enddo
         enddo
         var(k)='u3 '
         units(k)='m/s'
      enddo
      print*,'u3'
      call write_laps(bgtime,bgvalid,outdir,ext,
     .                nx_laps,ny_laps,nz_laps,nz_laps,var,
     .                ip,lvl_coord,units,comment,uw,istatus)
      if (istatus .ne. 1) then
          print *,'error writing interpolated data to laps database.'
          return
      endif

      do k=1,nz_laps                 ! v-component wind
         do j=1,ny_laps
         do i=1,nx_laps
            if(vw(i,j,k) .ge. missingflag .and. warncnt.lt.100) 
     +              then
               print*,'missingflag at ',i,j,k,' in vw'
               warncnt=warncnt+1
            endif
         enddo
         enddo
         var(k)='v3 '
      enddo
      print*,'v3'
      call write_laps(bgtime,bgvalid,outdir,ext,
     .                nx_laps,ny_laps,nz_laps,nz_laps,var,
     .                ip,lvl_coord,units,comment,vw,istatus)
      if (istatus .ne. 1) then
          print *,'error writing interpolated data to laps database.'
          return
      endif
c
c if laps is background then we'll allow writing missing ww data.
c --------------------------------------------------------------
      if(cmodel.eq.'laps')then
         do k=1,nz_laps                 ! w-component wind
            do j=1,ny_laps
            do i=1,nx_laps
             if(ww(i,j,k) .ge. missingflag .and. warncnt.lt.50)
     +              then
               print*,'missingflag at ',i,j,k,' in ww'
               warncnt=warncnt+1
             endif

!            convert from omega to w (approximately - not yet coded)            

            enddo
            enddo
            var(k)='w3 '
            units(k)='m/s'
         enddo


      else

         do k=1,nz_laps                 ! w-component wind
            do j=1,ny_laps
            do i=1,nx_laps
             if(ww(i,j,k) .ge. missingflag)then
               if(warncnt.lt.100)
     +              then
                 print*,'missingflag at ',i,j,k,' in ww'
                 warncnt=warncnt+1
               endif
             else
!              convert from omega to w (approximately)            
               scale_height = 8000.
               ww(i,j,k) = - (ww(i,j,k) / pr(i,j,k)) * scale_height

             endif

            enddo
            enddo
            var(k)='w3 '
            units(k)='m/s'
         enddo
      endif

      print*,'w3'
      call write_laps(bgtime,bgvalid,outdir,ext,
     .                nx_laps,ny_laps,nz_laps,nz_laps,var,
     .                ip,lvl_coord,units,comment,ww,istatus)
      if (istatus .ne. 1) then
          print *,'error writing interpolated data to laps database.'
          return
      endif

      return
      end
c
c -------------------------------------------------------------------------------
c
      subroutine write_lgb(nx_laps,ny_laps,bgtime,bgvalid,cmodel
     .,missflag,uw_sfc,vw_sfc,tp_sfc,t_sfc,qsfc,pr_sfc,mslp,td_sfc
     .,rp_sfc,pcp_sfc,cw_sfc,istatus)

      implicit none

      integer nx_laps,ny_laps
      integer bgtime,bgvalid
      integer istatus,ic,len_dir
      integer i,j
      integer warncnt

      real      missflag
      real      tp_sfc(nx_laps,ny_laps),
     .          t_sfc (nx_laps,ny_laps),
     .          td_sfc(nx_laps,ny_laps),
     .          qsfc(nx_laps,ny_laps),
     .          uw_sfc(nx_laps,ny_laps),
     .          vw_sfc(nx_laps,ny_laps),
     .          pr_sfc(nx_laps,ny_laps),
     .          rp_sfc(nx_laps,ny_laps),
     .          pcp_sfc(nx_laps,ny_laps),
     .          cw_sfc(nx_laps,ny_laps),
     .          mslp(nx_laps,ny_laps)

      character*256 outdir
      character*31  ext
      character*3   var
      character*4   lvl_coord
      character*10  units
      character*125 comment
      character*(*) cmodel

      ext = 'lgb'
      call get_directory(ext,outdir,len_dir)
      call s_len(cmodel,ic)

c        print *,'writing to ',outdir(1:len_dir),fname//af(nf)(3:4)
c    &,'00.',ext
c        call write_laps(bgtime,bgvalid,outdir,ext
c    .           ,nx_laps,ny_laps,nsfc_fields,nsfc_fields
c    .           ,var,ip,lvl_coord
c    .           ,units,comment,sfcgrid,istatus)
c        if(istatus.eq.1)then
c           goto 88
c        else
c           print*,'error writing interpolated data to laps lgb'
c           return
c        endif

      comment=cmodel(1:ic)//' interpolated to laps surface.'

c --- u
      warncnt=0
      do i=1,nx_laps
      do j=1,ny_laps
         if(uw_sfc(i,j).ge.missflag.and.warncnt.lt.100)
     +              then
            print*,'missing data at ',i,j,' in uw_sfc'
            warncnt=warncnt+1
         endif
      enddo
      enddo
      units='m/s'
      var='usf'
      print*,'usf'
      call write_laps(bgtime,bgvalid,outdir,ext
     .           ,nx_laps,ny_laps,1,1,var,0,lvl_coord
     .           ,units,comment,uw_sfc,istatus)
      if (istatus .ne. 1) then
         print*,'error writing interpolated data to laps lgb - u'
         return
      endif
c --- v
      warncnt=0
      do i=1,nx_laps
      do j=1,ny_laps
         if(vw_sfc(i,j).ge.missflag.and.warncnt.lt.100)
     +              then
            print*,'missing data at ',i,j,' in vw_sfc'
            warncnt=warncnt+1
         endif
      enddo
      enddo
      var='vsf'
      print*,'vsf'
      call write_laps(bgtime,bgvalid,outdir,ext
     .           ,nx_laps,ny_laps,1,1,var,0,lvl_coord
     .           ,units,comment,vw_sfc,istatus)
      if (istatus .ne. 1) then
         print*,'error writing interpolated data to laps lgb - v'
         return
      endif
c --- t
      warncnt=0
      do i=1,nx_laps
      do j=1,ny_laps
         if(tp_sfc(i,j).ge.missflag.and.warncnt.lt.100)
     +              then
            print*,'missing data at ',i,j,' in tp_sfc'
            warncnt=warncnt+1
         endif
      enddo
      enddo
      units='k'
      var='tsf'
      print*,'tsf'
      comment=cmodel(1:ic)//' 2-m temp fg'
      call write_laps(bgtime,bgvalid,outdir,ext
     .           ,nx_laps,ny_laps,1,1,var,0,lvl_coord
     .           ,units,comment,tp_sfc,istatus)
      if (istatus .ne. 1) then
         print*,'error writing interpolated data to laps lgb - t'
         return
      endif
c --- t at sfc
      warncnt=0
      do i=1,nx_laps
      do j=1,ny_laps
         if(t_sfc(i,j).ge.missflag.and.warncnt.lt.100)
     +              then
            print*,'missing data at ',i,j,' in t_sfc'
            warncnt=warncnt+1
         endif
      enddo
      enddo
      units='k'
      var='tgd'
      print*,'tgd'
      comment=cmodel(1:ic)//' ground temp fg'
      call write_laps(bgtime,bgvalid,outdir,ext
     .           ,nx_laps,ny_laps,1,1,var,0,lvl_coord
     .           ,units,comment,t_sfc,istatus)
      if (istatus .ne. 1) then
         print*,'error writing interpolated data to laps lgb - tsfc'
         return
      endif
c --- q
      warncnt=0
      do i=1,nx_laps
      do j=1,ny_laps
         if(qsfc(i,j).ge.missflag.and.warncnt.lt.100)
     +              then
            print*,'missing data at ',i,j,' in qsfc'
            warncnt=warncnt+1
         endif
      enddo
      enddo
      units='kg/kg'
      var='rsf'
      print*,'rsf'
      call write_laps(bgtime,bgvalid,outdir,ext
     .           ,nx_laps,ny_laps,1,1,var,0,lvl_coord
     .           ,units,comment,qsfc,istatus)
      if (istatus .ne. 1) then
         print*,'error writing interpolated data to laps lgb - sh'
         return
      endif
c --- p sfc
      warncnt=0
      do i=1,nx_laps
      do j=1,ny_laps
         if(pr_sfc(i,j).ge.missflag.and.warncnt.lt.100)
     +              then
            print*,'missing data at ',i,j,' in pr_sfc'
            warncnt=warncnt+1
         endif
      enddo
      enddo
      units='pa'
      var='psf'
      print*,'psf'
      call write_laps(bgtime,bgvalid,outdir,ext
     .           ,nx_laps,ny_laps,1,1,var,0,lvl_coord
     .           ,units,comment,pr_sfc,istatus)
      if (istatus .ne. 1) then
         print*,'error writing interpolated data to laps lgb - psf'
         return
      endif
c --- mslp
      warncnt=0
      do i=1,nx_laps
      do j=1,ny_laps
         if(mslp(i,j).ge.missflag.and.warncnt.lt.100)
     +              then
            print*,'missing data at ',i,j,' in mslp'
            warncnt=warncnt+1
         endif
      enddo
      enddo
      var='slp'
      print*,'slp'
      call write_laps(bgtime,bgvalid,outdir,ext
     .           ,nx_laps,ny_laps,1,1,var,0,lvl_coord
     .           ,units,comment,mslp,istatus)
      if (istatus .ne. 1) then
         print*,'error writing interpolated data to laps lgb - slp'
         return
      endif
c --- td
      warncnt=0
      do i=1,nx_laps
      do j=1,ny_laps
         if(td_sfc(i,j).ge.missflag.and.warncnt.lt.100)
     +              then
            print*,'missing data at ',i,j,' in td_sfc'
            warncnt=warncnt+1
         endif
      enddo
      enddo
      units='k'
      var='dsf'
      print*,'dsf'
      call write_laps(bgtime,bgvalid,outdir,ext
     .           ,nx_laps,ny_laps,1,1,var,0,lvl_coord
     .           ,units,comment,td_sfc,istatus)
      if (istatus .ne. 1) then
         print*,'error writing interpolated data to laps lgb - dsf'
         return
      endif
c --- reduced pressure
      warncnt=0
      do i=1,nx_laps
      do j=1,ny_laps
         if(rp_sfc(i,j).ge.missflag.and.warncnt.lt.100)
     +              then
            print*,'missing data at ',i,j,' in reduced press'
            warncnt=warncnt+1
         endif
      enddo
      enddo
      var='p'
      units='pa'
      print*,'p'
      call write_laps(bgtime,bgvalid,outdir,ext
     .           ,nx_laps,ny_laps,1,1,var,0,lvl_coord
     .           ,units,comment,rp_sfc,istatus)
      if (istatus .ne. 1) then
         print*,'error writing interpolated data to laps lgb - p'
         return
      endif

c --- precipitation 
      warncnt=0
      do i=1,nx_laps
      do j=1,ny_laps
         if(pcp_sfc(i,j).ge.missflag.and.warncnt.lt.100) ! modified rp_sfc -> pcp_sfc by wei-ting (130312)
     +              then
            print*,'missing data at ',i,j,' in precip - pcp '
            warncnt=warncnt+1
         endif
      enddo
      enddo
      var='r01'
      units='m' 
      print*,'r01'
      call write_laps(bgtime,bgvalid,outdir,ext
     .           ,nx_laps,ny_laps,1,1,var,0,lvl_coord
     .           ,units,comment,pcp_sfc,istatus)
      if (istatus .ne. 1) then
         print*,'error writing interpolated data to laps lgb - r01'
         return
      endif

c --- cloud water
      warncnt=0
      do i=1,nx_laps
      do j=1,ny_laps
         if(cw_sfc(i,j).ge.missflag.and.warncnt.lt.100)
     +              then
            print*,'missing data at ',i,j,' in cw_sfc'
            warncnt=warncnt+1
         endif
      enddo
      enddo
      var='cw'
      units='m' 
      print*,'cw'
      call write_laps(bgtime,bgvalid,outdir,ext
     .           ,nx_laps,ny_laps,1,1,var,0,lvl_coord
     .           ,units,comment,cw_sfc,istatus)
      if (istatus .ne. 1) then
         print*,
     .      'warning: not writing interpolated data to laps lgb - cw'
         istatus = 1 ! keep the program running normally
         return
      endif

      return
      end
