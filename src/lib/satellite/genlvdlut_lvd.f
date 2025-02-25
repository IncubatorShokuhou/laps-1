cdis    forecast systems laboratory
cdis    noaa/oar/erl/fsl
cdis    325 broadway
cdis    boulder, co     80303
cdis 
cdis    forecast research division
cdis    local analysis and prediction branch
cdis    laps 
cdis 
cdis    this software and its documentation are in the public domain and 
cdis    are furnished "as is."  the united states government, its 
cdis    instrumentalities, officers, employees, and agents make no 
cdis    warranty, express or implied, as to the usefulness of the software 
cdis    and documentation for any purpose.  they assume no responsibility 
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis    
cdis    permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  all modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making 
cdis    the modifications.  if significant modifications or enhancements 
cdis    are made to this software, the fsl software policy manager  
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
      subroutine genlvdlut_lvd(nx_l,ny_l,maxch,lat,lon,it,js,
     &gri,grj,istatus)
c
      implicit none

      integer nx_l,ny_l

      real    lat(nx_l,ny_l)
      real    lon(nx_l,ny_l)
      real    gri(nx_l,ny_l,maxch) 
      real    grj(nx_l,ny_l,maxch) 

      integer istatus
      integer it,js,lc
      integer maxch
      integer jxd                  !amount of domain expansion needed for mapping
c
c dimensions for lat/lon
c

      include 'satellite_dims_lvd.inc'
      include 'satellite_common_lvd.inc'
c
c ======================== start ==============================
c acquiring laps latitude and longitude arrays.
c -------------------------------------------------------------------
c     write(6,*)'get laps lat/lon grid'
c     call get_laps_domain(nx_l,ny_l,'nest7grid',lat,lon,topo,istatus)
c     if(istatus.ne.1)then
c        write(6,*)'error - unable to get lat/lon data'
c        goto 900 
c     end if
c-----------------------------------------------------------------
c
c
c input values "js" (representing the satellite id) and "it" (representing
c the format type for that satellite) direct this subroutine to
c compute the appropriate mapping nav arrays.
c
      write(6,*)'subroutine genlvdlut_lvd: ',it,js,c_sat_types(it,js)

      istatus = -1

      if(c_sat_types(it,js).eq.'cdf'.or.
     &   c_sat_types(it,js).eq.'wfo'.or.
     &   c_sat_types(it,js).eq.'twn'.or.
     &   c_sat_types(it,js).eq.'gr2'.or.
     &   c_sat_types(it,js).eq.'gnp')then

         do lc=1,maxchannel

            if(ichannels(lc,it,js).eq.1)then

             write(6,59)c_sat_id(js),c_sat_types(it,js),
     &c_channel_types(lc,it,js)
59           format(/,'generate mapping arrays: ',a6,"/",a3,"/",a3)

             if(c_sat_types(it,js).eq.'gnp' .or.
     &          c_sat_types(it,js).eq.'gr2')then
                call gen_lut_fx(js,it,lc,nx_l,ny_l,
     &                lat,lon,gri(1,1,lc),grj(1,1,lc),istatus)
             else
                call gen_lut_lambert(js,it,lc,nx_l,ny_l,
     &                lat,lon,gri(1,1,lc),grj(1,1,lc),istatus)
             endif

             if(istatus.eq.1)then
              write(6,*)'finished '
             elseif(istatus.eq.0)then
              write(*,*)'mapping arrays already generated'
             else
              write(6,*)'error! mapping arrays not generated ',
     &c_sat_id(js),'/',c_sat_types(it,js),'/',c_channel_types(lc,it,js)
             endif

            else
             write(6,49)c_sat_id(js),c_sat_types(it,js),
     &                  c_channel_types(lc,it,js)
            endif

         enddo

      elseif(c_sat_types(it,js).eq.'gvr'.or.
     &       c_sat_types(it,js).eq.'gwc')then

         do lc=1,maxchannel

            if(ichannels(lc,it,js).eq.1)then

             write(6,59)c_sat_id(js),c_sat_types(it,js),
     & c_channel_types(lc,it,js)
             call gen_gvarimage_lut(js,it,lc,nx_l,ny_l,lat,lon,
     & gri(1,1,lc),grj(1,1,lc),istatus)
             if(istatus.eq.1)then
              write(6,*)'mapping array successfully generated'
             elseif(istatus.eq.0)then
              write(*,*)'ir mapping already generated'
             else
              write(6,*)'ir mapping not generated ',
     &c_channel_types(lc,it,js)
             endif
            else
             write(6,49)c_sat_id(js),c_sat_types(it,js),
     &c_channel_types(lc,it,js)
49           format(3x,a6,"/",a3,"/",a3,1x
     1             ,'inactive in satellite namelist (via ichannels-1)')       

            endif

         enddo

c         elseif(c_sat_types(it,js).eq.'asc')then
c
c          write(6,*)'generate luts for ascii data'
c           do lc=1,maxchannel
c              lc=ichannels(lc,it,js)
c              call get_ascii_dimensions(path_to_data_cdf,
c    &c_channel_types(lch,i,j),nelem,nlines,istatus)
c              call gen_ascii_lut(path_to_data_cdf,
c    &csatid,c_sat_types(i,j),c_channel_types(lch,i,j),
c    &nelem,nlines,nx_l,ny_l,lat,lon,istatus)
c              if(istatus.eq.1)then
c                 write(6,*)'lut generated'
c              else
c                 write(6,*)'lut not generated ',c_channel_types(lch,i)
c              endif
c           enddo

      elseif(c_sat_types(it,js).eq.'hko')then

            print*,'compute mapping arrays for merc ll - hko'
            do lc=1,maxchannel

               if(ichannels(lc,it,js).eq.1)then

                  write(6,59)c_sat_id(js),c_sat_types(it,js),
     &c_channel_types(lc,it,js)

                  call gen_lut_mercator(js,it,lc,nx_l,ny_l,
     & lat,lon,gri(1,1,lc),grj(1,1,lc),istatus)

                  if(istatus.eq.1)then
                     write(6,*)'lut generated'
                  elseif(istatus.eq.0)then
                     write(*,*)'ir lut already generated'
                  else
                     write(6,*)'error! lut not generated ',
     &c_sat_id(js),'/',c_sat_types(it,js),'/',c_channel_types(lc,it,js)
                  endif

               else
                  write(6,49)c_sat_id(js),c_sat_types(it,js),
     &c_channel_types(lc,it,js)
               endif

            enddo

      elseif(c_sat_types(it,js).eq.'nll' .or.
     &       c_sat_types(it,js).eq.'jma'      )then

            print*,'compute mapping arrays for ',c_sat_types(it,js)
            do lc=1,maxchannel

               if(ichannels(lc,it,js).eq.1)then

                  write(6,59)c_sat_id(js),c_sat_types(it,js),
     &c_channel_types(lc,it,js)

                  call gen_lut_ll(js,it,lc,nx_l,ny_l,
     & lat,lon,gri(1,1,lc),grj(1,1,lc),istatus)

                  if(istatus.eq.1)then
                     write(6,*)'lut generated'
                  elseif(istatus.eq.0)then
                     write(*,*)'ir lut already generated'
                  else
                     write(6,*)'error! lut not generated ',
     &c_sat_id(js),'/',c_sat_types(it,js),'/',c_channel_types(lc,it,js)
                  endif

               else
                  write(6,49)c_sat_id(js),c_sat_types(it,js),
     &c_channel_types(lc,it,js)
               endif

            enddo

c polar orbiter/netcdf/polar stereographic projected data
c ---------------------------------------------------------
      elseif(c_sat_types(it,js).eq.'ncp')then

        print*,'compute mapping arrays for polar stereo - ncp'
        do lc=1,maxchannel
           if(ichannels(lc,it,js).eq.1)then

              write(6,59)c_sat_id(js),c_sat_types(it,js),
     &c_channel_types(lc,it,js)

              call gen_lut_ps(js,it,lc,nx_l,ny_l,
     & lat,lon,gri(1,1,lc),grj(1,1,lc),istatus)

              if(istatus.eq.1)then
                 write(6,*)'mapping arrays successfully generated'
              else
                write(6,*)'error! lut not generated ',
     &c_sat_id(js),'/',c_sat_types(it,js),'/',c_channel_types(lc,it,js)
              endif

           else
              write(6,49)c_sat_id(js),c_sat_types(it,js),
     &c_channel_types(lc,it,js)
           endif

        enddo

      elseif(c_sat_types(it,js).ne.'     ')then

        write(6,*)'unknown satellite data type! '
        write(6,*)'check static/satellite_lvd.nl'

      endif

c        elseif(c_sat_types(it,js).ne.'   ')then
c         write(6,47)c_sat_id(js), c_sat_types(it,js)
c47     format(2x,a6,"/",a3,' not on in satellite namelist')
c        endif
c       enddo
c      else
c       write(*,*)c_sat_id(js),' not on in satellite namelist'
c      endif
c     enddo

900   return 
      end
