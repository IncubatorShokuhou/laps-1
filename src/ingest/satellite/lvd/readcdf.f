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
      subroutine readcdf(csat_id,csat_type,chtype,
     1record,n_elem,n_lines,r4_image,scale_img,
     1latitude,longitude,
     1la1,lo1,dx,dy,latin,
     1lov, ivalidtime , ncid, istatus)

c     written by dan birkenheuer february 1994
c     j smart   5/95    included code to deal with gvar /public files
c     j smart   6/96    included code to deal with wfo satellite netcdf files
c     j smart   3/97    included subroutine rdblock_line_elem to read only the
c			appropriate sector of visible satellite data. argument
c                       list lets subroutine know whether logical or i*2 data type.
c     j smart   3/97    converted _vis routine to be the only netcdf reader. works
c			for vis, ir, wv and sounder. block reading allows this.
c     j smart   4/97    adapted code further to work with gvar. netcdf headers for
c                       raw gvar are different than for fsl-conus (ie., no la1, lo1,
c                       lov, or latin).
c     j smart   5/07    added capability to read fmi (finish met inst) data
c=====================================================================================
      implicit none

      include 'netcdf.inc'

      integer n_elem,n_lines,record
      integer rcode
c
      real      r4_image(n_elem,n_lines)
      real      latitude(n_elem,n_lines)
      real      longitude(n_elem,n_lines)
      integer  ib
      integer       bi (2)
      equivalence (ib,bi (1) )
      integer istatus
      integer   idx
      integer   id
      integer   nx
      integer   ny
      integer   nz
      real      la1     
      real      lo1     
      real      dx      
      real      dy      
      real      latin   
      real      lov     
      real      polat
      real      la2,lo2,level
      real      scale_img,offset_img ! o
      real arg1,arg2
      integer start(10)
      integer count(10)
      integer varid,ncid
      integer center_id, process_id,
     +     wmo_sat_id(record),ivalidtime,nav
      integer imax,jmax,kmax,kdim
      integer channel_fcinv

      integer lineres, elemres
      common /cdfdata/ lineres, elemres

      double precision reftime(record), valtime(record)
      character*30  c_valtime
      character*30  c_lov
      character*30  c_xres
      character*30  c_yres
      character*3   csat_type
      character*3   chtype
      character*6   csat_id
      character*132 origin_name
      character*132 x_dim
      character*132 y_dim
      character*132 earth_shape
      character*132 wavelength(record)
      character*132 grid_name
      character*132 process_name
      character*132 grid_type
      character*132 channel_comment_
      character*132 asctime
      character*10 varname
c---------------------------------------------------------
c   code

      print*,'subroutine readcdf...'
      print*,'n_elem/n_lines = ',n_elem,n_lines

      istatus = 0 ! bad istatus

      call ncpopt(0)

      if(csat_type.eq.'rll')then ! read lat/lon
         rcode = nf_inq_varid(ncid,'latitude',varid)
         if(rcode.ne.nf_noerr) then
            print *, nf_strerror(rcode)
            print *,'in var latitude'
         else
            write(6,*)
            write(6,*)'calling rdblock_line_elem - latitude'

            call rdblock_line_elem(csat_id,csat_type,chtype,
     &          ncid,varid,n_elem,n_lines,latitude,istatus)

            if(istatus .ne. 1)then
               write(6,*)'error in rdblock_line_elem'
               return
            endif
         endif
         call clean_nan2(latitude,n_elem,n_lines,istatus)
         arg1 = minval(latitude)
         arg2 = maxval(latitude)
         write(6,*)'readcdf latitude range:       ',arg1,arg2
         if(arg1 .eq. 0. .and. arg2 .eq. 0.)then
            write(6,*)' error in latitude range'
            istatus = 0
            return
         endif

         rcode = nf_inq_varid(ncid,'longitude',varid)
         if(rcode.ne.nf_noerr) then
            print *, nf_strerror(rcode)
            print *,'in var longitude'
         else
            write(6,*)
            write(6,*)'calling rdblock_line_elem - longitude'

            call rdblock_line_elem(csat_id,csat_type,chtype,
     &          ncid,varid,n_elem,n_lines,longitude,istatus)

            if(istatus .ne. 1)then
               write(6,*)'error in rdblock_line_elem'
               return
            endif
         endif
         call clean_nan2(longitude,n_elem,n_lines,istatus)
         write(6,*)'readcdf longitude range:      '
     1          ,minval(longitude),maxval(longitude)

      elseif(csat_type.eq.'nll' .or. csat_type.eq.'jma')then ! read 1d lat/lon
         write(6,*)' read 1d lat/lon arrays (under construction)'

      endif

      if(csat_type.eq.'ncp')then
         varname = 'channel'
      elseif(csat_type.eq.'nll')then
         varname = 'data'
      elseif(csat_type.eq.'jma')then ! determine varname from chtype
         if(chtype .eq. '10p')then
            varname = 'tbb_13'
         elseif(chtype .eq. 'vis')then
            varname = 'albedo_03'
         endif
      else
         varname = 'image'
      endif

      rcode = nf_inq_varid(ncid,trim(varname),varid)
      if(rcode.ne.nf_noerr) then
         print *, nf_strerror(rcode)
         print *,'in var ',trim(varname)
      else
         print *,'varid obtained for ',trim(varname)
      endif
c
c    statements to fill image                          
c
      write(6,*)
      write(6,*)'calling rdblock_line_elem - image read sub'

      call rdblock_line_elem(csat_id,csat_type,chtype,
     &ncid,varid,n_elem,n_lines,r4_image,istatus)

      if(istatus .ne. 1)then
         write(6,*)'error in rdblock_line_elem'
         return
      endif

      if(csat_type.eq.'rll' .or. csat_type.eq.'gnp' 
     1                      .or. csat_type.eq.'jma')then 

!         read scaling attribute
          rcode=nf_get_att_real(ncid,varid,'scale_factor'
     1                              ,scale_img)
          if(rcode.ne.nf_noerr) then
             write(6,*)'error reading image scaling attribute'
             scale_img = .01
             write(6,*)' use default value for image scale ',scale_img
          else
             write(6,*)' successfully read image scale ',scale_img
          endif

!         read offset attribute
          rcode=nf_get_att_real(ncid,varid,'add_offset'
     1                              ,offset_img)
          if(rcode.ne.nf_noerr) then
             write(6,*)'error reading image offset attribute'
             offset_img = 0.
             write(6,*)' use default value for image offset ',offset_img
          else
             write(6,*)' successfully read image offset ',offset_img
          endif

          r4_image(:,:) = r4_image(:,:) * scale_img + offset_img

      endif

      if(csat_type.eq.'rll')then 
          rcode = nf_inq_varid(ncid,'lineres',varid)
          if(rcode.ne.nf_noerr) then
             print *, nf_strerror(rcode)
             print *,'in var lineres'
          endif
          rcode = nf_get_var_int(ncid,varid,lineres)
          if(rcode.ne.nf_noerr) then
             print *, nf_strerror(rcode)
             print *,'in var lineres'
          endif
          if(rcode.ne.nf_noerr) then
             write(6,*)'error reading lineres'
          else
             write(6,*)' successfully read lineres ',lineres
          endif

          rcode = nf_inq_varid(ncid,'elemres',varid)
          if(rcode.ne.nf_noerr) then
             print *, nf_strerror(rcode)
             print *,'in var elemres'
          endif
          rcode = nf_get_var_int(ncid,varid,elemres)
          if(rcode.ne.nf_noerr) then
             print *, nf_strerror(rcode)
             print *,'in var elemres'
          endif
          if(rcode.ne.nf_noerr) then
             write(6,*)'error reading elemres'
          else
             write(6,*)' successfully read elemres ',elemres
          endif
c
      endif

      istatus = 0
c
c jsmart: 6-4-97.  wfo satellite netcdf files changed rather dramatically
c                  such that no header info exists and we must now jump over
c                  the statements to read the header info. 
c
      if(csat_type.eq.'cdf')then
         call read_netcdf_sat_head(ncid,record,
     + nx, ny, center_id,process_id, wmo_sat_id,dx,dy,
     + la1, latin, lo1, lov, reftime, valtime, earth_shape,
     + grid_name, grid_type, origin_name, process_name,
     + wavelength, x_dim, y_dim, istatus)

         ivalidtime = int(valtime(1))
      endif

      write(6,*)'readcdf image range: '
     1          ,minval(r4_image),maxval(r4_image)

      istatus = 1  ! ok!

      return
      end
