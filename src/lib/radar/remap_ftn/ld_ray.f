cdis   
cdis    open source license/disclaimer, forecast systems laboratory
cdis    noaa/oar/fsl, 325 broadway boulder, co 80305
cdis    
cdis    this software is distributed under the open source definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    in particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - all modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - if significant modifications or enhancements are made to this
cdis    software, the fsl software policy manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    this software and its documentation are in the public domain
cdis    and are furnished "as is."  the authors, the united states
cdis    government, its instrumentalities, officers, employees, and
cdis    agents make no warranty, express or implied, as to the usefulness
cdis    of the software and documentation for any purpose.  they assume
cdis    no responsibility (1) for the use of the software and
cdis    documentation; or (2) to provide technical support to users.
cdis   
cdis
cdis
cdis   
cdis
      subroutine ld_ray(i_ray,                  ! input
     :                  ngates_remap,           ! input
     :                  r_missing_data,         ! input
!    :                  i_missing_data,         ! input
     :                  velocity,               ! output
     :                  reflect,                ! output
     :                  istatus)                ! output
c
c     purpose:
c       transform a ray in the common block to output arrays that correspond
c       to the lookup tables. in the output arrays, the reflectivity and
c       velocity have been remapped to the same gate spacing.
c
c     changed buffer storage so vel, ref stored as integer*2.
c     this program converts reflectivity and velocity to real
c     keith brewster, caps
c
      implicit none
c
c     input variables
c
      integer i_ray
      integer ngates_remap
      real r_missing_data
c
c     output variables
c
      real velocity(ngates_remap)
      real reflect(ngates_remap)
      integer istatus
c
c     include files
c
      include 'remap_dims.inc'
      include 'remap_buffer.cmn'
      include 'remap_constants.dat'
c
c     misc. local variables
c
      integer i,igate_88d
      real ratio_ref, ratio_vel
c
c     make sure the array dimensions are consistent
c
      if (ngates_remap .ne. 1840) then
        write(6,805)
  805   format
     :(' ld_ray is hard wired for 1840 gates, error in ld_ray')
        write(6,810) ngates_remap
  810   format(' ngates_remap = ',i12)
        istatus = 0
        return
      end if

      if (n_vel_gates_cmn .ne. 920 .and. i_ray .eq. 1) then
        write(6,*)' note: ld_ray using non-standard # of vel gates'
      end if

!     if (n_vel_gates_cmn .ne. max_vel_gates) then
!       write(6,*)
!    1       ' error: ld_ray is expecting n_vel_gates_cmn=max_vel_gates'       
!       write(6,820) n_vel_gates_cmn,max_vel_gates
! 820   format(' n_vel_gates_cmn/max_vel_gates = ',2i12)
!       istatus = 0
!       return
!     end if

      if (n_ref_gates_cmn .ne. 460 .and. i_ray .eq. 1) then
        write(6,*)' note: ld_ray using non-standard # of ref gates'
      end if

!     if (n_ref_gates_cmn .ne. max_ref_gates) then
!       write(6,*)
!    1       ' error: ld_ray is expecting n_ref_gates_cmn=max_ref_gates'       
!       write(6,830) n_ref_gates_cmn,max_ref_gates
! 830   format(' n_ref_gates_cmn/max_ref_gates = ',2i12)
!       istatus = 0
!       return
!     end if
c
c     fill velocity and reflectivity ray.
c      output gates are 250 m spacing,
c             1840 are used for both reflecivity and velocity.
c
c      for a typical nexrad radar....
c     
c      input gates are as follows:
c             460 gates, 1000m spacing for reflectivity,
c             920 gates,  250m for velocity.
c      reflectivity is remapped by gate replication (factor of 4).
c      velocity is remapped by a 1 to 1 mapping for the 920 gates,
c        output gates 921-1840 are filled by the missing data value.
c
c      note that other input gate spacings may not be fully supported yet.
c

c     velocity

      ratio_vel = gsp_vel_m_cmn / gate_spacing_m
      velocity = r_missing_data           ! initialize array all at once

      do i = 1,n_vel_gates_cmn

        if (b_vel_cmn(i,i_ray) .ne. b_missing_cmn) then
          velocity(i) = 0.5*float(b_vel_cmn(i,i_ray))
!       else 
!         velocity(i) = r_missing_data
        end if

      enddo

c     reflectivity

      ratio_ref = gsp_ref_m_cmn / gate_spacing_m
      reflect = r_missing_data            ! initialize array all at once

      do igate_88d = 1,n_ref_gates_cmn

        i = nint(float(igate_88d) * ratio_ref)

        if(i .le. ngates_remap .and. i .ge. 1)then
            if (b_ref_cmn(igate_88d,i_ray) .ne. b_missing_cmn) then
              reflect(i) = 0.5*float(b_ref_cmn(igate_88d,i_ray))
!           else
!             reflect(i) = r_missing_data
            end if

        endif ! i

      enddo
c
      if(i_ray .eq. 1)then
          write(6,*)' ld_ray: ratio_ref = ',ratio_ref
     1                     ,' ratio_vel = ',ratio_vel
      endif

      istatus = 1
      return
      end
