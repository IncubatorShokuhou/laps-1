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
        function i4_to_byte(i4_in)

cdoc    converts an integer to a byte (character) variable

        character i4_to_byte

        integer i4_in,i4

        character barg(4)
        equivalence (barg,i4)

        i4 = i4_in
        i4_to_byte = barg(4)

        return

        end

        function byte_to_i4(b_in)

cdoc    converts a byte (character) to an integer variable

        integer i4,byte_to_i4
        character b_in

        character barg(4)
        equivalence (barg,i4)

        i4 = 0

        barg(4) = b_in
#ifndef hpux
        if(i4 .ge. 128)i4 = i4 - 256
#endif

        byte_to_i4 = i4

!       write(6,*)i4,byte_to_i4

        return
        end

#ifdef use_trigd
      module trigd
         public
         private pi_180
         interface sind
            module procedure rsind, dsind
         end interface
         interface cosd
            module procedure rcosd, dcosd
         end interface
         interface tand
            module procedure rtand, dtand
         end interface
         interface asind
            module procedure rasind, dasind
         end interface
         interface acosd
            module procedure racosd, dacosd
         end interface
         interface atand
            module procedure ratand, datand
         end interface

         interface atan2d
            module procedure ratan2d, datan2d
         end interface
         contains
            double precision function pi_180()
              pi_180 = 2.d0*acos(0.d0)/180.d0
            return 
            end function pi_180

            real function rsind(val)
            real, intent(in) :: val
              rsind = sin(val*pi_180())
            return
            end function rsind

            real function rasind(val)
            real, intent(in) :: val
              rasind = asin(val)/pi_180()
            return
            end function rasind

            real function rcosd(val)
            real, intent(in) :: val
              rcosd = cos(val*pi_180())
            return
            end function rcosd

            real function racosd(val)
            real, intent(in) :: val
              racosd = acos(val)/pi_180()
            return
            end function racosd


            real function rtand(val)
            real, intent(in) :: val
              rtand = tan(val*pi_180())
            return
            end function rtand

            real function ratand(val)
            real, intent(in) :: val
              ratand = atan(val)/pi_180()
            return
            end function ratand

            real function ratan2d(val1,val2)
            real, intent(in) ::  val1,val2
              ratan2d = atan2(val1,val2)/pi_180()
            return
            end function ratan2d

            double precision function dsind(val)
            double precision, intent(in) :: val
              dsind = sin(val*pi_180())
            end function dsind

            double precision function dasind(val)
            double precision, intent(in) :: val
              dasind = asin(val)/pi_180()
            return
            end function dasind

            double precision function dcosd(val)
            double precision, intent(in) :: val
              dcosd = cos(val*pi_180())
            end function dcosd

            double precision function dacosd(val)
            double precision, intent(in) :: val
              dacosd = acos(val)/pi_180()
            return
            end function dacosd

            double precision function dtand(val)
            double precision, intent(in) :: val
              dtand = tan(val*pi_180())
            end function dtand
            
            double precision function datand(val)
            double precision, intent(in) :: val
              datand = atan(val)/pi_180()
            return
            end function datand


            double precision function datan2d(val1,val2)
            double precision, intent(in) ::  val1,val2
              datan2d = atan2(val1,val2)/pi_180()
            return
            end function datan2d
      end module

cdoc  the above functions operate with the use_trigd compiler directive.
cdoc  this is normally set by configure via 'trigd.inc' and supplies degree 
cdoc  based trig functions if they are not intrinsic.

#endif
      subroutine open_append(lun,file,status,istat)

cdoc  opens a file for appending (or a new file if non-existent)

      integer istat
      character*(*) file,status
      logical exist

      istat=1

#if defined(hpux) && !defined(f90)
      open(lun,file=file,status=status,access='append',err=998)
#elif defined(irix) && !defined(f90)
      open(lun,file=file,status=status,access='append',err=998)
#else
      inquire(file=file, exist=exist)
      if (exist) then
        open(lun, file=file, status="old", position="append"
     1                                   , action="write")
      else
        open(lun, file=file, status="new", action="write")
      end if
#endif

      return
 998  istat=0
      return
      end

#if defined(hpux) && !defined(f90)
      real function transfer(in,out)
      integer in
      real out
      print*, 'transfer is not supported by this compiler'
      print*, 'please report this error to laps_bugs@fsl.noaa.gov'
      stop
      transfer=0
      return
      end
#endif

cdoc below are some taiwan fgge routines

cdoc  below are dummy routines for other platforms

      subroutine read_fa(lun,filename                   ! i
     .               ,nx,ny,nz                          ! i
     .               ,r_missing_data                    ! i
     .               ,pressures_pa                      ! o
     .               ,ht,tp,rh,uw,vw                    ! o
     .               ,mslp                              ! o
     .               ,istatus)                          ! o
      write(6,*)' warning: routine read_fa not supported '
     1         ,'on this platform'
      istatus = 1 ! failure
      return
      end

      subroutine read_fa_nf ( lun, filename,                   ! i
     .                        nx, ny, nz,                      ! i
     .                        r_missing_data,                  ! i
     .                        pressures_pa,                    ! o
     .                        ht, tp, rh, uw, vw, ww,          ! o
     .                        pss, tps, rhs, uws, vws,         ! o
     .                        mslp,                            ! o
     .                        istatus )                        ! o
      write(6,*)' warning: routine read_fa_nf not supported '
     1         ,'on this platform'
      istatus = 1 ! failure
      return
      end

