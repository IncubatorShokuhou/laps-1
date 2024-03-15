       subroutine w3tagb(prog,kyr,jd,lf,org)
c$$$   subprogram documentation block
c
c subprogram: w3tagb        operational job identifier
c   prgmmr: farley          org: np11          date: 1998-03-17
c
c abstract: prints identifying information for operational
c   codes. called at the beginning of a code, w3tagb prints
c   the program name, the year and julian day of its
c   compilation, and the responsible organization. on a 2nd
c   line it prints the starting date-time. called at the
c   end of a job, entry routine, w3tage prints a line with the
c   ending date-time and a 2nd line stating the program name 
c   and that it has ended.
c
c program history log:
c   85-10-29  j.newell
c   89-10-20  r.e.jones   convert to cray cft77 fortran
c   91-03-01  r.e.jones   add machine name to ending line
c   92-12-02  r.e.jones   add start-ending time-date
c   93-11-16  r.e.jones   add day of year, day of week, and julian day
c                         number. 
c   97-12-24  m.farley    print statements modified for 4-digit yr 
c   98-03-17  m.farley    replaced datimx with calls to w3locdat/w3doxdat 
c   99-01-29  b. vuong    converted to ibm rs/6000 sp
c
c   99-06-17  a. spruill  adjusted the size of program name to accommodate
c                         the 20 character name convention on the ibm sp.
c 1999-08-24  gilbert     added call to start() in w3tagb and a call
c                         to summary() in w3tage to print out a 
c                         resource summary list for the program using
c                         w3tags.
c
c usage:  call w3tagb(prog, kyr, jd, lf, org)
c         call w3tage(prog)
c
c   input variables:
c     names  interface description of variables and types
c     ------ --------- -----------------------------------------------
c     prog   arg list  program name   character*1
c     kyr    arg list  year of compilation   integer
c     jd     arg list  julian day of compilation   integer
c     lf     arg list  hundreths of julian day of compilation
c                      integer     (range is 0 to 99 inclusive)
c     org    arg list  organization code (such as wd42)
c                      character*1
c
c   output variables:
c     names  interface description of variables and types
c     ----------------------------------------------------------------
c     ddate  print     year and julian day (nearest hundreth)
c            file      of compilation  real
c
c   subprograms called: clock, date
c
c   remarks: full word used in order to have at least
c            seven decimal digits accuracy for value of ddate.
c            subprogram clock and date may differ for each type
c            computer. you may have to change them for another
c            type of computer.
c
c attributes:
c   language: fortran 90
c
c$$$
c
         character *(*) prog,org
         character * 3 jmon(12)
         character * 3 dayw(7)
c
         integer       idat(8), jdow, jdoy, jday
c
         save
c
         data  dayw/'sun','mon','tue','wen','thu','fri','sat'/
         data  jmon  /'jan','feb','mar','apr','may','jun',
     &                'jul','aug','sep','oct','nov','dec'/
c 
         call start()

         dyr   = kyr
         dyr   = 1.0e+03 * dyr
         djd   = jd
         dlf   = lf
         dlf   = 1.0e-02 * dlf
         ddate = dyr + djd + dlf
         print 600
  600    format(//,10('* . * . '))
         print 601, prog, ddate, org
  601    format(5x,'program ',a,' has begun. compiled ',f10.2,
     &   5x, 'org: ',a)
c
         call w3locdat(idat)
         call w3doxdat(idat,jdow,jdoy,jday)
         print 602, jmon(idat(2)),idat(3),idat(1),idat(5),idat(6),
     &   idat(7),idat(8),jdoy,dayw(jdow),jday
  602    format(5x,'starting date-time  ',a3,1x,i2.2,',',
     &   i4.4,2x,2(i2.2,':'),i2.2,'.',i3.3,2x,i3,2x,a3,2x,i8,//)
         return
c
         entry w3tage(prog)
c
         call w3locdat(idat)
         call w3doxdat(idat,jdow,jdoy,jday)
         print 603, jmon(idat(2)),idat(3),idat(1),idat(5),idat(6),
     &   idat(7),idat(8),jdoy,dayw(jdow),jday
  603    format(//,5x,'ending date-time    ',a3,1x,i2.2,',',
     &   i4.4,2x,2(i2.2,':'),i2.2,'.',i3.3,2x,i3,2x,a3,2x,i8)
         print 604, prog
  604    format(5x,'program ',a,' has ended.  ibm rs/6000 sp')
c 604    format(5x,'program ',a,' has ended.  cray j916/2048')
c 604    format(5x,'program ',a,' has ended.  cray y-mp el2/256')
         print 605
  605    format(10('* . * . '))

         call summary()
c
         return
         end
