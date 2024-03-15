      subroutine timpr(kfilds,kfildo,ititle) 
c 
c        january   1994   chambers, glahn   tdl  mos-2000
c        september 2002   glahn   added commas in format 110
c 
c        purpose 
c            a date/time stamping function which will write the date and
c            time and a user message of up to 20 characters on unit 
c            no. kfildo.  under the d compiler option, message will 
c            also be written to the current console, unit no. kfilds
c            (only) when kfilds ne kfildo. 
c 
c        data set use 
c            kfilds - unit number for output to current console.  (output) 
c            kfildo - output (print) file unit number.  (output) 
c 
c        variables 
c 
c            input 
c              kfilds = unit number for output to current console. 
c              kfildo = output (print) file unit number. 
c              ititle = user message.  (character*20) 
c 
c            internal 
c              dtatme = date and time for printing.  (character*24)
c 
c        nonsystem subroutines called 
c            none. 
c 
      character*24 dtatme
      character*20 ititle
c 
      call fdate(dtatme) 
c     if(kfilds.ne.kfildo)write(kfilds,110)ititle,dtatme
      write(kfildo,110)ititle,dtatme
 110  format(' ',a20,'   date/time:  ',a24) 
      return 
      end 
