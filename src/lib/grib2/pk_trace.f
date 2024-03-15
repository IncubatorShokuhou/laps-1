      subroutine pk_trace(kfildo,jer,ndier,kier,nerr,isevere)
c
c        june    2000   lawrence   error posting and tracing for grib2
c        january 2001   glahn      comments; added kier=kier-1; 
c                                  changed ier to jer
c 
c        purpose 
c            inserts an error code and its severity into jer( , ).
c            inserts 999 when jer( , ) is full.  in that case,  
c            two error codes are lost.  since the program proceeds,
c            more may be lost before program completion.
c
c        data set use 
c           kfildo - unit number for output (print) file. (output) 
c
c        variables 
c             kfildo = unit number for output (print) file.  (input)
c           jer(j,k) = array of errors (j=1,ndier) (k=1,2)
c                      (input/output)
c              ndier = first dimension of jer( , ). the maximum number
c                      of errors allowed to be reported in jer( , ).
c                      (input)
c               kier = number of values in jer( , ) upon exit.
c                      (input/output)
c               nerr = value to insert into jer(kier,1).  (input)
c            isevere = the severity level of the diagnostic to
c                      insert into jer(kier,2).  (input)
c                      valid severity levels are:
c                      0 = not an error
c                      1 = warning
c                      2 = fatal
c
c       local variables
c               none
c
c        non system subroutines called 
c           none
c
      dimension jer(ndier,2)
c
      kier=kier+1
c
      if(kier.le.ndier)then
         jer(kier,1)=nerr
         jer(kier,2)=isevere
      else
         kier=kier-1
         jer(kier,1)=999
         jer(kier,2)=2
      endif
c
      return
      end
