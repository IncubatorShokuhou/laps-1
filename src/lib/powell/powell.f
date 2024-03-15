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

      subroutine powell(p,xi,n,np,ftol,iter,fret,func,print_switch)
      parameter (nmax=40,itmax=5)
      external func
      integer print_switch
      real func                 ! funciton type
      dimension p(np),xi(np,np),pt(nmax),ptt(nmax),xit(nmax)
      fret=func(p)
      if(fret.eq.0.0) then      !notify
         if (print_switch .eq. 1) then
            write(6,*)'powell:fret = 0.0'
         endif
      endif
      do 11 j=1,n
         pt(j)=p(j)
 11   continue
      iter=0
 1    iter=iter+1
      fp=fret
      ibig=0
      del=0.
      do 13 i=1,n
         do 12 j=1,n
            xit(j)=xi(j,i)
 12      continue
         fptt=fret
         call linmin(p,xit,n,fret)
         if(abs(fptt-fret).gt.del)then
            del=abs(fptt-fret)
            ibig=i
         endif
 13   continue
      if(2.*abs(fp-fret).le.ftol*(abs(fp)+abs(fret)))then
c         write(6,*) 'powell difference less than ftol'
c         write(6,*) fp, fret, ftol,'fp, fret,ftol'
         return
      endif
c     if(iter.eq.itmax) pause 'powell exceeding maximum iterations.'
      if(iter.eq.itmax) then
         if (print_switch .eq. 1) then
            write(6,*) 'powell exceeding maximum iterations.'
         endif
         return
      endif
      do 14 j=1,n
         ptt(j)=2.*p(j)-pt(j)
         xit(j)=p(j)-pt(j)
         pt(j)=p(j)
 14   continue
      fptt=func(ptt)
      if(fptt.ge.fp)go to 1
      t=2.*(fp-2.*fret+fptt)*(fp-fret-del)**2-del*(fp-fptt)**2
      if(t.ge.0.)go to 1
      call linmin(p,xit,n,fret)
      do 15 j=1,n
         xi(j,ibig)=xit(j)
 15   continue
      go to 1
      return
      end
      
