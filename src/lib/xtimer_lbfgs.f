c                                                                                      
c  l-bfgs-b is released under the “new bsd license” (aka “modified bsd license”        
c  or “3-clause license”)                                                              
c  please read attached file license.txt                                               
c
c  yuanfu xie (2013) changes double precision to real for saving memory
c                                        
      subroutine timer_lbfgs(ttime)
      real ttime
c
      real temp
c
c     this routine computes cpu time in double precision; it makes use of 
c     the intrinsic f90 cpu_time therefore a conversion type is
c     needed.
c
c           j.l morales  departamento de matematicas, 
c                        instituto tecnologico autonomo de mexico
c                        mexico d.f.
c
c           j.l nocedal  department of electrical engineering and
c                        computer science.
c                        northwestern university. evanston, il. usa
c
c           yuanfu xie   change this routine name to timer_lbfgs from timer
c                        to avoid confusion with laps/intel compiler routine.
c                        december 9, 2013
c                         
c                        january 21, 2011
c
      temp = ttime
      call cpu_time(temp)
      ttime = temp

      return

      end
      
