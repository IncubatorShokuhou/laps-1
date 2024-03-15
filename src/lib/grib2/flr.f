      real function flr(value)
c
      if(value.ge.0)then
         flr=aint(value)
      else
c
         if(amod(value,1.0).eq.0)then
            flr=value
         else
            flr=aint(value)-1 
         endif
c
      endif
c
      return
      end
