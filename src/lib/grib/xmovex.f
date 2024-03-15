       subroutine xmovex(out,in,ibytes)
c
c      this subroutine may not be needed, its was in
c      assembler language to move data, it ran about three
c      times faster than a fortan do loop, it was used to
c      make sure the data to be unpacked was on a word boundary,
c      this may not be needed on some brands of computers.
c
       character * 1 out(*)
       character * 1 in(*)
c
       integer       ibytes
c
       do 100 i = 1,ibytes
         out(i) = in(i)
  100  continue  
c
       return    
       end

