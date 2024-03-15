       
         subroutine madis_qc_r(var,dd,level_qc,badflag)

         real var
         character*1 dd
         real badflag

         if(level_qc .ge. 1)then

             if(dd .eq. 'x')then ! failed level 1 qc
                 var = badflag
             endif

             if(level_qc .ge. 2)then
                 if(dd .eq. 'q')then ! failed level 2 or 3 qc
                     var = badflag
                 endif
             endif             

         endif

         if(dd .eq. 'b')then ! subjective qc (reject list)
             var = badflag
         endif
 
         return
         end
                 

         subroutine madis_qc_b(var,iqc,ibmask,idebug,badflag)

         integer ibmask(8)

         logical l_bit(8), btest, l_pass

         l_pass = .true.

         if(iqc .eq. -99)then ! missing value of qc word (bit flags)
             return
         endif

         if(var .eq. badflag)then ! already missing or flagged as bad
             return
         endif

         do i = 1,8
             l_bit(i) = btest(iqc,i-1)
         enddo

         do i = 1,8
             if(l_bit(i) .and. ibmask(i) .eq. 1)then
                 if(idebug .eq. 1)then
                     write(6,*)' failed madis_qc_b test ',iqc,i,var
                 endif
                 var = badflag
                 l_pass = .false.
             endif
         enddo

         if(l_pass .and. idebug .eq. 1)then
             write(6,*)' passed madis_qc_b tests ',iqc,var
         endif

         return
         end          
