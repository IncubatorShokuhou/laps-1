
       program ingest_aircraft

!      driver for aircraft ingest (pin intermediate file)

!      steve albers      may-1999       original version

       write(6,*)
       write(6,*)' call ingest_pireps'
       call ingest_pireps(istatus)

       write(6,*)
       write(6,*)' call ingest_acars'
       call ingest_acars(istatus)

       write(6,*)
       write(6,*)' call ingest_wisdom'
       call ingest_wisdom(istatus)

       end

