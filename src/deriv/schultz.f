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
cdis
cdis   
cdis
!--------------------------------------------------------------------------
!cloud water stays supercooled well below freezing, but there's almost zero
!supercooled water at -20 (253k).  pristine crystals melt pretty fast above
!freezing.

        subroutine convc2p (maxrate, t, rate)
        implicit none
        real maxrate, t, rate
        real pwr
        data pwr/2./

        rate = 0.

        if (t .lt. 253.) then
         rate = maxrate
         return
        else if (t .lt. 267.) then
         rate = maxrate * ((267.-t)/(267.-253.))**pwr
         return
        else if (t .lt. 273.15) then
         rate = 0.
         return
        else if (t .lt. 278.) then
         rate = -maxrate * (t-273.15)/(278.-273.15)
         return
        else
         rate = -maxrate
         return
        end if

        return
        end

!--------------------------------------------------------------------------
!this includes both autoconversion and collision/coalescence (collection).
!obviously, the former is just a function of cloud water content, but the
!latter is also dependent on the presence of rain.  note that collection
!can cause the rate to exceed the "maxrate" if there are very high
!concentrations of rain.

        subroutine convc2r (maxrate, qc, qcmin, qr, rate)
        implicit none
        real maxrate, qc, qcmin, qr, rate

        rate = 0.

        if (qc .lt. qcmin) then
         return
        else if (qc .lt. .0015) then
         rate = maxrate * (qc-qcmin)/(.0015-qcmin) * (1.+qr/.002)
         return
        else
         rate = maxrate * (1. + qr/.002)
         return
        end if

        return
        end

!-----------------------------------------------------------------------
!this includes both autoconversion and aggregration.  the former is just
!a function of cloud ice content, but the latter is also dependent on the
!presence of snow.  note that collection can cause the rate to exceed the
!"maxrate" if there are very high concentrations of snow.

        subroutine convp2s (maxrate, qp, qpmin, qs, rate)
        implicit none
        real maxrate, qp, qpmin, qs, rate

        rate = 0.

        if (qp .lt. qpmin) then
         return
        else if (qp .lt. .0015) then
         rate = maxrate * (qp-qpmin)/(.0015-qpmin) * (1.+qs/.002)
         return
        else
         rate = maxrate * (1. + qs/.002)
         return
        end if

        return
        end

!--------------------------------------------------------------------------
!riming.  note that the rate can exceed the "maxrate" if there are very high
!concentrations of snow.

        subroutine convc2s (maxrate, qc, qcmin, qs, rate)
        implicit none
        real maxrate, qc, qcmin, qs, rate

        rate = 0.
        if (qs .lt. .0000001) return    ! no autoconversion

        if (qc .lt. qcmin) then
         return
        else if (qc .lt. .0015) then
         rate = maxrate * (qc-qcmin)/(.0015-qcmin) * (1.+qs/.002)
         return
        else
         rate = maxrate * (1. + qs/.002)
         return
        end if

        return
        end

!--------------------------------------------------------------------------
!riming of ice particles.  note that the rate can exceed the "maxrate" if
!there are very high concentrations of ice.

        subroutine convc2i (maxrate, qc, qcmin, qi, rate)
        implicit none
        real maxrate, qc, qcmin, qi, rate

        rate = 0.
        if (qi .lt. .0000001) return    ! no autoconversion

        if (qc .lt. qcmin) then
         return
        else if (qc .lt. .0015) then
         rate = maxrate * (qc-qcmin)/(.0015-qcmin) * (1.+qi/.002)
         return
        else
         rate = maxrate * (1. + qi/.002)
         return
        end if

        return
        end

!--------------------------------------------------------------------------
!melting only, since rain does not freeze into snow.

        subroutine convs2r (maxrate, t, rate)
        implicit none
        real maxrate, t, rate

        rate = 0.

        if (t .lt. 273.15) then
         return
        else if (t .lt. 283.) then
         rate = maxrate * (t-273.15)/(283.-273.15)
         return
        else
         rate = maxrate
         return
        end if

        return
        end

!--------------------------------------------------------------------------
!freezing and melting.  for now, this is the same algorithm as c2p.

        subroutine convr2i (maxrate, t, rate)
        implicit none
        real maxrate, t, rate
        real pwr
        data pwr/2./

        rate = 0.

        if (t .lt. 253.) then
         rate = maxrate
         return
        else if (t .lt. 267.) then
         rate = maxrate * ((267.-t)/(267.-253.))**pwr
         return
        else if (t .lt. 273.15) then
         rate = 0.
         return
        else if (t .lt. 278.) then
         rate = -maxrate * (t-273.15)/(278.-273.15)
         return
        else
         rate = -maxrate
         return
        end if

        return
        end

!--------------------------------------------------------------------------
!evaporation of rain; function of vapor deficit only.

        subroutine convr2v (maxrate, rv, rvsatliq, rate)
        implicit none
        real maxrate, rv, rvsatliq, rate

        rate = 0.
        if (rv .ge. rvsatliq) return

        if (rvsatliq .lt. .000001) return

        rate = maxrate * (rvsatliq-rv)/rvsatliq

        return
        end

!--------------------------------------------------------------------------
!evaporation of snow; function of vapor deficit only.

        subroutine convs2v (maxrate, rv, rvsatice, rate)
        implicit none
        real maxrate, rv, rvsatice, rate

        rate = 0.
        if (rv .ge. rvsatice) return

        if (rvsatice .lt. .000001) return

        rate = maxrate * (rvsatice-rv)/rvsatice

        return
        end

!--------------------------------------------------------------------------
!evaporation of precipitating ice (graupel, sleet, hail); function of vapor
!deficit only.  maybe the vapor deficit should be with respect to water...

        subroutine convi2v (maxrate, rv, rvsatice, rate)
        implicit none
        real maxrate, rv, rvsatice, rate

        rate = 0.
        if (rv .ge. rvsatice) return

        if (rvsatice .lt. .000001) return

        rate = maxrate * (rvsatice-rv)/rvsatice

        return
        end
