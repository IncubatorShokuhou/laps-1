
        subroutine moon_brwn(t,mx,my,mz)

c       t is the julian date (tdt)
c       mx, my, mz are the geocentric coordinates of the moon using the
c               mean equinox and ecliptic of date (au)

        include 'trigd.inc'

        implicit real*8 (a-z)

        include '.././include/astparms.for'

c       integer decg,rah,minpx,latg,longg

        data    ll0  /   270.434164 d0/,
     1          ll1  /481267.8831   d0/,
     1  ll2  /     -.001133 d0/,
     1  ll3  /      .0000019d0/

        data    w0  /   334.3295556d0/,
     1  w1  /  4069.0340333d0/,
     1  w2  /     -.010325 d0/,
     1  w3  /     -.0000125d0/

        data    n0  /   259.183275 d0/,
     1  n1  / -1934.1420083d0/,
     1  n2  /      .0020777d0/,
     1  n3  /      .0000022d0/

        data    llp0 /   279.6966777d0/,
     1  llp1 / 36000.768925 d0/,
     1  llp2 /      .0003027d0/

        data    wp0 /   281.2208444d0/,
     1  wp1 /     1.719175 d0/,
     1  wp2 /      .000452 d0/

        o = obliquity(t)

c       write(6,*)' o = ',o

        t1=(t-2415020.0d0)/36525.d0
        t2=t1**2
        t3=t1**3

        ll =  ll0 +  ll1*t1 +  ll2*t2 +  ll3*t3
        w  =   w0 +   w1*t1 +   w2*t2 +   w3*t3
        n  =   n0 +   n1*t1 +   n2*t2 +   n3*t3
        llp= llp0 + llp1*t1 + llp2*t2
        wp =  wp0 +  wp1*t1 +  wp2*t2

c       l      =  ll - w
c       lp     =  llp - wp
c       f      =  ll - n
c       d      =  ll - llp

c       l_deg  =  dmod(l+360.d4,360.d0)
c       lp_deg =  dmod(lp+360.d4,360.d0)
c       f_deg =   dmod(f+360.d4,360.d0)
c       d_deg  =  dmod(d+360.d4,360.d0)
c       ll_deg = dmod(ll,360.d0)
c       w_deg  = dmod( w,360.d0)
c       n_deg  = dmod( n+360.d4,360.d0)
c       wp_deg = dmod(wp,360.d0)


c       write(6,1)ll_deg,w_deg,n_deg,lp_deg,wp_deg,d_deg,l_deg,f_deg
c1      format(1x,'mean values: ll,w,n,lp,wp,d,l,f'/8f10.4/)

!       additive terms
        venus_term = 346.560 + 132.870 * t1 - .0091731 * t2

        nadd =  95.96 * sind(n)
     1      + 15.58 * sind(n -2.3 * (t1-18.5)+276.2)

        n = n + nadd/3600.

        lladd = + 14.27  * sind(venus_term)
     1  + 10.71  * sind(140.0   * (t1-18.5)+170.7)
     1  +  7.261 * sind(n)

        ll = ll + lladd/3600.

        wp = wp
     1  +  2.076 * sind(n) / 3600.

        l      =  ll - w
        lp     =  llp - wp
        f      =  ll - n
        d      =  ll - llp

c       l_deg  =  dmod(l+360.d4,360.d0)
c       lp_deg =  dmod(lp+360.d4,360.d0)
c       f_deg =   dmod(f+360.d4,360.d0)
c       d_deg  =  dmod(d+360.d4,360.d0)
c       ll_deg = dmod(ll,360.d0)
c       w_deg  = dmod( w,360.d0)
c       n_deg  = dmod( n+360.d4,360.d0)
c       wp_deg = dmod(wp,360.d0)


c       write(6,2)ll_deg,w_deg,n_deg,lp_deg,wp_deg,d_deg,l_deg,f_deg
c2      format(1x,'after additive terms: ll,w,n,lp,wp,d,l,f'/8f10.4/)



        l  = l  * rpd
        w  = w  * rpd
        n  = n  * rpd
        lp = lp * rpd
        wp = wp * rpd
        d  = d  * rpd
        f  = f  * rpd

        longterms = 0.
     1  +          13.902 *  sin(4.*d )
     1  +        2369.902 *  sin(2.*d )

        longterms = longterms
     1  +         191.953 *  sin(l + 2.*d )
     1  +       22639.500 *  sin(l        )
     1  -        4586.426 *  sin(l - 2.*d )
     1  -          38.428 *  sin(l - 4.*d )

        longterms = longterms
     1  -          24.420 *  sin(   lp + 2.*d)
     1  -         668.111 *  sin(   lp)
     1  -         165.145 *  sin(   lp -2.*d )
     1  -         125.154 *  sin(    d)

        longterms = longterms
     1  +          14.387 *  sin(2.*l + 2.*d )
     1  +         769.016 *  sin(2.*l )
     1  -         211.656 *  sin(2.*l - 2.*d )
     1  -          30.773 *  sin(2.*l - 4.*d )
     1  -         109.667 *  sin(l + lp)
     1  -         205.962 *  sin(l+lp - 2.*d )

        longterms = longterms
     1  +          14.577 *  sin(l - lp + 2.*d)
     1  +         147.693 *  sin(l - lp)
     1  +          28.475 *  sin(l - lp - 2.*d)

        longterms = longterms
     1  -           7.486 *  sin(2.*lp)
     1  -           8.096 *  sin(2.*lp - 2.*d)

        longterms = longterms
     1  -           5.471 *  sin(2.*f + 2.*d )
     1  -         411.608 *  sin(2.*f )
     1  -          55.173 *  sin(2.*f - 2.*d )

        longterms = longterms
     1  +          18.023 *  sin(lp +    d)

        longterms = longterms
     1  -           8.466 *  sin(l +    d)
     1  +          18.609 *  sin(l -    d)

        longterms = longterms
     1  +          36.124 *  sin(3.*l)
     1  -          13.193 *  sin(3.*l - 2.*d )

        longterms = longterms
     1  -           7.649 *  sin(2.*l + lp)
     1  -           8.627 *  sin(2.*l + lp -2.*d)
     1  +           9.703 *  sin(2.*l - lp)
     1  -           7.412 *  sin(l + 2.*lp - 2.*d)

        longterms = longterms
     1  -          45.099 *  sin(l + 2.*f)
     1  -           6.382 *  sin(l - 2.*f + 2.*d)
     1  +          39.532 *  sin(l - 2.*f)
     1  +           9.366 *  sin(l - 2.*f - 2.*d)

        longterms = longterms/3600.d0
        long = ll + longterms
!       0 correction mag 1.0017
!       long = long -.0072 ! correct longitude
        long = long +.0050 ! correct longitude (1.0120)
!       long = long +.0055 ! correct longitude (1.0118)
!       long = long +.0130 ! correct longitude (0.9965)
        long = dmod(long+360.d0,360.d0)

c       calculate latitude

        sterms = 0.
     1  -         112.79 *   sin(   d)
     1  +        2373.36 *   sin(2.*d)
     1  +          14.06 *   sin(4.*d)
     1  +         192.72 *   sin(l + 2.*d)
     1  +       22609.07 *   sin(l       )
     1  -        4578.13 *   sin(l - 2.*d)
     1  -          38.64 *   sin(l - 4.*d)

        sterms = sterms
     1  +          14.78 *   sin(2.*l + 2.*d)
     1  +         767.96 *   sin(2.*l       )
     1  -         152.53 *   sin(2.*l - 2.*d)
     1  -          34.07 *   sin(2.*l - 4.*d)

        sterms = sterms
     1  +          50.64 *   sin(3.*l       )
     1  -          16.40 *   sin(3.*l - 2.*d)

        sterms = sterms
     1  -          25.10 *   sin(lp + 2.*d)
     1  +          17.93 *   sin(lp +    d)
     1  -         126.98 *   sin(lp       )
     1  -         165.06 *   sin(lp - 2.*d)

        sterms = sterms
     1  -          16.35 *   sin(2.*lp - 2.*d)

        sterms = sterms
     1  -          11.75 *   sin(l+lp + 2.*d)
     1  -         115.18 *   sin(l+lp       )
     1  -         182.36 *   sin(l+lp - 2.*d)
     1  -           9.66 *   sin(l+lp - 4.*d)

        sterms = sterms
     1  -          23.59 *   sin(lp-l + 2.*d)
     1  -         138.76 *   sin(lp-l)
     1  -          31.70 *   sin(lp-l - 2.*d)

        sterms = sterms
     1  -          10.56 *   sin(2.*l +lp)
     1  -           7.59 *   sin(2.*l +lp - 2.*d)

        sterms = sterms
     1  +          11.67 *   sin(2.*l -lp)

        sterms = sterms
     1  -           6.12 *   sin(l + 2.*lp -2.*d)

        sterms = sterms
     1  -          52.14 *   sin(2.*f -2.*d)

        sterms = sterms
     1  -           9.52 *   sin(l + 2.*f -2.*d)

        nterms = 0.
     1  -         526.069 *  sin(f - 2.*d)
     1  -           3.352 *  sin(f - 4.*d)
     1  +          44.297 *  sin(f+l-2.*d)
     1  -          30.598 *  sin(f-l-2.*d)
     1  -          24.649 *  sin(f-2.*l)
     1  -          22.571 *  sin(f+lp-2.*d)
     1  +          20.599 *  sin(f-l)
     1  -           6.000 *  sin(f+l-4.*d)
     1  -           2.000 *  sin(f-2.*l-2.*d)
     1  +          10.985 *  sin(f-lp-2.*d)

        gamma1_c = -1.540 * cos(-l +lp -2.*d)
        omega1 = .0004664 * cos(n)
        omega2 = .0000754 * cos(n+ 275.05 - 2.30*t1)

        s = f + sterms/3600.d0*rpd

        lat = (18518.511+1.189+gamma1_c) * sin(s)
     1          - 6.241 * sin(3.*s) + nterms
        lat = lat * (1.-omega1-omega2)/3600.
!       lat = lat + .0026 ! correct latitude (1.0120)
!       lat = lat + .0027 ! correct latitude (1.0122)
!       lat = lat + .0029 ! correct latitude (1.0126)
!       lat = lat + .0033 ! correct latitude (1.0135)
!       lat = lat + .0037 ! correct latitude (1.0143)
        lat = lat + .0038 ! correct latitude (1.0144)

c       calculate distance

        sinepx = 0.
     1  +             .2607 * cos(4.*d)
     1  +           28.2333 * cos(2.*d)
     1  +         3422.7000

        sinepx = sinepx
     1  +            3.0861 * cos(l + 2.*d)
     1  +          186.5398 * cos(l       )
     1  +           34.3117 * cos(l - 2.*d)
     1  +             .6008 * cos(l - 4.*d)

        sinepx = sinepx
     1  -             .3000 * cos(lp+ 2.*d)
     1  -             .3997 * cos(lp      )
     1  +            1.9178 * cos(lp- 2.*d)

        sinepx = sinepx
     1  -             .9781 * cos(d)

        sinepx = sinepx
     1  +             .2833 * cos(2.*l + 2.*d)
     1  +           10.1657 * cos(2.*l       )
     1  -             .3039 * cos(2.*l - 2.*d)
     1  +             .3722 * cos(2.*l - 4.*d)

        sinepx = sinepx
     1  -             .9490 * cos(l+lp)
     1  +            1.4437 * cos(l+lp-2.*d)
     1  +             .2302 * cos(l-lp+2.*d)
     1  -             .2257 * cos(l-lp-2.*d)

        sinepx = sinepx
     1  -             .1052 * cos(2.*f - 2.*d)

        sinepx = sinepx
     1  -             .1093 * cos(l+d)

        sinepx = sinepx
     1  +             .1494 * cos(lp + d)

        sinepx = sinepx
     1  +            1.1528 * cos(l-lp)

        sinepx = sinepx
     1  +             .6215 * cos(3.*l)
     1  -             .1187 * cos(3.*l - 2.*d)

        sinepx = sinepx
     1  -             .1038 * cos(2.*l + lp)
     1  +             .1268 * cos(2.*l - lp)
     1  -             .0833 * cos(l + 2.*f - 2.*d)

        sinepx = sinepx
     1  -             .7136 * cos(l-2.*f)

c       minpx = idint(sinepx/60.d0)
c       secpx = sinepx - minpx*60.d0
c       write(6,20)long,lat,minpx,secpx
c20     format(1x,'long,lat,sinepx ',2f12.6,i4,f8.3)

c       latd = dabs(lat)
c       latg = idint(latd)
c       latm = (latd - latg) * 60.d0
c       latg = latg * lat/dabs(lat)
c       lats = (latm - int(latm)) * 60.d0

c       longd = dabs(long)
c       longg = idint(longd)
c       longm = (longd - longg) * 60.d0
c       longg = longg * long/dabs(long)
c       longs = (longm - int(longm)) * 60.d0

c       write(6,29)latg,int(latm),lats,longg,int(longm),longs
c29     format(/1x,'coords of date:   lat',
c       1       i4,i3,f6.2,'   long',i4,i3,f6.2/)

        long = long * rpd
        lat  = lat  * rpd
        sinepx=(sinepx/3600.)*rpd
        r     =   1./sinepx * 4.2635d-5

        dec  = asin (cos(o)*sin(lat)+sin(o)*cos(lat)*sin(long))
        ra   =-asin((sin(o)*sin(lat)-cos(o)*cos(lat)*sin(long))/cos(dec)
     1)
        if(cos(long).lt.0.)ra = pi-ra
        ra = dmod(ra+2.d0*pi,2.d0*pi)

c       decd = dabs(dec) / rpd
c       decg = idint(decd)
c       decm = (decd - decg) * 60.d0
c       decg = decg * dec/dabs(dec)
c       decs = (decm - int(decm)) * 60.d0

c       rad  = (ra / rpd)/15.d0
c       rah  = idint(rad)
c       ram  = (rad - rah)*60.d0
c       ras  = (ram - int(ram)) * 60.d0

c       write(6,30)decg,int(decm),decs,rah,int(ram),ras,r
c30     format(/1x,'coords of date:   dec',
c       1       i4,i3,f6.2,'   ra',i4,i3,f6.2,'  r',f15.6/)

        mx = r * cos(dec) * cos(ra)
        my = r * cos(dec) * sin(ra)
        mz = r * sin(dec)

        return
        end

