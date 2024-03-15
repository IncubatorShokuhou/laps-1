
        function obliquity(t)

c       t is the julian date (tdt)

        implicit real*8 (a-z)

        pi = 3.1415926535897932d0
        rpd = pi/180.d0

        tu = (t-2451545d0)/36525d0
        o = 23.439291d0 + tu *
     1  (-.0130042d0 + tu * (-.00000016d0 + tu * .000000504d0))

        obliquity = o * rpd

        return
        end

