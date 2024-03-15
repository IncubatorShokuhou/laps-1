module src_global

implicit none
include 'variablen.incf'

integer, parameter :: mpnt = 32
real*8 :: u(mpnt),v(mpnt),wt(mpnt),xkr(mpnt),dkr(mpnt)
complex*16 ::pjzhx(mrank1,mrank1,mpnt),&
             pjzjx(mrank1,mrank1,mpnt),&
             rbjz(mrank1,mpnt),&
             rbhx(mrank1,mpnt), &
             refrc, eps, tmt(m2rank,m2rank,mrank1), pvz, phy
real*8 :: rbjx(mrank1,mpnt), wvnm
integer :: kxysym
real :: scshp, deq, dmx, wcnxl, wcnzc, wcnlm, &
          temp, scmix, refre, refim
real :: caupp(2), &
    cdtyp(2),&
    caavr(2),&
    cadev(2),&
        calow(2)
integer :: canum(2), ndr
real ::  extm(4,4,mdr),&
     bmum(4,4,mdr), &
         stks(14,mdr)
real :: cdprb(mca2,mdr),&        !canting parameters
    cdpth(mca2,mdr),&
        cdpab(mca2,mdr),&
    cdpdd(mca2,mdr),&
    cdpab2(mca2,mdr),&
        cdpdd2(mca2,mdr),&
    cdpabd(mca2,mdr)
integer ::  cdpnum(mdr), n2rank, n2mode
complex*16 :: tmat(m2rank,m2rank)
!complex*16 :: aaa(m2rank,m2rank),&
!              b(m2rank,m2rank)
save pjzhx, pjzjx, rbjz, rbhx, rbjx, wvnm, refrc, eps, kxysym, tmat,  &
  u, v, wt, xkr, dkr, scshp, deq, dmx, wcnxl, wcnzc, wcnlm, refre, refim, &
       caupp, cdtyp, caavr, cadev, calow, temp, scmix,  cdprb, cdpth, &
        cdpab,&
    cdpdd,&
    cdpab2,&
        cdpdd2,&
    cdpabd,ndr,  &
        cdpnum, extm, bmum, stks,  tmt, n2rank, n2mode, pvz, phy

end module src_global
