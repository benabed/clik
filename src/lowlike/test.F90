!===========================================================
program test_likelihood

!===========================================================

use planck_likelihood
use planck_options

implicit none

real(8), dimension(:), allocatable :: cl_tt,cl_te,cl_ee,cl_bb
character(LEN=128)          :: filename
real(8)                     :: like(num_pl),like_tot,expected_like_tot
integer                     :: lun, l, i, olun
integer :: tt_npix, teeebb_npix

!---------------------------------------------------

print *,""
print *,"Planck low-like test program"
print *,"==================================="
print *,""
!---------------------------------------------------
ttmin=2
temin=2
ttmax=32
temax=32
use_wmap_pol  = .false.
allocate( cl_tt(ttmin:ttmax) )
allocate( cl_te(ttmin:ttmax) )
allocate( cl_ee(ttmin:ttmax) )
allocate( cl_bb(ttmin:ttmax) )

cl_bb = 0d0

!---------------------------------------------------
! read in Cls
!---------------------------------------------------
filename = trim(Planck_data_dir)//'test_cls_v4.dat'

write(*,*)"Reading in Cls from: ",trim(filename)
call get_free_lun( lun )
open(unit=lun,file=filename,action='read',status='old')

do l=2,ttmax
   read(lun,*)i,cl_tt(l),cl_ee(l),cl_bb(l),cl_te(l)
enddo

close(lun)

!---------------------------------------------------
! put in likelihood options
!---------------------------------------------------

use_lowl_TT          = .true.
use_lowl_pol         = .true.

!---------------------------------------------------
! get likelihoods
!---------------------------------------------------
like =0.d0
use_wmap_pol  = .false.

call planck_lowlike_init
!call planck_print_options
use_wmap_pol  = .false.

call planck_lowlike_compute(cl_tt,cl_te,cl_ee,cl_bb,like)

!do i=1,num_pl
!  write (*,*) like(i)
!enddo
like_tot = sum(like(1:num_pl))

!---------------------------------------------------
! write outputs
!---------------------------------------------------
  print 1
  print 2
  print 1

  print 4, 'low-l TTTT gibbs        ', 2*like(ttlowllike)
  print 3, 'TT/TE/EE/BB low-l chi2  ', 2*like(lowllike)
  print 4, 'TT/TE/EE/BB low-l det   ', 2*like(lowldet)
  print 1
  print 4, 'TOTAL -2ln(L)           ', 2*like_tot
  print 1

  if(use_wmap_pol) then
     expected_like_tot = 1638.592658d0
  else
    expected_like_tot =  1375.238684d0
  endif
 
 print '(A,F13.6)', "Expected -2ln(L)         = ", expected_like_tot
  print '(A,F13.6)', "      Difference         = ", 2*like_tot-expected_like_tot
  print 1
  print *, ""
  print *, "Differences on the order of O(0.001) are normal between platforms."
  print *, ""

  stop
1 format ('------------------------------------------------------------------')
2 format ('Breakdown of -2ln(L)')
3 format (A24,' = ',F13.6,' for ',I6,' pixels')
4 format (A24,' = ',F13.6)
 end program test_likelihood
