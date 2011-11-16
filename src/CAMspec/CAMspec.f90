module temp_like

  implicit none

  real*8 X(5000), c_inv(5000,5000)
  integer :: lmin1, lmax1, lmin2, lmax2, lmin3, lmax3,np1, np2, np3, nX
  real*8 :: sz_temp(0:5000)
  integer :: lmax_sz
  real*8, parameter :: &
       zA =   36.0143  ,&          ! A
       zalpha =    0.558049,&      ! alpha
       zB = 0.112619 , &            ! B
       zbeta = 8.94585 , &          ! beta
       zell_c = 2.75599e+3 ,&       ! ell_c
       zgamma = 2.55989e-2  ,&      ! gamma
       zdelta = 7.29273          ! delta


  logical :: needinit=.true.

contains


  subroutine like_init(like_file, sz_file)

    integer :: i, j, l

    character*100 like_file, sz_file
    !
    !   read likelihood file
    !
    open(48, file=like_file, form='unformatted', status='unknown')
    read(48)  lmin1, lmax1, lmin2, lmax2, lmin3, lmax3, np1, np2, np3, nX
    if(nX.gt.5000) then
       print*, ' you need to increase the sizes of X and c_inv', nX
       stop
    end if
    read(48) (X(i), i = 1, nX)
    read(48) 
    read(48) ((c_inv(i, j), j = 1, nX), i = 1,  nX)
    close(48)

    open(48, file=sz_file, form='unformatted', status='unknown')
    read(48) lmax_sz
    read(48) (sz_temp(l), l = 0, lmax_sz)
    close(48)

    needinit=.false.

    return
  end subroutine like_init



  subroutine calc_like(zlike,  cell_cmb, A_ps_143, A_ps_217, A_cib_143, A_cib_217, A_sz,  &
       r_ps, r_cib, cal1, cal2)

    integer :: i, j, l, ipause
    real*8, dimension(:),  allocatable ::  X_theory, X_f, X_data, Y 
    real*8, dimension(0:) :: cell_cmb
    real*8 zlike, A_ps_143, A_ps_217, A_cib_143, A_cib_217, A_sz, r_ps, r_cib, &
         cal1, cal2
    real*8 zell, zGF, zCIB
    real*8 ztemp

    if (needinit) then
       print*, 'like_init should have been called before attempting to call calc_like.'
       stop
    end if
    allocate(X_theory(1:nX))      
    allocate(X_data(1:nX))      
    allocate(X_f(1:nX))      
    allocate(Y(1:nX))

    !
    !   143 foreground
    !
    do l = lmin1, lmax1
       zell = dfloat(l)
       zGF = zA*(100.d+00/zell)**zalpha +  &
            zB*(zell/1000.d+00)**zbeta/(1. + (zell/zell_c)**zgamma)**zdelta
       zGF = zGF/dfloat(l*(l+1))
       zGF = cal1*cal2*zGF*3.44*3.44/(3.44*3.44 - 1.)
       zCIB = A_cib_143*(dfloat(l)/3000.)**(0.7)/dfloat(l*(l+1))
       X_f(l - lmin1 + 1) = A_ps_143*1.d-6 + zCIB + A_sz*sz_temp(l) + zGF/3.44/3.44
       X_data(l - lmin1 + 1) = cal1*X(l - lmin1 + 1)
       X_theory(l-lmin1+1) = cell_cmb(l)
    end do

    if(lmax2.ne.0.and.lmax3.ne.0) then
       !
       !   217 foreground
       !
       do l = lmin2, lmax2
          zell = dfloat(l)
          zGF = zA*(100.d+00/zell)**zalpha +    &
               zB*(zell/1000.d+00)**zbeta/(1. + (zell/zell_c)**zgamma)**zdelta
          zGF = zGF/dfloat(l*(l+1))
          zGF = cal1*cal2*zGF*3.44*3.44/(3.44*3.44 - 1.)
          zCIB = A_cib_217*(dfloat(l)/3000.)**(0.7)/dfloat(l*(l+1))
          X_f(l - lmin2 + 1 + np1 ) = A_ps_217*1.d-6 + zCIB + zGF
          X_data(l - lmin2 + 1 + np1) = cal1*cal2*X(l - lmin2 + 1 + np1)
          X_theory(l - lmin2+ 1+ np1) = cell_cmb(l)
       end do

       !
       !   143x217 foreground
       !
       do l = lmin3, lmax3
          zell = dfloat(l)
          zGF = zA*(100.d+00/zell)**zalpha +    &
               zB*(zell/1000.d+00)**zbeta/(1. + (zell/zell_c)**zgamma)**zdelta
          zGF = zGF/dfloat(l*(l+1))
          zGF = dsqrt(cal1*cal2)*zGF*3.44*3.44/(3.44*3.44 - 1.)
          zCIB = dsqrt(A_cib_143*A_cib_217)*(dfloat(l)/3000.)**(0.7)/dfloat(l*(l+1))
          X_f(l - lmin3 + 1 + np1 + np2 ) =   &
               r_ps*dsqrt(A_ps_143*A_ps_217)*1.d-6 + r_cib*zCIB  + zGF/3.44
          X_data(l - lmin3 + 1 + np1 + np2) = cal1*dsqrt(cal2)*X(l - lmin3 + 1 + np1 + np2)
          X_theory(l-lmin3 + 1 + np1 + np2) = cell_cmb(l)
       end do
    end if

    do i = 1, nX
       Y(i) = X_data(i) - X_theory(i) - X_f(i)
    end do

    zlike = 0.d+00
    do  j = 1, nX
       ztemp=0.d0
       do  i = 1, nX
          ztemp = ztemp + Y(i)*c_inv(i, j)
!          zlike=zlike+Y(i)*Y(j)*c_inv(i,j)
       end do
       zlike=zlike+ztemp*Y(j)
    end do

    deallocate(X_theory)
    deallocate(X_data)      
    deallocate(X_f)      
    deallocate(Y)

  end subroutine calc_like


end module temp_like