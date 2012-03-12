module temp_like

  implicit none

  real*8, dimension(:), allocatable :: X
  real*8,  dimension(:,:), allocatable :: c_inv
  integer :: Nspec,nX,num_ells,nXfromdiff
  integer, dimension(:), allocatable  :: lminX, lmaxX, np, npt
  real*8 :: sz_100_temp(0:5000),sz_143_temp(0:5000)
  integer :: lmax_sz

  logical :: needinit=.true.

contains


  subroutine like_init_frommem(iNspec, inX,ilminX,ilmaxX,inp,inpt,ilmax_sz, sz_100,sz_143,mc_inv,mX)
    integer, intent(in):: iNspec, iNx,ilmax_sz
    integer,dimension(1:iNspec),intent(in)::ilmaxX, ilminX,inp,inpt
    real(8),dimension(0:ilmax_sz),intent(in)::sz_100,sz_143
    real(8),dimension(1:iNx),intent(in)::mX
    real(8),dimension(1:iNx,1:iNx),intent(in)::mc_inv

    Nspec = iNspec
    nX = inX

    allocate(lminX(Nspec))
    allocate(lmaxX(Nspec))
    allocate(np(Nspec))
    allocate(npt(Nspec))
    allocate(X(nX))
    allocate(c_inv(nX,nX))
    
    lminX = ilminX(:Nspec)
    lmaxX = ilmaxX(:Nspec)
    np = inp(:Nspec)
    npt = inpt(:Nspec)
    
    X(:nX) = mX(:nX)
    c_inv(:nX,:nX) = mc_inv(:nX,:nX)
    
    lmax_sz = ilmax_sz
    if(lmax_sz>5000) then
      print*, ' you need to increase the sizes of sz_temp', lmax_sz
      stop
    endif      

    sz_100_temp(0:lmax_sz) = sz_100(0:lmax_sz)
    sz_143_temp(0:lmax_sz) = sz_143(0:lmax_sz)
    
    needinit=.false.

  end subroutine like_init_frommem  


  subroutine like_init(like_file, sz100_file, sz143_file)

    integer :: i, j, l
    character*100 like_file, sz100_file, sz143_file
    integer:: iNspec, iNx,ilmax_sz
    integer,dimension(:),allocatable::ilmaxX, ilminX,inp,inpt
    real(8),dimension(0:5000)::sz_100,sz_143
    real(8),dimension(:),allocatable::mX
    real(8),dimension(:,:),allocatable::mc_inv
    

    if(needinit==.false.) then
       return
    endif 


    open(48, file=like_file, form='unformatted', status='unknown')

    read(48) iNspec,inX
    allocate(lminX(Nspec))
    allocate(lmaxX(Nspec))
    allocate(np(Nspec))
    allocate(npt(Nspec))
    allocate(X(nX))
    allocate(c_inv(nX,nX))

    read(48) (ilminX(i), ilmaxX(i), inp(i), inpt(i), i = 1, iNspec)
    read(48) (mX(i), i=1, inX)
    read(48) 
    read(48) ((mc_inv(i, j), j = 1, inX), i = 1,  inX)
    close(48)

    open(48, file=sz100_file, form='unformatted', status='unknown')
    read(48) ilmax_sz
    read(48) (sz_100(l), l = 0, ilmax_sz)
    close(48)

    open(48, file=sz143_file, form='unformatted', status='unknown')
    read(48) ilmax_sz
    read(48) (sz_143(l), l = 0, ilmax_sz)
    close(48)

    call like_init_frommem(iNspec, inX,ilminX,ilmaxX,inp,inpt,ilmax_sz, sz_100,sz_143,mc_inv,mX)
    
    deallocate(ilminX)
    deallocate(ilmaxX)
    deallocate(inp)
    deallocate(inpt)
    deallocate(mX)
    deallocate(mc_inv)


    return
  end subroutine like_init



  subroutine calc_like(zlike,  cell_cmb, A_ps_100,  A_ps_143, A_ps_217, A_cib_143, A_cib_217, A_sz,  &
       r_ps, r_cib, cal0, cal1, cal2)

    integer :: i, j, l, ipause
    real*8, dimension(:),  allocatable ::  X_theory, X_f, X_data, Y 
    real*8, dimension(0:) :: cell_cmb
    real*8 zlike, A_ps_100, A_ps_143, A_ps_217, A_cib_143, A_cib_217, A_sz, r_ps, r_cib, &
         cal0, cal1, cal2
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



    if(Nspec.ne.4) then
       print*, 'Nspec inconsistent with foreground corrections in calc_like.'
       stop
    end if

!   100 foreground
!
      do l = lminX(1), lmaxX(1)
       zell = dfloat(l)
       X_f(l - lminX(1) + 1) = A_ps_100*1.d-6 + &
          A_sz*sz_100_temp(l)  
       X_data(l - lminX(1) + 1) = cal0*X(l - lminX(1) + 1)
       X_theory(l-lminX(1) + 1) = cell_cmb(l)
      end do

!   143 foreground
!
      do l = lminX(2), lmaxX(2)
       zell = dfloat(l)
       zCIB = A_cib_143*(dfloat(l)/3000.)**(0.7)/dfloat(l*(l+1))
       X_f(l - lminX(2) + npt(2)) = A_ps_143*1.d-6 + zCIB + &
          A_sz*sz_143_temp(l) 
       X_data(l - lminX(2) +npt(2)) = cal1*X(l - lminX(2) + npt(2))
       X_theory(l-lminX(2) + npt(2)) = cell_cmb(l)
      end do

!
!   217 foreground
!
      do l = lminX(3), lmaxX(3)
       zell = dfloat(l)
       zCIB = A_cib_217*(dfloat(l)/3000.)**(0.7)/dfloat(l*(l+1))
       X_f(l - lminX(3) + npt(3) ) = A_ps_217*1.d-6 + zCIB 
       X_data(l - lminX(3) + npt(3)) = cal2*X(l - lminX(3) + npt(3))
       X_theory(l-lminX(3) + npt(3)) = cell_cmb(l)
      end do

!!
!!   100x143 foreground
!!
!      do l = lminX(4), lmaxX(4)
!       zell = dfloat(l)
!       zCIB = A_cib_217*(dfloat(l)/3000.)**(0.7)/dfloat(l*(l+1))
!       X_f(l - lminX(4) + npt(4) ) = dsqrt(A_ps_100*A_ps_143)*1.d-6 &
!        +  A_sz*dsqrt(sz_143_temp(l)*sz_100_temp(l))
!       X_data(l - lminX(4) + npt(4)) = dsqrt(cal1*cal0)*X(l - lminX(4) + npt(4))
!       X_theory(l-lminX(4) + npt(4)) = cell_cmb(l)
!      end do


!
!   143x217 foreground
!
      do l = lminX(4), lmaxX(4)
       zell = dfloat(l)
       zCIB = dsqrt(A_cib_143*A_cib_217)*(dfloat(l)/3000.)**(0.7) &
         /dfloat(l*(l+1))
       X_f(l - lminX(4) + npt(4) ) = &
          r_ps*dsqrt(A_ps_143*A_ps_217)*1.d-6 + r_cib*zCIB  
       X_data(l - lminX(4) + npt(4)) =  dsqrt(cal1*cal2)*X(l - lminX(4) + npt(4))
       X_theory(l-lminX(4) + npt(4)) = cell_cmb(l)
      end do


    do i = 1, nX
       Y(i) = X_data(i) - X_theory(i) - X_f(i)
    end do

    zlike = 0.d+00
    do  j = 1, nX
       ztemp=0.d0
       do  i = 1, nX
          ztemp = ztemp + Y(i)*c_inv(i, j)
       end do
       zlike=zlike+ztemp*Y(j)
    end do

    zlike=zlike+((cal2/cal1-1.0056d0)/0.00105d0)**2 &
               +((cal0/cal1-1.0127d0)/0.0005d0)**2


    deallocate(X_theory)
    deallocate(X_data)      
    deallocate(X_f)      
    deallocate(Y)

  end subroutine calc_like


end module temp_like
