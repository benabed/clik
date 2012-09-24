module temp_like

  implicit none

  real*8, dimension(:), allocatable :: X
  real*8,  dimension(:,:), allocatable :: c_inv
  integer :: Nspec,nX,num_ells,nXfromdiff
  integer, dimension(:), allocatable  :: lminX, lmaxX, np, npt
  real*8 :: sz_143_temp(0:5000)
  real*8 :: ksz_temp(0:5000), tszxcib_temp(0:5000)
  integer :: lmax_sz
  real*8, dimension(:,:), allocatable :: beam_cov_inv
  real*8, dimension(:,:,:), allocatable :: beam_modes ! mode#, l, spec#
  integer :: num_modes_per_beam,beam_lmax,beam_Nspec,cov_dim
  
  logical :: needinit=.true.

  !real*8 :: bestlike
  !integer :: tempinstance


  !character*100 :: bestroot,bestnum,bestname

contains

  !subroutine like_init(like_file, sz143_file, tszxcib_file, ksz_file, beam_file)
!  !  integer :: i, j, l,dummy
!  !  character*100 like_file, sz143_file, ksz_file, tszxcib_file, beam_file
!  !  integer,dimension(:),allocatable::ilminX,ilmaxX,inp,inpt
!  !  real*8, dimension(:), allocatable :: iX
!  !  real*8,  dimension(:,:), allocatable ::ic_inv
!  !  integer :: iNspec,inX
!  !  real*8 :: iksz_temp(0:5000), itszxcib_temp(0:5000)
!  !  real*8 :: isz_143_temp(0:5000)
!  !  real*8, dimension(:,:), allocatable :: ibeam_cov_inv
!  !  real*8, dimension(:,:,:), allocatable :: ibeam_modes ! mode#, l, spec#
!  !  integer :: inum_modes_per_beam,ibeam_lmax,ibeam_Nspec,icov_dim,ilmax_sz
!  !  
!
!  !  ! cl_ksz_148_tbo.dat file is in D_l, format l D_l, from l=2 to 10000
!  !  ! tsz_x_cib_template.txt is is (l D_l), from l=2 to 9999, normalized to unity 
!  !  !    at l=3000 
!
!  !  if(needinit==.false.) then
!  !     return
!  !  endif
!
!  !  open(48, file=like_file, form='unformatted', status='unknown')
!
!  !  read(48) iNspec,inX
!  !  allocate(ilminX(iNspec))
!  !  allocate(ilmaxX(iNspec))
!  !  allocate(inp(iNspec))
!  !  allocate(inpt(iNspec))
!  !  allocate(iX(inX))
!  !  allocate(ic_inv(inX,inX))
!
!  !  read(48) (ilminX(i), ilmaxX(i), inp(i), inpt(i), i = 1, iNspec)
!  !  read(48) (iX(i), i=1, inX)
!  !  read(48) 
!  !  read(48) ((ic_inv(i, j), j = 1, inX), i = 1,  inX)
!  !  close(48)
!
!  !  !  open(48, file=sz100_file, form='unformatted', status='unknown')
!  !  !  read(48) lmax_sz
!  !  !  read(48) (sz_100_temp(l), l = 0, lmax_sz)
!  !  !  close(48)
!
!  !  ilmax_sz=5000
!
!  !  open(48, file=sz143_file, form='formatted', status='unknown')
!  !  do i=2,ilmax_sz
!  !     read(48,*) dummy,isz_143_temp(i)
!  !  enddo
!  !  close(48)
!
!  !  open(48, file=ksz_file, form='formatted',status='unknown')
!  !  do i=2,ilmax_sz
!  !     read(48,*) dummy,iksz_temp(i)
!  !  enddo
!  !  close(48)
!
!  !  open(48, file=tszxcib_file,form='formatted',status='unknown')
!  !  do i=2,ilmax_sz
!  !     read(48,*) dummy,itszxcib_temp(i)
!  !  enddo
!  !  close(48)
!
!  !  open(48, file=beam_file, form='unformatted', status='unknown')
!  !  read(48) ibeam_Nspec,inum_modes_per_beam,ibeam_lmax
!  !  if(ibeam_Nspec.ne.iNspec) stop 'Problem: beam_Nspec != Nspec'
!  !  allocate(ibeam_modes(inum_modes_per_beam,0:ibeam_lmax,ibeam_Nspec))
!  !  icov_dim=ibeam_Nspec*inum_modes_per_beam
!  !  allocate(ibeam_cov_inv(icov_dim,icov_dim))
!  !  read(48) (((ibeam_modes(i,l,j),j=1,iNspec),l=0,ibeam_lmax),i=1,inum_modes_per_beam)
!  !  read(48) ! skipping beam_cov
!  !  read(48) ((ibeam_cov_inv(i,j),j=1,icov_dim),i=1,icov_dim)
!  !  close(48)
!
!  !  call like_init_frommem(iNspec, inX,ilminX,ilmaxX,inp,inpt, ic_inv,iX,ilmax_sz,isz_143_temp,iksz_temp,itszxcib_temp,ibeam_Nspec,inum_modes_per_beam,ibeam_lmax,icov_dim,ibeam_cov_inv,ibeam_modes)
!
  !  deallocate(ilmaxX)
  !  deallocate(ilminX)
  !  deallocate(inp)
  !  deallocate(inpt)
  !  deallocate(iX)
  !  deallocate(ic_inv)
  !  deallocate(ibeam_modes)
  !  deallocate(ibeam_cov_inv)
  !end subroutine

  subroutine like_init(like_file, sz143_file, tszxcib_file, ksz_file, beam_file)

    real*8, dimension(:), allocatable :: X
    real*8,  dimension(:,:), allocatable :: c_inv
    integer :: Nspec,nX,num_ells,nXfromdiff
    integer, dimension(:), allocatable  :: lminX, lmaxX, np, npt
    real*8 :: sz_143_temp(0:5000)
    real*8 :: ksz_temp(0:5000), tszxcib_temp(0:5000)
    integer :: lmax_sz
    real*8, dimension(:,:), allocatable :: beam_cov_inv
    real*8, dimension(:,:,:), allocatable :: beam_modes ! mode#, l, spec#
    integer :: num_modes_per_beam,beam_lmax,beam_Nspec,cov_dim
  
    integer :: i, j, l,dummy

    character*100 like_file, sz143_file, ksz_file, tszxcib_file, beam_file

    ! cl_ksz_148_tbo.dat file is in D_l, format l D_l, from l=2 to 10000
    ! tsz_x_cib_template.txt is is (l D_l), from l=2 to 9999, normalized to unity 
    !    at l=3000 

    if(needinit==.false.) then
       return
    endif

    open(48, file=like_file, form='unformatted', status='unknown')

    read(48) Nspec,nX
    allocate(lminX(Nspec))
    allocate(lmaxX(Nspec))
    allocate(np(Nspec))
    allocate(npt(Nspec))
    allocate(X(nX))
    allocate(c_inv(nX,nX))

    read(48) (lminX(i), lmaxX(i), np(i), npt(i), i = 1, Nspec)
    read(48) (X(i), i=1, nX)
    read(48) 
    read(48) ((c_inv(i, j), j = 1, nX), i = 1,  nX)
    close(48)

    !  open(48, file=sz100_file, form='unformatted', status='unknown')
    !  read(48) lmax_sz
    !  read(48) (sz_100_temp(l), l = 0, lmax_sz)
    !  close(48)

    lmax_sz=5000

    open(48, file=sz143_file, form='formatted', status='unknown')
    do i=2,lmax_sz
       read(48,*) dummy,sz_143_temp(i)
    enddo
    close(48)

    open(48, file=ksz_file, form='formatted',status='unknown')
    do i=2,lmax_sz
       read(48,*) dummy,ksz_temp(i)
    enddo
    close(48)

    open(48, file=tszxcib_file,form='formatted',status='unknown')
    do i=2,lmax_sz
       read(48,*) dummy,tszxcib_temp(i)
    enddo
    close(48)

    open(48, file=beam_file, form='unformatted', status='unknown')
    read(48) beam_Nspec,num_modes_per_beam,beam_lmax
    if(beam_Nspec.ne.Nspec) stop 'Problem: beam_Nspec != Nspec'
    allocate(beam_modes(num_modes_per_beam,0:beam_lmax,beam_Nspec))
    cov_dim=beam_Nspec*num_modes_per_beam
    allocate(beam_cov_inv(cov_dim,cov_dim))
    read(48) (((beam_modes(i,l,j),j=1,Nspec),l=0,beam_lmax),i=1,num_modes_per_beam)
    read(48) ! skipping beam_cov
    read(48) ((beam_cov_inv(i,j),j=1,cov_dim),i=1,cov_dim)
    close(48)
    call like_init_frommem(Nspec, nX,lminX,lmaxX,np,npt, c_inv,X,lmax_sz,sz_143_temp,ksz_temp,tszxcib_temp,beam_Nspec,num_modes_per_beam,beam_lmax,cov_dim,beam_cov_inv,beam_modes)

    !bestlike=1.e30
    !tempinstance=MPIrank+1
!
!    !write(bestnum,*) tempinstance
!    !bestroot='/home/stg20/src/cosmomc/chains/best'
    !bestname=trim(bestroot)//'_'//trim(adjustl(bestnum))//'.txt'


    deallocate(lminX)
    deallocate(lmaxX)
    deallocate(np)
    deallocate(npt)
    deallocate(X)
    deallocate(c_inv)
    deallocate(beam_cov_inv)
    deallocate(beam_modes)

    return
  end subroutine like_init

  subroutine like_init_frommem(iNspec, inX,ilminX,ilmaxX,inp,inpt, ic_inv,iX,ilmax_sz,isz_143_temp,iksz_temp,itszxcib_temp,ibeam_Nspec,inum_modes_per_beam,ibeam_lmax,icov_dim,ibeam_cov_inv,ibeam_modes)
    integer,intent(in)::iNspec,inX,inum_modes_per_beam,ibeam_lmax,ibeam_Nspec,icov_dim,ilmax_sz
    integer,dimension(:)::ilminX,ilmaxX,inp,inpt
    real*8, dimension(:) :: iX
    real*8,  dimension(:,:) ::ic_inv
    real*8,intent(in) :: iksz_temp(:), itszxcib_temp(:)
    real*8,intent(in) :: isz_143_temp(:)
    real*8, dimension(:,:) :: ibeam_cov_inv
    real*8, dimension(:,:,:) :: ibeam_modes ! mode#, l, spec#
    
    
    Nspec = iNspec
    nX = inX

    allocate(lminX(Nspec))
    allocate(lmaxX(Nspec))
    allocate(np(Nspec))
    allocate(npt(Nspec))
    allocate(X(nX))
    allocate(c_inv(nX,nX))

    lminX = ilminX
    lmaxX = ilmaxX
    X = iX
    c_inv = ic_inv
    np = inp
    npt = inpt
    
    lmax_sz = ilmax_sz
    ksz_temp = iksz_temp
    sz_143_temp = isz_143_temp
    tszxcib_temp = itszxcib_temp

    beam_Nspec = ibeam_Nspec
    num_modes_per_beam = inum_modes_per_beam
    beam_lmax = ibeam_lmax
    cov_dim = icov_dim

    allocate(beam_modes(num_modes_per_beam,0:beam_lmax,beam_Nspec))
    allocate(beam_cov_inv(cov_dim,cov_dim))

    beam_modes = ibeam_modes
    beam_cov_inv = ibeam_cov_inv

    !bestlike=1.e30
    !tempinstance=MPIrank+1

    !write(bestnum,*) tempinstance
    !bestroot='/home/stg20/src/cosmomc/chains/best'
    !bestname=trim(bestroot)//'_'//trim(adjustl(bestnum))//'.txt'

    needinit=.false.

  end subroutine like_init_frommem



  subroutine calc_like(zlike,  cell_cmb, A_ps_100,  A_ps_143, A_ps_217, A_cib_143, A_cib_217, A_sz,  &
       r_ps, r_cib, xi, A_ksz, cal0, cal1, cal2, beam_coeffs)

    integer :: i, j, l, ipause,ii,jj
    real*8, dimension(:),  allocatable ::  X_theory, X_f, X_data, X_beam_corr_model, Y 
    real*8, dimension(0:) :: cell_cmb
    real*8 zlike, A_ps_100, A_ps_143, A_ps_217, A_cib_143, A_cib_217, A_sz, r_ps, r_cib, &
         cal0, cal1, cal2, xi, A_ksz
    real*8 zell, zGF, zCIB
    real*8 ztemp

    real*8, dimension(:,:) :: beam_coeffs

    real*8 :: beamfac

    integer :: ie1,ie2,if1,if2

    
    if (needinit) then
       print*, 'like_init should have been called before attempting to call calc_like.'
       stop
    end if
    allocate(X_theory(1:nX))      
    allocate(X_data(1:nX))      
    allocate(X_f(1:nX))      
    allocate(X_beam_corr_model(1:nX))      
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
            A_ksz*ksz_temp(l)/dfloat(l*(l+1))+ &
            A_sz*2.022d0*sz_143_temp(l)/dfloat(l*(l+1))
       X_data(l - lminX(1) + 1) = X(l - lminX(1) + 1)
       X_theory(l-lminX(1) + 1) = cell_cmb(l)
       X_beam_corr_model(l-lminX(1)+1) = &
            ( X_theory(l - lminX(1) + 1)+ X_f(l - lminX(1) + 1))* &
            corrected_beam(1,l)/cal0
    end do

    !   143 foreground
    !
    do l = lminX(2), lmaxX(2)
       zell = dfloat(l)
       zCIB = 1.134d0*A_cib_143*(dfloat(l)/3000.)**(0.8)/dfloat(l*(l+1))
       X_f(l - lminX(2) + npt(2)) = A_ps_143*1.d-6 + zCIB + &
            A_ksz*ksz_temp(l)/dfloat(l*(l+1))+&
            A_sz*0.95d0*sz_143_temp(l)/dfloat(l*(l+1)) + &
            -2.0*sqrt(1.134d0*A_cib_143*0.95d0*A_sz*4.796)*xi*tszxcib_temp(l)/dfloat(l*(l+1))
       X_data(l - lminX(2) +npt(2)) = X(l - lminX(2) + npt(2))
       X_theory(l-lminX(2) + npt(2)) = cell_cmb(l) 
       X_beam_corr_model(l-lminX(2)+npt(2)) = &
            ( X_theory(l - lminX(2) + npt(2))+ X_f(l - lminX(2) + npt(2)))* &
            corrected_beam(2,l)/cal1
    end do

    !
    !   217 foreground
    !
    do l = lminX(3), lmaxX(3)
       zell = dfloat(l)
       zCIB = 1.33d0*A_cib_217*(dfloat(l)/3000.)**(0.8)/dfloat(l*(l+1))
       X_f(l - lminX(3) + npt(3) ) = A_ps_217*1.d-6 + zCIB &
            + A_ksz*ksz_temp(l)/dfloat(l*(l+1))   
       X_data(l - lminX(3) + npt(3)) = X(l - lminX(3) + npt(3))
       X_theory(l-lminX(3) + npt(3)) = cell_cmb(l)
       X_beam_corr_model(l-lminX(3)+npt(3)) = &
            ( X_theory(l - lminX(3) + npt(3))+ X_f(l - lminX(3) + npt(3)))* &
            corrected_beam(3,l)/cal2
    end do


    !
    !   143x217 foreground
    !
    do l = lminX(4), lmaxX(4)
       zell = dfloat(l)
       zCIB = 1.23d0*dsqrt(A_cib_143*A_cib_217)*(dfloat(l)/3000.)**(0.8) &
            /dfloat(l*(l+1))
       X_f(l - lminX(4) + npt(4) ) = &
            r_ps*dsqrt(A_ps_143*A_ps_217)*1.d-6 + r_cib*zCIB &
            +A_ksz*ksz_temp(l)/dfloat(l*(l+1))  &
            -sqrt(1.33d0*A_cib_217*0.95d0*A_sz*4.796)*xi*tszxcib_temp(l)/dfloat(l*(l+1))  
       X_data(l - lminX(4) + npt(4)) =  X(l - lminX(4) + npt(4))
       X_theory(l-lminX(4) + npt(4)) = cell_cmb(l)
       X_beam_corr_model(l-lminX(4)+npt(4)) = &
            ( X_theory(l - lminX(4) + npt(4))+ X_f(l - lminX(4) + npt(4)))* &
            corrected_beam(4,l)/dsqrt(cal1*cal2)
    end do

    do i = 1, nX
       Y(i) = X_data(i) - X_beam_corr_model(i)
    end do

    zlike = 0.d+00
    do  j = 1, nX
       ztemp=0.d0
       do  i = 1, nX
          ztemp = ztemp + Y(i)*c_inv(i, j)
       end do
       zlike=zlike+ztemp*Y(j)
    end do


    do if2=1,beam_Nspec
       do if1=1,beam_Nspec
          do ie2=1,num_modes_per_beam
             do ie1=1,num_modes_per_beam
                ii=ie1+num_modes_per_beam*(if1-1)
                jj=ie2+num_modes_per_beam*(if2-1)
                zlike=zlike+beam_coeffs(if1,ie1)*beam_cov_inv(ii,jj)*beam_coeffs(if2,ie2)
             enddo
          enddo
       enddo
    enddo


    zlike=zlike+((cal2/cal1-1.0056d0)/0.00105d0)**2 &
         +((cal0/cal1-1.0127d0)/0.0005d0)**2

    ! WARNING; weakening prior...
    !    zlike=zlike+((cal2/cal1-1.0056d0)/0.0063d0)**2 &
    !               +((cal0/cal1-1.0127d0)/0.003d0)**2


    ! correction to improve old way of handling cals...
    !zlike=zlike &
    !-2.d0*(lmaxX(1)-lminX(1)+1)*log(cal0) &
    !-2.d0*(lmaxX(2)-lminX(2)+1)*log(cal1) &
    !-2.d0*(lmaxX(3)-lminX(3)+1)*log(cal2) &
    !-2.d0*(lmaxX(4)-lminX(4)+1)*.5d0*(log(cal1)+log(cal2)) 



!!    if(zlike.lt.bestlike) then
!!       open(49,file=bestname,form='formatted',status='unknown')
!!       !do if1=1,beam_Nspec
!!       !write(49,101) (beam_coeffs(if1,ie1),ie1=1,num_modes_per_beam)
!!       !enddo
!!101    format(5(e12.4)) 
!!       do l = 0,2500
!!          write(49,*) l, l*(l+1)*cell_cmb(l)
!!       enddo
!!       close(49)
!!       bestlike=zlike
!!    endif

    deallocate(X_theory)
    deallocate(X_data)     
    deallocate(X_f)      
    deallocate(X_beam_corr_model)
    deallocate(Y)


  contains

    real*8 function corrected_beam(spec_num,l)
      integer, intent(in) :: spec_num,l
      integer :: i
      corrected_beam=1.d0
      do i=1,num_modes_per_beam   
         corrected_beam=corrected_beam+beam_coeffs(spec_num,i)*beam_modes(i,l,spec_num)
      enddo
    end function corrected_beam

  end subroutine calc_like



end module temp_like
