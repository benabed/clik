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

    if(needinit .eqv. .false.) then
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
    real(8)::renorm    
    
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
    renorm=1.d0/ksz_temp(3000)
    ksz_temp=ksz_temp*renorm
    
    sz_143_temp = isz_143_temp
    renorm=1.d0/sz_143_temp(3000)
    sz_143_temp=sz_143_temp*renorm

    tszxcib_temp = itszxcib_temp
    renorm=1.d0/tszxcib_temp(3000)
    tszxcib_temp=tszxcib_temp*renorm

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


  subroutine calc_like(zlike,  cell_cmb, freq_params)
    real(8), intent(in)  :: freq_params(:)
    real(8), dimension(0:) :: cell_cmb
    integer :: i, j, l, ii,jj
    real(8), dimension(:),  allocatable, save ::  X_theory, X_f, X_data, X_beam_corr_model, Y
    real(8) zlike, A_ps_100, A_ps_143, A_ps_217, A_cib_143, A_cib_217, A_sz, r_ps, r_cib, &
    cal0, cal1, cal2, xi, A_ksz, ncib
    real(8) zell, zCIB
    real(8) ztemp
    real(8) beam_coeffs(Nspec,num_modes_per_beam)
    integer :: ie1,ie2,if1,if2
    integer num_non_beam
    !    real(8) atime
    real(8) Lfac
    real(8), parameter :: sz_bandpass100_nom143 = 2.022d0
    real(8), parameter :: cib_bandpass143_nom143 = 1.134d0
    real(8), parameter :: sz_bandpass143_nom143 = 0.95d0
    real(8), parameter :: cib_bandpass217_nom217 = 1.33d0

    if (.not. allocated(lminX)) then
        print*, 'like_init should have been called before attempting to call calc_like.'
        stop
    end if
    if (.not. allocated(X_theory)) then
        allocate(X_theory(1:nX))
        allocate(X_data(1:nX))
        allocate(X_f(1:nX))
        allocate(X_beam_corr_model(1:nX))
        allocate(Y(1:nX))
    end if

    ! atime = MPI_Wtime()

    if(Nspec.ne.4) then
        print*, 'Nspec inconsistent with foreground corrections in calc_like.'
        stop
    end if

    num_non_beam = 14
    if (size(freq_params) < num_non_beam +  beam_Nspec*num_modes_per_beam) stop 'CAMspec: not enough parameters'
    A_ps_100=freq_params(1)
    A_ps_143 = freq_params(2)
    A_ps_217 = freq_params(3)
    A_cib_143=freq_params(4)
    A_cib_217=freq_params(5) 
    A_sz =freq_params(6)  !143
    r_ps = freq_params(7)
    r_cib=freq_params(8)
    ncib = freq_params(9)
    cal0=freq_params(10)
    cal1=freq_params(11) 
    cal2=freq_params(12)
    xi=freq_params(13)
    A_ksz=freq_params(14)
    do ii=1,beam_Nspec
        do jj=1,num_modes_per_beam
            beam_coeffs(ii,jj)=freq_params(num_non_beam+jj+num_modes_per_beam*(ii-1))
        enddo
    enddo

    
    !   100 foreground
    !
    do l = lminX(1), lmaxX(1)
        zell = real(l,8)
        Lfac = real(l*(l+1),8)
        X_f(l - lminX(1) + 1) = A_ps_100*1.d-6/9.d0 + &
        A_ksz*ksz_temp(l)/Lfac+ &
        A_sz*sz_bandpass100_nom143*sz_143_temp(l)/Lfac
        X_data(l - lminX(1) + 1) = X(l - lminX(1) + 1)
        X_theory(l-lminX(1) + 1) = cell_cmb(l)
        X_beam_corr_model(l-lminX(1)+1) = &
        ( X_theory(l - lminX(1) + 1)+ X_f(l - lminX(1) + 1))* &
        corrected_beam(1,l)/cal0
    end do

    !   143 foreground
    !
    do l = lminX(2), lmaxX(2)
        zell = real(l,8)
        Lfac = real(l*(l+1),8)
        zCIB = cib_bandpass143_nom143*A_cib_143*(zell/3000)**(ncib)/Lfac
        X_f(l - lminX(2) + npt(2)) = A_ps_143*1.d-6/9.d0 + zCIB + &
        A_ksz*ksz_temp(l)/Lfac + A_sz*sz_bandpass143_nom143*sz_143_temp(l)/Lfac  &
        -2.0*sqrt(cib_bandpass143_nom143*A_cib_143*sz_bandpass143_nom143*A_sz)*xi*tszxcib_temp(l)/Lfac
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
        zell = real(l,8)
        Lfac = real(l*(l+1),8)
        zCIB = cib_bandpass217_nom217*A_cib_217*(zell/3000)**(ncib)/Lfac
        X_f(l - lminX(3) + npt(3) ) = A_ps_217*1.d-6/9.d0 + zCIB &
        + A_ksz*ksz_temp(l)/Lfac
        X_data(l - lminX(3) + npt(3)) = X(l - lminX(3) + npt(3))
        X_theory(l-lminX(3) + npt(3)) = cell_cmb(l)
        X_beam_corr_model(l-lminX(3)+npt(3)) = &
        ( X_theory(l - lminX(3) + npt(3))+ X_f(l - lminX(3) + npt(3)))* &
        corrected_beam(3,l)/cal2
    end do

    !   143x217 foreground
    !
    do l = lminX(4), lmaxX(4)
        zell = real(l,8)
        Lfac = real(l*(l+1),8)
        zCIB = dsqrt(cib_bandpass143_nom143*A_cib_143*cib_bandpass217_nom217*A_cib_217)*(zell/3000)**(ncib) &
        /Lfac
        X_f(l - lminX(4) + npt(4) ) = &
        r_ps*dsqrt(A_ps_143*A_ps_217)*1.d-6/9.d0 + r_cib*zCIB &
        +A_ksz*ksz_temp(l)/Lfac  &
        -sqrt(cib_bandpass217_nom217*A_cib_217*sz_bandpass143_nom143*A_sz)*xi*tszxcib_temp(l)/Lfac
        X_data(l - lminX(4) + npt(4)) =  X(l - lminX(4) + npt(4))
        X_theory(l-lminX(4) + npt(4)) = cell_cmb(l)
        X_beam_corr_model(l-lminX(4)+npt(4)) = &
        ( X_theory(l - lminX(4) + npt(4))+ X_f(l - lminX(4) + npt(4)))* &
        corrected_beam(4,l)/dsqrt(cal1*cal2)
    end do


    do i = 1, nX
        Y(i) = X_data(i) - X_beam_corr_model(i)
    end do

    zlike = 0
    !$OMP parallel do private(j,i,ztemp) reduction(+:zlike) schedule(static,16)
    do  j = 1, nX
        ztemp= dot_product(Y(j+1:nX), c_inv(j+1:nX, j))
        zlike=zlike+ (ztemp*2 +c_inv(j, j)*Y(j))*Y(j)
    end do
    !    zlike = 0
    !    do  j = 1, nX
    !       ztemp= 0
    !       do  i = 1, nX
    !          ztemp = ztemp + Y(i)*c_inv(i, j)
    !       end do
    !       zlike=zlike+ztemp*Y(j)
    !    end do
    !   zlike = CAMSpec_Quad(c_inv, Y)

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

    zlike=zlike+((cal2/cal1-0.9966d0)/0.0015d0)**2  +((cal0/cal1-1.0006d0)/0.0004d0)**2

    ! print *,'CAMspec time:',  MPI_Wtime() - atime

    ! WARNING; weakening prior...
    !    zlike=zlike+((cal2/cal1-1.0056d0)/0.0063d0)**2 &
    !               +((cal0/cal1-1.0127d0)/0.003d0)**2


    ! correction to improve old way of handling cals...
    !zlike=zlike &
    !-2.d0*(lmaxX(1)-lminX(1)+1)*log(cal0) &
    !-2.d0*(lmaxX(2)-lminX(2)+1)*log(cal1) &
    !-2.d0*(lmaxX(3)-lminX(3)+1)*log(cal2) &
    !-2.d0*(lmaxX(4)-lminX(4)+1)*.5d0*(log(cal1)+log(cal2))


    !    if(writebest) then
    !    if(zlike.lt.bestlike) then
    !       open(49,file=bestname,form='formatted',status='unknown')
    !       !do if1=1,beam_Nspec
    !       !write(49,101) (beam_coeffs(if1,ie1),ie1=1,num_modes_per_beam)
    !       !enddo
    !101    format(5(e12.4))
    !       do l = 0,2500
    !          write(49,*) l, l*(l+1)*cell_cmb(l)
    !       enddo
    !       close(49)
    !       bestlike=zlike
    !    endif
    !    endif

    !if(storeall) then
    !countnum=countnum+1
    !write(storenumstring,*) countnum
    !storename=trim(storeroot)//'_'//trim(adjustl(bestnum))//'_' &
    !     //trim(adjustl(storenumstring))//'.txt'
    !   open(49,file=storename,form='formatted',status='unknown')
    !   do l=0,2500
    !   write(49,*) l,l*(l+1)*cell_cmb(l)
    !   enddo
    !   write(49,*)
    !   write(49,*) 'A_ps_100=',A_ps_100
    !   write(49,*) 'A_ps_143=',A_ps_143
    !   write(49,*) 'A_ps_217=',A_ps_217
    !   write(49,*) 'A_cib_143=',A_cib_143
    !   write(49,*) 'A_cib_217=',A_cib_217
    !   write(49,*) 'A_sz=',A_sz
    !   write(49,*) 'r_ps=',r_ps
    !   write(49,*) 'r_cib=',r_cib
    !   write(49,*) 'xi=',xi
    !   write(49,*) 'A_ksz=',A_ksz
    !   write(49,*) 'cal0=',cal0
    !   write(49,*) 'cal1=',cal1
    !   write(49,*) 'cal2=',cal2
    !   do i=1,4
    !   do j=1,num_modes_per_beam
    !      write(49,*) 'beam_coeffs(',i,',',j,')=',beam_coeffs(i,j)
    !      enddo
    !      enddo
    !      write(49,*)
    !      write(49,*) 'zlike(=chi^2)=',zlike
    !
    !
    !endif
    contains

    real(8) function corrected_beam(spec_num,l)
    integer, intent(in) :: spec_num,l
    integer :: i

    corrected_beam=1.d0
    do i=1,num_modes_per_beam
        corrected_beam=corrected_beam+beam_coeffs(spec_num,i)*beam_modes(i,l,spec_num)
    enddo
    end function corrected_beam

    end subroutine calc_like

end module temp_like
