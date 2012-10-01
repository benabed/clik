! ===========================================================================
MODULE spt_keisler_likelihood
! Parameters are defined in options module

  USE highell_options
  
  REAL(8) ::  btt_dat(0:bmax0_k-1),btt_var(0:bmax0_k-1)
  REAL(8) ::  fisher(1:bmax0_k,1:bmax0_k),bval(0:bmax0_k-1)
  REAL(8), dimension (:), allocatable :: cl_tsz,cl_ksz,cl_szcib,cl_c,cl_src
  REAL(8), dimension(:,:), allocatable :: win_func
  
  PRIVATE
  public :: spt_keisler_likelihood_init
  public :: spt_keisler_likelihood_compute
  public :: get_free_lun
  contains
  
  !=================================================================================================================
  SUBROUTINE spt_keisler_likelihood_init
    
    IMPLICIT NONE
    
    INTEGER  :: i,j,lun,il
    REAL(8)  :: dummy,ii
    CHARACTER(LEN=240) :: ttfilename, bblfilename, covfilename
    LOGICAL  :: good
    
    allocate(cl_tsz(2:tt_lmax_k),cl_ksz(2:tt_lmax_k),cl_szcib(2:tt_lmax_k))
    allocate(cl_c(2:tt_lmax_k),cl_src(2:tt_lmax_k))

    !-----------------------------------------------------------------------------------------
    ! set file names
    !-----------------------------------------------------------------------------------------
    
    ttfilename  = trim(SPT_data_dir)//'spt_lowell/spt_kspectrum_150x150.dat'
    covfilename = trim(SPT_data_dir)//'spt_lowell/inverse_short.dat'
    bblfilename = trim(SPT_data_dir)//'spt_lowell/BblMean_150x150.dat'

    !-----------------------------------------------------------------------------------------
    ! load spectrum, covariance and band power window functions 
    !-----------------------------------------------------------------------------------------
    
    inquire(file=ttfilename,exist = good)
    if(.not.good)then
       write(*,*) 'cant find', trim(ttfilename), trim(SPT_data_dir)
       stop
    endif
    call get_free_lun( lun )
    open(unit=lun,file=ttfilename,form='formatted',status='unknown',action='read')    
    do i=0,bmax0_k-1
       read(lun,*) bval(i),btt_dat(i), btt_var(i)
    enddo
    close(lun)

    call get_free_lun( lun )
    open(unit=lun,file=covfilename,form='formatted',status='unknown',action='read')
    do i=1,bmax0_k
       read(lun,*) fisher(i,1:bmax0_k)
    enddo
    close(lun)

    inquire (file=bblfilename,exist = good)
    if(.not.good)then
       write(*,*) 'cant find', trim(bblfilename), trim(SPT_data_dir)
       stop
    endif

    call get_free_lun( lun )
    open(unit=lun,file=bblfilename,form='formatted',status='unknown',action='read')
    allocate(win_func(0:bmax0_k-1,1:tt_lmax_k))
    do il = 2, tt_lmax_k
        read(lun,*) ii, (win_func(i,il), i=0,bmax0_k-1) 
        enddo
    close(lun)

   !-------------------------------------------------------------------------------------------
   ! Reading foreground components    
   !-------------------------------------------------------------------------------------------

    call get_free_lun( lun )
    open(unit=lun,file=trim(data_dir)//'Fg/tsz_143_eps0.50.dat',form='formatted',status='unknown')
    do il=2,tt_lmax_k
       read(lun,*) dummy,cl_tsz(il)
       cl_tsz(il) = cl_tsz(il)
    enddo
    close(lun) 

    call get_free_lun(lun)
    open(unit=lun,file=trim(data_dir)//'Fg/cl_ksz_148_trac.dat',form='formatted',status='unknown')
    do il=2,tt_lmax_k
       read(lun,*) dummy,cl_ksz(il)
       cl_ksz(il) = cl_ksz(il)
    enddo
    close(lun)
    
    do il=2,tt_lmax_k
       cl_c(il)=(il/3000.d0)**0.7
    enddo

    call get_free_lun(lun)
    open(unit=lun,file=trim(data_dir)//'Fg/sz_x_cib_template.dat',form='formatted',status='unknown')
    do il=2,tt_lmax_k
       read(lun,*) dummy,cl_szcib(il)
    enddo
    close(lun)

    
  END SUBROUTINE spt_keisler_likelihood_init
  !=================================================================================================================

  ! ===========================================================================================================================
  SUBROUTINE spt_keisler_likelihood_compute(cltt,amp_tsz,amp_ksz,xi,aps150,acib150,cal_2,like_sptk)
    
    IMPLICIT NONE
    REAL(8), intent(in) :: cltt(2:*), amp_tsz,amp_ksz,xi,aps150,acib150,cal_2
    REAL(8), intent(out) :: like_sptk
    INTEGER :: lun,il,i   
    REAL(8) :: cltt_temp(2:tt_lmax_k)
    REAL(8) :: btt_th(0:bmax0_k-1)
    REAL(8) :: diffs(bmax0_k,1),tmp(bmax0_k,1),diffs2(1,bmax0_k),chi2(1,1)
    REAL(8) :: fcibp2, f2_sz,f2_synch,f2_dust, beta_c
    REAL(8) :: planckratiod2,fluxtempd2,fluxtemps2

    fcibp2    =146.d0

    f2_sz     =152.9d0
    f2_synch  =150.2d0
    f2_dust   =153.8d0

    planckratiod2=planckfunctionratio(f2_dust,fcibp2)
    fluxtempd2=flux2tempratio(f2_dust,fcibp2)
    fluxtemps2=flux2tempratio(f2_synch,fcibp2)

    beta_c = 2.20


    !----------------------------------------------------------------
    ! Calculate theory
    !----------------------------------------------------------------

    do il=2,tt_lmax_k
       cl_src(il) = (aps150+9.2/9.0)*(il*(il+1)/1000**2)+acib150*cl_c(il)*(f2_dust/fcibp2)**(2.0*beta_c)*(planckratiod2*fluxtempd2)**2.0 &
                    -2.0*sqrt(acib150*amp_tsz*4.796*0.839)*xi*cl_szcib(il)*(f2_dust/fcibp2)**beta_c*(planckratiod2*fluxtempd2)
       cltt_temp(il) =cltt(il)+cl_src(il)+0.839*amp_tsz*cl_tsz(il)+amp_ksz*cl_ksz(il)
       cltt_temp(il) =cltt_temp(il)/(cal_2**2.0)
    enddo
    btt_th(0:bmax0_k-1)=MATMUL(win_func(0:bmax0_k-1,2:tt_lmax_k),cltt_temp(2:tt_lmax_k))
    

    !--------------------------------------------------------------
    ! chi2 calculation
    !--------------------------------------------------------------

    like_sptk = 0.d0

    do i = 1, bmax0_k
       diffs(i,1) = btt_dat(i-1) - btt_th(i-1)
       diffs2(1,i) = diffs(i,1)
    enddo 
 
    tmp(:,:) = matmul(fisher(:,:),diffs(:,:))
    chi2(:,:) = matmul(diffs2(:,:),tmp(:,:))

    like_sptk = like_sptk+chi2(1,1)/2.0

   10  continue
    
  end SUBROUTINE spt_keisler_likelihood_compute  
  ! ===========================================================================================================================

  !===============================================================================
  function planckfunctionratio(freq,feff)
  ! All the constants cancel out when rescaling the Planck function 
  ! to the effective 150 GHz frequency
  ! B(nu)/B(nu0)
  IMPLICIT NONE

  REAL(8) planckfunctionratio,freq,feff
  REAL(8), parameter ::  t_eff = 9.7, k_b = 1.3806503e-23, h_pl = 6.626068e-34
  REAL(8) :: nu,nu0,xx,xx0

  nu=freq*1.e9
  nu0=feff*1.e9
  xx=h_pl*nu/(k_b*t_eff)
  xx0=h_pl*nu0/(k_b*t_eff)
  planckfunctionratio = (nu/nu0)**3.0*(exp(xx0)-1)/(exp(xx)-1)

  END function
  !===============================================================================

  function flux2tempratio(freq,feff)
  ! rescaled to 150 GHz
  IMPLICIT NONE

  REAL(8) flux2tempratio,freq,feff
  REAL(8), parameter ::  t_cmb = 2.725, k_b = 1.3806503e-23, h_pl = 6.626068e-34
  REAL(8) :: nu,nu0,xx,xx0

  nu=freq*1.e9
  nu0=feff*1.e9
  xx=h_pl*nu/(k_b*t_cmb)
  xx0=h_pl*nu0/(k_b*t_cmb)
  flux2tempratio = (nu0/nu)**4.0*(exp(xx0)/exp(xx))*((exp(xx)-1)/(exp(xx0)-1))**2.0

  END function
  !================================================================================

  !===============================================================================
  subroutine get_free_lun( lun )

  implicit none
  integer, intent(out) :: lun
  integer, save :: last_lun = 19
  logical :: used
  lun = last_lun

  do
    inquire( unit=lun, opened=used )
    if ( .not. used ) exit
        lun = lun + 1
  end do

  last_lun = lun
  end subroutine
  !================================================================================

END MODULE spt_keisler_likelihood

