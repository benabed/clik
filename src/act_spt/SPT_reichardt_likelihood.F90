!===========================================================================
MODULE spt_reichardt_likelihood ! Parameters are defined in options module

  USE highell_options
  
  logical :: initialise_spt=.true.
  REAL(8), dimension(:,:) :: btt_dat(nspec_r,0:bmax0_r-1)
  REAL(8), dimension(:,:,:) :: btt_var(nspec_r,nspec_r,0:bmax0_r-1)
  REAL(8) ::  bval(nspec_r,0:bmax0_r-1), fisher(1:90,1:90)
  REAL(8), dimension (:), allocatable :: cl_tsz,cl_ksz,cl_szcib,cl_src,cl_c
  REAL(8) :: win_func(nspec_r,0:14,1:10000)
  
  PRIVATE
  public :: spt_reichardt_likelihood_init
  public :: spt_reichardt_likelihood_compute
  public :: get_free_lun

  contains
  
  !====================================================================================================================================
  SUBROUTINE spt_reichardt_likelihood_init
    
    IMPLICIT NONE
    
    INTEGER  :: i,j,lun,il
    REAL(8)  :: dummy
    CHARACTER(LEN=240) :: ttfilename(nspec_r),clfilename(nspec_r),winfilename(nspec_r),covfilename(90,90)
    LOGICAL  :: good

    allocate(cl_tsz(2:tt_lmax), cl_ksz(2:tt_lmax), cl_szcib(2:tt_lmax),cl_src(2:tt_lmax))
    allocate(cl_c(2:tt_lmax))

#ifdef TIMING
       call spt_timing_start( 'spt_reichardt_likelihood_init' )
#endif
    
    !----------------------------------------------------------------------------------
    ! set file names
    !----------------------------------------------------------------------------------
    
    ttfilename(1) = trim(SPT_data_dir)//'spt_high/spt_bandpower_95x95.dat'
    ttfilename(2) = trim(SPT_data_dir)//'spt_high/spt_bandpower_95x150.dat'
    ttfilename(3) = trim(SPT_data_dir)//'spt_high/spt_bandpower_95x220.dat'
    ttfilename(4) = trim(SPT_data_dir)//'spt_high/spt_bandpower_150x150.dat'
    ttfilename(5) = trim(SPT_data_dir)//'spt_high/spt_bandpower_150x220.dat'
    ttfilename(6) = trim(SPT_data_dir)//'spt_high/spt_bandpower_220x220.dat'
    covfilename   = trim(SPT_data_dir)//'spt_high/inverse_full.txt'
   
    winfilename(1) = trim(SPT_data_dir)//'spt_high/BblMean_95x95_new.dat'
    winfilename(2) = trim(SPT_data_dir)//'spt_high/BblMean_95x150_new.dat'
    winfilename(3) = trim(SPT_data_dir)//'spt_high/BblMean_95x220_new.dat'  
    winfilename(4) = trim(SPT_data_dir)//'spt_high/BblMean_150x150_new.dat'
    winfilename(5) = trim(SPT_data_dir)//'spt_high/BblMean_150x220_new.dat'
    winfilename(6) = trim(SPT_data_dir)//'spt_high/BblMean_220x220_new.dat'

    !-----------------------------------------------------------------------------------
    ! load TT data 
    !-----------------------------------------------------------------------------------
    do j=1,nspec_r

       inquire(file=ttfilename(j),exist = good)
       if(.not.good)then
          write(*,*) 'cant find', trim(ttfilename(j)), trim(SPT_data_dir)
          stop
       endif
       call get_free_lun( lun )
       open(unit=lun,file=ttfilename(j),form='formatted',status='unknown',action='read')
       do i=0,bmax0_r-1
          read(lun,*) bval(j,i),btt_dat(j,i), btt_var(j,j,i)
       enddo
       close(lun)
 
       inquire (file=winfilename(j),exist = good)
       if(.not.good)then
          write(*,*) 'cant find', trim(winfilename(j)), trim(SPT_data_dir)
          stop
       endif
       call get_free_lun( lun )
       open(unit=lun,file=winfilename(j),form='formatted',status='unknown',action='read')
       win_func(j,0:bmax0_r-1,1:10000)=0.d0 
       do il = 2, tt_lmax
          read(lun,*) i, (win_func(j,i,il), i=0,bmax0_r-1)       
       enddo
       close(lun) 
   enddo
   
   !-------------------------------------------------------------------------------------
   !Read spt inverse covariance matrix 
   !-------------------------------------------------------------------------------------

    call get_free_lun( lun )
    open(unit=lun,file=covfilename,form='formatted',status='unknown',action='read')
    do i=1,90
       read(lun,*) fisher(i,1:90)
    enddo
    close(lun)

   !------------------------------------------------------------------------------------
   ! Reading foreground components    
   !------------------------------------------------------------------------------------

    initialise_spt = .false.

    call get_free_lun(lun)
    open(unit=lun,file=trim(data_dir)//'Fg/tsz_143_eps0.50.dat',form='formatted',status='unknown')
    do il=2,tt_lmax
       read(lun,*) dummy,cl_tsz(il)
       cl_tsz(il) = cl_tsz(il)
    enddo
    close(lun)

    call get_free_lun(lun)
    open(unit=lun,file=trim(data_dir)//'Fg/cl_ksz_148_trac.dat',form='formatted',status='unknown')
    cl_ksz(2:tt_lmax) = 0.d0
    do il=2,tt_lmax-1
       read(lun,*) dummy,cl_ksz(il)
       cl_ksz(il) = cl_ksz(il)
    enddo
    close(lun)

    do il=2,tt_lmax
       cl_c(il)=(il/3000.d0)**(0.7)
    enddo

    call get_free_lun(lun)
    open(unit=lun,file=trim(data_dir)//'Fg/sz_x_cib_template.dat',form='formatted',status='unknown')
    cl_szcib(2:tt_lmax) = 0.d0
    do il=2,tt_lmax-1
       read(lun,*) dummy,cl_szcib(il)
    enddo
    close(lun)

#ifdef TIMING
      call spt_timing_end()
#endif

  END SUBROUTINE spt_reichardt_likelihood_init
 ! ======================================================================================================================================


 ! ==========================================================================================================================================
  SUBROUTINE spt_reichardt_likelihood_compute(cltt,amp_tsz,amp_ksz,xi,aps95,aps150,aps220,acib150,acib220,rps0,rps1,rps,rcib,cal_1,cal_2,cal_3,like_sptr)

    IMPLICIT NONE
    REAL(8), intent(in) :: cltt(2:*), amp_tsz,amp_ksz,xi,aps95,aps150,aps220,acib150,acib220,rps0,rps1,rps,rcib,cal_1,cal_2,cal_3
    REAL(8), intent(out) :: like_sptr
    INTEGER :: lun,il,i,j,k,jj
    REAL(8) :: cltt_temp(2:tt_lmax)
    REAL(8) :: btt_th(nspec_r,0:bmax0_r-1)
    REAL(8) :: diffs(90,1),tmp(90,1),diffs2(1,90),dlnlike, chi2(1,1)
    REAL(8) :: fcal_j
    REAL(8) :: beta_c,f1_sz,f1_synch,f1_dust,f2_sz,f2_synch,f2_dust,f3_sz,f3_synch,f3_dust
    REAL(8) :: fcibp2, fcibp3, planckratiod2,fluxtempd2,fluxtemps2, planckratiod3,fluxtempd3,fluxtemps3 

    f1_sz     =97.6d0
    f1_synch  =95.3d0
    f1_dust   =97.9d0
    f2_sz     =152.9d0
    f2_synch  =150.2d0
    f2_dust   =153.8d0
    f3_sz     =218.1d0
    f3_synch  =214.1d0
    f3_dust   =219.6d0


    fcibp2  = 146.d0
    fcibp3  = 226.d0

    planckratiod2=planckfunctionratio(f2_dust,fcibp2)
    fluxtempd2=flux2tempratio(f2_dust,fcibp2)
    fluxtemps2=flux2tempratio(f2_synch,fcibp2)

    planckratiod3=planckfunctionratio(f3_dust,fcibp3)
    fluxtempd3=flux2tempratio(f3_dust,fcibp3)
    fluxtemps3=flux2tempratio(f3_synch,fcibp3)

    beta_c = 2.20
    !----------------------------------------------------------------
    ! Calculate theory
    !----------------------------------------------------------------
    
    do j=1,nspec_r
       cltt_temp(2:tt_lmax)=0.d0
       do il=2,tt_lmax
          if(j==1) then
             cl_src(il) = aps95*(il*(il+1)/1000**2.0)
             cltt_temp(il) =cltt(il)+cl_src(il)+2.234*amp_tsz*cl_tsz(il)+amp_ksz*cl_ksz(il)
           else if (j==2) then
             cl_src(il) = rps0*sqrt(aps95*aps150)*(il*(il+1)/1000**2.0)
             cltt_temp(il) = cltt(il)+cl_src(il)+sqrt(2.234*0.839)*amp_tsz*cl_tsz(il)+amp_ksz*cl_ksz(il)-sqrt(acib150*amp_tsz*4.796*2.234)*xi*cl_szcib(il)*(f2_dust/fcibp2)**beta_c*(planckratiod2*fluxtempd2)
           else if (j==3) then
             cl_src(il) = rps1*sqrt(aps95*aps220)*(il*(il+1)/1000**2.0)
             cltt_temp(il) = cltt(il)+cl_src(il)+amp_ksz*cl_ksz(il)-sqrt(acib220*amp_tsz*4.796*2.234)*xi*cl_szcib(il)*(f3_dust/fcibp3)**beta_c*(planckratiod3*fluxtempd3)
           else if (j==4) then
             cl_src(il) = aps150*(il*(il+1)/1000**2.0)+acib150*cl_c(il)*(f2_dust/fcibp2)**(2.0*beta_c)*(planckratiod2*fluxtempd2)**2.0 &
                          -2.0*sqrt(acib150*amp_tsz*4.796*0.839)*xi*cl_szcib(il)*(f2_dust/fcibp2)**beta_c*(planckratiod2*fluxtempd2)
             cltt_temp(il) =cltt(il)+cl_src(il)+0.839*amp_tsz*cl_tsz(il)+amp_ksz*cl_ksz(il)
          else if (j==5) then
             cl_src(il) = rps*sqrt(aps150*aps220)*(il*(il+1)/1000**2.0)+rcib*sqrt(acib150*acib220)*cl_c(il)*(f2_dust/fcibp2)**beta_c*(planckratiod2*fluxtempd2)*(f3_dust/fcibp3)**beta_c*(planckratiod3*fluxtempd3)
             cltt_temp(il) =cltt(il)+cl_src(il)+amp_ksz*cl_ksz(il)-sqrt(acib220*amp_tsz*4.796*0.839)*xi*cl_szcib(il)*(f3_dust/fcibp3)**beta_c*(planckratiod3*fluxtempd3)
          else if(j ==6) then
             cl_src(il) = aps220*(il*(il+1)/1000**2)+acib220*cl_c(il)*(f3_dust/fcibp3)**(2.0*beta_c)*(planckratiod3*fluxtempd3)**2.0
             cltt_temp(il) =cltt(il)+cl_src(il)+amp_ksz*cl_ksz(il)
          endif
       enddo
      btt_th(j,0:bmax0_r-1)=MATMUL(win_func(j,0:bmax0_r-1,2:tt_lmax),cltt_temp(2:tt_lmax))
    enddo

    !--------------------------------------------------------------
    ! Calibrate theory
    !--------------------------------------------------------------

    do j=1,nspec_r
       if(j ==1 ) fcal_j = cal_1*cal_1
       if(j ==2 ) fcal_j = cal_1*cal_2
       if(j ==3 ) fcal_j = cal_1*cal_3
       if(j ==4 ) fcal_j = cal_2*cal_2
       if(j ==5 ) fcal_j = cal_2*cal_3
       if(j ==6 ) fcal_j = cal_3*cal_3
       btt_th(j,0:bmax0_r-1) = btt_th(j,0:bmax0_r-1)/fcal_j
    enddo

    !--------------------------------------------------------------
    ! chi2 calculation
    !--------------------------------------------------------------

    like_sptr = 0.d0

    do i =0,bmax0_r-1
          diffs(i+1,1)  = btt_dat(1,i) - btt_th(1,i)
          diffs(i+16,1) = btt_dat(2,i) - btt_th(2,i)
          diffs(i+31,1) = btt_dat(3,i) - btt_th(3,i)
          diffs(i+46,1) = btt_dat(4,i) - btt_th(4,i)
          diffs(i+61,1) = btt_dat(5,i) - btt_th(5,i)
          diffs(i+76,1) = btt_dat(6,i) - btt_th(6,i)
    enddo

    do i =1,90
       diffs2(1,i) = diffs(i,1)
    enddo


    tmp(:,:) = matmul(fisher(:,:),diffs(:,:))
    chi2(:,:) = matmul(diffs2(:,:),tmp(:,:)) 

    like_sptr = like_sptr+chi2(1,1)/2.0

   10  continue
    
#ifdef TIMING
       call spt_timing_end()
#endif

  end SUBROUTINE spt_reichardt_likelihood_compute
 ! ==========================================================================================================================================

  !===============================================================================
  function sz_func(freq)

  IMPLICIT NONE
  REAL(8) sz_func,freq
  REAL(8), parameter ::  t_cmb = 2.725, k_b = 1.3806503e-23, h_pl = 6.626068e-34
  REAL(8) :: nu, xx

  nu=freq*1.e9
  xx=h_pl*nu/(k_b*t_cmb)
  sz_func=(2.-(xx/2.)/tanh(xx/2.))

  end function

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

END MODULE spt_reichardt_likelihood

