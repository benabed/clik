!============================================================================
MODULE act_south_likelihood 
! Parameters are defined in ACT_options module
! ===========================================================================

  USE highell_options
  
  logical :: initialise_act=.true.
  REAL(8), dimension(:,:) :: btt_dat1(nsp11_s,0:nbin11-1),btt_dat12(nsp12_s,0:nbin12-1),btt_dat2(nsp22_s,0:nbin22-1)
  REAL(8), dimension(:,:) :: bval1(nsp11_s,0:nbin11-1),bval12(nsp12_s,0:nbin12-1),bval2(nsp22_s,0:nbin22-1)
  REAL(8) ::  inverse(1:datap_s,1:datap_s)
  REAL(8), dimension (:), allocatable :: cl_tsz,cl_ksz,cl_szcib,cl_src,cl_c,cl_p
  REAL(8) :: win_func(nspec,0:tbin-1,1:10000)
  
  PRIVATE
  public :: act_south_likelihood_init
  public :: act_south_likelihood_compute
  public :: get_free_lun

contains
  
  ! ============================================================================
  SUBROUTINE act_south_likelihood_init
  ! ============================================================================
    
    IMPLICIT NONE
    
    INTEGER  :: i,j,lun,il
    REAL(8)  :: dummy
    CHARACTER(LEN=240) :: ttfilename(nspec_s), winfilename(nspec), invcovfilename
    LOGICAL  :: good

    
    allocate(cl_tsz(2:tt_lmax),cl_ksz(2:tt_lmax),cl_szcib(2:tt_lmax),cl_src(2:tt_lmax))    
    allocate(cl_c(2:tt_lmax),cl_p(2:tt_lmax))

#ifdef TIMING
    call act_south_timing_start('act_likelihood_init')
#endif
    
    !------------------------------------------------
    ! set file names
    !------------------------------------------------
    
    ttfilename(1) = trim(ACT_data_dir)//'south/spectrum_148x148_season2sxseason2s.dat'
    ttfilename(2) = trim(ACT_data_dir)//'south/spectrum_148x148_season2sxseason3s.dat'
    ttfilename(3) = trim(ACT_data_dir)//'south/spectrum_148x148_season2sxseason4s.dat'
    ttfilename(4) = trim(ACT_data_dir)//'south/spectrum_148x148_season3sxseason3s.dat'
    ttfilename(5) = trim(ACT_data_dir)//'south/spectrum_148x148_season3sxseason4s.dat'
    ttfilename(6) = trim(ACT_data_dir)//'south/spectrum_148x148_season4sxseason4s.dat'
    ttfilename(7) = trim(ACT_data_dir)//'south/spectrum_148x220_season2sxseason2s.dat'
    ttfilename(8) = trim(ACT_data_dir)//'south/spectrum_148x220_season2sxseason3s.dat'
    ttfilename(9) = trim(ACT_data_dir)//'south/spectrum_148x220_season2sxseason4s.dat'
    ttfilename(10) = trim(ACT_data_dir)//'south/spectrum_148x220_season3sxseason2s.dat'
    ttfilename(11) = trim(ACT_data_dir)//'south/spectrum_148x220_season3sxseason3s.dat'
    ttfilename(12) = trim(ACT_data_dir)//'south/spectrum_148x220_season3sxseason4s.dat'
    ttfilename(13) = trim(ACT_data_dir)//'south/spectrum_148x220_season4sxseason2s.dat'
    ttfilename(14) = trim(ACT_data_dir)//'south/spectrum_148x220_season4sxseason3s.dat'
    ttfilename(15) = trim(ACT_data_dir)//'south/spectrum_148x220_season4sxseason4s.dat'
    ttfilename(16) = trim(ACT_data_dir)//'south/spectrum_220x220_season2sxseason2s.dat'
    ttfilename(17) = trim(ACT_data_dir)//'south/spectrum_220x220_season2sxseason3s.dat'
    ttfilename(18) = trim(ACT_data_dir)//'south/spectrum_220x220_season2sxseason4s.dat'
    ttfilename(19) = trim(ACT_data_dir)//'south/spectrum_220x220_season3sxseason3s.dat'
    ttfilename(20) = trim(ACT_data_dir)//'south/spectrum_220x220_season3sxseason4s.dat'
    ttfilename(21) = trim(ACT_data_dir)//'south/spectrum_220x220_season4sxseason4s.dat'

    invcovfilename= trim(ACT_data_dir)//'south/Inverse_South.dat'
   
    winfilename(1) = trim(ACT_data_dir)//'south/BblMean_148x148.dat'
    winfilename(2) = trim(ACT_data_dir)//'south/BblMean_148x220.dat'
    winfilename(3) = trim(ACT_data_dir)//'south/BblMean_220x220.dat'
    
    !----------------------------------------------
    ! load TT data 
    !----------------------------------------------
    do j=1,nspec_s

       inquire(file=ttfilename(j),exist = good)
       if(.not.good)then
          write(*,*) 'cant find', trim(ttfilename(j)), trim(ACT_data_dir)
          stop
       endif
       call get_free_lun( lun )
       open(unit=lun,file=ttfilename(j),form='formatted',status='unknown',action='read')
       if (j .le. nsp11_s) then
          do i=0,nbin11-1
             read(lun,*) bval1(j,i),btt_dat1(j,i)
          enddo
       else if (j .ge. (nsp11_s+1) .and. j .le. (nsp11_s+nsp12_s)) then  
          do i=0,nbin12-1
             read(lun,*) bval12(j-nsp11_s,i),btt_dat12(j-nsp11_s,i)
          enddo
       else 
          do i=0,nbin22-1
             read(lun,*) bval2(j-(nsp11_s+nsp12_s),i),btt_dat2(j-(nsp11_s+nsp12_s),i)
          enddo
       end if
       close(lun)
     enddo

    !----------------------------------------------
    ! read windows for theory 
    !----------------------------------------------
    do j=1,nspec
       inquire (file=winfilename(j),exist = good)
       if(.not.good)then
          write(*,*) 'cant find', trim(winfilename(j)), trim(ACT_data_dir)
          stop
       endif
       call get_free_lun( lun )
       open(unit=lun,file=winfilename(j),form='formatted',status='unknown',action='read')
       win_func(j,0:tbin-1,1:10000)=0.d0
       do il = 2, tt_lmax
          read(lun,*) i, (win_func(j,i,il), i=0,tbin-1)       
       enddo
       close(lun) 
    enddo
 

   !------------------------------------------------- 
   !Read inverse covariance matrix 
   !-------------------------------------------------

    call get_free_lun( lun )
    open(unit=lun,file=invcovfilename,form='formatted',status='unknown',action='read')
    do i=1,datap_s
       read(lun,*) inverse(i,1:datap_s)
    enddo
    close(lun)

   !-------------------------------------------------
   ! Reading foreground components    
   !-------------------------------------------------

    initialise_act = .false.

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
    call act_south_timing_end()
#endif
 

  END SUBROUTINE act_south_likelihood_init
  
 ! ===================================================================================================================================
  SUBROUTINE act_south_likelihood_compute(cltt,amp_tsz,amp_ksz,xi,aps148,aps217,acib150,acib220,rps,rcib,cas1,cas2,like_acts)
 ! ===================================================================================================================================

    IMPLICIT NONE
    REAL(8), intent(in) :: cltt(2:*), amp_tsz,amp_ksz,xi,aps148,aps217,acib150,acib220,rps,rcib,cas1,cas2
    REAL(8), intent(out) :: like_acts
    INTEGER :: lun,il,i,j,k
    REAL(8) :: cltt_temp(2:tt_lmax)
    REAL(8) :: btt_th(nspec,0:tbin-1)
    REAL(8) :: diffs(datap_s,1),tmp(datap_s,1),diffs2(1,datap_s),chi2(1,1)
    REAL(8) :: fcal_j,beta_c
    REAL(8) :: f1_sz,f1_synch,f1_dust,f2_sz,f2_synch,f2_dust,fcibp2, fcibp3
    REAL(8) :: planckratiod1,planckratiod2,fluxtempd1,fluxtempd2,fluxtemps1,fluxtemps2
 
    f1_sz     =146.9d0
    f1_synch  =147.6d0 
    f1_dust   =149.7d0
    f2_sz     =220.2d0
    f2_synch  =217.6d0
    f2_dust   =219.6d0   

    fcibp2  = 146.d0
    fcibp3  = 226.d0

 
    planckratiod1=planckfunctionratio(f1_dust,fcibp2)
    fluxtempd1=flux2tempratio(f1_dust,fcibp2)
    fluxtemps1=flux2tempratio(f1_synch,fcibp2)

    planckratiod2=planckfunctionratio(f2_dust,fcibp3)
    fluxtempd2=flux2tempratio(f2_dust,fcibp3)
    fluxtemps2=flux2tempratio(f2_synch,fcibp3)

    beta_c = 2.20

    !----------------------------------------------------------------
    ! Calculate theory
    !----------------------------------------------------------------
    
    do j=1,nspec
       cltt_temp(2:tt_lmax)=0.d0

       do il=2,tt_lmax
          if(j==1) then
             cl_src(il) = aps148*(il*(il+1)/1000**2.0)+acib150*cl_c(il)*(f1_dust/fcibp2)**(2.0*beta_c)*(planckratiod1*fluxtempd1)**2.0 &
                          -2.0*sqrt(acib150*amp_tsz*4.796*0.839)*xi*cl_szcib(il)*(f1_dust/fcibp2)**beta_c*(planckratiod1*fluxtempd1)
             cltt_temp(il) =cltt(il)+cl_src(il)+0.839*amp_tsz*cl_tsz(il)+amp_ksz*cl_ksz(il)
          else if (j==2) then
             cl_src(il) = rps*sqrt(aps148*aps217)*(il*(il+1)/1000**2.0)+rcib*sqrt(acib150*acib220)*cl_c(il)*(f1_dust/fcibp2)**beta_c*(planckratiod1*fluxtempd1)*(f2_dust/fcibp3)**beta_c*(planckratiod2*fluxtempd2)
             cltt_temp(il) =cltt(il)+cl_src(il)+amp_ksz*cl_ksz(il)-sqrt(acib220*amp_tsz*4.796*0.839)*xi*cl_szcib(il)*(f2_dust/fcibp3)**beta_c*(planckratiod2*fluxtempd2)
          else if(j ==3) then
             cl_src(il) = aps217*(il*(il+1)/1000**2)+acib220*cl_c(il)*(f2_dust/fcibp3)**(2.0*beta_c)*(planckratiod2*fluxtempd2)**2.0
             cltt_temp(il) =cltt(il)+cl_src(il)+amp_ksz*cl_ksz(il)
          endif
          cltt_temp(il) = cltt_temp(il)/((dble(il)*(dble(il)+1.0))/(2*PI))
      enddo
      btt_th(j,0:tbin-1)=MATMUL(win_func(j,0:tbin-1,2:tt_lmax),cltt_temp(2:tt_lmax))
   enddo

    !--------------------------------------------------------------
    ! Calibrate theory
    !--------------------------------------------------------------

    do j=1,nspec
       if(j ==1 ) fcal_j = cas1*cas1
       if(j ==2 ) fcal_j = cas1*cas2
       if(j ==3 ) fcal_j = cas2*cas2
       btt_th(j,0:tbin-1) = btt_th(j,0:tbin-1)/fcal_j
    enddo

    !--------------------------------------------------------------
    ! chi2 calculation
    !--------------------------------------------------------------

    like_acts = 0.d0

    diffs(datap_s,1) = 0.d0
    diffs2(1,datap_s) = 0.d0

    do i =0,nbin11-1
       do j =0,nsp11_s-1
          diffs(i+1+j*nbin11,1) = btt_dat1(j+1,i) - btt_th(1,i+4)
       enddo
    enddo 
    do i =0,nbin12-1
       do j=0,nsp12_s
          diffs(i+1+nsp11_s*nbin11+j*nbin12,1) = btt_dat12(j+1,i) - btt_th(2,i+14)
       enddo
    enddo
    do i=0,nbin22-1
       do j=0,nsp22_s-1
          diffs(i+1+nsp11_s*nbin11+nsp12_s*nbin12+j*nbin22,1) = btt_dat2(j+1,i) - btt_th(3,i+14)
       enddo 
   enddo

    do i =1,datap_s
       diffs2(1,i) = diffs(i,1)
    enddo

    tmp(:,:) = matmul(inverse(:,:),diffs(:,:))
    chi2(:,:) = matmul(diffs2(:,:),tmp(:,:))

    like_acts = like_acts+chi2(1,1)/2.0

10  continue
    
#ifdef TIMING
    call act_south_timing_end()
#endif

  end SUBROUTINE act_south_likelihood_compute
    
  !================================================================================
  function sz_func(freq)

  IMPLICIT NONE
  REAL(8) sz_func,freq
  REAL(8), parameter ::  t_cmb = 2.725, k_b = 1.3806503e-23, h_pl = 6.626068e-34
  REAL(8) :: nu, xx

  nu=freq*1.e9
  xx=h_pl*nu/(k_b*t_cmb)
  sz_func=(2.-(xx/2.)/tanh(xx/2.))

  end function

  !================================================================================

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
  !=================================================================================

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
  !=================================================================================

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

END MODULE act_south_likelihood
