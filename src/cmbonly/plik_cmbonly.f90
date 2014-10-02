module Plik_CMBonly 

  implicit none
  private
  integer, parameter :: campc = KIND(1.d0)

  character(LEN=500),public :: data_dir 
  character(LEN=*), parameter, public :: plik_like='Plik_TTTEEE_v13_cmbonly_like'

  !Possible combinations: TT only, TE only, EE only, TT+TE+EE
  logical :: use_tt  = .true.
  logical :: use_ee  = .true.
  logical :: use_te  = .true.

  integer, public :: tt_lmax  = 3000 
  integer, parameter,public :: plmax  = 2508
  integer, parameter,public :: plmin  = 50
  integer, parameter,public :: nbin   = 601 ! total bins   
  integer, public, parameter :: nspectt = 4
  integer, public, parameter :: nspecte = 6
  integer, public, parameter :: nspecee = 6
  integer :: nptt(nspectt), npte(nspecte), npee(nspecee) ! number of bins
  integer, public :: nbintt, nbinte, nbinee

  !-------------------------------------------------------
  real(campc), parameter :: PI    = 3.14159265358979323846264d0
  real(campc), dimension(:), allocatable ::  bval,X_data,X_sig,diff_vec
  real(campc), dimension(:,:), allocatable :: covmat, fisher
  real(campc), dimension(:), allocatable :: blmin,blmax,bin_w
  real(campc), dimension(:), allocatable :: bl,bm
 
  public like_init_cmbonly, calc_like_cmbonly
  contains
  
  ! ===========================================================================
  subroutine like_init_cmbonly
 
  integer  :: i,j,lun,il,info,dum,k
  character(LEN=1024) :: like_file, cov_file, blmin_file, blmax_file, binw_file
  logical  :: good
    
  print *, 'Initializing Planck likelihood, version '//plik_like
 
  nptt(1:nspectt)=211 !same numbers of bins for 100, 143, 100x143, 217
  npte(1:nspecte)=195 !same numbers of bins for 100, 100x143, 100x217, 143,100x143, 217
  npee(1:nspecee)=195 !same numbers of bins for 100, 100x143, 100x217, 143,100x143, 217
  nbintt = maxval(nptt) !in future they might cut some spectra
  nbinte = maxval(npte)
  nbinee = maxval(npee)
    
  like_file = trim(data_dir)//'cl_cmb_plik_v9.0.dat'
  cov_file  = trim(data_dir)//'c_matrix_plik_v9.0.dat'
  blmin_file = trim(data_dir)//'blmin.dat'
  blmax_file = trim(data_dir)//'blmax.dat'
  binw_file = trim(data_dir)//'bweight.dat'!!Bbl_150_spt_v2.dat'
 
  allocate(bval(nbin),X_data(nbin),X_sig(nbin))
  allocate(covmat(nbin,nbin))
   
  inquire(file=like_file, exist=good)
  if(.not.good)then
     write(*,*) 'file not found', trim(like_file), trim(data_dir)
     stop
  endif
   
  call get_free_lun(lun)
  open(unit=lun,file=like_file,form='formatted',status='unknown',action='read')
  do i=1,nbin !read Planck
     read(lun,*) bval(i),X_data(i),X_sig(i)
  enddo
  close(lun)

  inquire(file=cov_file, exist=good)
  if(.not.good)then
     write(*,*) 'file not found', trim(cov_file), trim(data_dir)
     stop
  endif
  call get_free_lun(lun)
  open(unit=lun,file=cov_file,form='unformatted',status='old')
      read(lun) covmat
  close(lun)
  do i=1,nbin
     do j=i+1,nbin
        covmat(j,i)=covmat(i,j) !covmat is now full matrix
     enddo
  enddo

  allocate(blmin(nbintt),blmax(nbintt),bin_w(0:plmax),bl(0:plmax),bm(nbintt))

  call get_free_lun(lun)
  open(unit=lun,file=blmin_file,form='formatted',status='old')
  do i=1,nbintt
       read(lun,*) blmin(i) 
  end do
  close(lun)

  call get_free_lun(lun)
  open(unit=lun,file=blmax_file,form='formatted',status='old')
  do i=1,nbintt
       read(lun,*) blmax(i)
  end do
  close(lun)

  call get_free_lun(lun)
  open(unit=lun,file=binw_file,form='formatted',status='old')
  do i=0,plmax
       read(lun,*) bin_w(i)
  end do
  close(lun)

  do i=0,plmax
     bl(i) = i*(i+1)/2./PI
  end do

  do i=1,nbintt !binned ell
     bm(i) = sum(bl(blmin(i)+plmin:blmax(i)+plmin)*bin_w(blmin(i):blmax(i)))
  end do

  end subroutine like_init_cmbonly
  
  ! ===========================================================================
  subroutine calc_like_cmbonly(plike,cell_tt,cell_te,cell_ee)

  real(campc), dimension(2:) :: cell_tt,cell_ee,cell_te
  real(campc) :: cl_tt(nbintt),cl_te(nbinte),cl_ee(nbinee)
  real(campc) plike
  integer :: bin_no,lun,il,i,j,info
  real(campc), allocatable, save ::  Y(:), X_model(:)
  real(campc), allocatable :: ptemp(:)

  if (.not. allocated(Y)) then
     allocate(X_model(nbin),Y(nbin))
     X_model = 0
     Y = 0
  end if

  do i=1,nbintt
     cl_tt(i) = sum(cell_tt(blmin(i)+plmin:blmax(i)+plmin)*bin_w(blmin(i):blmax(i)))
  end do
  do i=1,nbinte
     cl_te(i) = sum(cell_te(blmin(i)+plmin:blmax(i)+plmin)*bin_w(blmin(i):blmax(i)))
  end do
  do i=1,nbinee
     cl_ee(i) = sum(cell_ee(blmin(i)+plmin:blmax(i)+plmin)*bin_w(blmin(i):blmax(i)))
  end do

  cl_tt(1:nbintt)=cl_tt(1:nbintt)/bm(1:nbintt)
  cl_te(1:nbinte)=cl_te(1:nbinte)/bm(1:nbinte)
  cl_ee(1:nbinee)=cl_ee(1:nbinee)/bm(1:nbinee)

  X_model(1:nbintt) = cl_tt(1:nbintt) !TT
  X_model(nbintt+1:nbintt+nbinte) = cl_te(1:nbinte) !TE
  X_model(nbintt+nbinte+1:nbintt+nbinte+nbinee) = cl_ee(1:nbinee) !EE
 
!Start basic chisq 
  Y = X_data - X_model

  !Only TT
  if((use_tt .eqv. .true.) .and. (use_te .eqv. .false.) .and. (use_ee .eqv. .false.)) then
       bin_no=nbintt
       allocate(fisher(bin_no,bin_no))
       allocate(diff_vec(bin_no),ptemp(bin_no))
       diff_vec(:)=Y(1:bin_no)
       fisher(:,:)=covmat(1:bin_no,1:bin_no)
  !Only TE
  else if((use_tt .eqv. .false.) .and. (use_te .eqv. .true.) .and. (use_ee .eqv. .false.)) then 
       bin_no=nbinte
       allocate(fisher(bin_no,bin_no))
       allocate(diff_vec(bin_no),ptemp(bin_no))
       diff_vec(:)=Y(nbintt+1:nbintt+bin_no)
       fisher(:,:)=covmat(nbintt+1:nbintt+bin_no,nbintt+1:nbintt+bin_no)
  !Only EE
  else if((use_tt .eqv. .false.) .and. (use_te .eqv. .false.) .and. (use_ee .eqv. .true.)) then
       bin_no=nbinee
       allocate(fisher(bin_no,bin_no))
       allocate(diff_vec(bin_no),ptemp(bin_no))
       diff_vec(:)=Y(nbintt+nbinte+1:nbintt+nbinte+bin_no)
       fisher(:,:)=covmat(nbintt+nbinte+1:nbintt+nbinte+bin_no,nbintt+nbinte+1:nbintt+nbinte+bin_no)
  else if ((use_tt .eqv. .true.) .and. (use_te .eqv. .true.) .and. (use_ee .eqv. .true.)) then
       bin_no=nbin
       allocate(fisher(nbin,nbin))
       allocate(diff_vec(nbin),ptemp(nbin))
       diff_vec(:)=Y(:)
       fisher(:,:)=covmat(:,:)
    else
       write(*,*) 'Fail: no possible options chosen'
    endif

    !Invert covmat
    call dpotrf('U',bin_no,fisher,bin_no,info)
    if(info.ne.0)then
       print*, ' info in dpotrf =', info
       stop
    endif
    call dpotri('U',bin_no,fisher,bin_no,info)
    if(info.ne.0)then
       print*, ' info in dpotri =', info
       stop
    endif
    do i=1,bin_no
       do j=i,bin_no
          fisher(j,i)=fisher(i,j)
       enddo
    enddo

 
    ptemp=matmul(fisher,diff_vec)
    plike=sum(ptemp*diff_vec)
 
    deallocate(X_model,Y)
    deallocate(fisher,diff_vec,ptemp)
    
 end subroutine calc_like_cmbonly
  
 subroutine get_free_lun(lun)

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
 end subroutine get_free_lun

end module Plik_CMBonly
