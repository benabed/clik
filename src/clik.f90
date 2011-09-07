module clik

  integer(kind=4), parameter :: PN_SIZE=256
  integer(kind=4), parameter :: MAX_NUMNAMES=1000
  
  type clik_object
    integer(kind=8)::ptr
  end type clik_object


contains
  
  subroutine clik_init(clikid,hdffilepath)
    
    ! Input
    type(clik_object), intent(out) :: clikid
    character(len=*), intent(in) :: hdffilepath
    ! Local
    integer(kind=4) :: fpathlen
    fpathlen = len_trim(hdffilepath)
    ! Call wrapping routine
    call fortran_clik_init(clikid%ptr,trim(hdffilepath),fpathlen)
    
  end subroutine clik_init
  
  subroutine clik_get_has_cl(clikid,has_cl)

    ! Input
    type(clik_object), intent(in) :: clikid
    integer(kind=4), dimension(6), intent(out) :: has_cl
    ! Call wrapping routine
    call fortran_clik_get_has_cl(clikid%ptr,has_cl)

  end subroutine clik_get_has_cl

  integer(kind=4) function clik_get_extra_parameter_names(clikid,names)

    ! Input
    type(clik_object), intent(in) :: clikid
    character(len=256), dimension(:), pointer :: names
    integer(kind=4) :: numnames

    ! Local
    character(len=PN_SIZE*MAX_NUMNAMES) :: buf_names

    ! Get number of extra parameters
    call fortran_clik_get_extra_parameter_number(clikid%ptr,numnames)
    if (numnames > 0) then
       ! Allocate character buffer
       allocate(names(numnames))
       ! Get the numnames, and copy into the names array
       call fortran_clik_get_extra_parameter_names(clikid%ptr,buf_names)
				do i=1,numnames
          names(i) = buf_names(1+(i-1)*PN_SIZE:(i)*PN_SIZE)
       enddo
    endif

    clik_get_extra_parameter_names=numnames

  end function clik_get_extra_parameter_names

  subroutine clik_get_lmax(clikid,lmax)

    ! Input
    type(clik_object), intent(in) :: clikid
    integer(kind=4), dimension(6), intent(out) :: lmax
    
    call fortran_get_lmax(clikid%ptr,lmax)

  end subroutine clik_get_lmax

  real(kind=8) function clik_compute(clikid,cl_and_pars)

    !Input
    type(clik_object), intent(in) :: clikid
    real(kind=8), dimension(:) :: cl_and_pars
    ! Local 
    real(kind=8) :: lkl
    call fortran_clik_compute(clikid%ptr,cl_and_pars,lkl)

    clik_compute = lkl
    return

  end function clik_compute

  subroutine clik_cleanup(clikid)

    !Input
    type(clik_object), intent(in) :: clikid
    call fortran_clik_cleanup(clikid%ptr)

  end subroutine clik_cleanup
  
  subroutine clik_fill_cl(clikid,cl_and_pars,clTT,clEE,clBB,clTE,clTB,clEB,flag_cl_start_at_two,flag_is_llp1_over_2pi,extrapars)

    type(clik_object), intent(in) :: clikid
    real(kind=8), dimension(:), pointer,intent(out) :: cl_and_pars
    real(kind=8), dimension(:), intent(in) :: clTT,clEE,clBB,clTE,clTB,clEB,extrapars
    logical,intent(in) ::  flag_cl_start_at_two,flag_is_llp1_over_2pi
    integer :: l, offl, cur,cli
    real(kind=8):: llp1_over_2pi,over_2pi
    integer(kind=4), dimension(6) :: lmax
    integer(kind=4) :: nextra, ntot, clmax

    over_2pi = 1./(6.283185307179586476925286766559005768394) 
    
    if (flag_cl_start_at_two) then
      offl = 2
    else
      offl = 0
    end if

    call clik_get_lmax(clikid,lmax)
    if (size(clTT)+offl<lmax(1)+1) then
      stop "not enough data in clTT"
    end if
    if (size(clEE)+offl<lmax(2)+1) then
      stop "not enough data in clTT"
    end if
    if (size(clBB)+offl<lmax(3)+1) then
      stop "not enough data in clBB"
    end if
    if (size(clTE)+offl<lmax(4)+1) then
      stop "not enough data in clTE"
    end if
    if (size(clTB)+offl<lmax(5)+1) then
      stop "not enough data in clTB"
    end if
    if (size(clEB)+offl<lmax(6)+1) then
      stop "not enough data in clEB"
    end if
    
    call fortran_clik_get_extra_parameter_number(clikid%ptr,nextra)
    
    if (size(extrapars)<nextra) then
      stop "not enough data in extrapars"
    end if
    
    ntot = 6 + nextra
    do cli=1,6
      ntot = ntot+lmax(cli)

    end do

    allocate(cl_and_pars(ntot))

    curl = 1
    
    !TT
    clmax = lmax(1)
    if (clmax .ne. -1) THEN
      do l = 0,offl
        cl_and_pars(curl) = 0
        curl = curl + 1
      end do
      do l = offl,clmax
        llp1_over_2pi = 1.
        if (flag_is_llp1_over_2pi) THEN
          llp1_over_2pi = l*(l+1.) * over_2pi
        end if
        cl_and_pars(curl) = clTT(l+1-offl) / llp1_over_2pi 
        curl = curl + 1
      end do
    end if
    
    !EE
    clmax = lmax(2)
    if (clmax .ne. -1) THEN
      do l = 0,offl
        cl_and_pars(curl) = 0
        curl = curl + 1
      end do
      do l = offl,clmax
        llp1_over_2pi = 1.
        if (flag_is_llp1_over_2pi) THEN
          llp1_over_2pi = l*(l+1.) * over_2pi
        end if
        cl_and_pars(curl) = clEE(l+1-offl) / llp1_over_2pi 
        curl = curl + 1
      end do
    end if

    !BB
    clmax = lmax(3)
    if (clmax .ne. -1) THEN
      do l = 0,offl
        cl_and_pars(curl) = 0
        curl = curl + 1
      end do
      do l = offl,clmax
        llp1_over_2pi = 1.
        if (flag_is_llp1_over_2pi) THEN
          llp1_over_2pi = l*(l+1.) * over_2pi
        end if
        cl_and_pars(curl) = clBB(l+1-offl) / llp1_over_2pi 
        curl = curl + 1
      end do
    end if

    !TE
    clmax = lmax(4)
    if (clmax .ne. -1) THEN
      do l = 0,offl
        cl_and_pars(curl) = 0
        curl = curl + 1
      end do
      do l = offl,clmax
        llp1_over_2pi = 1.
        if (flag_is_llp1_over_2pi) THEN
          llp1_over_2pi = l*(l+1.) * over_2pi
        end if
        cl_and_pars(curl) = clTE(l+1-offl) / llp1_over_2pi 
        curl = curl + 1
      end do
    end if

    !TB
    clmax = lmax(5)
    if (clmax .ne. -1) THEN
      do l = 0,offl
        cl_and_pars(curl) = 0
        curl = curl + 1
      end do
      do l = offl,clmax
        llp1_over_2pi = 1.
        if (flag_is_llp1_over_2pi) THEN
          llp1_over_2pi = l*(l+1.) * over_2pi
        end if
        cl_and_pars(curl) = clTB(l+1-offl) / llp1_over_2pi 
        curl = curl + 1
      end do
    end if

    !EB
    clmax = lmax(6)
    if (clmax .ne. -1) THEN
      do l = 0,offl
        cl_and_pars(curl) = 0
        curl = curl + 1
      end do
      do l = offl,clmax
        llp1_over_2pi = 1.
        if (flag_is_llp1_over_2pi) THEN
          llp1_over_2pi = l*(l+1.) * over_2pi
        end if
        cl_and_pars(curl) = clEB(l+1-offl) / llp1_over_2pi 
        curl = curl + 1
      end do
    end if

    ! copy extra parameters
    cl_and_pars(curl:ntot) = extrapars(:nextra)

  end subroutine  clik_fill_cl  


end module clik
