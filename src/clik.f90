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
          names(i) = buf_names(1+i*PN_SIZE:(i+1)*PN_SIZE)
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
  

end module clik
