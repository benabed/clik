program clik_example_f90 

  use clik
  implicit none
  type(clik_object) :: pself
  character(len=128) :: hdffilename, clfilename, clnames(6)
  integer(kind=4), dimension(6) :: has_cl,lmax
  character(len=256), dimension(:), pointer :: names
  integer(kind=4) :: numnames, counter, nl, i, j, l, nargc
  real(kind=8), dimension(:), allocatable :: cl_and_pars, cls
  real(kind=8) :: lkl
  
  nargc = iargc()
  if (nargc < 2) then
     print*,'usage: clik_example_f90 hdffile clfile1 [clfile2 ...]'
     stop
  end if
  call getarg(1,hdffilename)

  call clik_init(pself,hdffilename)

  call clik_get_has_cl(pself,has_cl)

  call clik_get_lmax(pself,lmax)
  
  ! Log what was read
  clnames(1)='TT'
  clnames(2)='EE'
  clnames(3)='BB'
  clnames(4)='TE'
  clnames(5)='TB'
  clnames(6)='EB'
  print*,'Likelihood use Cl'
  do i=1,6
     if (has_cl(i)==1) then
        print*,'  ',trim(clnames(i)),' from l=0 to l=',lmax(i),' (incl.)'
     endif
  enddo

  ! See if we have extra parameters
  numnames=clik_get_extra_parameter_names(pself,names)
  print*,'Number of extra parameters: ',numnames
  do i=1,numnames
	  print *,'  ',trim(names(i))
  enddo

  ! Total number of multipoles to read

  nl = numnames ! Place for parameters values
  do i=1,6
     nl = nl + lmax(i)+1
  enddo
  print*,'parameter vector has ',nl,' elements'

  ! Fill cls
  do j=2,nargc
    call getarg(j,clfilename)  
    open(unit=100,file=clfilename,form='formatted')
    allocate(cl_and_pars(nl))
    counter=1
    do i=1,6
	    if (has_cl(i)==1) then
	      do l=0,lmax(i)
	        read(100,*),cl_and_pars(counter)
          counter = counter + 1
        enddo
      endif
    enddo

		do i=1,numnames
			read(100,*),cl_and_pars(counter)
      counter = counter + 1
		enddo
			
    lkl = clik_compute(pself,cl_and_pars)
    print*,'Log likelihood for this file ',trim(clfilename),' :',lkl
    close(unit=100)
    deallocate(cl_and_pars)
  enddo

  ! Free stuff
  if (numnames > 0) then
     deallocate(names)
  endif

end program clik_example_f90
