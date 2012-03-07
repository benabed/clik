MODULE CAMSPEC_EXTRA

    IMPLICIT NONE

    INTEGER:: BOK = 0
    real(8), dimension(:), allocatable :: cltt
    INTEGER::lmin,lmax
    INTEGER::xcase,npar
END MODULE CAMSPEC_EXTRA



SUBROUTINE CAMSPEC_EXTRA_ONLY_ONE(MOK)
    USE CAMSPEC_EXTRA
    INTEGER,INTENT(OUT)::MOK
    MOK = BOK
    BOK = 1
END SUBROUTINE  CAMSPEC_EXTRA_ONLY_ONE

SUBROUTINE CAMSPEC_EXTRA_FREE()
    USE CAMSPEC_EXTRA
    DEALLOCATE(cltt)
    BOK =0
END SUBROUTINE  CAMSPEC_EXTRA_FREE

SUBROUTINE CAMSPEC_EXTRA_INIT(iNspec, inX,ilminX,ilmaxX,inp,inpt,ilmax_sz, sz_100,sz_143,mc_inv,mX)
    USE CAMSPEC_EXTRA
    USE temp_like
    implicit none
    integer, intent(in):: iNspec, iNx,ilmax_sz
    integer,dimension(1:iNspec),intent(in)::ilmaxX, ilminX,inp,inpt
    real(8),dimension(0:ilmax_sz),intent(in)::sz_100,sz_143
    real(8),dimension(1:iNx),intent(in)::mX
    real(8),dimension(1:iNx,1:iNx),intent(in)::mc_inv

    integer::i
    
    call like_init_frommem(iNspec, inX,ilminX,ilmaxX,inp,inpt,ilmax_sz, sz_100,sz_143,mc_inv,mX)
    
    lmax = ilmaxX(1)
    DO i=2,iNspec
        IF (lmax<ilmaxX(i)) lmax = ilmaxX(i)
    ENDDO
    lmin = ilminX(1)
    DO i=2,iNspec
        IF (lmin>ilminX(i)) lmin = ilminX(i)
    ENDDO
    
    xcase = 1
    !IF (mlmax2.ne.0) xcase = 1

    ALLOCATE(cltt(0:lmax+1))
    npar = lmax+1-lmin+11
    
END SUBROUTINE CAMSPEC_EXTRA_INIT

SUBROUTINE CAMSPEC_EXTRA_GETCASE(xc)
    USE CAMSPEC_EXTRA
    INTEGER::xc

    xc = xcase
END SUBROUTINE CAMSPEC_EXTRA_GETCASE

SUBROUTINE CAMSPEC_EXTRA_LKL(LKL,CL)
    USE CAMSPEC_EXTRA
    use temp_like
    implicit none

    REAL(8),INTENT(OUT)::LKL
    REAL(8),DIMENSION(0:npar-1)::CL
    real(8)::zlike, A_ps_100, A_ps_143, A_ps_217, A_cib_143, A_cib_217, A_sz, r_ps, r_cib, &
         cal0, cal1, cal2
    INTEGER::l
    real(8)::tlkl

    cltt(:lmin-1) = 0
    DO l=lmin,lmax
        ! camspec expects cl/2pi !!! argl !
        cltt(l)=CL(l-lmin)/2./3.14159265358979323846264338328
    ENDDO

    A_ps_100 = CL(lmax+1-lmin + 0)
    A_ps_143 = CL(lmax+1-lmin + 1)
    A_ps_217 = CL(lmax+1-lmin + 2)
    A_cib_143 = CL(lmax+1-lmin + 3)
    A_cib_217 = CL(lmax+1-lmin + 4)
    A_sz      = CL(lmax+1-lmin + 5)
    r_ps      = CL(lmax+1-lmin + 6)
    r_cib     = CL(lmax+1-lmin + 7) 
    cal0      = CL(lmax+1-lmin + 8) 
    cal1      = CL(lmax+1-lmin + 9) 
    cal2      = CL(lmax+1-lmin + 10)     

    call calc_like(tlkl,  cltt, A_ps_100, A_ps_143, A_ps_217, A_cib_143, A_cib_217, A_sz, r_ps, r_cib, cal0,cal1, cal2)
    ! lkl is -2loglike clik returns loglik
    !print *,tlkl
    lkl = -tlkl/2.

END SUBROUTINE CAMSPEC_EXTRA_LKL

