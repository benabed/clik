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

SUBROUTINE CAMSPEC_EXTRA_INIT(mlmin1, mlmax1, mlmin2, mlmax2, mlmin3, mlmax3, mnp1, mnp2, mnp3, mnX,mX,mc_inv,mlmax_sz,msz_temp)
    USE CAMSPEC_EXTRA
    USE temp_like
    implicit none
    integer,intent(in) :: mlmin1, mlmax1, mlmin2, mlmax2, mlmin3, mlmax3,mnp1, mnp2, mnp3, mnX,mlmax_sz
    real(8),dimension(1:mnX)::mX
    real(8),dimension(0:mlmax_sz),intent(in)::msz_temp
    real(8),dimension(1:mnX,1:mnX),intent(in)::mc_inv
    integer::i
    
    call like_init_frommem(mlmin1, mlmax1, mlmin2, mlmax2, mlmin3, mlmax3, mnp1, mnp2, mnp3, mnX,mX,mc_inv,mlmax_sz,msz_temp)

    lmax = mlmax1
    IF (lmax<mlmax2) lmax = mlmax2
    IF (lmax<mlmax3) lmax = mlmax3
    lmin = mlmin1
    IF (lmin>mlmin2 .and. mlmax2.ne.0) lmin = mlmin2
    IF (lmin>mlmin3 .and. mlmax3.ne.0) lmin = mlmin3
    

    xcase = 0
    IF (mlmax2.ne.0) xcase = 1

    ALLOCATE(cltt(0:lmax+1))
    npar = lmax+1-lmin+5+4*xcase
            
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
    REAL(8)::A_ps_143, A_ps_217, A_cib_143, A_cib_217, A_sz, r_ps, r_cib, cal1, cal2
    INTEGER::l
    real(8)::tlkl

    cltt(:lmin-1) = 0
    DO l=lmin,lmax
        ! camspec expects cl/2pi !!! argl !
        cltt(l)=CL(l-lmin)/2./3.14159265358979323846264338328
    ENDDO

    IF (xcase==0) THEN
        A_ps_143  = CL(lmax+1-lmin + 0) 
        A_ps_217  = 0
        A_cib_143 = CL(lmax+1-lmin + 1)
        A_cib_217 = 0
        A_sz      = CL(lmax+1-lmin + 2)
        r_ps      = 0
        r_cib     = 0
        cal1      = CL(lmax+1-lmin + 3) 
        cal2      = CL(lmax+1-lmin + 4) 
    ELSE
        A_ps_143  = CL(lmax+1-lmin + 0) 
        A_ps_217  = CL(lmax+1-lmin + 1)
        A_cib_143 = CL(lmax+1-lmin + 2)
        A_cib_217 = CL(lmax+1-lmin + 3)
        A_sz      = CL(lmax+1-lmin + 4)
        r_ps      = CL(lmax+1-lmin + 5)
        r_cib     = CL(lmax+1-lmin + 6) 
        cal1      = CL(lmax+1-lmin + 7) 
        cal2      = CL(lmax+1-lmin + 8) 
    END IF
    call calc_like(tlkl,  cltt, A_ps_143, A_ps_217, A_cib_143, A_cib_217, A_sz, r_ps, r_cib, cal1, cal2)
    ! lkl is -2loglike clik returns loglik
    !print *,tlkl
    lkl = -tlkl/2.

END SUBROUTINE CAMSPEC_EXTRA_LKL

