MODULE CAMSPEC_EXTRA

    IMPLICIT NONE

    INTEGER:: BOK = 0
    real(8), dimension(:), allocatable :: cltt
    INTEGER::lmax
    INTEGER::xcase
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

SUBROUTINE CAMSPEC_EXTRA_INIT(LIKEFILE,SZFILE)
    USE CAMSPEC_EXTRA
    USE temp_like
    CHARACTER(LEN=2048),INTENT(IN)::LIKEFILE,SZFILE
    CHARACTER(LEN=2048)::DLIKEFILE,DSZFILE

    DLIKEFILE = LIKEFILE
    DLIKEFILE(len_trim(LIKEFILE):len_trim(LIKEFILE)) = ' '   
    DLIKEFILE = TRIM(ADJUSTL(DLIKEFILE))
    DSZFILE = SZFILE
    DSZFILE(len_trim(SZFILE):len_trim(SZFILE)) = ' '     
    DSZFILE = TRIM(ADJUSTL(DSZFILE))

    call like_init(DLIKEFILE, DSZFILE)

    lmax = lmax1
    IF (lmax<lmax2) lmax = lmax2
    IF (lmax<lmax3) lmax = lmax3
    lmin = lmin1
    IF (lmin>lmin2 .and. lmax2.ne.0) lmin = lmin2
    IF (lmin>lmin3 .and. lmax3.ne.0) lmin = lmin3
    

    xcase = 0
    IF (lmax2.ne.0) xcase = 1

    ALLOCATE(cltt(0:lmax+1))

            
END SUBROUTINE CAMSPEC_EXTRA_INIT

SUBROUTINE CAMSPEC_EXTRA_GETCASE(xc)
    USE CAMSPEC_EXTRA
    INTEGER::xc

    xc = xcase
END SUBROUTINE CAMSPEC_EXTRA_GETCASE

SUBROUTINE CAMSPEC_EXTRA_LKL(LKL,CL)
    USE CAMSPEC_EXTRA
    use temp_like
    REAL(8),INTENT(OUT)::LKL
    REAL(8),DIMENSION(0:)::CL
    REAL(8)::A_ps_143, A_ps_217, A_cib_143, A_cib_217, A_sz, r_ps, r_cib, cal1, cal2
    INTEGER::l

    cltt(:lmin-1) = 0
    DO l=lmin,lmax
        cltt(l)=CL(l-lmin)
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

    call calc_like(lkl,  cltt, A_ps_143, A_ps_217, A_cib_143, A_cib_217, A_sz, r_ps, r_cib, cal1, cal2)
    ! lkl is -2loglike clik returns loglik
    lkl = -lkl/2.

END SUBROUTINE CAMSPEC_EXTRA_LKL

