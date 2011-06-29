!
!  BoPix, a pixel base likelihood code developped at IASF Bologna by A. De Rosa, F. Finelli, A. Gruppuso (2011)
!
   MODULE BOPIX

    PRIVATE
     
    PUBLIC BOPIX_LIKELIHOOD

    CONTAINS
    
      SUBROUTINE BOPIX_LIKELIHOOD(CL,MENOLOGLIK_S)

      USE LIBRARY
      USE PARAMETER_MODULE

      IMPLICIT NONE

      REAL(4),DIMENSION(:,:)            :: CL
      REAL(4)                           :: MENOLOGLIK_S  
      REAL(8),DIMENSION(0:5,0:200)      :: CL1
      INTEGER (4)                       :: I,J,K
      REAL (8)                          :: MENOLOGLIK 
      REAL (8)                          :: LOG_DOT_EIGEN 
      REAL(8)                           :: DDOT
      INTEGER(4)                        :: this_recl
      EXTERNAL DDOT

!===============================================================================================================================
!      WRITE(6,*)'CALCOLO PUNTO DELLA LIKELIHOOD - BOPIX'
!      TOT_WALL_TIME = ETIME(ELAPSED)      
!      WRITE(6,*)'TEMPO TRASCORSO ',TOT_WALL_TIME, 'USER ', ELAPSED(1), 'SYSTEM ', ELAPSED(2)

	WRITE (*,*) "TT in ",CL(2:5,1)
	WRITE (*,*) "EE in ",CL(2:5,3)
	WRITE (*,*) "BB in ",CL(2:5,4)
	WRITE (*,*) "TE in ",CL(2:5,2)
	WRITE (*,*) "TT fid",BOPIX_CL(0,2:5)
	WRITE (*,*) "EE fid",BOPIX_CL(2,2:5)
	WRITE (*,*) "BB fid",BOPIX_CL(3,2:5)
	WRITE (*,*) "TE fid",BOPIX_CL(1,2:5)
	
  SELECT CASE (BOPIX_CL_FILE)


   CASE (1)
!      WRITE (6,*) 'TT l=2,32,33', BOPIX_CL(0,2), BOPIX_CL(0,32), BOPIX_CL(0,33)
!      WRITE (6,*) 'EE l=2,32,33', BOPIX_CL(2,2), BOPIX_CL(2,32), BOPIX_CL(2,33)
!      WRITE (6,*) 'BB l=2,32,33', BOPIX_CL(3,2), BOPIX_CL(3,32), BOPIX_CL(3,33)
!      WRITE (6,*) 'TE l=2,32,33', BOPIX_CL(1,2), BOPIX_CL(1,32), BOPIX_CL(1,33)
      CL1(0:5,0)                   = 0.D0
      CL1(0:5,1)                   = 0.D0
      CL1(5,:)                     = 0.D0
      CL1(4,:)                     = 0.D0
      CL1(3,2:BOPIX_CL_LMAX)       = DBLE(CL(2:BOPIX_CL_LMAX,2))
      CL1(3,BOPIX_CL_LMAX+1:200)   = DBLE(BOPIX_CL(1,BOPIX_CL_LMAX+1:200))
!      CL1(3,2:23)                   = DBLE(CL(2:23,2))
!      CL1(3,24:200)   = DBLE(BOPIX_CL(2,24:200))
!B modes in the computation
!      CL1(2,2:BOPIX_CL_LMAX)       = DBLE(CL(2:BOPIX_CL_LMAX,4))
!      CL1(2,BOPIX_CL_LMAX+1:200)   = DBLE(BOPIX_CL(3,BOPIX_CL_LMAX+1:200))
!no B modes in the computation
      CL1(2,:)                   = 0.D0
      CL1(1,2:BOPIX_CL_LMAX)       = DBLE(CL(2:BOPIX_CL_LMAX,3))
      CL1(1,BOPIX_CL_LMAX+1:200)   = DBLE(BOPIX_CL(2,BOPIX_CL_LMAX+1:200))
!      CL1(1,2:23)                  = DBLE(CL(2:23,3))
!      CL1(1,24:200)                = DBLE(BOPIX_CL(1,24:200))
      CL1(0,2:BOPIX_CL_LMAX)       = DBLE(CL(2:BOPIX_CL_LMAX,1))
      CL1(0,BOPIX_CL_LMAX+1:200)   = DBLE(BOPIX_CL(0,BOPIX_CL_LMAX+1:200))
!      WRITE (6,*) 'TT l=2,32,33', CL1(0,2), CL1(0,32), CL1(0,33)
!      WRITE (6,*) 'EE l=2,32,33', CL1(1,2), CL1(1,32), CL1(1,33)
!      WRITE (6,*) 'BB l=2,32,33', CL1(2,2), CL1(2,32), CL1(2,33)
!      WRITE (6,*) 'TE l=2,32,33', CL1(3,2), CL1(3,32), CL1(3,33)


   CASE DEFAULT
      CL1(0:5,0)                   = 0.D0
      CL1(0:5,1)                   = 0.D0
      CL1(5,:)                     = 0.D0
      CL1(4,:)                     = 0.D0
      CL1(3,2:200)                 = DBLE(CL(2:200,2))
      CL1(2,:)                     = 0.D0
      CL1(1,2:200)                 = DBLE(CL(2:200,3))
      CL1(0,2:200)                 = DBLE(CL(2:200,1))
   END SELECT

	WRITE (*,*) "TT out ",CL1(0,2:5)
	WRITE (*,*) "EE out ",CL1(1,2:5)
	WRITE (*,*) "BB out ",CL1(2,2:5)
	WRITE (*,*) "TE out ",CL1(3,2:5)

!      WRITE (6,*) 'TT l=2,35', CL1(0,2), CL1(0,35) 

      DO K=0,KMAX
         CL1(K,BOPIX_LMIN:BOPIX_LMAX_COV)=CL1(K,BOPIX_LMIN:BOPIX_LMAX_COV)*BEAMWQ(BOPIX_LMIN:BOPIX_LMAX_COV,K)
      END DO

        CALL COV(CL1(0,0))

        COVSIG=COVSIG+COVNOISE

         CALL DPOTRF( 'U', N, COVSIG(1,1), N, INFO ) 
       IF (INFO .NE. 0) WRITE (6,*) 'ESITO DELLA DECOMPOSIZIONE DI CHOLESKY SU SIGNAL+NOISE(0=OK)', INFO

         LOG_DOT_EIGEN = 0.D0
!$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:LOG_DOT_EIGEN)
       DO I=1,N
         LOG_DOT_EIGEN=LOG_DOT_EIGEN + LOG(COVSIG(I,I))
       END DO
!$OMP END PARALLEL DO

       CONTR_RES = MAPOUT_ALL
       
          CALL DPOTRS( 'U', N, 1,COVSIG(1,1),N,CONTR_RES,N,INFO)
        IF (INFO .NE. 0) WRITE (6,*) 'CALCOLO DELLA SOLUZIONE AX=B (0=OK)', INFO

         MENOLOGLIK = DDOT( N,MAPOUT_ALL,1,CONTR_RES,1)/2.D0 
         MENOLOGLIK = MENOLOGLIK + LOG_DOT_EIGEN 
         MENOLOGLIK = MENO_LOGLIK_FACTOR * MENOLOGLIK
         MENOLOGLIK_S = SNGL(MENOLOGLIK) 

 END SUBROUTINE BOPIX_LIKELIHOOD
 
 END MODULE BOPIX

