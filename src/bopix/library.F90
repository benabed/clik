MODULE LIBRARY

      USE PARAMETER_MODULE

  IMPLICIT NONE

  PUBLIC COV

   
 CONTAINS

 SUBROUTINE COV(CL_COV)

      USE PARAMETER_MODULE
      
      IMPLICIT NONE

      REAL(8),DIMENSION(0:5,0:200)  :: CL_COV
      INTEGER(4)                    :: K, L, I, J, K_INDEX
      REAL(8)                       :: POL,POL1,POL2,CO1,CO2,CO3,CO4,CO5,CO6,DEN,COVA00,COVA11,COVA22,COVA12,COVA01,&
                                       COVA02, LD, COVA11_A, COVA22_A, COVA12_A, COVA00_1, COVA11_D, COVA22_D
      REAL(8)                       :: CONST, VAL, CONST_1B_1C, POL_1, POL_2, SQRT24M1
      INTEGER(4)                    :: VAL1, INDI, INDJ
      INTEGER (4),PARAMETER         :: LL = 2
      REAL(8),PARAMETER             :: BOPIX_PI=3.14159265358979323846264338328D0
      REAL(8),PARAMETER             :: PI4M1=1.D0/4.D0/BOPIX_PI, PIM1=1.D0/BOPIX_PI
      REAL(8),DIMENSION(BOPIX_LMIN:BOPIX_LMAX_COV) :: SQRT_DEN, COVA00_CONST, CONST1B_1, CONST1B_2, DEN_LD
!      EXTERNAL PLGNDR_D1_0
!--------------------------------------------------------------------------------------------

      SQRT24M1=1.D0/SQRT(24.D0)

      CALL OMP_SET_NUM_THREADS(BOPIX_NTHREADS)
!$OMP PARALLEL DEFAULT(SHARED)

!$OMP SECTIONS

!$OMP SECTION
         COVA00_1 = 0.D0
           DO L = BOPIX_LMIN,BOPIX_LMAX_COV
             COVA00_1 = COVA00_1 + (2.D0*DBLE(L)+1.D0)*PI4M1*PLGNDR_D1_0_(L,0,1.D0)*CL_COV(0,L)
           ENDDO
!$OMP SECTION           
           DO L = BOPIX_LMIN,BOPIX_LMAX_COV
             COVA00_CONST(L) = (2.D0*DBLE(L)+1.D0)*PI4M1*CL_COV(0,L)
           ENDDO
!$OMP SECTION
           DO L  = BOPIX_LMIN+1,BOPIX_LMAX_COV
              LD = DBLE(L)    
              SQRT_DEN(L)  = 2.D0/SQRT((LD-1.D0)*LD*(LD+1.D0)*(LD+2.D0))
           END DO
!$OMP SECTION
           DO L  = BOPIX_LMIN+1,BOPIX_LMAX_COV
              LD = DBLE(L)
              DEN_LD(L)    = (LD-1.D0)*LD*(LD+1.D0)*(LD+2.D0)
           END DO
!$OMP SECTION
           DO L  = BOPIX_LMIN+1,BOPIX_LMAX_COV
              LD = DBLE(L)
              CONST1B_1(L) = (2.D0*LD+1.D0)*PI4M1
           END DO
!$OMP SECTION
           DO L  = BOPIX_LMIN+1,BOPIX_LMAX_COV
              LD = DBLE(L)
              CONST1B_2(L) = (LD*(LD-1.D0)*0.5D0)
           END DO
!$OMP SECTION
           COVA11_A=0.D0
           COVA12_A=0.D0 
           COVA22_A=0.D0
            DO L=BOPIX_LMIN,BOPIX_LMAX_COV
             LD = DBLE(L)
             CO3 =(2.D0*LD+1.D0)*PI4M1*((-1.D0)**L)
             COVA11_A = COVA11_A + 0.5D0*CO3*(CL_COV(1,L)-CL_COV(2,L))
             COVA22_A = COVA22_A + 0.5D0*CO3*(CL_COV(2,L)-CL_COV(1,L))
             COVA12_A = COVA12_A + CO3*CL_COV(5,L)
            ENDDO
!$OMP SECTION
           COVA11_D = 0.D0
!           COVA22_D = 0.D0
            DO L=BOPIX_LMIN,BOPIX_LMAX_COV
             COVA11_D = COVA11_D + (2.D0*DBLE(L)+1.D0)*PI4M1*0.5D0*(CL_COV(1,L)+CL_COV(2,L))
!             COVA22_D = COVA22_D + (2.D0*DBLE(L)+1.D0)/4.D0/BOPIX_PI*0.5D0*(CL_COV(2,L)+CL_COV(1,L))
            ENDDO
!$OMP END SECTIONS

!======================================================

!$OMP DO PRIVATE(J,INDI,INDJ,VAL,VAL1,L,COVA00) SCHEDULE(STATIC)
     DO J=1,NMAT_T
!      DO I=1,LDMAT_T
      DO I=1,J
       INDI = NPIXVECT_T(I)
       INDJ = NPIXVECT_T(J)
       IF (INDI .NE. INDJ) THEN    
                     VAL = PIJ_T(I,J)
                     VAL1 = INT(VAL -1.D0 -1.D-10)
                       SELECT CASE (VAL1)
                         CASE (-2)
!1A
                          COVA00=0.D0
                          DO L=BOPIX_LMIN,BOPIX_LMAX_COV
                            COVA00 = COVA00 + COVA00_CONST(L)*PLGNDR_0_T(L,I,J)
                          ENDDO
                          COVSIG(I,J) = COVA00    
                       CASE DEFAULT
!1A
                          COVA00 = 2.5D0*PI4M1*(3.D0*VAL**2-1)*CL_COV(0,2)
                            DO L=BOPIX_LMIN+1,BOPIX_LMAX_COV
                               COVA00 = COVA00 + COVA00_CONST(L)*PLGNDR_0_T(L,I,J)
                            ENDDO
                          COVSIG(I,J) = COVA00
		   END SELECT
        ELSE
         COVSIG(I,J) = COVA00_1
      END IF
     END DO
    END DO
!$OMP END DO

!======================================================

!$OMP DO PRIVATE(I,INDI,INDJ,VAL,VAL1,DEN,POL1,POL2,CONST_1B_1C,COVA01,COVA02,POL_1,L,LD,POL_2,CONST)&
!$OMP SCHEDULE(STATIC)
    DO J = 1,NMAT_TP
     DO I = 1,LDMAT_TP
       INDI = NPIXVECT_T(I)
       INDJ = NPIXVECT_P(J)   
       
       IF (INDI .NE. INDJ) THEN
        VAL = PIJ_TP(I,J)
        VAL1 = INT(VAL -1.D0 -1.D-10)
                       SELECT CASE (VAL1)
                         CASE (-2)
!1B 1C
                          COVSIG(I,J+NMAT_T)         = 0.0D0
                          COVSIG(I,J+NMAT_T+NMAT_TP) = 0.0D0
                       CASE DEFAULT
!1B-1C
                         DEN = 1/PIJ2_TP(I,J)
                         POL1 = VAL
                         POL2 = PLGNDR_0_TP(2,I,J)
                         CONST_1B_1C = - 1.25D0*PIM1*(2.D0*VAL*POL1*DEN-(2.D0*DEN+1.D0)*POL2)*2.D0*SQRT24M1
                         COVA01 = CONST_1B_1C*CL_COV(3,2)
                         COVA02 = CONST_1B_1C*CL_COV(4,2)
                         POL_1 = 0.5D0*(3.D0*VAL**2-1)
                          DO L=BOPIX_LMIN+1,BOPIX_LMAX_COV
                            LD = DBLE(L)
                            POL_2 = PLGNDR_0_TP(L,I,J)
                            CONST = CONST1B_1(L)*(LD*VAL*POL_1*DEN-(LD*DEN+CONST1B_2(L))*POL_2)*SQRT_DEN(L)
                            COVA01 = COVA01 - CONST*CL_COV(3,L)
                            COVA02 = COVA02 - CONST*CL_COV(4,L)
                            POL_1 = POL_2
                          ENDDO
                         COVSIG(I,J+NMAT_T)         = COVA01*COSGAMMA_TP(I,J)+COVA02*SINGAMMA_TP(I,J)
                         COVSIG(I,J+NMAT_T+NMAT_TP) = COVA02*COSGAMMA_TP(I,J)-COVA01*SINGAMMA_TP(I,J)
		   END SELECT
	ELSE
         COVSIG(I,J+NMAT_T)         = 0.0D0
         COVSIG(I,J+NMAT_T+NMAT_TP) = 0.0D0
    END IF	 
 
    END DO 
  END DO
!$OMP END DO

!======================================================

!$OMP DO PRIVATE(I,INDI,INDJ,VAL,VAL1,DEN,POL,POL1,POL2,CO1,CO2,CO3,CO4,CO5,CO6,COVA11,COVA22,COVA12,LD,L)&
!$OMP SCHEDULE(STATIC)
    DO J=1,NMAT_P
     DO I=1,LDMAT_P
       INDI = NPIXVECT_P(I)
       INDJ = NPIXVECT_P(J)
       IF (INDI .NE. INDJ) THEN
          VAL  = PIJ_P(I,J)
          VAL1 = INT(VAL -1.D0 -1.D-10)
                       SELECT CASE (VAL1)
                         CASE (-2)
!2B 2C 3C
   COVSIG(I+LDMAT_T,J+NMAT_T)                =  COVA12_A*(COSG_SINA_P(I,J)+SING_COSA_P(I,J))&
                                                    +COVA11_A*COSG_COSA_P(I,J) + COVA22_A*SING_SINA_P(I,J)
   COVSIG(I+LDMAT_T,J+NMAT_T+NMAT_P)         =  COVA12_A*(COSG_COSA_P(I,J)-SING_COSA_P(I,J))&
                                                    +COVA22_A*COSG_SINA_P(I,J) - COVA11_A*SING_COSA_P(I,J)
   COVSIG(I+LDMAT_T+LDMAT_P,J+NMAT_T+NMAT_P) = -COVA12_A*(COSG_SINA_P(I,J)+SING_COSA_P(I,J))&
                                                    +COVA11_A*SING_SINA_P(I,J) + COVA22_A*COSG_COSA_P(I,J)

                       CASE DEFAULT
                         DEN = 1.D0/PIJ2_P(I,J)
!2B 2C 3C
                        POL = 3*PIJ2_P(I,J)
                        CO5 = (((2.D0*DEN-1.D0)*POL)/12.D0)
                        CO6 = -VAL*POL/6.D0*DEN
                        COVA11 = 1.25D0*PIM1*(CO5*CL_COV(1,2)-CO6*CL_COV(2,2))
                        COVA22 = 1.25D0*PIM1*(CO5*CL_COV(2,2)-CO6*CL_COV(1,2))
                        COVA12 = (CO5+CO6)*CL_COV(5,2)
!-----------------
                         LD = 3.D0
                         L  = 3
                         POL1 = POL
                         POL2 = VAL*5.D0*POL
                         CO3 = 1.D0/120.D0
                         CO4 = 7.D0*PI4M1
                         CO1 = ((5.D0*VAL*POL1*DEN)-(-DEN+3.D0)*POL2)*2.D0*CO3
                         CO2 = (5.D0*POL1-2.D0*VAL*POL2)*4.D0*CO3*DEN                      
                           COVA11 = COVA11 + CO4*(CO1*CL_COV(1,L) - CO2*CL_COV(2,L))
                           COVA22 = COVA22 + CO4*(CO1*CL_COV(2,L) - CO2*CL_COV(1,L))
                           COVA12 = COVA12 + CO4*(CO1 + CO2)*CL_COV(5,L)
!------------------
                         LD = 4.D0
                         L  = 4
                         POL1 = POL2
                         POL2 = PLGNDR_2_P(L,I,J)
                         CO3 = 1.D0/360.D0
                         CO4 = 9.D0*PI4M1
                         CO1 = (VAL*POL1*DEN - POL2)*12.D0*CO3
                         CO2 = (2.D0*POL1-VAL*POL2)*12.D0*CO3*DEN                      
                           COVA11 = COVA11 + CO4*(CO1*CL_COV(1,L) - CO2*CL_COV(2,L))
                           COVA22 = COVA22 + CO4*(CO1*CL_COV(2,L) - CO2*CL_COV(1,L))
                           COVA12 = COVA12 + CO4*(CO1 + CO2)*CL_COV(5,L)
!-------------------
                       DO L=5,BOPIX_LMAX_COV
                         LD = DBLE(L)
                         POL1 = POL2
                         POL2 = PLGNDR_2_P(L,I,J)
                         CO3 = 1.D0/DEN_LD(L)
                         CO4 = (2.D0*LD+1.D0)*PI4M1
                         CO1 = (((LD+2.D0)*VAL*POL1*DEN)-((LD-4.D0)*DEN+CONST1B_2(L))*POL2)*2.D0*CO3
                         CO2 = ((LD+2.D0)*POL1-(LD-1.D0)*VAL*POL2)*4.D0*CO3*DEN                      
                           COVA11 = COVA11 + CO4*(CO1*CL_COV(1,L) - CO2*CL_COV(2,L))
                           COVA22 = COVA22 + CO4*(CO1*CL_COV(2,L) - CO2*CL_COV(1,L))
                           COVA12 = COVA12 + CO4*(CO1 + CO2)*CL_COV(5,L)
                       ENDDO
                
         COVSIG(I+LDMAT_T,J+NMAT_T)                =  COVA12*(SING_COSA_P(I,J)+COSG_SINA_P(I,J))&
                                                     +COVA22*SING_SINA_P(I,J)+COVA11*COSG_COSA_P(I,J)
         COVSIG(I+LDMAT_T,J+NMAT_T+NMAT_P)         =  COVA12*(COSG_COSA_P(I,J)-SING_SINA_P(I,J))&
                                                     +COVA22*COSG_SINA_P(I,J)-COVA11*SING_COSA_P(I,J)
         COVSIG(I+LDMAT_T+LDMAT_P,J+NMAT_T+NMAT_P) = -COVA12*(COSG_SINA_P(I,J)+SING_COSA_P(I,J))&
                                                     +COVA22*COSG_COSA_P(I,J)+COVA11*SING_SINA_P(I,J)

          END SELECT

      ELSE
         COVSIG(I+LDMAT_T,J+NMAT_T)                = COVA11_D
         COVSIG(I+LDMAT_T,J+NMAT_T+NMAT_P)         = 0.0D0	  
         COVSIG(I+LDMAT_T+LDMAT_P,J+NMAT_T+NMAT_P) = COVA11_D
      END IF
    END DO
    END DO
!$OMP END DO

!$OMP END PARALLEL
!--------------------    
    END SUBROUTINE COV   

!================================================================================================================================
    REAL*8 FUNCTION PLGNDR_D1_0_(L,M,X)
        IMPLICIT NONE
        INTEGER(4), INTENT(IN) :: L,M
        REAL(8), INTENT(IN) :: X
!        REAL(8),INTENT(OUT) :: PLGNDR_D1_0
        INTEGER(4) :: LL
        REAL(8) :: PLL,PMM,PMMP1,SOMX2,XX, XX1, LLM1

        XX=X
        PMM=1.0D0
        PMMP1=X
        DO LL=2,L
         LLM1 = 1.0D0/DBLE(LL)
         PLL=(XX*(2.D0-LLM1)*PMMP1-(1.D0-LLM1)*PMM)
         PMM=PMMP1
         PMMP1=PLL
        END DO
     PLGNDR_D1_0_=PLL

   END FUNCTION PLGNDR_D1_0_
!==================================================================================================================================

  END MODULE
