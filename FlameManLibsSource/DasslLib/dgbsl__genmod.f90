        !COMPILER-GENERATED INTERFACE MODULE: Wed May 18 11:05:29 2011
        MODULE DGBSL__genmod
          INTERFACE 
            SUBROUTINE DGBSL(ABD,LDA,N,ML,MU,IPVT,B,JOB)
              INTEGER(KIND=4) :: LDA
              REAL(KIND=8) :: ABD(LDA,1)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: ML
              INTEGER(KIND=4) :: MU
              INTEGER(KIND=4) :: IPVT(1)
              REAL(KIND=8) :: B(1)
              INTEGER(KIND=4) :: JOB
            END SUBROUTINE DGBSL
          END INTERFACE 
        END MODULE DGBSL__genmod
