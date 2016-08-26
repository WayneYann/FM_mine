        !COMPILER-GENERATED INTERFACE MODULE: Thu Nov 11 12:25:31 2010
        MODULE CKDUP__genmod
          INTERFACE 
            SUBROUTINE CKDUP(I,MAXSP,NS,NR,NU,NUNK,NFAL,IFAL,KFAL,ISAME)
              INTEGER(KIND=4) :: MAXSP
              INTEGER(KIND=4) :: I
              INTEGER(KIND=4) :: NS(*)
              INTEGER(KIND=4) :: NR(*)
              INTEGER(KIND=4) :: NU(MAXSP,*)
              INTEGER(KIND=4) :: NUNK(MAXSP,*)
              INTEGER(KIND=4) :: NFAL
              INTEGER(KIND=4) :: IFAL(*)
              INTEGER(KIND=4) :: KFAL(*)
              INTEGER(KIND=4) :: ISAME
            END SUBROUTINE CKDUP
          END INTERFACE 
        END MODULE CKDUP__genmod
