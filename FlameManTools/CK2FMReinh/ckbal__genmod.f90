        !COMPILER-GENERATED INTERFACE MODULE: Thu Nov 11 12:25:31 2010
        MODULE CKBAL__genmod
          INTERFACE 
            SUBROUTINE CKBAL(MXSPEC,KSPEC,KCOEF,MDIM,MM,KCHRG,KNCF,IERR)
              INTEGER(KIND=4) :: MDIM
              INTEGER(KIND=4) :: MXSPEC
              INTEGER(KIND=4) :: KSPEC(*)
              INTEGER(KIND=4) :: KCOEF(*)
              INTEGER(KIND=4) :: MM
              INTEGER(KIND=4) :: KCHRG(*)
              INTEGER(KIND=4) :: KNCF(MDIM,*)
              LOGICAL(KIND=4) :: IERR
            END SUBROUTINE CKBAL
          END INTERFACE 
        END MODULE CKBAL__genmod
