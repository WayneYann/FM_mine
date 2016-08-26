        !COMPILER-GENERATED INTERFACE MODULE: Thu Nov 11 12:25:31 2010
        MODULE CKRBAL__genmod
          INTERFACE 
            SUBROUTINE CKRBAL(MXSPEC,KSPEC,RCOEF,MDIM,MM,KCHRG,KNCF,    &
     &CKMIN,IERR)
              INTEGER(KIND=4) :: MDIM
              INTEGER(KIND=4) :: MXSPEC
              INTEGER(KIND=4) :: KSPEC(*)
              REAL(KIND=8) :: RCOEF(*)
              INTEGER(KIND=4) :: MM
              INTEGER(KIND=4) :: KCHRG(*)
              INTEGER(KIND=4) :: KNCF(MDIM,*)
              REAL(KIND=8) :: CKMIN
              LOGICAL(KIND=4) :: IERR
            END SUBROUTINE CKRBAL
          END INTERFACE 
        END MODULE CKRBAL__genmod
