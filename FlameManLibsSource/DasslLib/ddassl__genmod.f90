        !COMPILER-GENERATED INTERFACE MODULE: Wed May 18 11:05:28 2011
        MODULE DDASSL__genmod
          INTERFACE 
            SUBROUTINE DDASSL(RES,NEQ,T,Y,YPRIME,TOUT,INFO,RTOL,ATOL,   &
     &IDID,RWORK,LRW,IWORK,LIW,RPAR,IPAR,JAC)
              EXTERNAL RES
              INTEGER(KIND=4) :: NEQ
              REAL(KIND=8) :: T
              REAL(KIND=8) :: Y(*)
              REAL(KIND=8) :: YPRIME(*)
              REAL(KIND=8) :: TOUT
              INTEGER(KIND=4) :: INFO(15)
              REAL(KIND=8) :: RTOL(*)
              REAL(KIND=8) :: ATOL(*)
              INTEGER(KIND=4) :: IDID
              REAL(KIND=8) :: RWORK(*)
              INTEGER(KIND=4) :: LRW
              INTEGER(KIND=4) :: IWORK(*)
              INTEGER(KIND=4) :: LIW
              REAL(KIND=8) :: RPAR(*)
              INTEGER(KIND=4) :: IPAR(*)
              EXTERNAL JAC
            END SUBROUTINE DDASSL
          END INTERFACE 
        END MODULE DDASSL__genmod
