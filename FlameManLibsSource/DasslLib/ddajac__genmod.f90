        !COMPILER-GENERATED INTERFACE MODULE: Wed May 18 11:05:28 2011
        MODULE DDAJAC__genmod
          INTERFACE 
            SUBROUTINE DDAJAC(NEQ,X,Y,YPRIME,DELTA,CJ,H,IER,WT,E,WM,IWM,&
     &RES,IRES,UROUND,JAC,RPAR,IPAR,NTEMP)
              INTEGER(KIND=4) :: NEQ
              REAL(KIND=8) :: X
              REAL(KIND=8) :: Y(*)
              REAL(KIND=8) :: YPRIME(*)
              REAL(KIND=8) :: DELTA(*)
              REAL(KIND=8) :: CJ
              REAL(KIND=8) :: H
              INTEGER(KIND=4) :: IER
              REAL(KIND=8) :: WT(*)
              REAL(KIND=8) :: E(*)
              REAL(KIND=8) :: WM(*)
              INTEGER(KIND=4) :: IWM(*)
              EXTERNAL RES
              INTEGER(KIND=4) :: IRES
              REAL(KIND=8) :: UROUND
              EXTERNAL JAC
              REAL(KIND=8) :: RPAR(*)
              INTEGER(KIND=4) :: IPAR(*)
              INTEGER(KIND=4) :: NTEMP
            END SUBROUTINE DDAJAC
          END INTERFACE 
        END MODULE DDAJAC__genmod
