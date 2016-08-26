        !COMPILER-GENERATED INTERFACE MODULE: Wed May 18 11:05:28 2011
        MODULE DDAINI__genmod
          INTERFACE 
            SUBROUTINE DDAINI(X,Y,YPRIME,NEQ,RES,JAC,H,WT,IDID,RPAR,IPAR&
     &,PHI,DELTA,E,WM,IWM,HMIN,UROUND,NONNEG,NTEMP)
              INTEGER(KIND=4) :: NEQ
              REAL(KIND=8) :: X
              REAL(KIND=8) :: Y(*)
              REAL(KIND=8) :: YPRIME(*)
              EXTERNAL RES
              EXTERNAL JAC
              REAL(KIND=8) :: H
              REAL(KIND=8) :: WT(*)
              INTEGER(KIND=4) :: IDID
              REAL(KIND=8) :: RPAR(*)
              INTEGER(KIND=4) :: IPAR(*)
              REAL(KIND=8) :: PHI(NEQ,*)
              REAL(KIND=8) :: DELTA(*)
              REAL(KIND=8) :: E(*)
              REAL(KIND=8) :: WM(*)
              INTEGER(KIND=4) :: IWM(*)
              REAL(KIND=8) :: HMIN
              REAL(KIND=8) :: UROUND
              INTEGER(KIND=4) :: NONNEG
              INTEGER(KIND=4) :: NTEMP
            END SUBROUTINE DDAINI
          END INTERFACE 
        END MODULE DDAINI__genmod
