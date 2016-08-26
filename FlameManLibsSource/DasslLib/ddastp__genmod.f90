        !COMPILER-GENERATED INTERFACE MODULE: Wed May 18 11:05:28 2011
        MODULE DDASTP__genmod
          INTERFACE 
            SUBROUTINE DDASTP(X,Y,YPRIME,NEQ,RES,JAC,H,WT,JSTART,IDID,  &
     &RPAR,IPAR,PHI,DELTA,E,WM,IWM,ALPHA,BETA,GAMMA,PSI,SIGMA,CJ,CJOLD, &
     &HOLD,S,HMIN,UROUND,IPHASE,JCALC,K,KOLD,NS,NONNEG,NTEMP)
              INTEGER(KIND=4) :: NEQ
              REAL(KIND=8) :: X
              REAL(KIND=8) :: Y(*)
              REAL(KIND=8) :: YPRIME(*)
              EXTERNAL RES
              EXTERNAL JAC
              REAL(KIND=8) :: H
              REAL(KIND=8) :: WT(*)
              INTEGER(KIND=4) :: JSTART
              INTEGER(KIND=4) :: IDID
              REAL(KIND=8) :: RPAR(*)
              INTEGER(KIND=4) :: IPAR(*)
              REAL(KIND=8) :: PHI(NEQ,*)
              REAL(KIND=8) :: DELTA(*)
              REAL(KIND=8) :: E(*)
              REAL(KIND=8) :: WM(*)
              INTEGER(KIND=4) :: IWM(*)
              REAL(KIND=8) :: ALPHA(*)
              REAL(KIND=8) :: BETA(*)
              REAL(KIND=8) :: GAMMA(*)
              REAL(KIND=8) :: PSI(*)
              REAL(KIND=8) :: SIGMA(*)
              REAL(KIND=8) :: CJ
              REAL(KIND=8) :: CJOLD
              REAL(KIND=8) :: HOLD
              REAL(KIND=8) :: S
              REAL(KIND=8) :: HMIN
              REAL(KIND=8) :: UROUND
              INTEGER(KIND=4) :: IPHASE
              INTEGER(KIND=4) :: JCALC
              INTEGER(KIND=4) :: K
              INTEGER(KIND=4) :: KOLD
              INTEGER(KIND=4) :: NS
              INTEGER(KIND=4) :: NONNEG
              INTEGER(KIND=4) :: NTEMP
            END SUBROUTINE DDASTP
          END INTERFACE 
        END MODULE DDASTP__genmod
