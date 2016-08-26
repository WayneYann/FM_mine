        !COMPILER-GENERATED INTERFACE MODULE: Thu Nov 11 12:25:31 2010
        MODULE CKPRNT__genmod
          INTERFACE 
            SUBROUTINE CKPRNT(MDIM,MAXTP,MM,ENAME,KK,KNAME,WTM,KPHSE,   &
     &KCHRG,NT,T,TLO,TMID,THI,KNCF,ITHRM,LOUT,KERR)
              INTEGER(KIND=4) :: MAXTP
              INTEGER(KIND=4) :: MDIM
              INTEGER(KIND=4) :: MM
              CHARACTER(*) :: ENAME(*)
              INTEGER(KIND=4) :: KK
              CHARACTER(*) :: KNAME(*)
              REAL(KIND=8) :: WTM(*)
              INTEGER(KIND=4) :: KPHSE(*)
              INTEGER(KIND=4) :: KCHRG(*)
              INTEGER(KIND=4) :: NT(*)
              REAL(KIND=8) :: T(MAXTP,*)
              REAL(KIND=8) :: TLO
              REAL(KIND=8) :: TMID
              REAL(KIND=8) :: THI
              INTEGER(KIND=4) :: KNCF(MDIM,*)
              LOGICAL(KIND=4) :: ITHRM(*)
              INTEGER(KIND=4) :: LOUT
              LOGICAL(KIND=4) :: KERR
            END SUBROUTINE CKPRNT
          END INTERFACE 
        END MODULE CKPRNT__genmod
