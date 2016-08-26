        !COMPILER-GENERATED INTERFACE MODULE: Thu Nov 11 12:25:31 2010
        MODULE CKTHRM__genmod
          INTERFACE 
            SUBROUTINE CKTHRM(LUNIT,MDIM,ENAME,MM,AWT,KNAME,KK,KNCF,    &
     &KPHSE,KCHRG,WTM,MAXTP,NT,NTR,TLO,TMID,THI,T,NPCP2,A,ITHRM,KERR,   &
     &LOUT,ISTR)
              INTEGER(KIND=4) :: NPCP2
              INTEGER(KIND=4) :: NTR
              INTEGER(KIND=4) :: MAXTP
              INTEGER(KIND=4) :: MDIM
              INTEGER(KIND=4) :: LUNIT
              CHARACTER(*) :: ENAME(*)
              INTEGER(KIND=4) :: MM
              REAL(KIND=8) :: AWT(*)
              CHARACTER(*) :: KNAME(*)
              INTEGER(KIND=4) :: KK
              INTEGER(KIND=4) :: KNCF(MDIM,*)
              INTEGER(KIND=4) :: KPHSE(*)
              INTEGER(KIND=4) :: KCHRG(*)
              REAL(KIND=8) :: WTM(*)
              INTEGER(KIND=4) :: NT(*)
              REAL(KIND=8) :: TLO
              REAL(KIND=8) :: TMID
              REAL(KIND=8) :: THI
              REAL(KIND=8) :: T(MAXTP,*)
              REAL(KIND=8) :: A(NPCP2,NTR,*)
              LOGICAL(KIND=4) :: ITHRM(*)
              LOGICAL(KIND=4) :: KERR
              INTEGER(KIND=4) :: LOUT
              CHARACTER(LEN=80) :: ISTR
            END SUBROUTINE CKTHRM
          END INTERFACE 
        END MODULE CKTHRM__genmod
