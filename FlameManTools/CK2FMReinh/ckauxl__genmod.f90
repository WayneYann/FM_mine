        !COMPILER-GENERATED INTERFACE MODULE: Thu Nov 11 12:25:31 2010
        MODULE CKAUXL__genmod
          INTERFACE 
            SUBROUTINE CKAUXL(SUB,NSUB,II,KK,KNAME,LOUT,MAXSP,NPAR,NSPEC&
     &,NTHB,ITHB,NTBS,MAXTB,NKTB,AIK,NFAL,IFAL,IDUP,NFAR,PFAL,IFOP,NLAN,&
     &ILAN,NLAR,PLAN,NREV,IREV,RPAR,NRLT,IRLT,RLAN,NWL,IWL,WL,KERR,NORD,&
     &IORD,MAXORD,KORD,RORD,NUNK,NU,NRNU,IRNU,RNU)
              INTEGER(KIND=4) :: MAXORD
              INTEGER(KIND=4) :: NLAR
              INTEGER(KIND=4) :: NFAR
              INTEGER(KIND=4) :: MAXTB
              INTEGER(KIND=4) :: NPAR
              INTEGER(KIND=4) :: MAXSP
              CHARACTER(*) :: SUB(*)
              INTEGER(KIND=4) :: NSUB
              INTEGER(KIND=4) :: II
              INTEGER(KIND=4) :: KK
              CHARACTER(*) :: KNAME(*)
              INTEGER(KIND=4) :: LOUT
              INTEGER(KIND=4) :: NSPEC(*)
              INTEGER(KIND=4) :: NTHB
              INTEGER(KIND=4) :: ITHB(*)
              INTEGER(KIND=4) :: NTBS(*)
              INTEGER(KIND=4) :: NKTB(MAXTB,*)
              REAL(KIND=8) :: AIK(MAXTB,*)
              INTEGER(KIND=4) :: NFAL
              INTEGER(KIND=4) :: IFAL(*)
              INTEGER(KIND=4) :: IDUP(*)
              REAL(KIND=8) :: PFAL(NFAR,*)
              INTEGER(KIND=4) :: IFOP(*)
              INTEGER(KIND=4) :: NLAN
              INTEGER(KIND=4) :: ILAN(*)
              REAL(KIND=8) :: PLAN(NLAR,*)
              INTEGER(KIND=4) :: NREV
              INTEGER(KIND=4) :: IREV(*)
              REAL(KIND=8) :: RPAR(NPAR,*)
              INTEGER(KIND=4) :: NRLT
              INTEGER(KIND=4) :: IRLT(*)
              REAL(KIND=8) :: RLAN(NLAR,*)
              INTEGER(KIND=4) :: NWL
              INTEGER(KIND=4) :: IWL(*)
              REAL(KIND=8) :: WL(*)
              LOGICAL(KIND=4) :: KERR
              INTEGER(KIND=4) :: NORD
              INTEGER(KIND=4) :: IORD(*)
              INTEGER(KIND=4) :: KORD(MAXORD,*)
              REAL(KIND=8) :: RORD(MAXORD,*)
              INTEGER(KIND=4) :: NUNK(MAXSP,*)
              INTEGER(KIND=4) :: NU(MAXSP,*)
              INTEGER(KIND=4) :: NRNU
              INTEGER(KIND=4) :: IRNU(*)
              REAL(KIND=8) :: RNU(MAXSP,*)
            END SUBROUTINE CKAUXL
          END INTERFACE 
        END MODULE CKAUXL__genmod
