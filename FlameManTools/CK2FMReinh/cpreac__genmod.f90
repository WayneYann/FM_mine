        !COMPILER-GENERATED INTERFACE MODULE: Thu Nov 11 12:25:31 2010
        MODULE CPREAC__genmod
          INTERFACE 
            SUBROUTINE CPREAC(II,MAXSP,NSPEC,NPAR,PAR,RPAR,AUNITS,EUNITS&
     &,NREAC,NUNK,NU,KCHRG,MDIM,MM,KNCF,IDUP,NFAL,IFAL,KFAL,NFAR,PFAL,  &
     &IFOP,NREV,IREV,NTHB,ITHB,NLAN,ILAN,NRLT,IRLT,KERR,LOUT,NRNU,IRNU, &
     &RNU,CKMIN)
              INTEGER(KIND=4) :: NFAR
              INTEGER(KIND=4) :: MDIM
              INTEGER(KIND=4) :: NPAR
              INTEGER(KIND=4) :: MAXSP
              INTEGER(KIND=4) :: II
              INTEGER(KIND=4) :: NSPEC(*)
              REAL(KIND=8) :: PAR(NPAR,*)
              REAL(KIND=8) :: RPAR(NPAR,*)
              CHARACTER(*) :: AUNITS
              CHARACTER(*) :: EUNITS
              INTEGER(KIND=4) :: NREAC(*)
              INTEGER(KIND=4) :: NUNK(MAXSP,*)
              INTEGER(KIND=4) :: NU(MAXSP,*)
              INTEGER(KIND=4) :: KCHRG(*)
              INTEGER(KIND=4) :: MM
              INTEGER(KIND=4) :: KNCF(MDIM,*)
              INTEGER(KIND=4) :: IDUP(*)
              INTEGER(KIND=4) :: NFAL
              INTEGER(KIND=4) :: IFAL(*)
              INTEGER(KIND=4) :: KFAL(*)
              REAL(KIND=8) :: PFAL(NFAR,*)
              INTEGER(KIND=4) :: IFOP(*)
              INTEGER(KIND=4) :: NREV
              INTEGER(KIND=4) :: IREV(*)
              INTEGER(KIND=4) :: NTHB
              INTEGER(KIND=4) :: ITHB(*)
              INTEGER(KIND=4) :: NLAN
              INTEGER(KIND=4) :: ILAN(*)
              INTEGER(KIND=4) :: NRLT
              INTEGER(KIND=4) :: IRLT(*)
              LOGICAL(KIND=4) :: KERR
              INTEGER(KIND=4) :: LOUT
              INTEGER(KIND=4) :: NRNU
              INTEGER(KIND=4) :: IRNU(*)
              REAL(KIND=8) :: RNU(MAXSP,*)
              REAL(KIND=8) :: CKMIN
            END SUBROUTINE CPREAC
          END INTERFACE 
        END MODULE CPREAC__genmod
