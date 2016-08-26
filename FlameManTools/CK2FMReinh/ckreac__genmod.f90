        !COMPILER-GENERATED INTERFACE MODULE: Thu Nov 11 12:25:31 2010
        MODULE CKREAC__genmod
          INTERFACE 
            SUBROUTINE CKREAC(LINE,II,KK,KNAME,LOUT,MAXSP,NSPEC,NREAC,  &
     &NUNK,NU,NPAR,PAR,NTHB,ITHB,NFAL,IFAL,KFAL,NWL,IWL,WL,NRNU,IRNU,RNU&
     &,KERR)
              INTEGER(KIND=4) :: NPAR
              INTEGER(KIND=4) :: MAXSP
              CHARACTER(*) :: LINE
              INTEGER(KIND=4) :: II
              INTEGER(KIND=4) :: KK
              CHARACTER(*) :: KNAME(*)
              INTEGER(KIND=4) :: LOUT
              INTEGER(KIND=4) :: NSPEC(*)
              INTEGER(KIND=4) :: NREAC(*)
              INTEGER(KIND=4) :: NUNK(MAXSP,*)
              INTEGER(KIND=4) :: NU(MAXSP,*)
              REAL(KIND=8) :: PAR(NPAR,*)
              INTEGER(KIND=4) :: NTHB
              INTEGER(KIND=4) :: ITHB(*)
              INTEGER(KIND=4) :: NFAL
              INTEGER(KIND=4) :: IFAL(*)
              INTEGER(KIND=4) :: KFAL(*)
              INTEGER(KIND=4) :: NWL
              INTEGER(KIND=4) :: IWL(*)
              REAL(KIND=8) :: WL(*)
              INTEGER(KIND=4) :: NRNU
              INTEGER(KIND=4) :: IRNU(*)
              REAL(KIND=8) :: RNU(MAXSP,*)
              LOGICAL(KIND=4) :: KERR
            END SUBROUTINE CKREAC
          END INTERFACE 
        END MODULE CKREAC__genmod
