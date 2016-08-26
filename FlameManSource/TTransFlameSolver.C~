#include "FlameMaster.h"
#include "TTransFlameSolver.h"
#include "dassl.h"
#include "MapMan.h"
#include "Interrupt.h"
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
// Tau = 1.0e-3 is too large. Tau = 1.0e-4 seems to be sufficient
// Kappa = 1 seems to be ok

#define CVODE

#define DIFFCORR
#define MOLARDIFF
#define CONVECTION


#define UNSTLIBOUTDELTAT 20.0 //40 by default

#undef LAURENTLIB
#undef ENTHALPYLIB //off by default

#undef CHECKCONVECTION
#define CHIFACT 1.0
#define NEWBETAPDF

#undef READCHI //off by default
#undef TIMEDEPCHI //off by default
#undef ONECOLCHI

#undef READU
#undef READRADHEATLOSS

#define FINEPDF
#undef LEWISCHANGE
#define LEWISSWITCHTIME 5.6e-3
#define LEWISSWITCHPERIOD 0.5e-3

#undef EQUILIBRIUM

// currently grid movement is central MOVEZR is upwind
#undef NEWCONVECTION //off by default
#define INGOOPT
#undef SECONDORDUPZRCONV
#undef CENTRALZRCONV

#undef CENTRALGRID
#define COMPMONFUNC
#undef MOVEGRID //off by default
#define MONFROMTEMP

#ifdef MONFROMTEMP
#	undef MONFROMSPECIES
#else
#	define MONFROMSPECIES
#endif
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo

#undef MOVEZRIGHT

#undef TEMPFROMENTHALPY

#define EXACTCHI
#define NEWCHI
#undef TOTENT
#undef ROSSELAND
#undef NOSOOTRAD
#undef TESTNEWNUCLEATION

#define DELTAPROG 100

#undef SPARSEJAC

#ifdef MOVEGRID
#	undef SPARSEJAC
#endif

#undef LOGNORMALCHI

#define NEWPROPERTIES
#define ENTFLUX
#define HEATCAPGRAD

#undef UPWIND
#define CENTRALTEMP

#define SOOTCONVECTION
#define SIZEDEPDIFFUSION
#define FRACMOMFACT 0.0

#define UPWINDSOOT
#undef SECORDSOOTUPWIND
#undef NEWSOOTTERMDISC
#define NEWTHERMOPHOR

#define UPWINDTHERMOPHOR
#define PRONE
#undef FIRST
#undef SECOND

#ifdef NEWSOOTTERMDISC
#	undef UPWINDSOOT
#endif

#undef DEBUGRES
#undef DEBUGFTOCMAT
#undef DEBUGGRID
#define DEBUGINITIAL //on by default
#undef DEBUGGETSOL
#undef DEBUGBAG
#undef DEBUGSOOTSOURCE
#define DEBUGINITIALDENS //on by default

#define SEMIIMPLI
#undef VECTORTOL
#define FULLIMPLI

#undef DPDT

#ifdef FULLIMPLI
#	define SEMIIMPLI
#	undef VECTORTOL
#	undef DELTAPROG
#	define DELTAPROG 1
#endif

#define DELTAINEW

//Soot -- Mueller
#undef V
#define VS
#undef NOSOOTSOLVE

Double gasdev( int *idum );

void TTransFlameSolver::InitTTransFlameSolver( void )
{
	int	i;

#ifdef ENTFLUX
		fprintf( stderr, "ENTFLUX defined\n" );
#else
		fprintf( stderr, "ENTFLUX undefined\n" );
#endif
#ifdef HEATCAPGRAD
		fprintf( stderr, "HEATCAPGRAD defined\n" );
#else
		fprintf( stderr, "HEATCAPGRAD undefined\n" );
#endif
#ifdef DIFFCORR
		fprintf( stderr, "DIFFCORR defined\n" );
#else
		fprintf( stderr, "DIFFCORR undefined\n" );
#endif
#ifdef MOLARDIFF
		fprintf( stderr, "MOLARDIFF defined\n" );
#else
		fprintf( stderr, "MOLARDIFF undefined\n" );
#endif
#ifdef CONVECTION
		fprintf( stderr, "CONVECTION defined\n" );
#else
		fprintf( stderr, "CONVECTION undefined\n" );
#endif

	if ( fSoot ) {
		if ( fSoot->WithThermoPhoresis() ) {
			fprintf( stderr, "ThermoPhoresis defined\n" );
		}
		else {
			fprintf( stderr, "ThermoPhoresis undefined\n" );
		}

#ifdef FIRST
		fprintf( stderr, "FIRST defined\n" );
#else
		fprintf( stderr, "FIRST undefined\n" );
#endif

#ifdef SECOND
		fprintf( stderr, "SECOND defined\n" );
#else
		fprintf( stderr, "SECOND undefined\n" );
#endif

#ifdef MOVEGRID
		fprintf( stderr, "MOVEGRID defined\n" );
#else
		fprintf( stderr, "MOVEGRID undefined\n" );
#endif

#ifdef CONVECTION
		fprintf( stderr, "CONVECTION defined\n" );
#else
		fprintf( stderr, "CONVECTION undefined\n" );
#endif

#ifdef SOOTCONVECTION
		fprintf( stderr, "SOOTCONVECTION defined\n" );
#else
		fprintf( stderr, "SOOTCONVECTION undefined\n" );
#endif

#ifdef SECORDSOOTUPWIND
		fprintf( stderr, "SECORDSOOTUPWIND defined\n" );
#else
		fprintf( stderr, "SECORDSOOTUPWIND undefined\n" );
#endif

#ifdef SIZEDEPDIFFUSION
		fprintf( stderr, "SIZEDEPDIFFUSION defined\n" );
#else
		fprintf( stderr, "SIZEDEPDIFFUSION undefined\n" );
#endif

#ifdef UPWINDSOOT
		fprintf( stderr, "UPWINDSOOT defined\n" );
#else
		fprintf( stderr, "UPWINDSOOT undefined\n" );
#	ifdef NEWSOOTTERMDISC
		fprintf( stderr, "NEWSOOTTERMDISC defined\n" );
#	endif
#endif
	}
	
// set output file
	if ( !fInputData->fOutFileName || strcmp( fInputData->fOutFileName, "sterr" ) == 0 ) {
		fOutFilePtr = stderr;
		fprintf( stderr, "write output to 'stderr'\n" );
	}
	else if ( strcmp( fInputData->fOutFileName, "stdout" ) == 0 ) {
		fOutFilePtr = stdout;
		fprintf( stderr, "write output to 'stdout'\n" );
	}
	else {
		fOutFilePtr = GetOutfile( fInputData->fOutFileName, kText );
		fprintf( stderr, "write output to file '%s'\n\n"
			, GetOutfileName( fInputData->fOutFileName, kText ) );
	}

#ifdef LOGNORMALCHI
	SetRandomNumber();
#else
	fRandomNumber = 1.0;
#endif
	fFirstCall = TRUE;
	fMaxOrd = 5;
	fEquidistant = fInputData->fEquidistant;
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo	
	fKappa = fInputData->fKappa;
	fTauGrid = fInputData->fTauGrid;
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
	fPrintMolarFractions = fInputData->fPrintMolarFractions;
	if ( fInputData->fDeltaTMax > 0.0 ) {
		fDeltaTMax = fInputData->fDeltaTMax;
	}
	else {
		fDeltaTMax = 1.0e-4;
	}

	fNGridPoints = fInputData->fInitialGridPoints-2;

	fFirstTimeStep = 1.0e-18; //1.0e-18; //1.0e-9 by default

	int	nOfSpecies = fSpecies->GetNOfSpecies();
	int	nOfSpeciesIn = fSpecies->GetNSpeciesInSystem();

#ifdef VECTORTOL
	fRTol = New1DArray( fNOfEquations );
	fATol = New1DArray( fNOfEquations );
	for ( i = 0; i < fNOfEquations; ++i ) {
		fRTol[i] = ( fInputData->fTolRes > 0.0 ) ? fInputData->fTolRes : 1.0e-4;
		fATol[i] = ( fInputData->fTolDy > 0.0 ) ? fInputData->fTolDy : 1.0e-12;
	}
//	fRTol[GetFuelIndex()+fFirstSpecies] = 1.0e-6;
//	fATol[GetFuelIndex()+fFirstSpecies] = 1.0e-13;
#else
	fRTol = New1DArray( 1 );
	fATol = New1DArray( 1 );
	fRTol[0] = ( fInputData->fTolRes > 0.0 ) ? fInputData->fTolRes : 1.0e-4;
	fATol[0] = ( fInputData->fTolDy > 0.0 ) ? fInputData->fTolDy : 1.0e-12;
	cerr << "Tolerance: " << fATol[0] << endl;
#endif
	
// workspace for ddassl
#ifdef 	FULLIMPLI
#	ifdef LOGNORMALCHI
	fprintf( fOutFilePtr, "use gaussian distribution of scalar dissipation\n" );
#	endif
	fNOfSolver = 1;
#ifdef SPARSEJAC
	fML = fNOfEquations;
#else
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
#	ifdef SECONDORDUPZRCONV
	fML = 3 * fNOfEquations;
#	else
	fML = 2 * fNOfEquations;
#	endif
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
#endif
	int	NEQ = ( fNGridPoints + 2 ) * fNOfEquations;

	fNActualStep = new int*[fNOfSolver];
	fNActualOrd = new int*[fNOfSolver];
   	if ( !fNActualStep || !fNActualOrd ) {
		FatalError( "memory allocation of TTransFlameSolver failed" );
	}
	fTime = NewVector( fNOfSolver );
	fSolution = NewMatrix( fNOfEquations, fNGridPoints+2, kColumnPointers );
	fSolution->mat = &fSolution->mat[kNext];
	fSolPrime = NewMatrix( fNOfEquations, fNGridPoints+2, kColumnPointers );
	fSolPrime->mat = &fSolPrime->mat[kNext];

#ifdef CVODE
	fprintf(stderr, "\nUse CVODE solver\n\n");
	fNActualStep[0] = &fNActualStepCV;
	fNActualOrd[0] = &fNActualOrdCV;
	fActualTimeStepSize = &fActualTimeStepSizeCV;
#else
	fprintf(stderr, "\nUse DASSL solver\n\n");
	fLRW =  40 + ( fMaxOrd + 4 ) * NEQ + ( 2 * fML + fML + 1 ) * NEQ + 2 * ( NEQ / ( fML + fML + 1 ) + 1 );
	fLIW = 20 + NEQ;

	fDasslNEq = NEQ;
	fInfo = new IntVectorPtr[fNOfSolver];
	fRWork = new VectorPtr[fNOfSolver];
	fIWork = new IntVectorPtr[fNOfSolver];
	for ( i = 0; i < fNOfSolver; ++i ) {
		fInfo[i] = NewIntVector( 16 );
		fRWork[i] = NewVector( fLRW );
		fIWork[i] = NewIntVector( fLIW );
		fNActualStep[i] = &fIWork[i]->vec[kF11];
		fNActualOrd[i] = &fIWork[i]->vec[kF8];
	}
	InitDassl( fNOfSolver );
	if ( fInfo[0]->vec[kF10] == 1 ) {
		fprintf( fOutFilePtr, "###attention: clip negative concentrations\n" );
	}
	fActualTimeStepSize = &fRWork[0]->vec[kF7];
#endif
#else
	fNOfSolver = fNGridPoints;
	fLRW =  40 + ( fMaxOrd + 4 ) * fNOfEquations + fNOfEquations * fNOfEquations;
	fLIW =  20 + fNOfEquations;

	fDmDy = New2DArray( nOfSpeciesIn+1, nOfSpeciesIn );
	fDasslNEq = fNOfEquations;
#endif

	fMaxStepTaken = New1DArray( fNOfSolver );


	
// other workspace
	
	fSolTime = NewVector( fNGridPoints+2 );
	fSolMassFracs = NewMatrix( nOfSpecies, fNGridPoints+2, kColumnPointers );
	if ( fSoot ) {
		fSolSootMoments = NewMatrix( fSoot->GetNSootMoments(), fNGridPoints+2, kColumnPointers );
		fSolOldSootMoments = NewMatrix( fSoot->GetNSootMoments(), fNGridPoints+2, kColumnPointers );
		fSootMomentsWork = NewMatrix( nOfSpecies, fNGridPoints+2, kColumnPointers );

		fSolSootMoments->mat = &fSolSootMoments->mat[kNext];
		fSolOldSootMoments->mat = &fSolOldSootMoments->mat[kNext];
		fSootMomentsWork->mat = &fSootMomentsWork->mat[kNext];
	}
	else {
		fSolSootMoments = NULL;
		fSolOldSootMoments = NULL;
		fSootMomentsWork = NULL;
	}
	fSolTemp = NewVector( fNGridPoints+2 );
	fSolProg = NewVector( fNGridPoints+2 );
	fSolEnth = NewVector( fNGridPoints+2 );
	fSolOldTime = NewVector( fNGridPoints+2 );
	fSolOldMassFracs = NewMatrix( nOfSpecies, fNGridPoints+2, kColumnPointers );
	fSolOldTemp = NewVector( fNGridPoints+2 );
	fSolOldProg = NewVector( fNGridPoints+2 );
	fSolOldEnth = NewVector( fNGridPoints+2 );
	fSolGrid = NewVector( fNGridPoints+2 );
	fFDWCurr = New1DArray( fNGridPoints );
	fFDWPlus = New1DArray( fNGridPoints );
	fFDWMinus = New1DArray( fNGridPoints );
	fWCurr = New1DArray( fNGridPoints );
	fWPlus = New1DArray( fNGridPoints );
	fWMinus = New1DArray( fNGridPoints );
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
    fh = New1DArray( fNGridPoints+2 );
    fhm = New1DArray( fNGridPoints+2 );	
	fMonFct = New1DArray( fNGridPoints+2 );
	fgpDens = NewVector( fNGridPoints+2 );
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
	fMassFracsWork = NewMatrix( nOfSpecies, fNGridPoints+2, kColumnPointers );
	fTempWork = NewVector( fNGridPoints+2 );
	fProgWork = NewVector( fNGridPoints+2 );
	fEnthWork = NewVector( fNGridPoints+2 );
	fOutSolWork = NewVector( fNGridPoints+2 );
	fTotEntStart = NewVector( fNGridPoints+2 );
	fTotEntEnd = NewVector( fNGridPoints+2 );

	fHeatCpMix = NewVector( fNGridPoints+2 );
	fViscosity = NewVector( fNGridPoints+2 );
	fMolarMassMix = NewVector( fNGridPoints+2 );
	fDensity = NewVector( fNGridPoints+2 );
	fLambdaOverCpMix = NewVector( fNGridPoints+2 );

	fSpecHeatCp = NewMatrix( nOfSpeciesIn, fNGridPoints+2, kColumnPointers );
	fSpecEnthalpy = NewMatrix( nOfSpeciesIn, fNGridPoints+2, kColumnPointers );
	fProdRate = NewMatrix( nOfSpeciesIn, fNGridPoints+2, kColumnPointers );

	fDiffTermY = NewVector( nOfSpeciesIn );
	fDiffTermW = NewVector( nOfSpeciesIn );
	fDiffCorrY = NewVector( nOfSpeciesIn );
	fDiffCorrW = NewVector( nOfSpeciesIn );

	fSolTime->vec = &fSolTime->vec[kNext];
	fSolMassFracs->mat = &fSolMassFracs->mat[kNext];
	fSolTemp->vec = &fSolTemp->vec[kNext];
	fSolProg->vec = &fSolProg->vec[kNext];
	fSolEnth->vec = &fSolEnth->vec[kNext];
	fSolOldTime->vec = &fSolOldTime->vec[kNext];
	fSolOldMassFracs->mat = &fSolOldMassFracs->mat[kNext];
	fSolOldTemp->vec = &fSolOldTemp->vec[kNext];
	fSolOldProg->vec = &fSolOldProg->vec[kNext];
	fSolOldEnth->vec = &fSolOldEnth->vec[kNext];
	fSolGrid->vec = &fSolGrid->vec[kNext];
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
    fh = &fh[kNext];
    fhm = &fhm[kNext];
    fgpDens->vec = &fgpDens->vec[kNext];
    fMonFct = &fMonFct[kNext];
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
	fMassFracsWork->mat = &fMassFracsWork->mat[kNext];
	fTempWork->vec = &fTempWork->vec[kNext];
	fProgWork->vec = &fProgWork->vec[kNext];
	fEnthWork->vec = &fEnthWork->vec[kNext];
	fTotEntStart->vec = &fTotEntStart->vec[kNext];
	fTotEntEnd->vec = &fTotEntEnd->vec[kNext];

	fHeatCpMix->vec = &fHeatCpMix->vec[kNext];
	fViscosity->vec = &fViscosity->vec[kNext];
	fMolarMassMix->vec = &fMolarMassMix->vec[kNext];
	fDensity->vec = &fDensity->vec[kNext];
	fLambdaOverCpMix->vec = &fLambdaOverCpMix->vec[kNext];

	fSpecHeatCp->mat = &fSpecHeatCp->mat[kNext];
	fSpecEnthalpy->mat = &fSpecEnthalpy->mat[kNext];
	fProdRate->mat = &fProdRate->mat[kNext];

	fZRin = NewVector( 10 );
	for ( i = 0; i < fZRin->len; ++i ) {
		fZRin->vec[i] = 1.0;
	}

	fMaxVals = NewVector( fNOfEquations );

//	SetWeights();

	if ( fSoot ) {
		fSoot->SetMomentsOffset( fSootMoments );
	}

	// set variable names
	fVariableNames = new String[fNOfEquations];

	fVariableNames[fTemperature] = new char[2];
	strcpy( fVariableNames[fTemperature], "T" );

	fVariableNames[fProg] = new char[5];
	strcpy( fVariableNames[fProg], "prog" );

	fVariableNames[fEnth] = new char[5];
	strcpy( fVariableNames[fEnth], "enth" );

	fVariableNames[fGrid] = new char[5];
	strcpy( fVariableNames[fGrid], "Grid" );

	for ( i = 0; i < nOfSpeciesIn; ++i ) {
		fVariableNames[fFirstSpecies + i] = new char[strlen( fSpecies->GetNames()[i] ) + 1];
		strcpy( fVariableNames[fFirstSpecies + i], fSpecies->GetNames()[i] );
	}
	if (fSoot)
	{
	  int offset = fSoot->GetOffsetSootMoments();
	  for (i = 0; i < fSoot->GetNSootMoments(); ++i)
	  {
	    int j = fSoot->Geti(i);
	    int k = fSoot->Getj(i);
#ifdef V
	    fVariableNames[offset + i] = new char[8];
	    sprintf(fVariableNames[offset + i], "M%d", j);
#endif
#ifdef VS
	    fVariableNames[offset + i] = new char[11];
	    sprintf(fVariableNames[offset + i], "M%d-%d", j, k);
#endif
	  }
	}
	
	FILE *fpNames = GetOutfile( "VarNames", kText );
	for ( i = 0; i < fNOfEquations; ++i ) {
		fprintf( fpNames, "%s\n", fVariableNames[i] );
	}
	fclose( fpNames );

	fpNames = GetOutfile( "WanNames", kText );
	for ( i = 0; i < fNOfEquations; ++i ) {
		fprintf( fpNames, "      VANAMES(%d)='%-20s'\n", i+1, fVariableNames[i] );
	}
	fclose( fpNames );

	CompLewisNumbers( fInputData->fLewisNumberFile );
	fflush( fOutFilePtr );

	char	buff[128];
	char	*theFile;
	int		count = 0;
	FILE	*fpCA = NULL;
	Flag	isNew = FALSE;
	do {
		sprintf( buff, "CA%d.in", count );
		theFile = GetOutfileName( buff, kNone );
		fpCA = fopen( theFile, "r" );
		if ( fpCA ) {
			fclose( fpCA );
			++count;
		}
		else {
			isNew = TRUE;
		}
	} while( !isNew && count < 50 );
	
	if ( !isNew ) {
		fprintf( stderr, "cannot open CA.in file ( count = %d )\n", count );
		exit(2);
	}
	else {
		fCAinFile = fopen( theFile, "w" );
		fprintf( fCAinFile, "RPM = -1\n" );
		fprintf( fCAinFile, "VarsIn = 6\n" );
		fprintf( fCAinFile, "Time(s)	Pressure(Pa)	TOx(K)	TFuel(K)	Sci(1/s)	ZR\n" );
		fflush( fCAinFile );
	}

#ifdef READCHI
#	ifdef TIMEDEPCHI
	fprintf( stderr, "### use scalar dissipation rate from input file\n" );
	Double	*ZIn = New1DArray( 1000 );
	Double	*TimeIn = New1DArray( 1000 );
	char 	dummy[128];
	int		conv, j;

// read Z
	FILE *fpZIn = fopen( "zi.dat", "r" );
	if ( !fpZIn ) {
		fprintf( stderr, "#error: couldn't open input file 'zi.dat'\n" );
		exit( 2 );
	}
	for ( j = 0; j < 1000; ++j ) {
		conv = fscanf( fpZIn, "%lg", &ZIn[j] );
		//fprintf( stderr, "Zin[%d] = %g\n", j, ZIn[j] );
		if ( !conv || conv == EOF ) {
			break;
		}
	}
	if ( j == 0 ) {
		fprintf( stderr, "### error reading file 'zi.dat'\n" );
		exit( 2 );
	}
	fZIn = NewVector( j );

	for ( j = 0; j < fZIn->len; ++j ) {
		fZIn->vec[j] = ZIn[j];
	}
	fclose(fpZIn);

// read Z
	FILE *fpTimeIn = fopen( "Time.dat", "r" );
	if ( !fpTimeIn ) {
		fprintf( stderr, "#error: couldn't open input file 'Time.dat'\n" );
		exit( 2 );
	}
	for ( j = 0; j < 1000; ++j ) {
		conv = fscanf( fpTimeIn, "%lg", &TimeIn[j] );
		//fprintf( stderr, "Timein[%d] = %g\n", j, TimeIn[j] );
		if ( !conv || conv == EOF ) {
			break;
		}
	}
	if ( j == 0 ) {
		fprintf( stderr, "### error reading file 'Time.dat'\n" );
		exit( 2 );
	}
	fTimeIn = NewVector( j );

	for ( j = 0; j < fTimeIn->len; ++j ) {
		fTimeIn->vec[j] = TimeIn[j];
	}
	fclose(fpTimeIn);

// read Chi
	fChiIn = NewMatrix( fZIn->len, fTimeIn->len, kColumnPointers );
	Double	**chiIn = fChiIn->mat;
	FILE *fpChiCount = fopen( "chi.dat", "r" );
	if ( !fpChiCount ) {
		fprintf( stderr, "#error: couldn't open input file 'chi.dat'\n" );
		exit( 2 );
	}

#		ifdef ONECOLCHI
	for ( i = 0; i < fTimeIn->len; ++i ) {
		for ( j = 0; j < fZIn->len; ++j ) {
			conv = fscanf( fpChiCount, "%lg", &chiIn[i][j] );
			if ( !conv || conv == EOF ) {
				fprintf( stderr, "something wrong at i=%d j=%d\n", i, j );
				fprintf( stderr, "iLen actually %d\tjLen actually %d\n", fTimeIn->len, fZIn->len );
				break;
			}
		}
	}
#		else
	for ( j = 0; j < fZIn->len; ++j ) {
		for ( i = 0; i < fTimeIn->len; ++i ) {
			conv = fscanf( fpChiCount, "%lg", &chiIn[i][j] );
			if ( !conv || conv == EOF ) {
				fprintf( stderr, "something wrong at i=%d j=%d\n", i, j );
				break;
			}
		}
	}
#		endif
	fclose( fpChiCount );

/* extrapolate chi to ZMean + 2 sqrt(ZVar)*/
	// first read ZR
	VectorPtr	ZRinVec = NewVector( fTimeIn->len );
	Double		*zrin = ZRinVec->vec;
	FILE 		*fpZRin = fopen( "ZR.in", "r" );
	if ( !fpZRin ) {
		fprintf( stderr, "#error: couldn't open input file 'ZR.in' -> don't extrapolate\n" );
		for ( i = 0; i < fTimeIn->len; ++i ) {
			zrin[i] = 1.0;
		}
	}
	else {
		for ( i = 0; i < fTimeIn->len; ++i ) {
			conv = fscanf( fpZRin, "%lg", &zrin[i] );
			zrin[i] = MIN( 1.0, zrin[i] );
			if ( !conv || conv == EOF ) {
				fprintf( stderr, "something wrong at i=%d j=%d\n", i );
				break;
			}
		}
	}
	
	int zrind, zc;
	for ( i = 0; i < fTimeIn->len; ++i ) {
		// get last non-zero point
		j = fZIn->len-2;
		while ( j > 0 && chiIn[i][j] < 1.0e-8 ) --j;
		// now j points at last nonzero scalar diss rate
		// go to next if zr is not larger than first zero chi mixture fraction
		fprintf( stderr, "do i = %d\tzr = %g\tzl = %g\tzl+1 = %g\n", i, zrin[i], fZIn->vec[j], fZIn->vec[j+1] );
		if ( j == 0 || fZIn->vec[j+1] >= zrin[i]-1.0e-8  ) {
			continue;
		}
		// then get zrpoint 
		zrind = j+1;
		while ( fZIn->vec[zrind] < zrin[i]-1.0e-8 ) ++zrind;
		// now zrind points at first z larger than zr
		for ( zc = j+1; zc < zrind; ++zc ) {
			chiIn[i][zc] = chiIn[i][j] * ( 1.0 - ( fZIn->vec[zc] - fZIn->vec[j] ) / ( zrin[i] - fZIn->vec[j] ) );
		}
	}
	
	

	FILE	*fpZR = GetOutfile( "ZR", kData );
	fprintf( fpZR, "*\ntime\tZR\n" );

	for ( i = 0; i < fTimeIn->len; ++i ) {
		j = fZIn->len-2;
		while ( chiIn[i][j] < 1.0e-8 ) --j;
		fprintf( fpZR, "%g\t%g\n", fTimeIn->vec[i], fZIn->vec[j+1] );
	}
	fclose( fpZR );

	FILE	*fp = GetOutfile( "ChiIn", kData );
	fprintf( fp, "*\nZ" );
	for ( i = 0; i < fTimeIn->len; ++i ) {
		fprintf( fp, "\tt%g", fTimeIn->vec[i] );
	}
	fprintf( fp, "\n" );

	for ( j = 0; j < fZIn->len; ++j ) {
		fprintf( fp, "%g", fZIn->vec[j] );
		for ( i = 0; i < fTimeIn->len; ++i ) {
			fprintf( fp, "\t%g", chiIn[i][j] );
		}
		fprintf( fp, "\n" );
	}
	fclose( fp );
	Free1DArray( ZIn );
	Free1DArray( TimeIn );
#	else
	Double	*ZIn = New1DArray( 1000 );
	Double	*chiIn = New1DArray( 1000 );
	char 	dummy[128];
	int		conv, j;

	FILE *fpChiCount = fopen( "Chi.tout", "r" );
	if ( !fpChiCount ) {
		fprintf( stderr, "#error: couldn't open input file 'Chi.tout'\n" );
		exit( 2 );
	}

	fscanf( fpChiCount, "%s%s", dummy, dummy );
	for ( j = 0; j < 1000; ++j ) {
		conv = fscanf( fpChiCount, "%lg%lg", &ZIn[j], &chiIn[j] );
		if ( !conv || conv == EOF ) {
			break;
		}
	}
	fclose( fpChiCount );

	fZCount = NewVector( j );
	fChiCount = NewVector( j );
	for ( j = 0; j < fZCount->len; ++j ) {
		fZCount->vec[j] = ZIn[j];
		fChiCount->vec[j] = chiIn[j];
	}
/*	fprintf( stderr, "*\nZ\tChi\n" );
	for ( j = 0; j < fZCount->len; ++j ) {
		fprintf( stderr, "%g\t%g\n", fZCount->vec[j], fChiCount->vec[j] );
	}*/

#	endif
#endif

#ifdef READU
#	ifndef READCHI
	fprintf( stderr, "###error: can't use READU without READCHI\n" );
	exit( 2 );
#	endif
#	ifndef TIMEDEPCHI
	fprintf( stderr, "###error: can't use READU without TIMEDEPCHI\n" );
	exit( 2 );
#	endif
	fprintf( stderr, "### use axial velocity from input file\n" );

// read U
	fUstOverU = NewMatrix( fZIn->len, fTimeIn->len, kColumnPointers );
	Double	**UIn = fUstOverU->mat;
	FILE *fpUCount = fopen( "Ucond.dat", "r" );
	if ( !fpUCount ) {
		fprintf( stderr, "#error: couldn't open input file 'Ucond.dat'\n" );
		exit( 2 );
	}

	for ( i = 0; i < fTimeIn->len; ++i ) {
		for ( j = 0; j < fZIn->len; ++j ) {
			conv = fscanf( fpUCount, "%lg", &UIn[i][j] );
			if ( !conv || conv == EOF ) {
				fprintf( stderr, "something wrong with UIn at i=%d j=%d\n", i, j );
				fprintf( stderr, "iLen actually %d\tjLen actually %d\n", fTimeIn->len, fZIn->len );
				break;
			}
			if ( UIn[i][j] <= 0.0 ) {
				if ( j == 0 ) {
					if ( i > 0 ) {
						UIn[i][j] = UIn[i-1][j];
					}
					else {
						// check next value
						++j;
						conv = fscanf( fpUCount, "%lg", &UIn[i][j] );
						if ( !conv || conv == EOF ) {
							fprintf( stderr, "something wrong with UIn at i=%d j=%d\n", i, j );
							fprintf( stderr, "iLen actually %d\tjLen actually %d\n", fTimeIn->len, fZIn->len );
							break;
						}
						if ( UIn[i][j] <= 0.0 ) {
							fprintf( stderr, "###error: something wrong with Uin 1\n" );
							exit( 2 );
						}
						else {
							UIn[i][j-1] = UIn[i][j];
						}
					}
				}
				else {
					UIn[i][j] = UIn[i][j-1];
				}
			}
		}
	}
	fclose( fpUCount );

	fp = GetOutfile( "UIn", kData );
	fprintf( fp, "*\nZ" );
	for ( i = 0; i < fTimeIn->len; ++i ) {
		fprintf( fp, "\tt%g", fTimeIn->vec[i] );
	}

	for ( j = 0; j < fZIn->len; ++j ) {
		fprintf( fp, "%g", fZIn->vec[j] );
		for ( i = 0; i < fTimeIn->len; ++i ) {
			fprintf( fp, "\t%g", UIn[i][j] );
		}
		fprintf( fp, "\n" );
	}
	fclose( fp );
	
#endif

#ifdef READRADHEATLOSS
#	ifndef READCHI
	fprintf( stderr, "###error: can't use READRADHEATLOSS without READCHI\n" );
	exit( 2 );
#	endif
#	ifndef TIMEDEPCHI
	fprintf( stderr, "###error: can't use READRADHEATLOSS without TIMEDEPCHI\n" );
	exit( 2 );
#	endif
	fprintf( stderr, "### use axial velocity from input file\n" );

// read U
	fDelQ = NewMatrix( fZIn->len, fTimeIn->len, kColumnPointers );
	Double	**delQ = fDelQ->mat;
	FILE *fpdelQ = fopen( "delQ.dat", "r" );
	if ( !fpdelQ ) {
		fprintf( stderr, "#error: couldn't open input file 'delQ.dat'\n" );
		exit( 2 );
	}

	for ( i = 0; i < fTimeIn->len; ++i ) {
		for ( j = 0; j < fZIn->len; ++j ) {
			conv = fscanf( fpdelQ, "%lg", &delQ[i][j] );
			if ( !conv || conv == EOF ) {
				fprintf( stderr, "something wrong with delQ at i=%d j=%d\n", i, j );
				fprintf( stderr, "iLen actually %d\tjLen actually %d\n", fTimeIn->len, fZIn->len );
				break;
			}
		}
	}
	fclose( fpdelQ );
	int jjNow;
	for ( i = 0; i < fTimeIn->len; ++i ) {
		for ( j = 1; j < fZIn->len-1; ++j ) {
			if ( fabs(delQ[i][j]) <= 1e-5 && fabs(delQ[i][j-1]) > 1e-5 ) {
				jjNow = j;
				while ( jjNow < fZIn->len-1 ) {
					if ( fabs(delQ[i][jjNow]) > 1e-5 ) {
						delQ[i][j] = Interpol( fZIn->vec[j], delQ[i][j-1], fZIn->vec[j-1], delQ[i][jjNow], fZIn->vec[jjNow] );
						break;
					}
					++jjNow;
				}
			}
		}
	}

	fp = GetOutfile( "delQ", kData );
	fprintf( fp, "*\nZ" );
	for ( i = 0; i < fTimeIn->len; ++i ) {
		fprintf( fp, "\tt%g", fTimeIn->vec[i] );
	}

	for ( j = 0; j < fZIn->len; ++j ) {
		fprintf( fp, "%g", fZIn->vec[j] );
		for ( i = 0; i < fTimeIn->len; ++i ) {
			fprintf( fp, "\t%g", delQ[i][j] );
		}
		fprintf( fp, "\n" );
	}
	fclose( fp );
	
#endif

#ifdef DELTAINEW
	fTempGSave = New1DArray( fNGridPoints+2 );
    fYGSave = New2DArray( fNGridPoints+2, nOfSpeciesIn );
    fDeltaI = New2DArray( fNGridPoints+2, nOfSpeciesIn );
	fG_ij = New3DArray( fNGridPoints+2, nOfSpeciesIn, nOfSpeciesIn );
	fTempGSave = &fTempGSave[kNext];
    fYGSave = &fYGSave[kNext];
    fDeltaI = &fDeltaI[kNext];
	fG_ij = &fG_ij[kNext];
#endif

	// Density Correction
	rhodot = NewVector(fNGridPoints+2);
	enthdot = NewVector(fNGridPoints+2);

}

TTransFlameSolver::~TTransFlameSolver( void )
{
#ifdef DELTAINEW
	fG_ij = &fG_ij[kPrev];
	fDeltaI = &fDeltaI[kPrev];
	fYGSave = &fYGSave[kPrev];
	fTempGSave = &fTempGSave[kPrev];
	Free3DArray( fG_ij );
	Free2DArray( fDeltaI );
	Free2DArray( fYGSave );
	Free1DArray( fTempGSave );
#endif

#ifdef READU
	DisposeMatrix( fUstOverU );
#endif

#ifdef READCHI
#	ifdef TIMEDEPCHI
	DisposeMatrix( fChiIn );
	DisposeVector( fZIn );
	DisposeVector( fTimeIn );
#	else
	DisposeVector( fChiCount );
	DisposeVector( fZCount );
#	endif
#endif

	int nOfSpeciesIn = fSpecies->GetNSpeciesInSystem();

	for ( int i = 0; i < nOfSpeciesIn+fVariablesWithoutSpecies; ++i ) {
		delete fVariableNames[i];
	}
	delete fVariableNames;
	
	fMonFct = &fMonFct[kPrev];
	fhm = &fhm[kPrev];
	fh = &fh[kPrev];
	fgpDens->vec = &fgpDens->vec[kPrev];
	fTempWork->vec = &fTempWork->vec[kPrev];
	fProgWork->vec = &fProgWork->vec[kPrev];
	fEnthWork->vec = &fEnthWork->vec[kPrev];
	fMassFracsWork->mat = &fMassFracsWork->mat[kPrev];
	fSolGrid->vec = &fSolGrid->vec[kPrev];
	fSolOldTemp->vec = &fSolOldTemp->vec[kPrev];
	fSolOldProg->vec = &fSolOldProg->vec[kPrev];
	fSolOldEnth->vec = &fSolOldEnth->vec[kPrev];
	fSolOldMassFracs->mat = &fSolOldMassFracs->mat[kPrev];
	fSolOldTime->vec = &fSolOldTime->vec[kPrev];
	fSolTemp->vec = &fSolTemp->vec[kPrev];
	fSolProg->vec = &fSolProg->vec[kPrev];
	fSolEnth->vec = &fSolEnth->vec[kPrev];
	fSolMassFracs->mat = &fSolMassFracs->mat[kPrev];
	fSolTime->vec = &fSolTime->vec[kPrev];
	fTotEntStart->vec = &fTotEntStart->vec[kPrev];
	fTotEntEnd->vec = &fTotEntEnd->vec[kPrev];

	fHeatCpMix->vec = &fHeatCpMix->vec[kPrev];
	fViscosity->vec = &fViscosity->vec[kPrev];
	fMolarMassMix->vec = &fMolarMassMix->vec[kPrev];
	fDensity->vec = &fDensity->vec[kPrev];
	fLambdaOverCpMix->vec = &fLambdaOverCpMix->vec[kPrev];

	fSpecHeatCp->mat = &fSpecHeatCp->mat[kPrev];
	fSpecEnthalpy->mat = &fSpecEnthalpy->mat[kPrev];
	fProdRate->mat = &fProdRate->mat[kPrev];

//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
	DisposeVector( fgpDens );
	Free1DArray( fMonFct );
	Free1DArray( fhm );
	Free1DArray( fh );
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
	DisposeVector( fOutSolWork );
	DisposeVector( fTempWork );
	DisposeVector( fProgWork );
	DisposeVector( fEnthWork );
	DisposeMatrix( fMassFracsWork );
	Free1DArray( fWMinus );
	Free1DArray( fWPlus );
	Free1DArray( fWCurr );
	Free1DArray( fFDWMinus );
	Free1DArray( fFDWPlus );
	Free1DArray( fFDWCurr );
	DisposeVector( fSolGrid );
	DisposeVector( fSolOldTemp );
	DisposeVector( fSolOldProg );
	DisposeVector( fSolOldEnth );
	DisposeMatrix( fSolOldMassFracs );
	DisposeVector( fSolOldTime );
	DisposeVector( fSolTemp );
	DisposeVector( fSolProg );
	DisposeVector( fSolEnth );
	DisposeVector( fTotEntStart );
	DisposeVector( fTotEntEnd );
	DisposeVector( fZRin );
	DisposeVector( fDiffCorrW );
	DisposeVector( fDiffCorrY );
	DisposeVector( fDiffTermW );
	DisposeVector( fDiffTermY );
	DisposeMatrix( fProdRate );
	DisposeMatrix( fSpecEnthalpy );
	DisposeVector( fLambdaOverCpMix );
	DisposeVector( fDensity );
	DisposeVector( fMolarMassMix );
	DisposeVector( fHeatCpMix );
	DisposeVector( fViscosity );
	if ( fSoot ) {
		fSootMomentsWork->mat = &fSootMomentsWork->mat[kPrev];
		fSolOldSootMoments->mat = &fSolOldSootMoments->mat[kPrev];
		fSolSootMoments->mat = &fSolSootMoments->mat[kPrev];
		DisposeMatrix( fSootMomentsWork );
		DisposeMatrix( fSolOldSootMoments );
		DisposeMatrix( fSolSootMoments );
	}
	DisposeMatrix( fSolMassFracs );
	DisposeVector( fSolTime );

	fSolPrime->mat = &fSolPrime->mat[kPrev];
	fSolution->mat = &fSolution->mat[kPrev];
	DisposeMatrix( fSolPrime );
	DisposeMatrix( fSolution );


#ifndef CVODE
	DisposeIntVector( fIWork[0] );
	DisposeVector( fRWork[0] );
	DisposeIntVector( fInfo[0] );	

	delete fIWork;
	delete fRWork;
	delete fInfo;
#endif

	Free1DArray( fMaxStepTaken );
	delete fNActualOrd;
	delete fNActualStep;

	DisposeVector( fTime );

	Free1DArray( fATol );
	Free1DArray( fRTol );

	DisposeVector(rhodot);
	DisposeVector(enthdot);
}

void TTransFlameSolver::SolutionToCVode( void )
{
	copy( ( fNGridPoints + 2 ) * fNOfEquations, fSolution->mat[kPrev], 1, fCYdata, 1 );
}

void TTransFlameSolver::SolutionFromCVode( void )
{
	copy( ( fNGridPoints + 2 ) * fNOfEquations, fCYdata, 1, fSolution->mat[kPrev], 1 );

	//Enforce non-negative species
	for (int k = -1; k < fNGridPoints; k++)
	  for (int i = 0; i < fSpecies->GetNOfSpecies(); i++)
	    fSolution->mat[k][fFirstSpecies+i] = max(fSolution->mat[k][fFirstSpecies+i],1.0e-60);
}

int TTransFlameSolver::InitCVODE( void )
{
	int flag;  // for checking return values from cvode functions
	int	NEQ = ( fNGridPoints + 2 ) * fNOfEquations;

  /* set the solution vector and initialise from solution YInit
     also assign pointer to first memory location in Y
     NOTE: use compiler directive here later for parallel */
	fCVY = N_VMake_Serial(NEQ, fSolution->mat[kPrev]);

	fCYdata = NV_DATA_S(fCVY);

  /* create the memory object for CVODE:
     The system is stiff due to chemistry, thus choose solver appropriately
       CV_BDF:    linear multistep method (Backward Differentiation Formulas) 
       CV_NEWTON: Newton iteration (modified newton iteration for banded) */
	fMem = CVodeCreate(CV_BDF, CV_NEWTON);

  /* allocate internal memory object:
     initialise time:tInit
     RHS function:   cvodeRHS
     CV_SS:          scalar relative and absolute tolerances
     reltol:         relative tolerance
     abstol:         absolute tolerance */
#ifdef SUNDIALS23
    flag = CVodeMalloc(fMem, ResTransFlameImpliSolverCV, fTStart, fCVY, CV_SS, fRTol[0], fATol);
#else

	CVodeInit(fMem, ResTransFlameImpliSolverCV, fTStart, fCVY);
      
	flag = CVodeSStolerances(fMem, fRTol[0], fATol[0]);
#endif

	if (flag != CV_SUCCESS) {
      if (flag == CV_MEM_NULL)
		printf("CVODE ERROR: CVodeMalloc requires CVodeCreate");
		printf("exiting: unable to allocate with CVodeMalloc\n");
		exit(1);
	}

  // set pointer to parameters (structs) needed by the RHS function
#ifdef SUNDIALS23
  flag = CVodeSetFdata(fMem, this);
#else
  flag = CVodeSetUserData(fMem, this);
#endif
  
  // defined dense matrix, set number of equations
  flag = CVBand(fMem, NEQ, fML, fML);
  
  // options
  CVodeSetInitStep(fMem, fFirstTimeStep );
  CVodeSetMaxStep(fMem, fDeltaTMax );
  
  return flag;
}

void TTransFlameSolver::SetInfo( int i )
{
	int		*info = fInfo[i]->vec;

	info[kF1] = 0;	// first call
#ifdef VECTORTOL
	info[kF2] = 1;	// RTOL, ATOL are scalars
#else
	info[kF2] = 0;	// RTOL, ATOL are scalars
#endif
	info[kF3] = 1;	// solution at TOUT, no intermediates
	info[kF4] = 1;	// integration can be carried out beyond TOUT ( and interpolated )
//	if ( fUseNumericalDM ) {
		info[kF5] = 0;	// numerical jacobian
//	}
//	else {
//		info[kF5] = 1;
//	}
#ifdef FULLIMPLI
	info[kF6] = 1;	// full ( not banded ) jacobian
#else
	info[kF6] = 0;	// full ( not banded ) jacobian
#endif
	info[kF7] = 1;	// no choice for maximum stepsize
	info[kF8] = 1;	// no choice for first stepsize
	if ( fMaxOrd == 5 ) {
		info[kF9] = 0;	// choose default for MAXORD ( = 5 )
	}
	else {
		info[kF9] = 1;	// don't choose default for MAXORD ( = 5 )
	}
//	info[kF10] = 1;	// solve without non negative constraint // this is clip
	info[kF10] = 0;	// solve without non negative constraint // this is don't clip
	info[kF11] = 1;	// initial values are consistent
}

void TTransFlameSolver::InitDassl( int nGridPoints )
{

	for ( int i = 0; i < nGridPoints; ++i ) {
		SetInfo( i );
		InitDasslWorkSpace( i );
	}
}

void TTransFlameSolver::InitDasslWorkSpace( int i )
{
	int		*info = fInfo[i]->vec;
	int		*iWork = fIWork[i]->vec;
	Double	*rWork = fRWork[i]->vec;

	if ( info[kF6] == 1 ) {
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
		iWork[kF1] = fML;
		iWork[kF2] = fML;		
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
	}
	if ( info[kF8] == 1 ) {
		rWork[kF3] = fFirstTimeStep;
	}
	if ( info[kF9] == 1 ) {
		iWork[kF3] = fMaxOrd;
	}
	if ( info[kF7] == 1 ) {
		rWork[kF2] = fDeltaTMax;
	}
}

#ifdef SEMIIMPLI
void TTransFlameSolver::SetMaxTimeStep( int kAct, Double /*t*/ )
{
	fRWork[kAct]->vec[kF2] = fDeltaTMax;
#else
void TTransFlameSolver::SetMaxTimeStep( int kAct, Double t )
{
	if ( fInfo[kAct]->vec[kF7] == 1 ) {
		Double	chi = GetDissRate( t, fSolGrid->vec[kAct] );
//		Double	chi = Interpol( t, fChiStart, fTStart, fChiEnd, fTEnd );
		if ( chi ) {
			Double	dZSquare = ( fSolGrid->vec[kAct+1] - fSolGrid->vec[kAct] ) *
							( fSolGrid->vec[kAct] - fSolGrid->vec[kAct-1] );
			fRWork[kAct]->vec[kF2] = 0.5 * dZSquare / chi;
		}
		else {
			fRWork[kAct]->vec[kF2] = fDeltaTMax;
		}
	}
//	fprintf( fOutFilePtr, "maxStep = %g\n", fRWork[kAct]->vec[kF2] );
#endif
#ifdef SEMIIMPLI
}
#else
}
#endif

void TTransFlameSolver::GetSpeciesData( char **newNames, MatrixPtr highMat
					, MatrixPtr lowMat, Double *molarMass )
{
	int		ind;
	Double	*theMolarMass = fSpecies->GetMolarMass()->vec;
	Double	**high = highMat->mat;
	Double	**low = lowMat->mat;
	Double	**theHigh = fSpecies->GetCoeffHigh()->mat;
	Double	**theLow = fSpecies->GetCoeffLow()->mat;
	
	for ( int i = 0; i < highMat->cols; ++i ) {
		if ( ( ind = fSpecies->FindSpecies( newNames[i] ) ) < 0 ) {
			fprintf( fOutFilePtr, "%s'%s'\n", "#error in function 'GetSpeciesData': no species ", newNames[i]  );
			exit( 2 );
		}
		copy( NOFCOEFFS, theLow[ind], 1, low[i], 1 );
		copy( NOFCOEFFS, theHigh[ind], 1, high[i], 1 );
		molarMass[i] = theMolarMass[ind];
	}
}

void GETSPECIESNAMES( int **object, char *names, int *nChars, int *len )
{
	getspeciesnames( object, names, nChars, len );
}

void getspeciesnames( int **object, char *names, int *nChars, int *len )
{
	TTransFlameSolverPtr flame = *( TTransFlameSolverPtr *)object;
	
	CToFortranCharArray( names, flame->GetSpecies()->GetNames(), *nChars
		, ( flame->GetSpecies()->GetNOfSpecies() > *len ) 
				? *len : flame->GetSpecies()->GetNOfSpecies() );
}

void GETFUEL( int **object, char *name, int *nChars )
{
	getfuel( object, name, nChars );
}

void getfuel( int **object, char *name, int *nChars )
{
	TTransFlameSolverPtr flame = *( TTransFlameSolverPtr *)object;
	
	CToFortranChar( name, flame->GetSpecies()->GetNames()[flame->GetFuelIndex()], *nChars );
}

int GETMAXSPECNAMELEN( int **object )
{
	return getmaxspecnamelen( object );
}

int getmaxspecnamelen( int **object )
{
	TTransFlameSolverPtr	flame = *( TTransFlameSolverPtr *)object;
	char					**names = flame->GetSpecies()->GetNames();
	int						nSpecies = flame->GetSpecies()->GetNOfSpecies();
	int 					maxlen = strlen( names[0] );
	
	for ( int i = 1; i < nSpecies; ++i ) {
		maxlen = maxint( strlen( names[i] ), maxlen );
	}

	return maxlen;
}

void GETNOFSPECIES( int **object, int *nSpecies, int *nSpeciesIn )
{
	getnofspecies( object, nSpecies, nSpeciesIn );
}

void getnofspecies( int **object, int *nSpecies, int *nSpeciesIn )
{
	TTransFlameSolverPtr flame = *( TTransFlameSolverPtr *)object;
	
	*nSpecies = flame->GetSpecies()->GetNOfSpecies();
	*nSpeciesIn = flame->GetSpecies()->GetNSpeciesInSystem();
}

void TTransFlameSolver::Initialize( Double timeStart
					, ConstStringArray names
					, Double **startSolution, int vars
					, Double *grid, int gridPointsA
					, Double pressureStart, Double scalarDissRateStart
					, Double firstTimeStep
					, Double ZRef, Double ZRStart )
{
	int 	k;
	Double	*t = fTime->vec;
	
	if ( firstTimeStep > 1.0e-13 ) {
		fprintf( fOutFilePtr, "initial timestep is %g s\n", firstTimeStep );
		fFirstTimeStep = firstTimeStep;
#ifndef CVODE
		for ( k = 0; k < fNOfSolver; ++k ) {
			fRWork[k]->vec[kF3] = firstTimeStep;
		}
#endif
	}
	fTStart = fTEnd = timeStart;
	if ( timeStart > 0.0 ) {
		fprintf( fOutFilePtr, "start at time %g s\n", timeStart );
	}
	for ( k = 0; k < fNOfSolver; ++k ) {
		t[k] = timeStart;
	}
	fPressStart = fPressEnd = pressureStart;
	fChiStart = fChiEnd = scalarDissRateStart;

	fZRStart = fZREnd = ZRStart;
	
	fZl = grid[0];
	fZr = grid[gridPointsA-1];

	if ( gridPointsA != fSolGrid->len ) {
		fprintf( stderr, "make grid, nGIn = %d\tnG = %d\n", gridPointsA, fSolGrid->len );
		MakeGrid( &fSolGrid->vec[kPrev], fSolGrid->len, fZl, fZr, fEquidistant );
	}
	else{
		fprintf( stderr, "use inputgrid, nGIn = %d\tnG = %d\n", gridPointsA, fSolGrid->len );
		Double	*newGrid = &fSolGrid->vec[kPrev];
		for ( k = 0; k < gridPointsA; ++k ) {
			newGrid[k] = grid[k];
		}
	}
	SetWeights();

	SetInitial( names, startSolution, grid, gridPointsA, vars );

#ifdef CVODE
	InitCVODE();
#endif
	
	fTempOxEnd = fSolution->mat[kPrev][fTemperature];
	fTempFuelEnd = fSolution->mat[fNGridPoints][fTemperature];

	fprintf( fCAinFile, "%g\t%g\t%g\t%g\t%g\t%g\n", fTStart, fPressStart
				, fTempOxEnd, fTempFuelEnd, fChiStart, fZRStart );
	fflush( fCAinFile );

	if ( ZRef < 0.0 ) {
		fZRef = GetZStoi();
		fprintf( fOutFilePtr, "ZStoi = %g\n", fZRef );
	}
	else {
		fZRef = ZRef;
	}
	if ( fZRef < grid[0] || fZRef > grid[gridPointsA-1] ) {
		fprintf( fOutFilePtr, "###error: ZRef = %g out of bounds (Zmin = %g Zmax = %g)\n"
				, fZRef, grid[0], grid[gridPointsA-1] );
		exit(2);
	}

#ifdef READU
	/* find Zstoi */

	int		i, j;
	Double	**UIn = fUstOverU->mat;
	for ( j = 0; j < fZIn->len; ++j ) {
		if ( fZIn->vec[j] > GetZStoi() ) break;
	}
	// now zstoi is between j-1 and j
	int zloc = j;
	Double	ustoi;

	for ( i = 0; i < fTimeIn->len; ++i ) {
		ustoi = Interpol( GetZStoi(), UIn[i][zloc-1], fZIn->vec[zloc-1], UIn[i][zloc-1], fZIn->vec[zloc] );
		for ( j = 0; j < fZIn->len; ++j ) {
			UIn[i][j] = ustoi / UIn[i][j];
		}
	}
	
	FILE	*fp = GetOutfile( "UstOverU", kData );
	fprintf( fp, "*\nZ" );
	for ( i = 0; i < fTimeIn->len; ++i ) {
		fprintf( fp, "\tt%g", fTimeIn->vec[i] );
	}

	for ( j = 0; j < fZIn->len; ++j ) {
		fprintf( fp, "%g", fZIn->vec[j] );
		for ( i = 0; i < fTimeIn->len; ++i ) {
			fprintf( fp, "\t%g", UIn[i][j] );
		}
		fprintf( fp, "\n" );
	}
	fclose( fp );
#endif

#ifdef EQUILIBRIUM
	int		nSpeciesIn = fSpecies->GetNSpeciesInSystem();
	int		nndd = 7;
	Double	*xmol_i = fSpecies->GetMolarMass()->vec;
	Double	p = GetPressure( GetCurrentTime() ) * 1.0e-5;
	Double	*c0;
	char	*symbol = new char[(nSpeciesIn+4)*20];
	Double	*ponalo = new Double[nSpeciesIn*nndd];
	Double	*ponahi = new Double[nSpeciesIn*nndd];
	Double	totent;

	Double	**lowVals = fSpecies->GetCoeffLow()->mat;
	Double	**highVals = fSpecies->GetCoeffHigh()->mat;
	MatrixPtr	lowMat = FortranToCMat( ponalo, nndd, nndd, nSpeciesIn, nSpeciesIn );
	MatrixPtr	highMat = FortranToCMat( ponahi, nndd, nndd, nSpeciesIn, nSpeciesIn );
	Double	**low = lowMat->mat;
	Double	**high = highMat->mat;
	for ( int i = 0; i < nSpeciesIn; ++i ) {
		copy( nndd, lowVals[i], 1, low[i], 1 );
		copy( nndd, highVals[i], 1, high[i], 1 );
	}

	CToFortranCharArray( symbol, fSpecies->GetNames(), 20, nSpeciesIn );

	for ( k = -1; k <= fNGridPoints; ++k ) {
//		c0 = fSolMassFracs->mat[k];
		c0 = &fSolution->mat[k][fFirstSpecies];
		totent = fTotEntStart->vec[k]*1.0e-3;
		ADIABFLAMETEMP( &fSolGrid->vec[k], &nSpeciesIn, c0, &fSolution->mat[k][fTemperature], &p, symbol, 
                                 xmol_i, &totent,ponalo, ponahi,
                                 &nndd,
                                 &fSolution->mat[k][fTemperature],c0 );
//		fprintf( stderr, "temp[%g] = %g\n", fSolGrid->vec[k], fSolTemp->vec[k] );
		SaveSolution( k, GetCurrentTime(), fSolution->mat[k] );
		SaveSolution( k, GetCurrentTime(), fSolution->mat[k] );
	}

	DisposeFToCMat( highMat );
	DisposeFToCMat( lowMat );
	fprintf( stderr, "write equilibrium data and exit\n" );
	WriteFlameletFile( NULL, NULL, "equi" );
	exit( 2 );
#endif
}

Double TTransFlameSolver::GetTotEnt( Double temp, Double *Y )
{
	int		i;
	int		nSpeciesIn = fSpecies->GetNSpeciesInSystem();
	Double	*h = fSpecies->GetEnthalpy()->vec;
	Double	totEnt;
	
	fSpecies->ComputeSpeciesProperties( temp );
	totEnt = 0.0;
	for ( i = 0; i < nSpeciesIn; ++i ) {
		totEnt += Y[i] * h[i];
	}
	
	return totEnt;
}


/*void TTransFlameSolver::SetTotEnthalpy( void )
{
	int		i, k;
	int		nSpeciesIn = fSpecies->GetNSpeciesInSystem();
	Double	*h = fSpecies->GetEnthalpy()->vec;
	Double	*totEnt = fTotEntStart->vec;
	Double	**YMat = fSolMassFracs->mat;
	Double	*Y;
	
	for ( k = -1; k <= fNGridPoints; ++k ) {
		fSpecies->ComputeSpeciesProperties( fSolTemp->vec[k] );
		Y = YMat[k];
		totEnt[k] = 0.0;
		for ( i = 0; i < nSpeciesIn; ++i ) {
			totEnt[k] += Y[i] * h[i];
		}
	}
	FILE	*fp = fopen( "TotEnt.dout", "w" );
	fprintf( fp, "*\nZ\tTotEnt\n" );
	
	for ( k = -1; k <= fNGridPoints; ++k ) {
		fprintf( fp, "%g\t%g\n", fSolGrid->vec[k], totEnt[k] );
	}

	fclose( fp );
}
*/

Flag TTransFlameSolver::Solve( Double timeEnd, Double pressureEnd, Double scalarDissRateEnd
					, Double temOxEnd, Double tempFuelEnd, int deltaStepsOut )
{
	return Solve( timeEnd, pressureEnd, scalarDissRateEnd, temOxEnd, tempFuelEnd, 1.0, deltaStepsOut );
}

Flag TTransFlameSolver::Solve( Double timeEnd, Double pressureEnd, Double scalarDissRateEnd
					, Double temOxEnd, Double tempFuelEnd, Double ZREnd, int deltaStepsOut )
{
	Flag	error = FALSE;
	Flag	leave = FALSE;
	int		counter = 0;
		
	fprintf( fCAinFile, "%g\t%g\t%g\t%g\t%g\t%g\n", timeEnd, pressureEnd
				, temOxEnd, tempFuelEnd, scalarDissRateEnd, ZREnd );
	fflush( fCAinFile );

	SetEndValues( timeEnd, pressureEnd, scalarDissRateEnd, temOxEnd, tempFuelEnd, ZREnd );

// check dissrate
//	PrintDissRate( 0.0 );

#ifdef LOGNORMALCHI
	SetRandomNumber();
#endif
	
#ifdef FULLIMPLI
/*	if ( fFirstCall ) {*/
/*		leave = FirstImpliStep();*/
/*<<		{*/
/*			char	countstring[128];*/
/*			FILE	*fp;*/
/*			if ( counter > 0 ) {*/
/*				sprintf( countstring, "%d", counter );*/
/*				fp = GetOutputFile( fTime->vec[0], NULL, countstring, TFlame::kText );*/
/*			}*/
/*			else {*/
/*				fp = GetOutputFile( fTime->vec[0], NULL, NULL, TFlame::kText );*/
/*			}*/
/*			fprintf( fOutFilePtr, "dump output\n" );*/
/*			WriteFlameletFile( fp, NULL, NULL );*/
/*			fclose( fp );*/
/*			counter++;*/
/*		}>>*/
/*	}*/
/*	error = leave;*/
/*	if ( fTime->vec[0] + MIN( ( fTEnd - fTStart ) * 1.0e-7, 1.0e-10 ) >= fTEnd ) {*/
/*		leave = TRUE;*/
/*	}*/
/*#	ifdef HP*/
/*		if ( gExit ) {*/
/*			DoExit();*/
/*		}*/
/*#	endif*/
/*	if ( leave ) {*/
/*		return error;*/
/*	}*/
/*	*/
	do {
		leave = OneImpliStep();
		if ( deltaStepsOut > 0 && *fNActualStep[0] % deltaStepsOut == 0 ) {
			char	countstring[128];
			FILE	*fp;
			if ( counter > 0 ) {
				sprintf( countstring, "%d", counter );
				fp = GetOutputFile( fTime->vec[0], NULL, countstring, TFlame::kText );
			}
			else {
				fp = GetOutputFile( fTime->vec[0], NULL, NULL, TFlame::kText );
			}
			fprintf( fOutFilePtr, "dump output\n" );
			WriteFlameletFile( fp, NULL, NULL );
			fclose( fp );

/*			if ( counter > 0 ) {*/
/*				sprintf( countstring, "%d", counter );*/
/*				fp = GetOutputFile( fTime->vec[0], NULL, countstring, TFlame::kData );*/
/*			}*/
/*			else {*/
/*				fp = GetOutputFile( fTime->vec[0], NULL, NULL, TFlame::kData );*/
/*			}*/
/**/
/*			MatrixPtr	theMat = NewMatrix( fNOfEquations, fNGridPoints+2, kColumnPointers );*/
/**/
/*			fSolPrime->mat = &fSolPrime->mat[kPrev];*/
/*			fSolution->mat = &fSolution->mat[kPrev];*/
/**/
/*			copy_mat( theMat, fSolPrime );*/
/*			div_mat( 1, theMat, fSolution );*/
/**/
/*			theMat->mat = &theMat->mat[kNext];*/
/*			fSolPrime->mat = &fSolPrime->mat[kNext];*/
/*			fSolution->mat = &fSolution->mat[kNext];*/
/**/
/*			PrintSolution( fp, fSolGrid->len-2, fNOfEquations, fSolGrid->vec, theMat->mat, fVariableNames );*/
/*			fclose( fp );*/
/*			theMat->mat = &theMat->mat[kPrev];*/
/*			DisposeMatrix( theMat );*/

			counter++;
		}
		error = leave;
		if ( fTime->vec[0] + MIN( ( fTEnd - fTStart ) * 1.0e-7, 1.0e-12 ) >= fTEnd ) {
			leave = TRUE;
			cerr << "dump libout" << NEWL;
			char   tl[128];
			Double         tNow = fSolution->mat[LocationOfMax( fNGridPoints, &fSolution->mat[0][fTemperature], fSolution->phys_rows )][fTemperature];
			sprintf( tl, "Lib_Chi%05g_T%04.0ft%05.2g", GetRefDissRate( fTime->vec[0] ), tNow, fTime->vec[0]*1000.0 );
			FILE *fp = GetOutfile( tl, kText );
			//                FILE *fp = GetOutputFile( fTime->vec[0], tl, NULL, TFlame::kText );                           
			WriteFlameletFile( fp, NULL, NULL );
			fclose( fp );

/*			char	countstring[128];*/
/*			FILE	*fp;*/
/*			if ( counter > 0 ) {*/
/*				sprintf( countstring, "%d", counter );*/
/*				fp = GetOutputFile( fTime->vec[0], NULL, countstring, TFlame::kText );*/
/*			}*/
/*			else {*/
/*				fp = GetOutputFile( fTime->vec[0], NULL, NULL, TFlame::kText );*/
/*			}*/
/*			fprintf( fOutFilePtr, "dump output\n" );*/
/*			WriteFlameletFile( fp, NULL, NULL );*/
/*			fclose( fp );*/
/*			counter++;*/
		}
#	ifdef HP
		if ( gExit ) {
			DoExit();
		}
#	endif
	} while ( !leave );
	
#else
/*	if ( fFirstCall ) {*/
/*		leave = FirstStep();*/
/*	}*/
	do {
		if ( ( fActualPoint = GetActualPoint( timeEnd ) ) >= 0 ) {
			leave = OneStep( fActualPoint );
			error = leave;
#	ifdef HP
			if ( gExit ) {
				DoExit();
			}
#	endif
		}
		else {
			leave = TRUE;
			error = FALSE;
		}
	} while ( !leave );
#endif

	return error;
}

void TTransFlameSolver::SetEndValues( Double timeEnd, Double pressureEnd, Double scalarDissRateEnd
					, Double temOxEnd, Double tempFuelEnd, Double ZREnd )
{
	int	k;
	fprintf( fOutFilePtr, "ZREnd = %g\n", ZREnd );
	if ( fTEnd >= timeEnd ) {
		fprintf( stderr, "###Error: tStart (%g) >= tEnd (%g)\n", fTEnd, timeEnd );
		exit( 2 );
	}
	fTStart = fTEnd;
	fPressStart = fPressEnd;
	fChiStart = fChiEnd;
	fTempOxStart = fTempOxEnd;
	fTempFuelStart = fTempFuelEnd;
	fZRStart = fZREnd;
	
//set ZREnd
	Double sum = 0.0;
	for ( int i = fZRin->len-1; i > 0 ; --i ) {
		fZRin->vec[i] = fZRin->vec[i-1];
		sum += fZRin->vec[i];
	}
//	if ( timeEnd < 4.0e-4 ) {
//		fZRin->vec[0] = 1.0;
//	}
//	else {
		fZRin->vec[0] = ZREnd;
//	}
	fZREnd = MIN( fZRStart, ( sum + fZRin->vec[0] ) / fZRin->len );

	Double	*totEntStart = fTotEntStart->vec;
	Double	*totEntEnd = fTotEntEnd->vec;
	Double	*Z = fSolGrid->vec;
	copy( fNGridPoints+2, &totEntEnd[kPrev], 1, &totEntStart[kPrev], 1 );
		
	fTEnd = timeEnd;
	fPressEnd = pressureEnd;
	fChiEnd = scalarDissRateEnd;

#ifdef LAURENTLIB
	totEntEnd[kPrev] = temOxEnd;
	totEntEnd[fNGridPoints] = tempFuelEnd;
	Double	*Y = fSolMassFracs->mat[kPrev];
	fTempOxEnd = GetTempOfEnthalpy( totEntEnd[kPrev], Y, fTempOxStart );
	Y = fSolMassFracs->mat[fNGridPoints];
	fTempFuelEnd = GetTempOfEnthalpy( totEntEnd[fNGridPoints], Y, fTempFuelStart );
#else
	fTempOxEnd = temOxEnd;
	fTempFuelEnd = tempFuelEnd;
	totEntEnd[kPrev] = GetTotEnt( fTempOxEnd, fSolMassFracs->mat[kPrev] );
	totEntEnd[fNGridPoints] = GetTotEnt( fTempFuelEnd, fSolMassFracs->mat[fNGridPoints] );
#endif
	for ( k = 0; k <= fNGridPoints-1; ++k ) {
		totEntEnd[k] = Interpol( Z[k], totEntEnd[kPrev], Z[kPrev], totEntEnd[fNGridPoints], Z[fNGridPoints] );
	}

/*	FILE	*fp = fopen( "TotEntEnd.dout", "w" );*/
/*	fprintf( fp, "*\nZ\tTotEnt\n" );*/
/*	*/
/*	for ( k = -1; k <= fNGridPoints; ++k ) {*/
/*		fprintf( fp, "%g\t%g\n", fSolGrid->vec[k], totEntEnd[k] );*/
/*	}*/
/**/
/*	fclose( fp );*/
	
	fprintf( fOutFilePtr, "solve flamelet\n" );
	fprintf( fOutFilePtr, "from t = %.4g ms with p = %.4g bar  chi = %.4g 1/s  Tox = %.4g K  Tfu = %.4g K  ZR = %.4g\n"
		, fTStart * 1.0e3, fPressStart / 1.0e5, fChiStart, fTempOxStart, fTempFuelStart, fZRStart );
	fprintf( fOutFilePtr, "  to t = %.4g ms with p = %.4g bar  chi = %.4g 1/s  Tox = %.4g K  Tfu = %.4g K  ZR = %.4g\n"
		, fTEnd * 1.0e3, fPressEnd / 1.0e5, fChiEnd, fTempOxEnd, fTempFuelEnd, fZREnd );

	fDTdtOx = ( fTempOxEnd - fTempOxStart ) / ( fTEnd - fTStart );
	fDTdtFu = ( fTempFuelEnd - fTempFuelStart ) / ( fTEnd - fTStart );

	fdDeltaZdt = ( fZREnd - fZRStart ) / ( fTEnd - fTStart );

#ifdef DPDT
	fDPdt = ( fPressEnd - fPressStart ) / ( fTEnd - fTStart );
#else
	fDPdt = 0.0;
#endif

/*#ifdef MOVEZRIGHT*/
/*	if ( fDPdt != 0.0 || fabs( fDTdtOx ) > 1.0e-10 || fabs( fDTdtFu ) > 1.0e-10 ) {*/
/*		fprintf( stderr, "fDPdt = %g\tfDTdtOx = %g\tfDTdtFu = %g\n", fDPdt, fDTdtOx, fDTdtFu );*/
/*		exit( 2 );*/
/*	}*/
/*#endif*/

	SetMolarMassOverRInf();

#ifdef CVODE
	CVodeSetStopTime(fMem, fTEnd);
#else
	for ( k = 0; k < fNOfSolver; ++k ) {
		if ( fInfo[k]->vec[kF4] == 1 ) {
			fRWork[k]->vec[kF1] = fTEnd;
		}
	}
#endif
}

void TTransFlameSolver::SetMolarMassOverRInf( void )
{
	fProperties->ComputeMixtureMolarMass( fWOverRInf, fSolMassFracs->mat[kPrev]
				, fSpecies->GetMolarMass()->vec, fSpecies->GetNSpeciesInSystem() );
	fWOverRInf /= RGAS;
}

void TTransFlameSolver::PrintDissRate( Double t )
{
	static int	counter = 0;		
	char	name[64];
	sprintf( name, "DissRate%d", counter++ );
	FILE	*fp = GetOutfile( name, kData );

	fprintf( fp, "*\nZ\tChi\n" );
	for ( int k = -1; k <= fNGridPoints; ++k ) 
	{
		fprintf( fp, "%g\t%g\n", fSolGrid->vec[k], GetDissRate( t, fSolGrid->vec[k] ) );
	}
	
	fclose( fp );
}

Double TTransFlameSolver::GetTempOfEnthalpy( Double ent, Double *Y, Double initialGuess )
{
	Double	*h = fSpecies->GetEnthalpy()->vec;
	Double	*cp = fSpecies->GetHeatCapacity()->vec;
	Double	temp = ( initialGuess > 0.0 ) ? initialGuess : 1000;
	Double	deltaT;
	Double	entSum, cpSum;
	int		i;
	int		nOfSpeciesIn = fSpecies->GetNSpeciesInSystem();
	int		count = 0;

	do {
		fSpecies->ComputeSpeciesProperties( temp );
		for ( i = 0, entSum = cpSum = 0; i < nOfSpeciesIn; ++i ) {
			cpSum += cp[i] * Y[i];
			entSum += h[i] * Y[i];
		}
		
		deltaT = -( entSum - ent ) / cpSum;
		temp += deltaT;
		
		count++;
		if ( ++count > 1000 ) {
			fprintf( stderr, "#Error: temperature iteration for h = %g not converged\n"
					, ent );
			exit( 2 );
		}
	} while ( fabs( deltaT / temp ) > 1.0e-3 );

	return temp;
}

void TTransFlameSolver::SetInitial( ConstStringArray names, Double **startSol, Double *grid, int gridPointsA, int vars )
{
	int			i, k, ind;
	int			nOfSpeciesIn = fSpecies->GetNSpeciesInSystem();
	MMDataBag	bag( vars );
	
// copy input to fSolution
	bag.Initialize();
	bag.SetOldInpedVar( grid, gridPointsA, 1, "xIn" );
	for ( i = 0; i < vars; ++i ) {
		bag.Insert( &startSol[0][i], gridPointsA, startSol[1] - startSol[0], names[i] );
	}
	bag.SetNewInpedVar( &fSolGrid->vec[kPrev], fSolGrid->len, 1, "xNew" );

	Double		*newY = &fOutSolWork->vec[kNext];
	Double		**sol = fSolution->mat;
	for ( i = 0; i < bag.NumElems(); ++i ) {
		ind = GetVariableIndex( bag[i].Name() );
		if ( ind >= 0 ) {
			bag[i].Map( fOutSolWork );
			for ( k = -1; k < fSolGrid->len-1; ++k ) {
				sol[k][ind] = newY[k];
			}
			if ( fSoot && ind >= fSootMoments 
						&& ind < fSootMoments + fSoot->GetNSootMoments() ) {
//				momentSet = TRUE;
			}
		}
		else {
			fprintf( fOutFilePtr, "#warning: no match for variable %s\n", bag[i].Name() );
		}
	}

//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
	for ( k = -1; k < fSolGrid->len-1; ++k ) {
        sol[k][fGrid] = fSolGrid->vec[k];
	}
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo

	Double	*totEntStart = fTotEntStart->vec;
	Double	*h = fSpecies->GetEnthalpy()->vec;
#ifdef TEMPFROMENTHALPY
// set Enthalpy
	Double	*solatk;

	// set total enthalpy
		// at k = -1
	totEntStart[kPrev] = 0.0;
	solatk = fSolution->mat[kPrev];
	fSpecies->ComputeSpeciesProperties( solatk[fTemperature] );
	for ( i = 0; i < nOfSpeciesIn; ++i ) {
		totEntStart[kPrev] += solatk[fFirstSpecies+i] * h[i];
	}
		// at k = -1
	totEntStart[fNGridPoints] = 0.0;
	solatk = fSolution->mat[fNGridPoints];
	fSpecies->ComputeSpeciesProperties( solatk[fTemperature] );
	for ( i = 0; i < nOfSpeciesIn; ++i ) {
		totEntStart[fNGridPoints] += solatk[fFirstSpecies+i] * h[i];
	}
		// linearly interpolate enthalpy
	for ( k = 0; k < fNGridPoints; ++k ) {
		solatk = fSolution->mat[k];
		totEntStart[k] = Interpol( fSolGrid->vec[k]
					, totEntStart[kPrev], fSolGrid->vec[kPrev]
					, totEntStart[fNGridPoints], fSolGrid->vec[fNGridPoints] );
		solatk[fTemperature] = GetTempOfEnthalpy( totEntStart[k], &solatk[fFirstSpecies], fSolution->mat[k-1][fTemperature] );
	}
	
#endif

#ifdef DEBUGINITIAL
	FILE	*fp = GetOutfile( "InitialSolOld", kData );
	PrintSolution( fp, fSolGrid->len-2, fNOfEquations, fSolGrid->vec, fSolution->mat, fVariableNames );
	fclose( fp );
#endif

	if ( !fEquidistant ) {
//		SetGrid();
	}

	if (fSoot)
	{
	  // solver works with log(M/rho), incoming moments already in log
	  int nSootMoments = fSoot->GetNSootMoments();
	  int momOff = fSoot->GetOffsetSootMoments();

	  double mixMolarMass;
	  double rho;
	  double minval, ii, jj;

	  double * solk;

	  for (k = -1; k <= fNGridPoints; ++k)
	  {
	    solk = sol[k];
	    fProperties->ComputeMixtureMolarMass( mixMolarMass, &sol[k][fFirstSpecies], fSpecies->GetMolarMass()->vec, nOfSpeciesIn );
	    rho = GetPressure(fTime->vec[0]) * mixMolarMass / (RGAS * solk[fTemperature]);
	    for (i = 0; i < nSootMoments; ++i)
	    {
	      // Make sure the incoming moments make sense
	      ii = double(GetSoot()->Geti(i))/6.0;
	      jj = double(GetSoot()->Getj(i))/6.0;
	      minval = pow(2.0*GetSoot()->nucl_nbrC2, ii+(2.0/3.0)*jj) * 1.0e-20;
	      if (i == fSoot->GetNSootMoments()-1) minval = 1.0e-20;
	      minval = log(minval/rho);
	      solk[momOff+i] = MAX(solk[momOff+i],minval);
	      //solk[momOff+i] /= rho;
	      //cerr << solk[momOff+i] << "\t" << minval << endl;
	      //solk[momOff+i] = log(MAX(solk[momOff+i], minval) / rho);
	    }
	  }
	}
	
#ifdef DEBUGINITIAL
	for ( k = -1; k < fSolGrid->len-1; ++k ) {
        fSolGrid->vec[k] = sol[k][fGrid];
	}
	fp = GetOutfile( "InitialSol", kData );
	PrintSolution( fp, fSolGrid->len-2, fNOfEquations, fSolGrid->vec, fSolution->mat, fVariableNames );
	fclose( fp );
#endif

// save initial solution to fSol
	int		nSpeciesIn = fSpecies->GetNSpeciesInSystem();
	Double	*Y;
	Double	*solVec;
	Double	pressure;
	Double	*totEntEnd = fTotEntEnd->vec;
	for ( k = -1; k <= fNGridPoints; ++k ) {
		solVec = fSolution->mat[k];
		pressure = GetPressure( fTime->vec[0] );
		Y = fSolMassFracs->mat[k];
		// put fSolution to fMassFracs
		SaveSolution( k, fTime->vec[0], fSolution->mat[k] );
		// update steady state concs

		// Don't send soot moments to thermoproperties
#ifdef DELTAINEW
		//UpdateThermoProps(k, Y, fSolTemp->vec[k], pressure, fProperties->GetDensityRef(), kDensFromPress, (fSoot)?&solVec[fSoot->GetOffsetSootMoments()]:NULL);
		UpdateThermoProps(k, Y, fSolTemp->vec[k], pressure, fProperties->GetDensityRef(), kDensFromPress, NULL);
#else 
		T0DFlame::UpdateThermoProps( Y, fSolTemp->vec[k], pressure
					, fProperties->GetDensityRef(), kDensFromPress
					, (fSoot)?&solVec[fSoot->GetOffsetSootMoments()]:NULL );
#endif
		fDensity->vec[k] = fProperties->GetDensity();
		// second call puts solution including steady state concs to  down to fMassFracs
		SaveSolution( k, fTime->vec[0], fSolution->mat[k] );
		// third call puts solution down to fSolOld
		SaveSolution( k, fTime->vec[0], fSolution->mat[k] );
		// set total enthalpy
		totEntStart[k] = 0.0;
		for ( i = 0; i < nSpeciesIn; ++i ) {
			totEntStart[k] += Y[i] * h[i];
		}
		totEntEnd[k] = totEntStart[k];
	}

	SetMaxVals();
	
/*	FILE	*fpEnt = fopen( "TotEnt.dout", "w" );*/
/*	fprintf( fpEnt, "*\nZ\tTotEnt\n" );*/
/*	*/
/*	for ( k = -1; k <= fNGridPoints; ++k )*/
/*	{*/
/*		fprintf( fpEnt, "%g\t%g\n", fSolGrid->vec[k], totEntEnd[k] );*/
/*	}*/
/**/
/*	fclose( fpEnt );*/
}

void TTransFlameSolver::SetGrid( void )
{
	int		i, k, countF = 0;
	Double	**sol = fSolution->mat;
	Double	*F = New1DArray( fNGridPoints+2 );
	Double	**solOld = New2DArray( fNGridPoints+2, fNOfEquations );
	Double	*Z = fSolGrid->vec;
	Double	deltaF, FCurr;

	F = &F[kNext];
	solOld = &solOld[kNext];

	for ( k = -1; k <= fNGridPoints; ++k ) {
		for ( i = 0; i < fNOfEquations; ++i ) {
			solOld[k][i] = sol[k][i];
		}
	}

	SetFDWeights( (Double *) *sol );
	ComputeGPDistributionWeights( (Double *) *sol );

// estimate weight fuction at Z=0
	fMonFct[kPrev] = fMonFct[0] - ( fMonFct[1] - fMonFct[0] ) / ( Z[1] - Z[0] ) * ( Z[0] - Z[kPrev] );


	TrapIntegrate( fNGridPoints+2, &fMonFct[kPrev], &Z[kPrev], &F[kPrev] );
	deltaF = F[fNGridPoints] / ( fNGridPoints + 1 );
	
	for ( k = 0; k < fNGridPoints; ++k ) {
		FCurr = ( k + 1 ) * deltaF;
		while ( countF < fNGridPoints && F[countF] < FCurr ) ++countF;
		for ( i = 0; i < fNOfEquations; ++i ) {
			sol[k][i] = Interpol( FCurr, solOld[countF-1][i], F[countF-1], solOld[countF][i], F[countF] );
		}
	}
//	sol[fNGridPoints][fGrid] = 1.0;

#ifdef DEBUGINITIAL
	FILE	*fp = GetOutfile( "F", kData );
	fprintf( fp, "*\nZ\tG\tF\tTOld\tZNew\tTNew\n" );
	for ( k = -1; k <= fNGridPoints; ++k )
	{
		fprintf( fp, "%g\t%g\t%g\t%g\t%g\t%g\n", fSolGrid->vec[k], fMonFct[k], F[k], solOld[k][fTemperature]
								, sol[k][fGrid], sol[k][fTemperature] );
	}
	fclose( fp );
#endif


	Free2DArray(&solOld[kPrev]);
	Free1DArray(&F[kPrev]);
}

void TTransFlameSolver::GetSolution( ConstStringArray names, Double **outSol, Double *grid, int gridPointsA, int vars
									, Double *density )
{
	int			i, k, ind;
	int			nOfSpecies = fSpecies->GetNOfSpecies();
	int			nSootMoments;
	MMDataBag	bag( nOfSpecies+fVariablesWithoutSpecies+( ( density ) ? 1 : 0 ) );
	Double		**theY = fMassFracsWork->mat;
	Double		**theMom;
	int			theYOff = fMassFracsWork->rows;
	int			theMomOff, momOff;
	Double		*theTemp = fTempWork->vec;
	double * theProg = fProgWork->vec;
	double * theEnth = fEnthWork->vec;
	char		**specNames = fSpecies->GetNames();
	VectorPtr	fOutSol = NewVector( gridPointsA );
	VectorPtr	ZGridVec = NewVector( fSolGrid->len );
	Double		*ZGrid = &ZGridVec->vec[kNext];
	Flag		*outSet = new Flag[vars];
	Flag		densitySet = FALSE;
	for ( i = 0; i < vars; ++i ) {
		outSet[i] = FALSE;
	}

	if ( fSoot ) {
		nSootMoments = fSoot->GetNSootMoments();
		momOff = fSoot->GetOffsetSootMoments();
		theMomOff = fSootMomentsWork->rows;
		theMom = fSootMomentsWork->mat;
	}
	
	Double	ZR = Interpol( GetCurrentTime(), fZRStart, fTStart, fZREnd, fTEnd );
	for ( i = -1; i < fSolGrid->len-1; ++i ) {
		ZGrid[i] = fSolGrid->vec[i] * ZR;
	}
	
// copy solution to bag
	bag.Initialize();
	bag.SetOldInpedVar( &ZGrid[kPrev], ZGridVec->len, 1, "xIn" );

	SetOutSolution();
	if ( density ) {
		int		nSpeciesIn = fSpecies->GetNSpeciesInSystem();
		Double	MM;
		Double	*rho = fDensity->vec;
		for ( k = -1; k <= fNGridPoints; ++k ) {
			fProperties->ComputeMixtureMolarMass( MM, theY[k], fSpecies->GetMolarMass()->vec
						, nSpeciesIn );
			rho[k] = GetPressure() * MM / ( RGAS * theTemp[k] );
		}
	}
	bag.Insert( &theTemp[-1], fNGridPoints+2, 1, fVariableNames[fTemperature] );
	bag.Insert( &theProg[-1], fNGridPoints+2, 1, fVariableNames[fProg] );
	bag.Insert( &theEnth[-1], fNGridPoints+2, 1, fVariableNames[fEnth] );
	for ( i = 0; i < nOfSpecies; ++i ) {
		bag.Insert( &theY[-1][i], fNGridPoints+2, theYOff, specNames[i] );
	}
	if ( fSoot ) {
		for ( i = 0; i < nSootMoments; ++i ) {
			bag.Insert( &theMom[-1][i], fNGridPoints+2, theMomOff, fVariableNames[momOff+i] );
		}
	}
	if ( density ) {
		bag.Insert( &fDensity->vec[-1], fNGridPoints+2, 1, "density" );
	}
	
	bag.SetNewInpedVar( grid, gridPointsA, 1, "xNew" );

#ifdef DEBUGBAG
	cout << bag;
#endif

	Double		**newSol = &outSol[kNext];
	Double		*outWork = &fOutSol->vec[kNext];
	for ( i = 0; i < bag.NumElems(); ++i ) {
		ind = GetVariableIndex( bag[i].Name(), names, vars );
		if ( ind >= 0 ) {
			bag[i].Map( fOutSol );
			for ( k = -1; k < gridPointsA-1; ++k ) {
				newSol[k][ind] = outWork[k];
			}
			outSet[ind] = TRUE;
		}
		else {
			if ( density && strcmp( bag[i].Name(), "density") == 0 ) {
				bag[i].Map( fOutSol );
				for ( k = -1; k < gridPointsA-1; ++k ) {
					density[k] = outWork[k];
				}
				densitySet = TRUE;
			}
//			cerr << "#warning: no match for variable " << bag[i].Name() << NEWL;
		}
	}
	
	for ( i = 0; i < vars; ++i ) {
		if ( outSet[i] == FALSE ) {
			fprintf( fOutFilePtr, "#warning from function 'GetSolution': no match for variable '%s'\n", names[i] );
//			cerr << "#warning from function 'GetSolution': no match for variable " 
//				<< names[i] << NEWL;
		}
	}
	if ( density && densitySet == FALSE ) {
		fprintf( fOutFilePtr, "#warning from function 'GetSolution': no match for variable '%s'\n", "density" );
	}
	
#ifdef DEBUGGETSOL
	FILE	*fp = GetOutfile( "NewSol", kData );
	PrintSolution( fp, gridPointsA-2, vars, &grid[kNext], newSol, fVariableNames );
	fclose( fp );
#endif
	
#ifdef DEBUGINITIALDENS
	if ( density ) {
		FILE	*fpdens = GetOutfile( "DensSol", kData );
		
		fprintf( fpdens, "*\nZ\tdensity\n" );
		for ( k = -1; k < gridPointsA-1; ++k ) {
			fprintf( fpdens, "%g\t%g\n", grid[k], density[k] );
		}
	
		fclose( fpdens );
	}
#endif

	// clean up
	delete outSet;
	DisposeVector( ZGridVec );
	DisposeVector( fOutSol );
}

void TTransFlameSolver::GetSolutionInfo( Double *timeStep )
{
  //*timeStep = fRWork[0]->vec[kF7];
  *timeStep = *fActualTimeStepSize;
}

void TTransFlameSolver::SetOutSolution( void )
{
	int			i, k;
	int			nOfSpecies = fSpecies->GetNOfSpecies();
	int			nSootMoments;
	Double		**theY = fMassFracsWork->mat;
	Double		*theTemp = fTempWork->vec;
	double * theProg = fProgWork->vec;
	double * theEnth = fEnthWork->vec;
	Double		*oTemp = fSolOldTemp->vec;
	Double		*nTemp = fSolTemp->vec;
	double * oProg = fSolOldProg->vec;
	double * nProg = fSolProg->vec;
	double * oEnth = fSolOldEnth->vec;
	double * nEnth = fSolEnth->vec;
	Double		**oY = fSolOldMassFracs->mat;
	Double		**nY = fSolMassFracs->mat;
	Double		*oTime = fSolOldTime->vec;
	Double		*nTime = fSolTime->vec;
	Double		currTime = GetCurrentTime();
	Double		**nMom;
	Double		**oMom;
	Double		**theMom;
	Double		oTimek, nTimek, *theYk, *oYk, *nYk;
	if ( fSoot ) {
		nSootMoments = fSoot->GetNSootMoments();
		nMom = fSolSootMoments->mat;
		oMom = fSolOldSootMoments->mat;
		theMom = fSootMomentsWork->mat;
	}

	for ( k = -1; k <= fNGridPoints; ++k )
	{
		oTimek = oTime[k];
		nTimek = nTime[k];
		theYk = theY[k];
		oYk = oY[k];
		nYk = nY[k];
		theTemp[k] = Interpol( currTime, oTemp[k], oTimek, nTemp[k], nTimek );
		theProg[k] = Interpol( currTime, oProg[k], oTimek, nProg[k], nTimek );
		theEnth[k] = Interpol( currTime, oEnth[k], oTimek, nEnth[k], nTimek );
		for ( i = 0; i < nOfSpecies; ++i ) {
			theYk[i] = Interpol( currTime, oYk[i], oTimek, nYk[i], nTimek );
		}

		if ( fSoot ) {
			for ( i = 0; i < nSootMoments; ++i ) {
			  theMom[k][i] = Interpol( currTime, oMom[k][i], oTimek, nMom[k][i], nTimek );
			}
		}
	}
}

void TTransFlameSolver::PrintSolution( FILE *fp, int nGPoints, int nEq, Double *grid
									, Double **sol, ConstStringArray names )
{
	int		i, k;
	char	format[128];
	int 	*len = new int[nEq];
	
	fprintf( fp, "*\n%-12s", "Z" );
	
	for ( i = 0; i < nEq; ++i ) {
		len[i] = maxint( strlen( names[i] ), 12 );
		sprintf( format, "%s%%-%ds", "\t", len[i] );
		fprintf( fp, format, names[i] );
	}
	
	for ( k = -1; k < nGPoints+1; ++k )
	{
		fprintf( fp, "\n%-12E", grid[k] );
		for ( i = 0; i < nEq; ++i ) {
			sprintf( format, "%s%%-%dE", "\t", len[i] );
            fprintf( fp, format, sol[k][i] );
		}
	}
	
	delete len;
}

void MAKEGRID( int **object, Double *grid, int *gridPoints
							, Double *left, Double *right, int *equidistant )
{
	makegrid( object, grid, gridPoints, left, right, equidistant );
}

void makegrid( int **object, Double *grid, int *gridPoints
							, Double *left, Double *right, int *equidistant )
{
	TTransFlameSolverPtr flame = *( TTransFlameSolverPtr *)object;

	flame->MakeGrid( grid, *gridPoints, *left, *right, ( *equidistant ) ? TRUE : FALSE );
}

void TTransFlameSolver::MakeGrid( VectorPtr theGrid, Double left, Double right, Flag equidistant )
{
	MakeGrid( theGrid->vec, theGrid->len, left, right, equidistant );
}

void TTransFlameSolver::MakeGrid( Double *grid, int len
							, Double /*left*/, Double /*right*/, Flag equidistant )
{
	Double	deltaZFine;
	
	if ( equidistant ) {
		deltaZFine = 0.0;
	}
	else {
	  deltaZFine = 0.2; //0.2
	}
	
	int			 	gridPoints = len - 2;
	Double			deltaZCoarse = 1.0 - deltaZFine;
	Double	deltaZ = ( deltaZCoarse + 2.0 * deltaZFine ) 
							/ ( ( Double ) ( gridPoints + 1 ) );
	Double	bound = deltaZFine - 0.5 * deltaZ + 1.0e-10;

	grid[0] = 0.0;
	grid[gridPoints+1] = 1.0;

	for ( int k = 1; k <= gridPoints; ++k ) {
		if ( grid[k-1] <= bound ) {
			grid[k] = grid[k-1] + 0.5 * deltaZ;
		}
		else {
			grid[k] = grid[k-1] + deltaZ;
		}
		if ( grid[k] >= 1.0 ) {
			fprintf( fOutFilePtr, "#error: invalid grid generated: Z[%d] = %g\n", k, grid[k] );
//			cerr << "#error: invalid grid generated: Z["<<k<<"] = " << grid[k] << NEWL;
		}
	}
#	ifdef DEBUGGRID
	static int	gridCount = 0;
	char		fName[128];
	sprintf( fName, "grid%2d", gridCount++ );
	FILE		*fp = GetOutfile( fName, kData );
	
	fprintf( fp, "*\nx\td\n" );
	for ( k = 0; k < gridPoints+2; ++k ) {
		fprintf( fp, "%g\t1.0\n", grid[k] );
	}
	
	fclose( fp );
#	endif
}

void TTransFlameSolver::SetWeights( void )
{
	Double	*grid = fSolGrid->vec;
	Double	h, hm, hnenn;
	
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
	fgpDens->vec[kPrev] = 1.0 / ( grid[fGrid] - fZl );
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
	for ( int k = 0; k < fNGridPoints; ++k )
	{
		h = grid[k+1] - grid[k];
		hm = grid[k] - grid[k-1];
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
		fh[k] = grid[k+1] - grid[k];
		fhm[k] = grid[k] - grid[k-1]; 
		fgpDens->vec[k] = 1.0 / h;
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
		hnenn = h * hm * ( h + hm );
		fFDWCurr[k] = ( h * h - hm * hm ) / hnenn;
		fFDWMinus[k] = - h * h / hnenn;
		fFDWPlus[k] = hm * hm / hnenn;
		fWCurr[k] = - 2.0 * ( h + hm ) / hnenn;
		fWMinus[k] = 2.0 * h / hnenn;
		fWPlus[k] = 2.0 * hm / hnenn;
	}
}

//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
void TTransFlameSolver::SetFDWeights( Double *y )
{
	int		gpOff;
	double  fhnenn;

	fgpDens->vec[kPrev] = 1.0 / ( y[fGrid] - y[fGrid-fNOfEquations] );
	for ( int k = 0; k < fNGridPoints; ++k ) {
		gpOff = k * fNOfEquations;
		fh[k] = y[gpOff+fGrid+fNOfEquations] - y[gpOff+fGrid];
		fhm[k] = y[gpOff+fGrid] - y[gpOff+fGrid-fNOfEquations];
		fh[kPrev] = fhm[kCurr];
		fhm[fNGridPoints] = fh[fNGridPoints-1];
		fhnenn = fh[k] * fhm[k] * ( fh[k] + fhm[k] );
		fFDWCurr[k] = ( fh[k] * fh[k] - fhm[k] * fhm[k] ) / fhnenn;
		fFDWMinus[k] = - fh[k] * fh[k] / fhnenn;
		fFDWPlus[k] = fhm[k] * fhm[k] / fhnenn;
		fWCurr[k] = - 2.0 * ( fh[k] + fhm[k] ) / fhnenn;
		fWMinus[k] = 2.0 * fh[k] / fhnenn;
		fWPlus[k] = 2.0 * fhm[k] / fhnenn;	
		fgpDens->vec[k] = 1.0 / fh[k];
	}
}
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo

//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
void TTransFlameSolver::ComputeGPDistributionWeights( Double *y )
{
	int			nSpeciesIn = fSpecies->GetNSpeciesInSystem();
	Double		dTdZ, dYdZ; 
//	Double		d2YdZ2, dTdZsq; 
	Double		*maxVals = fMaxVals->vec;
	int			j, k;
	int			gpOff, sootOff, nSootMoments;
	if ( fSoot ) {
		sootOff = fSoot->GetOffsetSootMoments();
		nSootMoments = fSoot->GetNSootMoments();
	}
	
	fMonFct[kPrev] = fMonFct[fNGridPoints] = 1.0;
	for ( k = 0; k < fNGridPoints; ++k ) {
		gpOff = k * fNOfEquations;
		fMonFct[k] = 1.0;

		dTdZ = ( y[gpOff+fTemperature+fNOfEquations] - y[gpOff+fTemperature] ) / fh[k];
/*		dTdZsq = dTdZ * dTdZ;*/
/*		if ( dTdZsq > maxdTdZsq ) {*/
/*			maxdTdZsq = dTdZsq;*/
/*		}*/
		fMonFct[k] += dTdZ * dTdZ / MAX( 1.0e-10, maxVals[fTemperature] * maxVals[fTemperature] );
		for ( j = 0; j < nSpeciesIn; ++j ) {
			dYdZ = ( y[gpOff+fNOfEquations+fFirstSpecies+j] - y[gpOff+fFirstSpecies+j] ) / fh[k];
			fMonFct[k] += dYdZ * dYdZ / MAX( 1.0e-10, maxVals[fFirstSpecies+j] * maxVals[fFirstSpecies+j] );
		}
		if ( fSoot ) {
			for ( j = 0; j < nSootMoments; ++j ) {
				dYdZ = ( y[gpOff+fNOfEquations+sootOff+j] - y[gpOff+sootOff+j] ) / fh[k];
					fMonFct[k] += dYdZ * dYdZ / MAX( 1.0e-10, maxVals[sootOff+j] * maxVals[sootOff+j] );
			}
		}
		fMonFct[k] = sqrt( fMonFct[k] );
		
		// smooth monitor function
/*		if ( k != 0 ) {*/
/*			fMonFct[k-1] = ( ValMinMin + ValMin + fMonFct[k] ) / 3.0;*/
/*		}*/
/*		else {*/
/*			fMonFct[k-1] = fMonFct[k];*/
/*			ValMin = fMonFct[k-1];*/
/*		}*/
/*		ValMinMin = ValMin;*/
/*		ValMin = fMonFct[k];*/
	}

/*#ifdef COMPMONFUNC*/
/*		d2YdZ2 =  fWMinus[k] * y[gpOff+fTemperature-fNOfEquations]*/
/*				+ fWCurr[k] * y[gpOff+fTemperature]*/
/*				+ fWPlus[k] * y[gpOff+fTemperature+fNOfEquations]/ */
/*		fMonFct[k] += d2YdZ2 * d2YdZ2 / MAX( 1.0e-10, dTdZ * dTdZ );*/
/*#endif*/
	
	fMonFct[fNGridPoints] = fMonFct[fNGridPoints-1];
}
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo

int TTransFlameSolver::GetActualPoint( Double tEnd )
{
	int		minPoint = 0;
	Double	*t = fSolTime->vec;
	Double	minValue = t[0];
	
	for ( int k = 1; k < fNGridPoints; ++k )
	{
		if ( t[k] < minValue ) {
			minPoint = k;
			minValue = t[k];
		}
	}
	
	if ( minValue < tEnd ) {
		return minPoint;
	}
	else {
		return -1;
	}
}

void TTransFlameSolver::SaveSolution( int k, Double t, Double *y )
{
	int		nOfSpecies = fSpecies->GetNOfSpecies();
	int		nOfSpeciesIn = fSpecies->GetNSpeciesInSystem();
	int		nSootMoments, momOff;
	Double	*YNew = fSolMassFracs->mat[k];
	Double	*YOld = fSolOldMassFracs->mat[k];
	Double	*momNew;
	Double	*momOld;
	
	if ( fSoot ) {
		nSootMoments = fSoot->GetNSootMoments();
		momOff = fSoot->GetOffsetSootMoments();
		momNew = fSolSootMoments->mat[k];
		momOld = fSolOldSootMoments->mat[k];
	}
	
	// save old solution
	fSolOldTime->vec[k] = fSolTime->vec[k];
	copy( nOfSpecies, YNew, 1, YOld, 1 );
	if ( fSoot ) {
		copy( nSootMoments, momNew, 1, momOld, 1 );
	}
	fSolOldTemp->vec[k] = fSolTemp->vec[k];
	fSolOldProg->vec[k] = fSolProg->vec[k];
	fSolOldEnth->vec[k] = fSolEnth->vec[k];

	// save new solution
	fSolTime->vec[k] = t;
	copy( nOfSpeciesIn, &y[fFirstSpecies], 1, YNew, 1 );

	if ( fSoot ) {
	  double mixMolarMass, rho;
	  fProperties->ComputeMixtureMolarMass( mixMolarMass, YNew, fSpecies->GetMolarMass()->vec, nOfSpeciesIn );
	  rho = GetPressure(t) * mixMolarMass / (RGAS * y[fTemperature]);
		
	  for ( int i = 0; i < nSootMoments; ++i ) {
	    momNew[i] = exp(y[momOff+i]) * rho;
		}
	}
	fSolTemp->vec[k] = y[fTemperature];
	fSolProg->vec[k] = y[fProg];
	fSolEnth->vec[k] = y[fEnth];
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
	fSolGrid->vec[k] = y[fGrid];
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
}

void TTransFlameSolver::SetMaxVals( void )
{
	int		nOfSpeciesIn = fSpecies->GetNSpeciesInSystem();
	int		nSootMoments, momOff;
	Double	**Y = fSolMassFracs->mat;
	Double	**mom;
	Double	*maxvals = fMaxVals->vec;
	
	if ( fSoot ) {
		nSootMoments = fSoot->GetNSootMoments();
		momOff = fSoot->GetOffsetSootMoments();
		mom = fSolSootMoments->mat;
	}
	
	// save new solution
	for ( int i = 0; i < nOfSpeciesIn; ++i ) {
		maxvals[fFirstSpecies+i] = Y[LocationOfMax( fNGridPoints+2, &Y[kPrev][i]
									, fSolMassFracs->phys_rows )-1][i];
	}
	
	if ( fSoot ) {		
		int loc;
		for ( int i = 0; i < nSootMoments; ++i ) {
			loc = LocationOfMax( fNGridPoints+2, &mom[kPrev][i]
									, fSolSootMoments->phys_rows )-1;
			maxvals[momOff+i] = mom[loc][i] / fDensity->vec[loc];
		}
	}

	maxvals[fTemperature] = fSolTemp->vec[LocationOfMax( fNGridPoints+2, &fSolTemp->vec[kPrev], 1 )-1];
	maxvals[fProg] = fSolProg->vec[LocationOfMax(fNGridPoints+2,&fSolProg->vec[kPrev],1)-1];
	maxvals[fEnth] = fSolEnth->vec[LocationOfMax(fNGridPoints+2,&fSolEnth->vec[kPrev],1)-1];
/*	fprintf( fOutFilePtr, "MaxTemp = %g\tOH_Max = %g\tMaxWeight = %g\n", maxvals[fTemperature]*/
/*						, maxvals[GetVariableIndex( "OH" )], fMonFct[LocationOfMax( fNGridPoints, fMonFct*/
/*							, 1 )] );*/
/*	fprintf( fOutFilePtr, "Z[%d]: %.2f\tZ[%d]: %.2f\n"*/
/*			, (int)((fNGridPoints+2)/4), fSolGrid->vec[(int)((fNGridPoints+2)/4)]*/
/*			, (int)((fNGridPoints+2)/2), fSolGrid->vec[(int)((fNGridPoints+2)/2)]  );*/
}

Flag TTransFlameSolver::FirstStep( void )
{
	Flag	leave;
	
	for ( int k = 0; k < fNGridPoints; ++k ) {
		fActualPoint = k;
		if ( leave = OneStep( k ) ) {
			break;
		}
		fInfo[k]->vec[kF11] = 0;
	}

	fFirstCall = FALSE;
	return leave;
}

/*#ifdef FULLIMPLI*/
/*Flag TTransFlameSolver::FirstImpliStep( void )*/
/*{*/
/*	Flag	leave;*/
/*	*/
/*	leave = OneImpliStep();*/
/*	fprintf( stderr, "info[11] = %d\n", fInfo[0]->vec[kF11] );*/
/*	fInfo[0]->vec[kF11] = 0;*/
/**/
/*	fFirstCall = FALSE;*/
/*	return leave;*/
/*}*/
/*#endif*/

Flag TTransFlameSolver::OneStep( int k )
{
	SetMaxTimeStep( k, fTime->vec[k] );
#ifdef UNICOS
	SDASSL( ::ResTransFlameSolver, &fDasslNEq, &fTime->vec[k], fSolution->mat[k]
			, fSolPrime->mat[k], &fTEnd, fInfo[k]->vec, fRTol, fATol, &fIdid
			, fRWork[k]->vec, &fLRW, fIWork[k]->vec, &fLIW, NULL, ( int * ) this
			, ::JacTransFlameSolver );
#else
	DDASSL( ::ResTransFlameSolver, &fDasslNEq, &fTime->vec[k], fSolution->mat[k]
			, fSolPrime->mat[k], &fTEnd, fInfo[k]->vec, fRTol, fATol, &fIdid
			, fRWork[k]->vec, &fLRW, fIWork[k]->vec, &fLIW, NULL, ( int * ) this
			, ::JacTransFlameSolver );
#endif

//		fprintf( fOutFilePtr, "gp = %3d    max stp = %8.2e    stp = %5d    ord = %1d    time = %.4g ms\ttemp = %6.1f K\n"
//			, k, fMaxStepTaken[k], *fNActualStep[k], *fNActualOrd[k], fTime->vec[k] * 1000.0, fSolution->mat[k][fTemperature] );
	fMaxStepTaken[k] = MAX( fMaxStepTaken[k], fRWork[k]->vec[kF7] );
	if ( *fNActualStep[k] % DELTAPROG == 0 ) {
		fprintf( fOutFilePtr, "gp = %3d    max stp = %8.2e    stp = %5d    ord = %1d    time = %.4g ms\ttemp = %6.1f K\n"
			, k, fMaxStepTaken[k], *fNActualStep[k], *fNActualOrd[k], fTime->vec[k] * 1000.0, fSolution->mat[k][fTemperature] );
		fflush( fOutFilePtr );
	}
//	cerr << "max stepsize at gridpoint " << k << " is " << fMaxStepTaken[k] << NEWL;

	if ( fIdid < 1 ) {
		fprintf( fOutFilePtr, "#error: ddassl error no. %d occured at gridpoint no. %d\n"
			, fIdid, k );
//		cerr << "#error: ddassl error no. " << fIdid << " occured at gridpoint no. "
//			<< k << NEWL;
		return TRUE;
	}
	else {
		if ( fRWork[k]->vec[kF3] > fTime->vec[k] ) {
			fprintf( fOutFilePtr, "#warning: integration carried out beyond tEnd\n" );
			fprintf( fOutFilePtr, "info[4] = %d\n", fInfo[k]->vec[kF4] );
			fprintf( fOutFilePtr, "fIdid = %d\n", fIdid );
//			cerr << "#warning: integration carried out beyond tEnd" << NEWL;
//			cerr << "info[4] = " << fInfo[k]->vec[kF4] << NEWL;
//			cerr << "fIdid = " << fIdid << NEWL;
		}
		SaveSolution( k, fTime->vec[k], fSolution->mat[k] );
		SetMaxVals();
		return FALSE;
	}
}

#ifdef FULLIMPLI
Flag TTransFlameSolver::OneImpliStep( void )
{
  static int   initrandomfile=0;
  static Double timrandomchange=0.0;
  const Double deltarandomchange = 1.0e-4;
  static FILE *chiout;

#ifdef ENTHALPYLIB
  static Double  tlastout = fSolution->mat[LocationOfMax(fNGridPoints, &fSolution->mat[0][fTemperature], fSolution->phys_rows)][fTemperature];//10000.0; //
#endif
  
    if (initrandomfile==0) {
        ++initrandomfile;
        chiout=GetOutfile( "RandomChi", kData );
        fprintf(chiout,"*\nChi\trand\tT\n");
	}

#ifdef CVODE
    int flag;
    long int nst;

#ifdef SUNDIALS23
    flag = CVode(fMem, fTEnd, fCVY, &fTime->vec[0], CV_ONE_STEP_TSTOP);
#else
    flag = CVode(fMem, fTEnd, fCVY, &fTime->vec[0], CV_ONE_STEP);
#endif
    SolutionFromCVode();

    CVodeGetNumSteps(fMem, &nst);
    fNActualStepCV = nst;
    CVodeGetLastStep(fMem, &fActualTimeStepSizeCV);
    CVodeGetLastOrder(fMem, &fNActualOrdCV);

	if ( *fNActualStep[0] % DELTAPROG == 0 ) {
		fprintf( fOutFilePtr, "stp = %5d    dt = %8.2e s    ord = %1d    time = %.4g ms\ttemp = %6.1f K\n"
			, *fNActualStep[0], *fActualTimeStepSize, *fNActualOrd[0], fTime->vec[0] * 1000.0
			, fSolution->mat[LocationOfMax( fNGridPoints, &fSolution->mat[0][fTemperature], fSolution->phys_rows )][fTemperature] );
		fflush( fOutFilePtr );
	}

	if (flag < 0) {
		cerr << "#error: cvode error no. " << flag << " occured" << NEWL;
		return TRUE;
	}
#else
#ifdef UNICOS
	SDASSL( ::ResTransFlameImpliSolver, &fDasslNEq, &fTime->vec[0], fSolution->mat[kPrev]
			, fSolPrime->mat[kPrev], &fTEnd, fInfo[0]->vec, fRTol, fATol, &fIdid
			, fRWork[0]->vec, &fLRW, fIWork[0]->vec, &fLIW, NULL, ( int * ) this
			, NULL );
#else
	DDASSL( ::ResTransFlameImpliSolver, &fDasslNEq, &fTime->vec[0], fSolution->mat[kPrev]
			, fSolPrime->mat[kPrev], &fTEnd, fInfo[0]->vec, fRTol, fATol, &fIdid
			, fRWork[0]->vec, &fLRW, fIWork[0]->vec, &fLIW, NULL, ( int * ) this
			, NULL );
#endif
	fMaxStepTaken[0] = MAX( fMaxStepTaken[0], fRWork[0]->vec[kF7] );
	if ( *fNActualStep[0] % DELTAPROG == 0 ) {
		fprintf( fOutFilePtr, "stp = %5d    dt = %8.2e s    ord = %1d    time = %.4g ms\ttemp = %6.1f K\n"
			, *fNActualStep[0], *fActualTimeStepSize, *fNActualOrd[0], fTime->vec[0] * 1000.0
			, fSolution->mat[LocationOfMax( fNGridPoints, &fSolution->mat[0][fTemperature], fSolution->phys_rows )][fTemperature] );
		fflush( fOutFilePtr );
	}

	if ( fIdid < 1 ) {
		fprintf( fOutFilePtr, "#error: dassl error no. %d occured at gridpoint no. %d\n"
			, fIdid, *fNActualStep[0] );
		fprintf( fOutFilePtr, "number of error test failures is %d\n", fIWork[0]->vec[kF14] );
		fprintf( fOutFilePtr, "number of convergence test failures is %d\n", fIWork[0]->vec[kF15] );
		fprintf( fOutFilePtr, "dassl suggested order %d and stepsize %g for the next step\n"
			,fIWork[0]->vec[kF7], fRWork[0]->vec[kF3] );
		return TRUE;
	}
	if ( fRWork[0]->vec[kF3] > fTime->vec[0] ) {
	  fprintf( fOutFilePtr, "#warning: integration carried out beyond tEnd\n" );
	  fprintf( fOutFilePtr, "info[4] = %d\n", fInfo[0]->vec[kF4] );
	  fprintf( fOutFilePtr, "fIdid = %d\n", fIdid );
	}
#endif


#if defined (applec) || defined (powerc)
	RotateCursor( 32 );
#endif


#ifdef LOGNORMALCHI
//	if ( ( int )( fTime->vec[0] / 1.0e-5 ) != ( int )( fTime->vec[0] - fRWork[0]->vec[kF7] / 5.0e-6 ) ) {
    if ( fTime->vec[0] > timrandomchange+deltarandomchange ) {
		fprintf( chiout, "%g\t%g\t%g\n", GetRefDissRate( fTime->vec[0] ), fRandomNumber, fSolution->mat[LocationOfMax( fNGridPoints, &fSolution->mat[0][fTemperature], fSolution->phys_rows )][fTemperature] );
        fflush( chiout );
		SetRandomNumber();
        timrandomchange=fTime->vec[0];
		fprintf( fOutFilePtr, "new random number: %g\n", fRandomNumber );
	}
#endif

	for ( int k = -1; k <= fNGridPoints; ++k ) {
	  SaveSolution( k, fTime->vec[0], fSolution->mat[k] );
	}

	SetMaxVals();
	
	if ( fSoot ) {
	  int maxloc0 = 0;
	  int maxloc1 = 0;
	  for ( int k = 1; k < fNGridPoints; ++k ) {
	    if ( fSolSootMoments->mat[k][0] > fSolSootMoments->mat[maxloc0][0] ) {
	      maxloc0 = k;
	    }
	    if ( fSolSootMoments->mat[k][1] > fSolSootMoments->mat[maxloc1][1] ) {
	      maxloc1 = k;
	    }
	  }
	  fprintf( fOutFilePtr, "M0_max at Z=%g is %g\tM1_max at Z=%g is %g\tfv_max = %g\n\n"
		   , fSolGrid->vec[maxloc0], fSolSootMoments->mat[maxloc0][0]
		   , fSolGrid->vec[maxloc1], fSolSootMoments->mat[maxloc1][1]
		   , fSolSootMoments->mat[maxloc1][1] * 24.0 / 1800.0 );
	}
		
#ifdef ENTHALPYLIB 
	
		Double         tNow = fSolution->mat[LocationOfMax( fNGridPoints, &fSolution->mat[0][fTemperature], fSolution->phys_rows )][fTemperature];
		if ( fabs(tNow-tlastout) > UNSTLIBOUTDELTAT ) { //
//		if ( CheckOutput( fTime->vec ) ) {
		  cerr << "dump libout" << NEWL; 
		  char   tl[128]; 
		  sprintf( tl, "Lib_Chi%05g_T%04.0f_t%05.2g", GetRefDissRate( fTime->vec[0] ), tNow, fTime->vec[0]*1000 ); 
		  FILE *fp = GetOutfile( tl, kText );
//		  FILE *fp = GetOutputFile( fTime->vec[0], tl, NULL, TFlame::kText ); 
		  WriteFlameletFile( fp, NULL, NULL ); 
		  fclose( fp );
		  tlastout = tNow; //
		}
#endif

		return FALSE;
}
#endif

Flag TTransFlameSolver::CheckOutput( Double *t )
{
//  T = 0; C = 1; W = 2
	const int		vars = 3;
	static Flag		init = FALSE;
	int				i,k;
	Flag			ord = FALSE;
	Double			max[vars];
	Double			pressure = GetPressure( *t );
	static Double	**valsnew, **valsold;
	Double			**nY = fSolMassFracs->mat;
	Double			*nTemp = fSolTemp->vec;
	Double			*prodRate = fSpecies->GetProductionRate()->vec;

	if ( !init ) {
		valsold = New2DArray( fNGridPoints, vars );
		valsnew = New2DArray( fNGridPoints, vars );
		init = TRUE;
	}

// compute W and C
	int	CO2 = fSpecies->FindSpecies( "CO2" );
	int	CO = fSpecies->FindSpecies( "CO" );
	int	H2O = fSpecies->FindSpecies( "H2O" );
	int	H2 = fSpecies->FindSpecies( "H2" );

	for ( int k = 0; k < fNGridPoints; ++k ) {
#ifdef DELTAINEW
		UpdateThermoProps( k, nY[k], nTemp[k], pressure, fProperties->GetDensityRef()
										, kDensFromPress, ( fSoot ) ? fSolSootMoments->mat[k] : NULL );
#else
		T0DFlame::UpdateThermoProps( nY[k], nTemp[k], pressure, fProperties->GetDensityRef()
										, kDensFromPress, ( fSoot ) ? fSolSootMoments->mat[k] : NULL );
#endif
		valsnew[k][0] = fSolution->mat[k][fTemperature];
		valsnew[k][1] = fSolution->mat[k][fFirstSpecies+CO2] + fSolution->mat[k][fFirstSpecies+CO] 
							+ fSolution->mat[k][fFirstSpecies+H2O] + fSolution->mat[k][fFirstSpecies+H2];
		valsnew[k][2] = prodRate[CO2] + prodRate[CO] + prodRate[H2O] + prodRate[H2];
	}

//	get maxes
	max[0] = fabs( fSolution->mat[LocationOfMax( fNGridPoints, &fSolution->mat[0][fTemperature], fSolution->phys_rows )][fTemperature] );
	max[1] = fabs( valsnew[LocationOfMax( fNGridPoints, &valsnew[0][1], vars )][1] );
	max[2] = fabs( valsnew[LocationOfMax( fNGridPoints, &valsnew[0][2], vars )][2] );
		
	for ( k = 0; k < fNGridPoints; ++k ) {
		for ( i = 0; i < vars; ++i ) {
			if ( fabs( ( valsnew[k][i] - valsold[k][i] ) / MAX( 1.0e-10, max[i] ) ) > 0.2 ) {
				fprintf( stderr, "valsnew[k][i] = %g  valsold[k][i] = %g, i = %d, k = %d, max = %g\n", valsnew[k][i], valsold[k][i], i, k, max[i]);
				ord = TRUE;
				break;
			}
		}
		if ( ord == TRUE ) break;
		if ( !ord && fabs( valsnew[k][0] - valsold[k][0] ) > 20.0 ) {
			fprintf( stderr, "Temp 20K up to %g at k = %d\n", valsnew[k][0], k );
			ord = TRUE;
			break;
		}
		if ( fabs( valsnew[k][1] - valsold[k][1] ) > 0.01 ) {
			fprintf( stderr, "ProgVar 0.01 up to %g at k = %d\n", valsnew[k][1], k );
			ord = TRUE;
			break;
		}

	}
	if ( ord == TRUE ) {
		for ( k = 0; k < fNGridPoints; ++k ) {
			for ( i = 0; i < vars; ++i ) {
				valsold[k][i] = valsnew[k][i];
			}
		}
	}
	
	return ord;
}

void TTransFlameSolver::PostIter( Double t )
{
	int		nGm1 = fNGridPoints-1;
	Double	**fSol = fSolution->mat;
	Double	slope = *fActualTimeStepSize * fdDeltaZdt / ( fSolGrid->vec[fNGridPoints] - fSolGrid->vec[nGm1] );
	
	fSolOldTime->vec[fNGridPoints] = fSolTime->vec[fNGridPoints];
	fSolTime->vec[fNGridPoints] = t;

	fSolOldTime->vec[kPrev] = fSolTime->vec[kPrev];
	fSolTime->vec[kPrev] = t;

#ifdef MOVEZRIGHT
	for ( int i = 0; i < fNOfEquations; ++i ) {
		fSol[i][fNGridPoints] = ( fSol[i][fNGridPoints] - slope * fSol[i][nGm1] )
								/ ( 1.0 - slope );
	}
#else
	fSol[fTemperature][fNGridPoints] = Interpol( t, fTempFuelStart, fTStart
									, fTempFuelEnd, fTEnd );
#endif
}

void ResTransFlameSolver( Double *T, Double *y, Double *yPrime, Double *delta
			, int *iRes, Double *rPar, int *iPar )
{
	( ( TTransFlameSolverPtr ) iPar )->ResTransFlameSolver( T, y, yPrime, delta
			, iRes, rPar, iPar );
}

void TTransFlameSolver::ResTransFlameSolver( Double *t, Double *y, Double *yPrime, Double *delta
			, int */*iRes*/, Double */*rPar*/, int */*iPar*/ )
{
	int						nSpeciesIn = fSpecies->GetNSpeciesInSystem();
	int						i, ieq, kAct = fActualPoint;
	Double					sum = 0.0;
	Double					tempMinus;
	Double					tempPlus;


	Double					*prodRate = fSpecies->GetProductionRate()->vec;
	Double					*enth = fSpecies->GetEnthalpy()->vec;
	Double					*molarMass = fSpecies->GetMolarMass()->vec;
	Double					*Le = fSpecies->GetLewisNumber()->vec;

	Double					temp = y[fTemperature];
	Double					*YF = &y[fFirstSpecies];

	Double					*oT = &fSolOldTemp->vec[kAct];
	Double					*oTime = &fSolOldTime->vec[kAct];
	Double					*nT = &fSolTemp->vec[kAct];
	Double					*nTime = &fSolTime->vec[kAct];
	Double					**nY = &fSolMassFracs->mat[kAct];

//	Double					**oMom;
	Double					**nMom;
	Double					*moments;
	int						sootOff;
	int						nSootMoments;
	if ( fSoot ) {
		nMom = &fSolSootMoments->mat[kAct];
//		oMom = &fSolOldSootMoments->mat[kAct];
		moments = fSoot->GetMoments()->vec;
		nSootMoments = fSoot->GetNSootMoments();
		sootOff = fSoot->GetOffsetSootMoments();
	}

	// copy nY to oY and then YF to nY
#ifndef SEMIIMPLI
	copy( nSpeciesIn, nY[kCurr], 1, oY[kCurr], 1 );
#endif
	copy( nSpeciesIn, YF, 1, nY[kCurr], 1 );

	Double					Z = fSolGrid->vec[kAct];
	Double					chi = GetDissRate( t[kCurr], Z );
	Double					pressure = GetPressure( t[kCurr] );

	
	if ( fSoot ) {
		fSoot->MOverRhoToM( &y[sootOff], moments, nSootMoments
				, nY[kCurr], temp, pressure
				, molarMass, nSpeciesIn, fProperties );
	}
	T0DFlame::UpdateThermoProps( nY[kCurr], temp, pressure, fProperties->GetDensityRef()
									, kDensFromPress, moments );

	Double				rho = fProperties->GetDensity();
	Double				heatCap = fProperties->GetMixHeatCapacity();
	Double				mixMolarMass = fProperties->GetMixMolarMass();

	for ( i = 0; i < nSpeciesIn; ++i ) {
		ieq = fFirstSpecies + i;
#ifdef SEMIIMPLI
		delta[ieq] = yPrime[ieq] - prodRate[i] / rho - chi / ( 2.0 * Le[i] ) 
			* ( fWMinus[kAct] * nY[kPrev][i]
			+ fWCurr[kAct] * YF[i]
			+ fWPlus[kAct] * nY[kNext][i] );
#else
		delta[ieq] = yPrime[ieq] - prodRate[i] / rho - chi / ( 2.0 * Le[i] ) 
			* ( fWMinus[kAct] * Interpol( nTime[kCurr], oY[kPrev][i], oTime[kPrev], nY[kPrev][i], nTime[kPrev] )
			+ fWCurr[kAct] * oY[kCurr][i]
			+ fWPlus[kAct] * Interpol( nTime[kCurr], oY[kNext][i], oTime[kNext], nY[kNext][i], nTime[kNext] ) );
#endif
		sum += enth[i] * prodRate[i];
	}
	
#ifdef SEMIIMPLI
//		Double	interTime = MIN( nTime[kPrev], nTime[kNext] );
//		tempMinus = Interpol( nTime[kCurr], oT[kPrev], oTime[kPrev], nT[kPrev], nTime[kPrev] );
//		tempPlus = Interpol( nTime[kCurr], oT[kNext], oTime[kNext], nT[kNext], nTime[kNext] );

	tempMinus = nT[kPrev];
	tempPlus = nT[kNext];

	if ( kAct == 0 ) {
		tempMinus = Interpol( t[kCurr], oT[kPrev], oTime[kPrev], nT[kPrev], nTime[kPrev] );
	}
	if ( kAct == fNGridPoints-1 ) {
		tempPlus = Interpol( t[kCurr], oT[kNext], oTime[kNext], nT[kNext], nTime[kNext] );
	}

	delta[fTemperature] = yPrime[fTemperature] 
			+ ( sum - fDPdt ) / ( rho * heatCap ) - 0.5 * chi
			* ( fWMinus[kAct] * tempMinus
			+ fWCurr[kAct] * temp
			+ fWPlus[kAct] * tempPlus )
			- GetTempSource( Z * GetZR() );
#else
	delta[fTemperature] = yPrime[fTemperature] 
			+ ( sum - fDPdt ) / ( rho * heatCap ) - 0.5 * chi
			* ( fWMinus[kAct] * Interpol( nTime[kCurr], oT[kPrev], oTime[kPrev], nT[kPrev], nTime[kPrev] )
			+ fWCurr[kAct] * nT[kCurr]
			+ fWPlus[kAct] * Interpol( nTime[kCurr], oT[kNext], oTime[kNext], nT[kNext], nTime[kNext] ) )
			- GetTempSource( Z * GetZR() );
#endif
	if ( fProperties->GetRadiation() ) {
		delta[fTemperature] -= ( fProperties->GetRadiation()->GetRadiation() ) / ( rho * heatCap );
	}
	
	if ( fSoot ) {
		// ATTENTION (hp)
		// calculation of soot with semi explicit solver currently doesn't
		// work, since y[sootOff] is M/rho and nMom is M
		fprintf( fOutFilePtr, "###error: calculation of soot with semi explicit solver currently doesn't work\n" );
		exit(2);

		Double	fracIndex;
		for ( i = 0; i < nSootMoments; ++i ) {
			fracIndex = i - 2.0 / 3.0;
			ieq = sootOff + i;
			delta[ieq] = yPrime[ieq] 
				- fSoot->NucleationNew( i, temp, fSoot->GetPAHMoments()->vec ) / rho 
				- chi / ( 2.0 * fSoot->GetLewis1() ) 
#ifdef SIZEDEPDIFFUSION
				* ( fWMinus[kAct] * fSoot->FracMom2( fracIndex, nMom[kPrev] )
				+ fWCurr[kAct] * fSoot->FracMom2( fracIndex, &y[sootOff] )
				+ fWPlus[kAct] * fSoot->FracMom2( fracIndex, nMom[kNext] ) );
#else
				* ( fWMinus[kAct] * nMom[kPrev][i]
				+ fWCurr[kAct] * y[sootOff+i]
				+ fWPlus[kAct] * nMom[kNext][i] );
#endif
//			fprintf( fOutFilePtr, "Nuc%d = %g\n", i
//					, fSoot->NucleationNew( i, tempCurr, fSoot->GetPAHMoments()->vec ) / rho);
		}
	}
//	fprintf( fOutFilePtr, "%s\t%20.10E\n", fVariableNames[fTemperature], temp );
#ifdef DEBUGRES
	for ( int j = 0; j < nSpeciesIn+fVariablesWithoutSpecies; ++j ) {
		fprintf( stderr, "%s\t%g\t%g\n", fVariableNames[j], y[j], delta[j] );
	}
	fprintf( stderr, "\n" );
	for ( j = nSpeciesIn; j < fSpecies->GetNOfSpecies(); ++j ) {
		fprintf( stderr, "SteadyState %s:\t%g\n", fSpecies->GetNames()[j], nY[kCurr][j] );
	}
	fprintf( stderr, "\n" );
#endif
	
}

#ifndef UNICOS
void ResTransFlameImpliSolver( Double *T, Double *y, Double *yPrime, Double *delta
			, int *iRes, Double *rPar, int *iPar )
{
	( ( TTransFlameSolverPtr ) iPar )->ResTransFlameImpliSolver( T, y, yPrime, delta
			, iRes, rPar, iPar );
}

void TTransFlameSolver::ResTransFlameImpliSolver(double * t, double * y, double * yPrime, double * delta, int * iRes, double * /*rPar*/, int * /*iPar*/)
{
  int k;
  int nSpeciesIn = fSpecies->GetNSpeciesInSystem();
  int i, ieq, gpOff, gpSpecOff;
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
  int gpGridOff;
  double * gpDens = fgpDens->vec;
  double kappa_ = fKappa * ( 1.0 + fKappa );
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo	
  double sum, sumYH;

  double * prodRate = fSpecies->GetProductionRate()->vec;
  double * enth = fSpecies->GetEnthalpy()->vec;
  double * molarMass = fSpecies->GetMolarMass()->vec;
#ifdef LEWISCHANGE
  double * LeOrig = fSpecies->GetLewisNumber()->vec;
  double * Le = new Double[nSpeciesIn];
  double LeFunc;
  double switchTime = LEWISSWITCHTIME;
  double switchPeriod = LEWISSWITCHPERIOD;
  const double Pi = 4.0 * atan(1.0);
  double omega = 2.0 * Pi / switchPeriod;
  double tStar = *t - (switchTime - 0.5 * switchPeriod);

  if (tStar >= 0.0)
  {
    if (tStar < switchPeriod)
    {
      LeFunc = 0.5 * (cos(omega * tStar) + 1.0);
      for (int iLe = 0; iLe < nSpeciesIn; ++iLe)
	Le[iLe] = LeFunc * (LeOrig[iLe] - 1.0) + 1.0;
    }
    else
      for (int iLe = 0; iLe < nSpeciesIn; ++iLe)
	Le[iLe] = 1.0;
  }
  else
    for (int iLe = 0; iLe < nSpeciesIn; ++iLe)
      Le[iLe] = LeOrig[iLe];
#else
  double * Le = fSpecies->GetLewisNumber()->vec;
#endif

  double ** nY = fSolMassFracs->mat;
  double * nTemp = fSolTemp->vec;
  double * nProg = fSolProg->vec;
  double * nEnth = fSolEnth->vec;

  double *Z = fSolGrid->vec;
  double chi, chiPrev, chiNext;
  double rhoChiCurr, rhoChiPrev, rhoChiNext, rhoChiCpOverLambda;
  double pressure = GetPressure(*t);
  double sumCpidYidZ;
  double * cpiOld = fSpecies->GetHeatCapacity()->vec;
  double * cpMix = fHeatCpMix->vec;
  double * mu = fViscosity->vec;
  double * mixMolarMass = fMolarMassMix->vec;
  double * lambdaOverCp = fLambdaOverCpMix->vec;
  double ** cpi = fSpecHeatCp->mat;
  double ** hi = fSpecEnthalpy->mat;
  double ** mDoti = fProdRate->mat;
  double * rho = fDensity->vec;
  double * diffTermY = fDiffTermY->vec;
  double diffCorrY;
  double sumdYdZ, sumdYdZOverLe;

#ifdef MOLARDIFF
  double diffCorrW;
  double * diffTermW = fDiffTermW->vec;
  double d2WdZ2, dWdZ;
  double diffWFact;
#endif
  double halfChi;
  double ** moments;
  double * MOverRhoPrev;
  double * MOverRhoCurr;
  double * MOverRhoNext;
  int sootOff;
  int nSootMoments;

  double * nrhodot = rhodot->vec;
  double * nenthdot = enthdot->vec;

  y = &y[fNOfEquations];
  yPrime = &yPrime[fNOfEquations];
  delta = &delta[fNOfEquations];

#ifdef READU			
  int iU, kU, ieqU;
  double UstOverU;
  for (kU = -1; kU <= fNGridPoints; ++kU)
  {
    ieqU = kU * fNOfEquations;
    UstOverU = GetU(*t, Z[kU]);
    for (iU = 0; iU < fNOfEquations; ++iU)
      yPrime[ieqU+iU] /= UstOverU;
  }
#endif

//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
  SetFDWeights(y);
  ComputeGPDistributionWeights(y);
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo

  if (fSoot)
  {
    moments = fSolSootMoments->mat;
    nSootMoments = fSoot->GetNSootMoments();
    sootOff = fSoot->GetOffsetSootMoments();
  }
//	fprintf( stderr, "now\n" );
	// copy YF to nY
//	Double	temp10, rho10, p10;
  for (k = -1; k <= fNGridPoints; ++k)
  {
    copy(nSpeciesIn, &y[k*fNOfEquations+fFirstSpecies], 1, nY[k], 1);
    nTemp[k] = y[k*fNOfEquations+fTemperature];
    nProg[k] = y[k*fNOfEquations+fProg];
    nEnth[k] = y[k*fNOfEquations+fEnth];
//		Z[k] = y[k*fNOfEquations+fGrid];
    if (fSoot)
    {
      if ( k > -1 ) MOverRhoPrev = &y[(k-1)*fNOfEquations+sootOff];
      if ( k < fNGridPoints ) MOverRhoNext = &y[(k+1)*fNOfEquations+sootOff];
      MOverRhoCurr = &y[k*fNOfEquations+sootOff];
      fSoot->MOverRhoToM(&y[k*fNOfEquations+sootOff], moments[k], nSootMoments, nY[k], nTemp[k], pressure, molarMass, nSpeciesIn, fProperties);
    }

    if (nTemp[k] < 100.0 || nTemp[k] > 5000.0)
    {
      fprintf(fOutFilePtr, "###Warning: Temp = %g @ gp = %d and Z = %g\n", nTemp[k], k, Z[k]);
      *iRes = -1;
      return;
    }
#ifdef DELTAINEW
    UpdateThermoProps(k, nY[k], nTemp[k], pressure, fProperties->GetDensityRef(), kDensFromPress, ( fSoot ) ? moments[k] : NULL , &nrhodot[k], &nenthdot[k]);
#else
    T0DFlame::UpdateThermoProps(nY[k], nTemp[k], pressure, fProperties->GetDensityRef(), kDensFromPress, ( fSoot ) ? moments[k] : NULL, &nrhodot[k], &nenthdot[k]);
#endif
    cpMix[k] = fProperties->GetMixHeatCapacity();
    rho[k] = fProperties->GetDensity();
    mu[k] = fProperties->GetMixViscosity();
    lambdaOverCp[k] = fProperties->GetMixConductivity() / fProperties->GetMixHeatCapacity();
    mixMolarMass[k] = fProperties->GetMixMolarMass();
    copy(nSpeciesIn, fSpecies->GetHeatCapacity()->vec, 1, cpi[k], 1);
    copy(nSpeciesIn, fSpecies->GetEnthalpy()->vec, 1, hi[k], 1);
    copy(nSpeciesIn, fSpecies->GetProductionRate()->vec, 1, mDoti[k], 1);
  }
	
  // Start Transport Equation Loop
  for (k = 0; k < fNGridPoints; ++k)
  {
    gpOff = k * fNOfEquations;
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
    gpGridOff = gpOff + fGrid;	
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
		
/*		chi = GetDissRate( *t, y[gpGridOff] );*/
/*		chiPrev = GetDissRate( *t, y[gpGridOff-fNOfEquations] );*/
/*		chiNext = GetDissRate( *t, y[gpGridOff+fNOfEquations] );*/

    //Simulation Chi: Conserved
    chi = GetDissRate(*t, Z[k]);
    chiPrev = GetDissRate(*t, Z[k-1]);
    chiNext = GetDissRate(*t, Z[k+1]);

    //Simulation Chi: Bilger
    //chi = GetDissRate(*t, ComputeZBilger(nY[k],nY[fNGridPoints],nY[kPrev]));
    //chiPrev = GetDissRate(*t, ComputeZBilger(nY[k-1],nY[fNGridPoints],nY[kPrev]));
    //chiNext = GetDissRate(*t, ComputeZBilger(nY[k+1],nY[fNGridPoints],nY[kPrev]));

    rhoChiCurr = rho[k] * chi;
    rhoChiPrev = rho[k-1] * chiPrev;
    rhoChiNext = rho[k+1] * chiNext;
    rhoChiCpOverLambda = rhoChiCurr / lambdaOverCp[k];
    halfChi = 0.5 * chi;
		
#ifdef MOLARDIFF
    d2WdZ2 = fWMinus[k] * mixMolarMass[k-1] + fWCurr[k] * mixMolarMass[k] + fWPlus[k] * mixMolarMass[k+1];
    dWdZ = fFDWMinus[k] * mixMolarMass[k-1] + fFDWCurr[k] * mixMolarMass[k] + fFDWPlus[k] * mixMolarMass[k+1];
    diffWFact = halfChi / mixMolarMass[k] * d2WdZ2;
    diffCorrW = 0.0;
#endif
    sumdYdZ = sumdYdZOverLe = 0.0;
    diffCorrY = 0.0;
    double dYdZ;
    for (i = 0; i < nSpeciesIn; ++i)
    {
      diffTermY[i] = halfChi / Le[i] * (fWMinus[k] * nY[k-1][i] + fWCurr[k] * nY[k][i]	+ fWPlus[k] * nY[k+1][i]);
      diffCorrY += diffTermY[i];
      dYdZ = fFDWMinus[k] * nY[k-1][i] + fFDWCurr[k] * nY[k][i]	+ fFDWPlus[k] * nY[k+1][i];
      sumdYdZ += dYdZ;
      sumdYdZOverLe += dYdZ / Le[i];
#ifdef MOLARDIFF
      diffTermW[i] = diffWFact * nY[k][i] / Le[i];
      diffCorrW += diffTermW[i];
#endif
    }
				
#ifdef NEWPROPERTIES
#else
#  ifdef DELTAINEW
    UpdateThermoProps(k, nY[k], nTemp[k], pressure, fProperties->GetDensityRef(), kDensFromPress, ( fSoot ) ? moments[k] : NULL);
#  else 
    T0DFlame::UpdateThermoProps(nY[k], nTemp[k], pressure, fProperties->GetDensityRef(), kDensFromPress, ( fSoot ) ? moments[k] : NULL);
#  endif

/*		cpi[k] = cpiOld;*/
/*		hi[k] = enth;*/
//		mDoti[k] = prodRate;
    copy(nSpeciesIn, fSpecies->GetHeatCapacity()->vec, 1, cpi[k], 1);
    copy(nSpeciesIn, fSpecies->GetEnthalpy()->vec, 1, hi[k], 1);
    copy(nSpeciesIn, fSpecies->GetProductionRate()->vec, 1, mDoti[k], 1);
#endif
    gpSpecOff = gpOff + fFirstSpecies;

    //double constConvVeloY = 0.25 / rho[k] * ((fFDWMinus[k] * rhoChiPrev + fFDWCurr[k] * rhoChiCurr + fFDWPlus[k] * rhoChiNext) + rhoChiCpOverLambda * (fFDWMinus[k] * lambdaOverCp[k-1]	+ fFDWCurr[k] * lambdaOverCp[k]	+ fFDWPlus[k] * lambdaOverCp[k+1]));
    double constConvVeloY = 0.25 / (rho[k] * lambdaOverCp[k]) * (fFDWMinus[k] * lambdaOverCp[k-1] * rhoChiPrev + fFDWCurr[k] * lambdaOverCp[k] * rhoChiCurr + fFDWPlus[k] * lambdaOverCp[k+1] * rhoChiNext);

    sum = 0.0;
    sumYH = 0.0;
    sumCpidYidZ = 0.0;

    // Species Mass Fraction Equations
    for (i = 0; i < nSpeciesIn; ++i)
    {
      ieq = gpSpecOff + i;
      delta[ieq] = yPrime[ieq] - mDoti[k][i] / rho[k] - diffTermY[i];
#ifdef MOLARDIFF
      delta[ieq] -= diffTermW[i];
#endif
#ifdef DIFFCORR
      delta[ieq] += nY[k][i] * diffCorrY;
#  ifdef MOLARDIFF
      delta[ieq] += nY[k][i] * diffCorrW;
#  endif
#endif

#ifdef CONVECTION
      // Mass Fraction Convection
      double convVeloY = (1.0 / Le[i] - 1.0) * constConvVeloY;
#  ifdef UPWIND
      if (convVeloY < 0.0)
	delta[ieq] -= convVeloY * (nY[k][i] - nY[k-1][i] ) / (fhm[k]);
      else
	delta[ieq] -= convVeloY * (nY[k+1][i] - nY[k][i] ) / (fh[k]);
#  else
      delta[ieq] -= convVeloY * (fFDWMinus[k] * nY[k-1][i] + fFDWCurr[k] * nY[k][i] + fFDWPlus[k] * nY[k+1][i]);
#  endif

#  ifdef MOLARDIFF
      // Molar Mass Convection
      double convVeloM = 0.0;
      convVeloM += 0.25 / (Le[i] * rho[k]) * (fFDWMinus[k] * rhoChiPrev * nY[k-1][i] / mixMolarMass[k-1] + fFDWCurr[k] * rhoChiCurr * nY[k][i] / mixMolarMass[k] + fFDWPlus[k] * rhoChiNext * nY[k+1][i] / mixMolarMass[k+1]);
      convVeloM += 0.25 / (Le[i] * rho[k]) * rhoChiCpOverLambda * (fFDWMinus[k] * lambdaOverCp[k-1] * nY[k-1][i] / mixMolarMass[k-1] + fFDWCurr[k] * lambdaOverCp[k] * nY[k][i] / mixMolarMass[k] + fFDWPlus[k] * lambdaOverCp[k+1] * nY[k+1][i] / mixMolarMass[k+1]);

#    ifdef UPWIND
      if (convVeloM > 0.0)
	delta[ieq] -= convVeloM * (mixMolarMass[k] - mixMolarMass[k-1]) / (fhm[k]);
      else
	delta[ieq] -= convVeloM * (mixMolarMass[k+1] - mixMolarMass[k]) / (fh[k]);
#    else
      delta[ieq] -= convVeloM * (fFDWMinus[k] * mixMolarMass[k-1] + fFDWCurr[k] * mixMolarMass[k] + fFDWPlus[k] * mixMolarMass[k+1]);
#    endif
#  endif
#  ifdef DIFFCORR
      // Mass Fraction Diffusivity Correction Convection
      delta[ieq] += 0.25 * sumdYdZOverLe / rho[k] * (fFDWMinus[k] * rhoChiPrev * nY[k-1][i] + fFDWCurr[k] * rhoChiCurr * nY[k][i] + fFDWPlus[k] * rhoChiNext * nY[k+1][i]);
      delta[ieq] += 0.25 * sumdYdZOverLe / rho[k] * + rhoChiCpOverLambda * (fFDWMinus[k] * lambdaOverCp[k-1] * nY[k-1][i] + fFDWCurr[k] * lambdaOverCp[k] * nY[k][i] + fFDWPlus[k] * lambdaOverCp[k+1] * nY[k+1][i]);
#    ifdef MOLARDIFF
      // Molar Mass Diffusivity Correction Convection
      double sumDiffCorrMM1 = 0.0;
      double sumDiffCorrMM2 = 0.0;
      double rhoChiYOverMPrev = rhoChiPrev * nY[k-1][i] / mixMolarMass[k-1];
      double rhoChiYOverMCurr = rhoChiCurr * nY[k  ][i] / mixMolarMass[k  ];
      double rhoChiYOverMNext = rhoChiNext * nY[k+1][i] / mixMolarMass[k+1];
      double lambdaOCpYOMPrev = lambdaOverCp[k-1] * nY[k-1][i] / mixMolarMass[k-1];
      double lambdaOCpYOMCurr = lambdaOverCp[k  ] * nY[k  ][i] / mixMolarMass[k  ];
      double lambdaOCpYOMNext = lambdaOverCp[k+1] * nY[k+1][i] / mixMolarMass[k+1];
      for (int j = 0; j < nSpeciesIn; ++j)
      {
	sumDiffCorrMM1 += (fFDWMinus[k] * rhoChiYOverMPrev * nY[k-1][j] + fFDWCurr[k] * rhoChiYOverMCurr * nY[k][j] + fFDWPlus[k] * rhoChiYOverMNext * nY[k+1][j]) / Le[j];
	sumDiffCorrMM2 += (fFDWMinus[k] * lambdaOCpYOMPrev * nY[k-1][j] + fFDWCurr[k] * lambdaOCpYOMCurr * nY[k][j] + fFDWPlus[k] * lambdaOCpYOMNext * nY[k+1][j]) / Le[j];
      }		

#      ifdef UPWIND
      fprintf( stderr, "upwind not completely implemented\n" );
      if (convVeloM > 0.0)
	delta[ieq] -= convVeloM * (mixMolarMass[k] - mixMolarMass[k-1]) / (fhm[k]);
      else
	delta[ieq] -= convVeloM * (mixMolarMass[k+1] - mixMolarMass[k]) / (fh[k]);
#      else
      delta[ieq] += 0.25 / rho[k] * (sumDiffCorrMM1 + rhoChiCpOverLambda * sumDiffCorrMM2) * dWdZ;
#      endif
#    endif
#  endif
#endif
      sum += hi[k][i] * mDoti[k][i];
      sumYH += hi[k][i] * nY[k][i];

      double cptmp = cpi[k][i];
#ifdef DIFFCORR
      cptmp -= cpMix[k];
#endif
      sumCpidYidZ += cptmp / Le[i] * (fFDWMinus[k] * nY[k-1][i] + fFDWCurr[k] * nY[k][i] + fFDWPlus[k]  * nY[k+1][i]);
#ifdef MOLARDIFF
      sumCpidYidZ += cptmp / Le[i] * (nY[k][i] / mixMolarMass[k] * dWdZ);
#endif

      // Soot Density Correction
      double convVelorhodot = (nrhodot[k] * Z[k] - ComputeZBilgerSource(mDoti[k],nY[fNGridPoints],nY[kPrev], nrhodot[k])) / rho[k];
#ifdef UPWIND
      if (convVeloYrhodot < 0.0)
	delta[ieq] -= convVelorhodot * (nY[k][i] - nY[k-1][i] ) / (fhm[k]);
      else
	delta[ieq] -= convVelorhodot * (nY[k+1][i] - nY[k][i] ) / (fh[k]);
#else
      delta[ieq] -= convVelorhodot * (fFDWMinus[k] * nY[k-1][i] + fFDWCurr[k] * nY[k][i] + fFDWPlus[k] * nY[k+1][i]);
#endif
      
      delta[ieq] += nrhodot[k] * nY[k][i] / rho[k];

//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
#ifndef NEWCONVECTION
#  ifdef CENTRALGRID
      delta[ieq] += - yPrime[gpGridOff] * (fFDWMinus[k] * nY[k-1][i] + fFDWCurr[k] * nY[k][i] + fFDWPlus[k] * nY[k+1][i]);
#  else
      if ((-yPrime[gpGridOff]) > 0.0)
	delta[ieq] += - (nY[k][i] - nY[k-1][i]) / fhm[k] * yPrime[gpGridOff];
      else
	delta[ieq] += - (nY[k+1][i] - nY[k][i]) / fh[k] * yPrime[gpGridOff];
#  endif
#endif
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo

    }
	
    // Energy Equation
    ieq = gpOff + fTemperature;
#ifdef TOTENT
//		Double	totentstart = Interpol( Z[k] , fTotEntStart->vec[kPrev], Z[kPrev] , fTotEntStart->vec[fNGridPoints], Z[fNGridPoints] );
//		Double	totentend = Interpol( Z[k] , fTotEntEnd->vec[kPrev], Z[kPrev] , fTotEntEnd->vec[fNGridPoints], Z[fNGridPoints] );
//		Double	totEntPresc = Interpol( *t, totentstart, fTStart, totentend, fTEnd );
//		delta[ieq] = totEntPresc - sumYH;
    delta[ieq] = GetTotEnt(k, Z, *t) - sumYH;
#else
    delta[ieq] = yPrime[ieq] + (sum - fDPdt) / (rho[k] * cpMix[k]) - 0.5 * chi * (fWMinus[k] * nTemp[k-1] + fWCurr[k] * nTemp[k] + fWPlus[k] * nTemp[k+1]) - GetTempSource(y[gpGridOff] * GetZR());
#  ifdef CENTRALTEMP
#    ifdef ENTFLUX
    delta[ieq] -= 0.5 * chi/ cpMix[k] * (fFDWMinus[k] * nTemp[k-1] +fFDWCurr[k] * nTemp[k] +fFDWPlus[k] * nTemp[k+1]) * sumCpidYidZ;
#    endif
#    ifdef HEATCAPGRAD
    delta[ieq] -= 0.5 * chi/ cpMix[k] * (fFDWMinus[k] * nTemp[k-1] +fFDWCurr[k] * nTemp[k] +fFDWPlus[k] * nTemp[k+1]) * (fFDWMinus[k] * cpMix[k-1] + fFDWCurr[k] * cpMix[k] + fFDWPlus[k] * cpMix[k+1]);
#    endif
#  else
#    ifdef ENTFLUX
    double convVeloEntFlux = 0.5 * chi/ cpMix[k] * sumCpidYidZ;
    if (convVeloEntFlux < 0.0)
      delta[ieq] -= convVeloEntFlux * (nTemp[k] - nTemp[k-1]) / (fhm[k]);
    else
      delta[ieq] -= convVeloEntFlux * (nTemp[k+1] - nTemp[k]) / (fh[k]);
#    endif
#    ifdef HEATCAPGRAD
    double convVelocpdT = 0.5 * chi/ cpMix[k] * (fFDWMinus[k] * cpMix[k-1] + fFDWCurr[k] * cpMix[k] + fFDWPlus[k] * cpMix[k+1]);
    if (convVelocpdT < 0.0)
      delta[ieq] -= convVelocpdT * (nTemp[k] - nTemp[k-1]) / (fhm[k]);
    else
      delta[ieq] -= convVelocpdT * (nTemp[k+1] - nTemp[k]) / (fh[k]);
#    endif
#  endif

    if (fProperties->GetRadiation())
    {
#  ifdef READRADHEATLOSS
      delta[ieq] += 1.00 * GetDelQ(*t, Z[k]) / (rho[k] * cpMix[k]);
#  else
      // wird zweimal aufgerufen, koennte aber in TFlame.cp auskommentiert werden
      fProperties->GetRadiation()->SetRadiation(nTemp[k], nY[k], molarMass, rho[k]); 
      delta[ieq] -= (fProperties->GetRadiation()->GetRadiation()) / (rho[k] * cpMix[k]);

      if (fSoot)
      {
#    ifdef NOSOOTRAD
#    else
#      ifdef ROSSELAND
	delta[ieq] += GetRosseRadiation(k, &nTemp[k], &moments[k], rhoChiPrev / lambdaOverCp[k-1], rhoChiCpOverLambda, rhoChiNext / lambdaOverCp[k+1]) / (rho[k] * cpMix[k]);
#      else
	delta[ieq] += fSoot->GetSootRadiation(nTemp[k], moments[k]) / (rho[k] * cpMix[k]);
#      endif
#    endif
      }
#  endif
    }

    // Soot Density Correction
    double convVelorhodot = (nrhodot[k] * Z[k] - ComputeZBilgerSource(mDoti[k],nY[fNGridPoints],nY[kPrev], nrhodot[k])) / rho[k];
#  ifndef CENTRALTEMP
    if (convVeloYrhodot < 0.0)
      delta[ieq] -= convVelorhodot * (nTemp[k] - nTemp[k-1]) / (fhm[k]);
    else
      delta[ieq] -= convVelorhodot * (nTemp[k+1] - nTemp[k]) / (fh[k]);
#  else
    delta[ieq] -= convVelorhodot * (fFDWMinus[k] * nTemp[k-1] + fFDWCurr[k] * nTemp[k] + fFDWPlus[k] * nTemp[k+1]);
#  endif
      
    delta[ieq] -= nenthdot[k] / (rho[k] * cpMix[k]);;

#endif
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
#ifndef NEWCONVECTION
#  ifdef CENTRALGRID
    delta[ieq] += - yPrime[gpGridOff] * (fFDWMinus[k] * nTemp[k-1] + fFDWCurr[k] * nTemp[k] + fFDWPlus[k] * nTemp[k+1]);
#  else
    if ((-yPrime[gpGridOff] ) > 0.0)
      delta[ieq] += - (nTemp[k] - nTemp[k-1]) / fhm[k] * yPrime[gpGridOff];
    else
      delta[ieq] += - (nTemp[k+1] - nTemp[k]) / fh[k] * yPrime[gpGridOff];
#  endif
#endif
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo

    // Progress Variable Equation
    ieq = gpOff + fProg;

    double diffTermProg = halfChi * (fWMinus[k]*nProg[k-1] + fWCurr[k]*nProg[k] + fWPlus[k]*nProg[k+1]);
    delta[ieq] = yPrime[ieq] - diffTermProg;

    int f_CO2 = fSpecies->FindSpecies("CO2");
    int f_CO = fSpecies->FindSpecies("CO");
    int f_H2O = fSpecies->FindSpecies("H2O");
    int f_H2 = fSpecies->FindSpecies("H2");
    double prog_src = mDoti[k][f_CO2] + mDoti[k][f_CO] + mDoti[k][f_H2O] + mDoti[k][f_H2];
    double C_MAX = ComputeCMAX(nY[k], nY[fNGridPoints], nY[kPrev]);
    delta[ieq] -= (prog_src) / C_MAX / rho[k];

    double convVeloCrhodot = (nrhodot[k]*Z[k] - ComputeZBilgerSource(mDoti[k],nY[fNGridPoints],nY[kPrev],nrhodot[k])) / rho[k];
#ifdef UPWIND
    if (convVeloCrhodot < 0.0)
      delta[ieq] -= convVeloCrhodot * (nProg[k]-nProg[k-1]) / fhm[k];
    else
      delta[ieq] -= convVeloCrhodot * (nProg[k+1]-nProg[k]) / fh[k];
#else
    delta[ieq] -= convVeloCrhodot * (fFDWMinus[k]*nProg[k-1] + fFDWCurr[k]*nProg[k] + fFDWPlus[k]*nProg[k+1]);
#endif

    delta[ieq] += nrhodot[k] * nProg[k] / rho[k];

//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
#ifndef NEWCONVECTION
#  ifdef CENTRALGRID
    delta[ieq] += - yPrime[gpGridOff] * (fFDWMinus[k] * nProg[k-1] + fFDWCurr[k] * nProg[k] + fFDWPlus[k] * nProg[k+1]);
#  else
    if ((-yPrime[gpGridOff]) > 0.0)
      delta[ieq] += - (nProg[k] - nProg[k-1]) / fhm[k] * yPrime[gpGridOff];
    else
      delta[ieq] += - (nProg[k+1] - nProg[k]) / fh[k] * yPrime[gpGridOff];
#  endif
#endif
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo

    // Unity Lewis Number Enthalpy Equation
    ieq = gpOff + fEnth;

    double diffTermEnth = halfChi * (fWMinus[k]*nEnth[k-1] + fWCurr[k]*nEnth[k] + fWPlus[k]*nEnth[k+1]);
    delta[ieq] = yPrime[ieq] - diffTermEnth;

    //delta[ieq] -= nenthdot[k] / rho[k];// / 1.0e20;

    double convVeloHrhodot = (nrhodot[k]*Z[k] - ComputeZBilgerSource(mDoti[k],nY[fNGridPoints],nY[kPrev],nrhodot[k])) / rho[k];
#ifdef UPWIND
    if (convVeloHrhodot < 0.0)
      delta[ieq] -= convVeloHrhodot * (nEnth[k]-nEnth[k-1]) / fhm[k];
    else
      delta[ieq] -= convVeloHrhodot * (nEnth[k+1]-nEnth[k]) / fh[k];
#else
    delta[ieq] -= convVeloHrhodot * (fFDWMinus[k]*nEnth[k-1] + fFDWCurr[k]*nEnth[k] + fFDWPlus[k]*nEnth[k+1]);
#endif

    //delta[ieq] += nrhodot[k] * nEnth[k] / rho[k];

    if (fProperties->GetRadiation())
    {
#  ifdef READRADHEATLOSS
      delta[ieq] += 1.00 * GetDelQ(*t, Z[k]) / rho[k];
#  else
      // wird zweimal aufgerufen, koennte aber in TFlame.cp auskommentiert werden
      fProperties->GetRadiation()->SetRadiation(nTemp[k], nY[k], molarMass, rho[k]); 
      delta[ieq] -= (fProperties->GetRadiation()->GetRadiation()) / rho[k];

      if (fSoot)
      {
#    ifdef NOSOOTRAD
#    else
#      ifdef ROSSELAND
	delta[ieq] += GetRosseRadiation(k, &nTemp[k], &moments[k], rhoChiPrev / lambdaOverCp[k-1], rhoChiCpOverLambda, rhoChiNext / lambdaOverCp[k+1]) / rho[k];
#      else
	delta[ieq] += fSoot->GetSootRadiation(nTemp[k], moments[k]) / rho[k];
#      endif
#    endif
      }
#  endif
    }

//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
#ifndef NEWCONVECTION
#  ifdef CENTRALGRID
    delta[ieq] += - yPrime[gpGridOff] * (fFDWMinus[k] * nEnth[k-1] + fFDWCurr[k] * nEnth[k] + fFDWPlus[k] * nEnth[k+1]);
#  else
    if ((-yPrime[gpGridOff]) > 0.0)
      delta[ieq] += - (nEnth[k] - nEnth[k-1]) / fhm[k] * yPrime[gpGridOff];
    else
      delta[ieq] += - (nEnth[k+1] - nEnth[k]) / fh[k] * yPrime[gpGridOff];
#  endif
#endif
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo

    // Soot Equations
    if (fSoot)
    {
      MOverRhoPrev = &y[(k-1)*fNOfEquations+sootOff];
      MOverRhoNext = &y[(k+1)*fNOfEquations+sootOff];
      MOverRhoCurr = &y[k*fNOfEquations+sootOff];

#ifdef NEWPROPERTIES
      fSoot->UpdateSoot(fReaction, fSpecies, moments[k], nTemp[k], nY[k], rho[k], mixMolarMass[k]);
#endif

      gpSpecOff = gpOff + sootOff;

      double * logPrev = &y[(k-1)*fNOfEquations+sootOff];
      double * logCurr = &y[k*fNOfEquations+sootOff];
      double * logNext = &y[(k+1)*fNOfEquations+sootOff];

      double dRhoChi_dZ = fFDWMinus[k] * rhoChiPrev + fFDWCurr[k] * rhoChiCurr + fFDWPlus[k] * rhoChiNext;

      double RhoDiffPrev = lambdaOverCp[k-1];
      double RhoDiff = lambdaOverCp[k];
      double RhoDiffNext = lambdaOverCp[k+1];
      double dRhoDiff_dZ = fFDWMinus[k] * RhoDiffPrev + fFDWCurr[k] * RhoDiff + fFDWPlus[k] * RhoDiffNext;
      
      double SootVel = 0.25 * (dRhoChi_dZ + rhoChiCurr / RhoDiff * dRhoDiff_dZ);

      double MuPrev = mu[k-1];
      double Mu = mu[k];
      double MuNext = mu[k+1];

      double dT_dZ = fFDWMinus[k] * nTemp[k-1] + fFDWCurr[k] * nTemp[k] + fFDWPlus[k] * nTemp[k+1];
      double d2T_dZ2 = fWMinus[k] * nTemp[k-1] + fWCurr[k] * nTemp[k] + fWPlus[k] * nTemp[k+1];
      
      double ThermVel = 0.5 * 0.55 * Mu/RhoDiff * rhoChiCurr / nTemp[k] * dT_dZ;

      double RhoChi_TPrev = rhoChiPrev / nTemp[k-1];
      double RhoChi_T = rhoChiCurr / nTemp[k];
      double RhoChi_TNext = rhoChiNext / nTemp[k+1];
      double dRhoChi_T_dZ = fFDWMinus[k] * RhoChi_TPrev + fFDWCurr[k] * RhoChi_T + fFDWPlus[k] * RhoChi_TNext;

      double dMu_dZ = fFDWMinus[k] * MuPrev + fFDWCurr[k] * Mu + fFDWPlus[k] * MuNext;

      SootVel = SootVel - Z[k] * nrhodot[k] + ComputeZBilgerSource(mDoti[k],nY[fNGridPoints],nY[kPrev], nrhodot[k]);
      if (fSoot->WithThermoPhoresis())
	SootVel = SootVel - ThermVel;

      double i1, i2;

      fSoot->rho = rho[k];
      fSoot->mu = mu[k];
      fSoot->Wmix = mixMolarMass[k];
      
      for (i = 0; i < nSootMoments; ++i)
      {
	ieq = gpSpecOff + i;

	delta[ieq] = yPrime[ieq];

#ifndef NOSOOTSOLVE
	// Differential Diffusion (Infinite Lewis Number): Convective Term
	if (SootVel > 0.0)
	  delta[ieq] += (SootVel / rho[k]) * (logCurr[i] - logPrev[i]) / fhm[k];
	else
	  delta[ieq] += (SootVel / rho[k]) * (logNext[i] - logCurr[i]) / fh[k];
	//delta[ieq] += (SootVel / rho[k]) * (fFDWMinus[k] * logPrev[i] + fFDWCurr[k] * logCurr[i] + fFDWPlus[k] * logNext[i]);

	// Density Correction Source Term
	delta[ieq] += nrhodot[k] / rho[k];

	// Thermophoresis Source Term
	if (fSoot->WithThermoPhoresis())
	{
	  delta[ieq] -= 0.5 * (0.55 * Mu/RhoDiff) * (dT_dZ*dRhoChi_T_dZ + RhoChi_T * d2T_dZ2) / rho[k];
	  delta[ieq] += 0.25 * (0.55 * Mu/RhoDiff) * dRhoChi_dZ / nTemp[k] * dT_dZ / rho[k];
	  delta[ieq] += 0.25 * (0.55 / nTemp[k] * Mu/RhoDiff/RhoDiff * dT_dZ * rhoChiCurr * dRhoDiff_dZ) / rho[k];
	  delta[ieq] -= 0.5 * 0.55 / nTemp[k] / RhoDiff * rhoChiCurr * dT_dZ * dMu_dZ / rho[k];
	}

	// Soot Source Terms
	i1 = double(fSoot->Geti(i))/6.0;
	i2 = double(fSoot->Getj(i))/6.0;

	//double fact = 0.0;
	//fact = 0.5*(1+tanh(1000.0*(*t-0.06)));
	if (i != (nSootMoments-1))
	{
	  if (fSoot->WithNucleation())
	    delta[ieq] -= fSoot->NucleationSource(i1, i2, nTemp[k]) / moments[k][i];
	  if (fSoot->WithCoagulation())
	    delta[ieq] -= fSoot->CoagulationSource(i1, i2, nTemp[k], moments[k]) / moments[k][i];
	  if (fSoot->WithCondensation())
	    delta[ieq] -= fSoot->CondensationSource(i1, i2, nTemp[k], moments[k]) / moments[k][i];
	  if (fSoot->WithSurfaceGrowth())
	    delta[ieq] -= fSoot->SurfaceGrowthSource(i1, i2, moments[k], nY[k], nTemp[k], rho[k], fSpecies->GetMolarMass()->vec) / moments[k][i];
	  if (fSoot->WithSurfaceOxidation())
	    delta[ieq] -= fSoot->OxidationSource(i1, i2, moments[k], nY[k], nTemp[k], rho[k], fSpecies->GetMolarMass()->vec) / moments[k][i];
	}
	else
	{
	  if (fSoot->WithNucleation())
	    delta[ieq] -= fSoot->NucleationSource(0.0, 0.0, nTemp[k]) / moments[k][i];
	  if (fSoot->WithCoagulation())
	    delta[ieq] -= fSoot->CoagulationSourceSmall(nTemp[k], moments[k]) / moments[k][i];
	  if (fSoot->WithCondensation())
	    delta[ieq] -= fSoot->CondensationSourceSmall(nTemp[k], moments[k]) / moments[k][i];
	  if (fSoot->WithSurfaceGrowth())
	    delta[ieq] -= fSoot->SurfaceGrowthSourceSmall(moments[k], nY[k], nTemp[k], rho[k], fSpecies->GetMolarMass()->vec) / moments[k][i];
	  if (fSoot->WithSurfaceOxidation())
	    delta[ieq] -= fSoot->OxidationSourceSmall(moments[k], nY[k], nTemp[k], rho[k], fSpecies->GetMolarMass()->vec) / moments[k][i];
	}
#endif



//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
#ifndef NEWCONVECTION
#  ifdef CENTRALGRID
	delta[ieq] += - yPrime[gpGridOff] * ( fFDWMinus[k] * y[ieq-fNOfEquations] + fFDWCurr[k] * y[ieq] + fFDWPlus[k] * y[ieq+fNOfEquations] );
#  else
	if ( ( -yPrime[gpGridOff] ) > 0.0 )
	  delta[ieq] += - ( y[ieq] - y[ieq-fNOfEquations] ) / fhm[k] * yPrime[gpGridOff];
	else
	  delta[ieq] += - ( y[ieq+fNOfEquations] - y[ieq] ) / fh[k] * yPrime[gpGridOff];
#  endif
#endif
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
      }
    }

// Grid
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
#ifdef MOVEGRID
		if ( k == 0 ) {
			delta[gpGridOff] = - 2.0 * yPrime[gpGridOff] + yPrime[gpGridOff+fNOfEquations]; 
		}
		else if ( k == 1 ) {
			delta[gpGridOff] = fTauGrid * (
				  ( - kappa_ * gpDens[k-1] * gpDens[k-1] / fMonFct[k]
			  	+ ( - kappa_ * gpDens[k-2] * gpDens[k-2] 
				  - ( 1.0 + 2.0 * kappa_ ) * gpDens[k-1] * gpDens[k-1] ) / fMonFct[k-1] ) * yPrime[gpGridOff-fNOfEquations]
			  + ( ( kappa_ * gpDens[k-1] * gpDens[k-1]
    			  + ( 1.0 + 2.0 * kappa_ ) * gpDens[k] * gpDens[k] ) / fMonFct[k] 
				+ ( kappa_ * gpDens[k] * gpDens[k] 
	  			+ ( 1.0 + 2.0 * kappa_ ) * gpDens[k-1] * gpDens[k-1] ) / fMonFct[k-1] ) * yPrime[gpGridOff]
			  + ( - kappa_ * gpDens[k] * gpDens[k] / fMonFct[k-1]
			  	+ ( - kappa_ * gpDens[k+1] * gpDens[k+1] 
				  - ( 1.0 + 2.0 * kappa_ ) * gpDens[k] * gpDens[k] ) / fMonFct[k] ) * yPrime[gpGridOff+fNOfEquations]
			  + ( kappa_ * gpDens[k+1] * gpDens[k+1] / fMonFct[k] ) * yPrime[gpGridOff+2*fNOfEquations]
			  )
			  + ( - kappa_ * gpDens[k-1] + ( 1.0 + 2.0 * kappa_ ) * gpDens[k] - kappa_ * gpDens[k+1] ) / fMonFct[k]
			  - ( - kappa_ * gpDens[k-2] + ( 1.0 + 2.0 * kappa_ ) * gpDens[k-1] - kappa_ * gpDens[k] ) / fMonFct[k-1];
		}	
		else if ( k == fNGridPoints-1 ) {
			delta[gpGridOff] = - 2.0 * yPrime[gpGridOff] + yPrime[gpGridOff-fNOfEquations];
		}
	 	else if ( k == fNGridPoints-2 ) {
			delta[gpGridOff] = fTauGrid * (
				( kappa_ * gpDens[k-2] * gpDens[k-2] / fMonFct[k-1] ) * yPrime[gpGridOff-2*fNOfEquations]
			  + ( - kappa_ * gpDens[k-1] * gpDens[k-1] / fMonFct[k]
			  	+ ( - kappa_ * gpDens[k-2] * gpDens[k-2] 
				  - ( 1.0 + 2.0 * kappa_ ) * gpDens[k-1] * gpDens[k-1] ) / fMonFct[k-1] ) * yPrime[gpGridOff-fNOfEquations]
			  + ( ( kappa_ * gpDens[k-1] * gpDens[k-1]
			      + ( 1.0 + 2.0 * kappa_ ) * gpDens[k] * gpDens[k] ) / fMonFct[k] 
				+ ( kappa_ * gpDens[k] * gpDens[k] 
				  + ( 1.0 + 2.0 * kappa_ ) * gpDens[k-1] * gpDens[k-1] ) / fMonFct[k-1] ) * yPrime[gpGridOff]
			  + ( - kappa_ * gpDens[k] * gpDens[k] / fMonFct[k-1]
			  	+ ( - kappa_ * gpDens[k+1] * gpDens[k+1] 
				  - ( 1.0 + 2.0 * kappa_ ) * gpDens[k] * gpDens[k] ) / fMonFct[k] ) * yPrime[gpGridOff+fNOfEquations]
			  )
			  + ( - kappa_ * gpDens[k-1] + ( 1.0 + 2.0 * kappa_ ) * gpDens[k] - kappa_ * gpDens[k+1] ) / fMonFct[k]
			  - ( - kappa_ * gpDens[k-2] + ( 1.0 + 2.0 * kappa_ ) * gpDens[k-1] - kappa_ * gpDens[k] ) / fMonFct[k-1];	
		}	
		else if ( k > 1 && k < fNGridPoints-2 ) {
			delta[gpGridOff] = fTauGrid * (
				( kappa_ * gpDens[k-2] * gpDens[k-2] / fMonFct[k-1] ) * yPrime[gpGridOff-2*fNOfEquations]
			  + ( - kappa_ * gpDens[k-1] * gpDens[k-1] / fMonFct[k]
			  	+ ( - kappa_ * gpDens[k-2] * gpDens[k-2] 
	    		  - ( 1.0 + 2.0 * kappa_ ) * gpDens[k-1] * gpDens[k-1] ) / fMonFct[k-1] ) * yPrime[gpGridOff-fNOfEquations]
  			+ ( ( kappa_ * gpDens[k-1] * gpDens[k-1]
  			    + ( 1.0 + 2.0 * kappa_ ) * gpDens[k] * gpDens[k] ) / fMonFct[k] 
				+ ( kappa_ * gpDens[k] * gpDens[k] 
			  + ( 1.0 + 2.0 * kappa_ ) * gpDens[k-1] * gpDens[k-1] ) / fMonFct[k-1] ) * yPrime[gpGridOff]
			  + ( - kappa_ * gpDens[k] * gpDens[k] / fMonFct[k-1]
			  	+ ( - kappa_ * gpDens[k+1] * gpDens[k+1] 
				  - ( 1.0 + 2.0 * kappa_ ) * gpDens[k] * gpDens[k] ) / fMonFct[k] ) * yPrime[gpGridOff+fNOfEquations]
			  + ( kappa_ * gpDens[k+1] * gpDens[k+1] / fMonFct[k] ) * yPrime[gpGridOff+2*fNOfEquations]
 			 )
 			 + ( - kappa_ * gpDens[k-1] + ( 1.0 + 2.0 * kappa_ ) * gpDens[k] - kappa_ * gpDens[k+1] ) / fMonFct[k]
			  - ( - kappa_ * gpDens[k-2] + ( 1.0 + 2.0 * kappa_ ) * gpDens[k-1] - kappa_ * gpDens[k] ) / fMonFct[k-1];	
  		}
#else
		delta[gpGridOff] = yPrime[gpGridOff];
#endif
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo			 
		//	fprintf( fOutFilePtr, "%s\t%20.10E\n", fVariableNames[fTemperature], temp );
#ifdef DEBUGRES
		fprintf( stderr, "gp = %d\n", k );
		for ( int j = 0; j < nSpeciesIn+fVariablesWithoutSpecies; ++j ) {
			fprintf( stderr, "%s\t%g\t%g\n", fVariableNames[j]
									, y[gpOff + j]
									, delta[gpOff + j] );
		}
		fprintf( stderr, "\n" );
/*		for ( j = nSpeciesIn; j < fSpecies->GetNOfSpecies(); ++j ) {*/
/*			fprintf( stderr, "SteadyState %s:\t%g\n", fSpecies->GetNames()[j], nY[kCurr][j] );*/
/*		}*/
/*		fprintf( stderr, "\n" );*/
#endif
	}
// boundary conditions
	for ( k = -1; k <= fNGridPoints; k+=fNGridPoints+1 ) {
		gpOff = k*fNOfEquations;
		gpSpecOff = gpOff + fFirstSpecies;
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
    	gpGridOff = gpOff + fGrid;	
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo

//	species
		sum = 0.0;
		sumYH = 0.0;
		for ( i = 0; i < nSpeciesIn; ++i ) {
			delta[gpSpecOff+i] = yPrime[gpSpecOff+i];
			delta[gpSpecOff+i] -= mDoti[k][i] / rho[k];
			sum += hi[k][i] * mDoti[k][i];
#	ifdef TOTENT
			sumYH += hi[k][i] * nY[k][i];
#	endif
		}

//	temperature
#	ifdef TOTENT
//		Double	totentstart = Interpol( Z[k] , fTotEntStart->vec[kPrev], Z[kPrev] , fTotEntStart->vec[fNGridPoints], Z[fNGridPoints] );
//		Double	totentend = Interpol( Z[k] , fTotEntEnd->vec[kPrev], Z[kPrev] , fTotEntEnd->vec[fNGridPoints], Z[fNGridPoints] );
//		Double	totEntPresc = Interpol( *t, totentstart, fTStart, totentend, fTEnd );
//		delta[gpOff+fTemperature] = totEntPresc - sumYH;
		delta[gpOff+fTemperature] = GetTotEnt( k, Z, *t ) - sumYH;
#	else
		delta[gpOff+fTemperature] = yPrime[gpOff+fTemperature];
		delta[gpOff+fTemperature] += ( sum - fDPdt ) / ( rho[k] * cpMix[k] ) 
				- GetTempSource( y[gpOff+fGrid] * GetZR() );
		if ( fProperties->GetRadiation() ) {
			// wird zweimal aufgerufen, koennte aber in TFlame.cp auskommentiert werden
			fProperties->GetRadiation()->SetRadiation( nTemp[k], nY[k], molarMass, rho[k] ); 
			delta[gpOff+fTemperature] -= ( fProperties->GetRadiation()->GetRadiation() ) 
												/ ( rho[k] * cpMix[k] );
			if ( fSoot ) {
#		ifdef NOSOOTRAD
#		else
			  delta[gpOff+fTemperature] += ( fSoot->GetSootRadiation( nTemp[k], moments[k] ) ) 
												/ ( rho[k] * cpMix[k] );
#		endif
			}
		}
#	endif

//	soot moments
		if ( fSoot ) {
			gpSpecOff = gpOff + sootOff;
			for ( i = 0; i < nSootMoments; ++i ) {
				ieq = gpSpecOff + i;
				delta[ieq] = yPrime[ieq];
				//fSoot->ComputePolymereConcs( nY[k], nTemp[k], rho[k], fSpecies->GetMolarMass()->vec
				//		, fSoot->GetPij()->mat, fSoot->GetSumPi()->vec, fSoot->GetPAHMoments()->vec
				//		, moments[k], fReaction );
				fSoot->UpdateSoot( fReaction, fSpecies, moments[k], nTemp[k], nY[k], rho[k]
									, mixMolarMass[k] );
				if ( fSoot->WithNucleation() ) {
					delta[ieq] -= fSoot->NucleationNew( i, nTemp[k]
							, fSoot->GetPAHMoments()->vec ) / rho[k];
				}
				if ( fSoot->WithCoagulation() ) {
					delta[ieq] -= fSoot->SourceCoagulationNew( i, nTemp[k]
										, moments[k] ) / rho[k];
				}
				if ( fSoot->WithCondensation() ) {
					delta[ieq] -= fSoot->SourceCondensationNew( i, nTemp[k]
									, fSoot->GetPAHMoments()->vec, moments[k],  nY[k], rho[k], molarMass ) 
									/ rho[k];
				}
				if ( fSoot->WithSurfaceGrowth() ) {
					delta[ieq] -= fSoot->SourceSurfGrowthNew( i, moments[k], nY[k], rho[k], molarMass ) 
									/ rho[k];
				}
				if ( fSoot->WithSurfaceOxidation() ) {
					delta[ieq] -= fSoot->SourceSootOxidationNew( i, moments[k], nY[k], rho[k], molarMass )
									/ rho[k];
				}
			}
		}
	
//	grid equation
		delta[gpOff+fGrid] = yPrime[gpOff+fGrid];
	}	
#ifdef MOVEZRIGHT
//	Double	ZR = Interpol( GetCurrentTime(), fZRStart, fTStart, fZREnd, fTEnd );
	Double	ZR = GetZR();
#	ifdef CENTRALZRCONV
// start of central part
	for ( k = 0; k < fNGridPoints; ++k ) {
		gpOff = k * fNOfEquations;
		for ( i = 0; i < fNOfEquations; ++i ) {
			if ( i != fGrid ) {
				ieq = gpOff+i;
#		ifdef NEWCONVECTION
				Double	theConvVelo = /*-y[gpOff+fGrid] / ZR * fdDeltaZdt*/ -yPrime[gpOff + fGrid];
				delta[ieq] += theConvVelo * ( 
						fFDWMinus[k] * y[ieq-fNOfEquations]
						+ fFDWCurr[k] * y[ieq]
						+ fFDWPlus[k] * y[ieq+fNOfEquations] );
#		else
				delta[ieq] -= y[gpOff+fGrid] / ZR * fdDeltaZdt
							* ( fFDWMinus[k] * y[ieq-fNOfEquations]
								+fFDWCurr[k] * y[ieq]
								+fFDWPlus[k] * y[ieq+fNOfEquations] );
#		endif
			}
		}
	}
	gpOff = fNGridPoints * fNOfEquations;
	for ( i = 0; i < fNOfEquations; ++i ) {
		if ( i != fGrid ) {
			ieq = gpOff+i;
#		ifdef NEWCONVECTION
			Double	theConvVelo = /*-y[gpOff+fGrid] / ZR * fdDeltaZdt */-yPrime[gpOff + fGrid];
			delta[ieq] += theConvVelo * ( y[ieq] - y[ieq-fNOfEquations] ) / fhm[k];
#		else
			delta[ieq] -= y[gpOff+fGrid] / ZR * fdDeltaZdt
						* ( y[ieq] - y[ieq-fNOfEquations] ) / ( fhm[k] );
#		endif
		}
	}
// end of central part
// start of upwind part
	for ( k = 0; k <= fNGridPoints; ++k ) {
		gpOff = k * fNOfEquations;
		for ( i = 0; i < fNOfEquations; ++i ) {
			if ( i != fGrid ) {
				ieq = gpOff+i;
#ifdef NEWCONVECTION
				Double	theConvVelo = -y[gpOff+fGrid] / ZR * fdDeltaZdt/* -yPrime[gpOff + fGrid]*/;
				if ( theConvVelo > 0.0 || k == fNGridPoints ) {
					delta[ieq] += theConvVelo * ( y[ieq] - y[ieq-fNOfEquations] ) / fhm[k];
				}
				else {
					delta[ieq] += theConvVelo * ( y[ieq+fNOfEquations] - y[ieq] ) / fh[k];
				}	
#else
				delta[ieq] -= y[gpOff+fGrid] / ZR * fdDeltaZdt
							* ( y[ieq] - y[ieq-fNOfEquations] ) / ( fhm[k] );
#endif
			}
		}
	}
// end of upwind part
#	else
#		ifdef SECONDORDUPZRCONV
/*	gpOff = 0;*/
/*	for ( i = 0; i < fNOfEquations; ++i ) {*/
/*		if ( i != fGrid ) {*/
/*			ieq = gpOff+i;*/
/*			delta[ieq] -= y[gpOff+fGrid] / ZR * fdDeltaZdt*/
/*						* ( y[ieq] - y[ieq-fNOfEquations] ) / ( fhm[0] );*/
/*		}*/
/*	}*/
/**/
/*	Double	h1, h2;*/
/*	for ( k = 1; k <= fNGridPoints; ++k ) {*/
/*		gpOff = k * fNOfEquations;*/
/*		h1 = fhm[k];*/
/*		h2 = fhm[k]+fhm[k-1];*/
/*		for ( i = 0; i < fNOfEquations; ++i ) {*/
/*			if ( i != fGrid ) {*/
/*				ieq = gpOff+i;*/
/*				delta[ieq] -= y[gpOff+fGrid] / ZR * fdDeltaZdt*/
/*							* ( ( h2 * h2 - h1 * h1 ) * y[ieq] */
/*								- h2 * h2 * y[ieq-fNOfEquations] */
/*								+ h1 * h1 * y[ieq-2*fNOfEquations] ) / ( h1 * h2 * ( h2 - h1 ) );*/
/*			}*/
/*		}*/
/*	}*/

	gpOff = 0;
	for ( i = 0; i < fNOfEquations; ++i ) {
		if ( i != fGrid ) {
			ieq = gpOff+i;
			delta[ieq] -= 1.0 / ZR * fdDeltaZdt
						* ( ( y[gpOff+fGrid] * y[ieq] - y[gpOff+fGrid-fNOfEquations] * y[ieq-fNOfEquations] ) 
									/ ( fhm[0] ) - y[ieq] );
		}
	}

	Double	h1, h2;
	for ( k = 1; k <= fNGridPoints; ++k ) {
		gpOff = k * fNOfEquations;
		h1 = fhm[k];
		h2 = fhm[k]+fhm[k-1];
		for ( i = 0; i < fNOfEquations; ++i ) {
			if ( i != fGrid ) {
				ieq = gpOff+i;
				delta[ieq] -= 1.0 / ZR * fdDeltaZdt
					* ( ( ( h2 * h2 - h1 * h1 ) * y[gpOff+fGrid] * y[ieq] 
						- h2 * h2 * y[gpOff+fGrid-fNOfEquations] * y[ieq-fNOfEquations] 
						+ h1 * h1 * y[gpOff+fGrid-2*fNOfEquations] * y[ieq-2*fNOfEquations] ) 
						/ ( h1 * h2 * ( h2 - h1 ) ) - y[ieq] );
			}
		}
	}

/*	Double	h1, h2, h1ph2, h1sq, h1ph2sq;*/
/*	Double	hnenn;*/
/*	Double	wCurr;*/
/*	Double	wPrev;*/
/*	Double	wPrevPrev;*/
/*	for ( k = 1; k <= fNGridPoints; ++k ) {*/
/*		gpOff = k * fNOfEquations;*/
/*		h1 = fhm[k];*/
/*		h2 = y[gpOff-fNOfEquations+fGrid] - y[gpOff-2*fNOfEquations+fGrid];*/
/*		if (h2 != fhm[k-1]) {*/
/*			fprintf( stderr, "h2 != fhm[k-1]@%d\n", k );*/
/*			exit(2);*/
/*		}*/
/*		h1ph2 = h1 + h2;*/
/*		hnenn = h1 * h2 * h1ph2;*/
/*		h1sq = h1 * h1;*/
/*		h1ph2sq = h1ph2 * h1ph2;*/
/*		wCurr = ( h1ph2sq - h1sq ) / hnenn;*/
/*		wPrev = -h1ph2sq / hnenn;*/
/*		wPrevPrev = h1sq / hnenn;*/
/*		*/
/*		for ( i = 0; i < fNOfEquations; ++i ) {*/
/*			if ( i != fGrid ) {*/
/*				ieq = gpOff+i;*/
/*				delta[ieq] -= y[gpOff+fGrid] / ZR * fdDeltaZdt*/
/*							* ( wCurr * y[ieq] */
/*							+ wPrev * y[ieq-fNOfEquations] */
/*							+ wPrevPrev * y[ieq-2*fNOfEquations] );*/
/*			}*/
/*		}*/
/*	}*/
#		else
	for ( k = 0; k <= fNGridPoints; ++k ) {
		gpOff = k * fNOfEquations;
		for ( i = 0; i < fNOfEquations; ++i ) {
			if ( i != fGrid ) {
				ieq = gpOff+i;
#ifdef NEWCONVECTION
				Double	theConvVelo = -y[gpOff+fGrid] / ZR * fdDeltaZdt -yPrime[gpOff + fGrid];
				if ( theConvVelo > 0.0 || k == fNGridPoints ) {
					delta[ieq] += theConvVelo * ( y[ieq] - y[ieq-fNOfEquations] ) / fhm[k];
				}
				else {
					delta[ieq] += theConvVelo * ( y[ieq+fNOfEquations] - y[ieq] ) / fh[k];
				}	
#else
				delta[ieq] -= y[gpOff+fGrid] / ZR * fdDeltaZdt
							* ( y[ieq] - y[ieq-fNOfEquations] ) / ( fhm[k] );
/*				delta[ieq] -= 1.0 / ZR * fdDeltaZdt*/
/*							* ( ( y[gpOff+fGrid] * y[ieq] - y[gpOff+fGrid-fNOfEquations] */
/*										* y[ieq-fNOfEquations] ) / ( fhm[k] ) - y[ieq] );*/

/*				delta[ieq] -= 1.0 / ZR * fdDeltaZdt*/
/*							* ( fFDWMinus[k] * y[gpOff+fGrid-fNOfEquations] * y[ieq-fNOfEquations]*/
/*								+fFDWCurr[k] * y[gpOff+fGrid] * y[ieq]*/
/*								+fFDWPlus[k] * y[gpOff+fGrid+fNOfEquations] * y[ieq+fNOfEquations] - y[ieq] );*/

/*				Double fHatPlus, fHatMinus;*/
/*				Double fPlus = y[gpOff+fGrid+fNOfEquations] * y[ieq+fNOfEquations];*/
/*				Double fCurr = y[gpOff+fGrid] * y[ieq];*/
/*				Double fMinus = y[gpOff+fGrid-fNOfEquations] * y[ieq-fNOfEquations];*/
/*				if ( k != 0 && ( fCurr - fMinus ) / fhm[k] < ( fPlus - fCurr ) / fh[k] ) {*/
/*					// negative curvature */
/*					Double fMinMin = y[gpOff+fGrid-fNOfEquations-fNOfEquations] * y[ieq-fNOfEquations-fNOfEquations];*/
/*					fHatPlus = fCurr + 0.5 * fh[k] / fhm[k] * ( fCurr - fMinus );*/
/*					fHatMinus = fMinus + 0.5 * fhm[k] / fhm[k-1] * ( fMinus - fMinMin );*/
/*				}*/
/*				else {*/
/*					// positive curvature */
/*					fHatPlus = 0.5 * ( fPlus + fCurr );*/
/*					fHatMinus = 0.5 * ( fCurr + fMinus );*/
/*				}*/
/*				delta[ieq] -= 1.0 / ZR * fdDeltaZdt*/
/*						* ( 2.0 * ( fHatPlus - fHatMinus ) / ( fh[k] + fhm[k] ) - y[ieq] );*/
#endif
			}
		}
	}
#		endif
#	endif
#endif

#ifdef READU			
	for ( kU = -1; kU <= fNGridPoints; ++kU ) {
		ieqU = kU * fNOfEquations;
		UstOverU = GetU( *t, Z[kU] );
		for ( iU = 0; iU < fNOfEquations; ++iU ) {
			yPrime[ieqU+iU] *= UstOverU;
			delta[ieqU+iU] *= UstOverU;
		}
	}
#endif

#ifdef LEWISCHANGE
	delete Le;
#endif
}
#endif


int ResTransFlameImpliSolverCV( realtype T, N_Vector u, N_Vector udot, void *data )
{
	TTransFlameSolverPtr flame = (TTransFlameSolverPtr) data;
	int		idum, iRes = 0;
	int		i, nEq = ( flame->fNGridPoints + 2 ) * flame->fNOfEquations;
	Double	ddum;
	
	Double	*y = NV_DATA_S(u);
	Double	*delta = NV_DATA_S(udot);

	Double	*yPrime = flame->fSolPrime->mat[kPrev];

	for ( i = 0; i < nEq; ++i ) {
		yPrime[i] = 0.0;
	}

	flame->ResTransFlameImpliSolver( &T, y, yPrime, delta, &idum, &ddum, &idum );

	for ( i = 0; i < nEq; ++i ) {
		delta[i] = -delta[i];
	}

	if (iRes < 0 ) {
		return 1; // means recoverable error
	}
	else {
		return 0; // successful
	}
}


void JacTransFlameSolver( Double *T, Double *y, Double *yPrime, Double *pd
			, Double *cj, Double *rPar, int *iPar )
{
	TTransFlameSolverPtr	flame = ( TTransFlameSolverPtr ) iPar;
	int						nEq = flame->fNOfEquations;
	MatrixPtr	pdNew = FortranToCMat( pd, nEq, nEq, nEq, nEq );
	( ( TTransFlameSolverPtr ) iPar )->JacTransFlameSolver( T, y, yPrime, pdNew->mat
			, *cj, rPar, iPar );
	DisposeFToCMat( pdNew );
}

void TTransFlameSolver::JacTransFlameSolver( Double */*T*/, Double *y, Double */*yPrime*/
			, Double **pd, Double cj, Double */*rPar*/, int */*iPar*/ )
{
	int						nSpeciesIn = fSpecies->GetNSpeciesInSystem();
	int						i, j, ieq, jeq, kAct = fActualPoint;
	Double					sum = 0.0;	

	Double					*prodRate = fSpecies->GetProductionRate()->vec;
	Double					*enth = fSpecies->GetEnthalpy()->vec;
	Double					*cp = fSpecies->GetHeatCapacity()->vec;
	Double					*molarMass = fSpecies->GetMolarMass()->vec;
	Double					*Le = fSpecies->GetLewisNumber()->vec;

	Double					temp = y[fTemperature];
	Double					*YF = &y[fFirstSpecies];

	Double					*nTime = &fSolTime->vec[kAct];
	Double					Z = fSolGrid->vec[kAct];
	Double					chi = GetDissRate( nTime[kCurr], Z );
	Double					pressure = GetPressure( nTime[kCurr] );
	Double					**dMdY = fDmDy;
	Double					*dMdT = fDmDy[nSpeciesIn];
							// dMdY[i][j] = dM_j/dY_i

	if ( fSoot ) {
		fprintf( fOutFilePtr, "###error: analytical evaluation of jacobian for soot not yet implemented\n" );
		exit( 2 );
	}

	T0DFlame::UpdateThermoProps( YF, temp, pressure, fProperties->GetDensityRef()
									, kDensFromPress, NULL );

	Double				rho = fProperties->GetDensity();
	Double				heatCap = fProperties->GetMixHeatCapacity();
	Double				mixMolarMass = fProperties->GetMixMolarMass();

	fReaction->FilldMdYOnePointAnal( dMdY, YF, fReaction->GetReactionRate()->vec
						, mixMolarMass, molarMass, fReaction->GetTBConc()->vec
						, prodRate, rho );
	fReaction->FilldMdTOnePointAnal( dMdT, temp, fReaction->GetReactionRate()->vec
						, fReaction->GetRateCoefficients()->vec, pressure, molarMass
						, fReaction->GetTBConc()->vec );

//	fill dG_i/dY_j
	for ( i = 0; i < nSpeciesIn; ++i ) {
		ieq = fFirstSpecies + i;
		pd[ieq][ieq] += cj;
#ifdef SEMIIMPLI
		pd[ieq][ieq] -= chi / ( 2.0 * Le[i] ) * fWCurr[kAct];
#endif
		pd[ieq][fTemperature] -= ( dMdT[i] + prodRate[i] / temp ) / rho;
		for ( j = 0; j < nSpeciesIn; ++j ) {
			jeq = fFirstSpecies + j;
			pd[ieq][jeq] -= ( dMdY[j][i] + prodRate[i] * mixMolarMass / molarMass[j] ) / rho;
		}
		sum += enth[i] * prodRate[i];
	}

	Double	oneOverRhoCp = 1.0 / ( rho * heatCap );
	
	pd[fTemperature][fTemperature] += cj;
#ifdef SEMIIMPLI
	pd[fTemperature][fTemperature] -= 0.5 * chi * fWCurr[kAct];
#endif
	pd[fTemperature][fTemperature] +=
			sum / temp * oneOverRhoCp;
	for ( i = 0; i < nSpeciesIn; ++i ) {
		pd[fTemperature][fTemperature] +=
				( enth[i] * dMdT[i] + cp[i] * prodRate[i] ) * oneOverRhoCp;
		pd[fTemperature][fFirstSpecies+i] +=
				sum * mixMolarMass / molarMass[i] * oneOverRhoCp;
		for ( j = 0; j < nSpeciesIn; ++j ) {
			jeq = fFirstSpecies + j;
			pd[fTemperature][jeq] +=
					enth[i] * dMdY[j][i] * oneOverRhoCp;
		}
	}

//	fprintf( fOutFilePtr, "%s\t%20.10E\n", fVariableNames[fTemperature], temp );
#ifdef DEBUGRES
	for ( j = 0; j < nSpeciesIn+fVariablesWithoutSpecies; ++j ) {
		fprintf( stderr, "%s\t%g\n", fVariableNames[j], pd[fTemperature][j] );
	}
	fprintf( stderr, "\n" );
#endif
	
}


Double TTransFlameSolver::GetPressure( Double t )
{
	return fPressStart + ( fPressEnd - fPressStart ) / MAX( 1e-30, fTEnd - fTStart ) 
					* ( t - fTStart );
//	return fPressStart + ( fPressEnd - fPressStart ) / ( fTEnd - fTStart ) 
//					* ( t - fTStart );
}

Double TTransFlameSolver::GetTempSource( Double Z )
{
#ifdef DPDT
	return 0.0;
#else
	return fDTdtOx + Z * ( fDTdtFu - fDTdtOx );
#endif
}

Double TTransFlameSolver::Interpol( Double t, Double valOld, Double tOld, Double valNew, Double tNew )
{
	if ( valNew == valOld ) {
		return valNew;
	}
	return valOld + ( valNew - valOld ) / ( tNew - tOld ) * ( t - tOld );
}

#ifdef READCHI
#	ifdef TIMEDEPCHI
Double TTransFlameSolver::GetDissRate( Double t, Double Z )
{
	int	i = 1;
	int	j = 1;
	int	zlen = fZIn->len;
	int	tlen = fTimeIn->len;
	Double	*zIn = fZIn->vec;
	Double	*tIn = fTimeIn->vec;
	Double	**chiIn = fChiIn->mat;
	Double	chiAtZt1, chiAtZt2;
//fprintf( stderr, "Z = %g\n", Z );
	while ( j < zlen-1 && Z > zIn[j] ) ++j;
	if ( j == zlen-1 && Z > zIn[j]+1.0e-8 ) {
		fprintf( stderr, "Z = %g out of range[%g,%g]\n", Z, zIn[0]
								, zIn[zlen-1] );
	}
	
	while ( i < tlen-1 && t > tIn[i] ) ++i;
	if ( i == tlen-1 && t > tIn[i] ) {
		fprintf( stderr, "t = %g out of range[%g,%g]: linear extrapolation\n", Z, tIn[0]
								, tIn[tlen-1] );
	}
	
	if ( zIn[j] == zIn[j-1] ) {
		fprintf( stderr, "something's wrong\n" );
	}
	else {
		chiAtZt1 = Interpol( Z, chiIn[i-1][j-1], zIn[j-1], chiIn[i-1][j], zIn[j] );
		chiAtZt2 = Interpol( Z, chiIn[i][j-1], zIn[j-1], chiIn[i][j], zIn[j] );
	}
	
	if ( tIn[i] == tIn[i-1] ) {
		fprintf( stderr, "something else's wrong\n" );
		return 0.0;
	}
	else {
		return CHIFACT * Interpol( t, chiAtZt1, tIn[i-1], chiAtZt2, tIn[i] );
	}
}
#	else
Double TTransFlameSolver::GetDissRate( Double /*t*/, Double z )
{
	int	i = 1;
	int	len = fZCount->len;
	Double	*zCount = fZCount->vec;
	Double	*chiCount = fChiCount->vec;

	while ( i < len && z > zCount[i] ) ++i;
	if ( i == len && z > zCount[len-1] ) {
		fprintf( stderr, "z = %g out of range[%g,%g]\n", z, zCount[0]
								, zCount[len-1] );
	}
	
	if ( zCount[i] == zCount[i-1] ) {
		return 0.0;
	}
	else {
		return chiCount[i-1] + ( chiCount[i] - chiCount[i-1] ) / ( zCount[i] - zCount[i-1] ) 
											* ( z - zCount[i-1] );
	}
}
#	endif
#else
#	ifdef EXACTCHI
#		ifdef MOVEZRIGHT
Double TTransFlameSolver::GetDissRate( Double t, Double z )
{
	if ( z < 1.0e-10 || z >= 1.0 ) {
		return 0.0;
	}

	Double	zetaRef = fZRef / Interpol( t, fZRStart, fTStart, fZREnd, fTEnd );

	return GetRefDissRate( t ) * z * z / ( zetaRef * zetaRef ) * log( z ) / log( zetaRef );
}
#		else
Double TTransFlameSolver::GetDissRate( Double t, Double z )
{
//   if (z<1.0e-10)
//     return 0.0;
//   else
//     return GetRefDissRate(t) * z*z*log(z)/fZRef/fZRef/log(fZRef);
	return GetRefDissRate( t ) * ExactChi( z ) / ExactChi( fZRef );
}
#		endif
#	else
Double TTransFlameSolver::GetDissRate( Double t, Double z )
{
	Double	zRefStar = ( fZRef - fZl ) / ( fZr - fZl );
	Double	zStar = ( z - fZl ) / ( fZr - fZl );
	Double	zRefNew = zRefStar - 0.5;
	Double	chi = GetRefDissRate( t );
	Double	zRefNew2 = zRefNew * zRefNew;
	Double	zRefNew4 = zRefNew2 * zRefNew2;
	Double	zRefNew6 = zRefNew4 * zRefNew2;
	Double	zRefNew8 = zRefNew6 * zRefNew2;
	Double	zNew = zStar - 0.5;
	Double	zNew2 = zNew * zNew;
	Double	zNew4 = zNew2 * zNew2;
	Double	zNew6 = zNew4 * zNew2;
	Double	zNew8 = zNew6 * zNew2;
	
#		ifdef NEWCHI
	if ( zStar > 0.5 ) {
		zStar = 1.0 - zStar;
	}
	Double	chiNow, fofZst;
	if ( zRefStar <= 0.0547 ) {
		fofZst = 8.0 * zRefStar * zRefStar;
	}
	else if ( zRefStar <= 0.2811 ) {
		fofZst = ( 0.8752 * zRefStar - 0.02394 );
	}
	else {
		fofZst = -2.0 * ( zRefStar - 0.5 ) * ( zRefStar - 0.5 ) + 0.318;
	}

	if ( zStar <= 0.0547 ) {
		chiNow = chi * 8.0 * zStar * zStar / fofZst;
	}
	else if ( zStar <= 0.2811 ) {
		chiNow = chi * ( 0.8752 * zStar - 0.02394 ) / fofZst;
	}
	else {
		chiNow = chi * ( -2.0 * ( zStar - 0.5 ) * ( zStar - 0.5 ) + 0.318 ) / fofZst;
	}
#		else
	Double	duedxst = chi / ( 1.00539 * ( 12.9041 - 
					82.123 * zRefNew2 +
                    115.29 * zRefNew4 -
                    201.898 * zRefNew6 +
                    912.136 * zRefNew8 ) );
	Double	chiNow = 1.00539       *  duedxst *
                    (12.9041       - 
                    82.123  * zNew2 +
                    115.29  * zNew4 -
                    201.898 * zNew6 +
                    912.136 * zNew8 );
				
#		endif					
	chiNow = MAX( chiNow, 0.0 );
	
	return chiNow;
}
#	endif
#endif

Double TTransFlameSolver::GetRefDissRate( Double t )
{
	return CHIFACT * fRandomNumber * Interpol( t, fChiStart, fTStart, fChiEnd, fTEnd );
}

Double TTransFlameSolver::GetDissRate( Double t, Double z, Double rho )
{
	Double	fRhoRef = 0.0;
	fprintf( stderr, "error in GetDissRate\n" );
	exit(2);
	return DissRateFact( t, rho ) / DissRateFact( t, fRhoRef ) * GetDissRate( t, z );
}

Double TTransFlameSolver::ExactChi( Double Z )
{
	int				i;
	const int		nSteps = 100;
	Double			twoZ = 2.0 * Z;
	Double			deltax = 4.0 / ( ( Double ) nSteps ), twoErfc;
	static Double	erFunc[4*nSteps];
	static Flag		init = FALSE;
	
	if ( !init ) {
		init = TRUE;
		for ( i = 0; i < 4 * nSteps; ++i ) {
			erFunc[i] = erfc( i*deltax );
		}
	}
	
	if ( Z > 0.5 ) {
		twoZ = 2.0 * ( 1.0 - Z );
	}
	if ( Z < 1.0e-7 ) {
		return 0.0;
	}
	
	i = 0;
	while ( erFunc[i] > twoZ ) ++i;
	
	twoErfc = Interpol( twoZ, (i-1)*deltax, erFunc[i-1], i*deltax, erFunc[i] );
	
	return exp( -2.0 * twoErfc * twoErfc );
}

Double TTransFlameSolver::DissRateFact( Double t, Double rho )
{
	Double rhoInf = fWOverRInf * GetPressure( t ) / Interpol( t, fSolOldTemp->vec[kPrev], fSolOldTime->vec[kPrev]
									, fSolTemp->vec[kPrev], fSolTime->vec[kPrev] );
	Double	densRatio = sqrt( rhoInf / rho );

	return 1.5 * ( densRatio + 1 ) * ( densRatio + 1 ) / ( 2.0 * densRatio + 1 );
}

Double TTransFlameSolver::GetZStoi( void )
{ 
	int				i;
	int				nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();
	int				*indexFuel = GetFuelIndexVec()->vec;
	int				indexOx;
	TInputDataPtr	input = GetInputData();
	Double			zStoi;
	Double			nuOx = fInputData->fOxIndex;
	Double			nuFuel;
	Double			*molarMass = fSpecies->GetMolarMass()->vec;
	Double			nu = 0.0;
	Double			sumYFuel = 0.0;
	Double			*Y2 = &fSolution->mat[kPrev][fFirstSpecies];
	Double			*Y1 = &fSolution->mat[fNGridPoints][fFirstSpecies];
	char			**names = fSpecies->GetNames();
	
	indexOx = fInputData->fOxIndex;
	
	nuOx = GetNu( fInputData->fGlobalReaction, names[indexOx] );

	if ( nuOx == -1 ) {
		fprintf( fOutFilePtr, "###error: oxidizer '%s' is not in the global reaction\n", names[indexOx] );
		exit(2);
	}
	
	for ( i = 0; i < GetNFuels(); ++i ) {
		nuFuel = GetNu( input->fGlobalReaction, names[indexFuel[i]] );
		if ( nuFuel == -1 ) {
			fprintf( fOutFilePtr, "###error: fuel '%s' is not in the global reaction\n", names[indexFuel[i]] );
			exit(2);
		}
		nu += nuFuel * molarMass[indexFuel[i]];
		sumYFuel += Y1[indexFuel[i]];
	}

	nu = nuOx * molarMass[indexOx] / nu;

	zStoi = 1.0 / ( 1.0 + nu * sumYFuel / Y2[indexOx] - Y1[indexOx] / Y2[indexOx] );

	if ( zStoi < 1.0e-10 || zStoi > 0.9 ) {
		fprintf( fOutFilePtr, "#warning: the value of ZStoi = %g indicates an error in the inital conditons\n"
				, zStoi );
	}
	else {
		fprintf( fOutFilePtr, "ZStoi = %g\n"
				, zStoi );
	}
	
	return zStoi;
}

Double TTransFlameSolver::GetDelQ( Double t, Double Z )
{
	int	i = 1;
	int	j = 1;
	int	zlen = fZIn->len;
	int	tlen = fTimeIn->len;
	Double	*zIn = fZIn->vec;
	Double	*tIn = fTimeIn->vec;
	Double	**delQIn = fDelQ->mat;
	Double	delQAtZt1, delQAtZt2;
//fprintf( stderr, "Z = %g\n", Z );
	while ( j < zlen-1 && Z > zIn[j] ) ++j;
	if ( j == zlen-1 && Z > zIn[j]+1.0e-8 ) {
		fprintf( stderr, "Z = %g out of range[%g,%g]\n", Z, zIn[0]
								, zIn[zlen-1] );
	}
	
	while ( i < tlen-1 && t > tIn[i] ) ++i;
	if ( i == tlen-1 && t > tIn[i] ) {
		fprintf( stderr, "t = %g out of range[%g,%g]: linear extrapolation\n", Z, tIn[0]
								, tIn[tlen-1] );
	}
	
	if ( zIn[j] == zIn[j-1] ) {
		fprintf( stderr, "something's wrong\n" );
	}
	else {
		delQAtZt1 = Interpol( Z, delQIn[i-1][j-1], zIn[j-1], delQIn[i-1][j], zIn[j] );
		delQAtZt2 = Interpol( Z, delQIn[i][j-1], zIn[j-1], delQIn[i][j], zIn[j] );
	}
	
	if ( tIn[i] == tIn[i-1] ) {
		fprintf( stderr, "something else's wrong\n" );
		return 0.0;
	}
	else {
		return Interpol( t, delQAtZt1, tIn[i-1], delQAtZt2, tIn[i] );
	}
}

int	TTransFlameSolver::GetVariableIndex( const char *name )
{
	for ( int i = 0; i < fNOfEquations; ++ i ) {
		if ( strcmp( name, fVariableNames[i] ) == 0 ) {
			return i;
		}
	}

	return -1;
}

int	TTransFlameSolver::GetVariableIndex( const char *name, ConstStringArray array, int len )
{
	for ( int i = 0; i < len; ++ i ) {
		if ( strcmp( name, array[i] ) == 0 ) {
			return i;
		}
	}

	return -1;
}

void WRITEFLAMELETFILE( int **object )
{
	writeflameletfile( object );
}

void writeflameletfile( int **object )
{
	TTransFlameSolverPtr flame = *( TTransFlameSolverPtr *)object;
	
	flame->WriteFlameletFile( NULL, NULL, NULL );
}

#undef WRITEGLOBALPRODRATE

void TTransFlameSolver::WriteFlameletFile( FILE *fp, char *head, char *tail )
{
	T0DPropertiesPtr	props = GetProperties();
	T0DSpeciesPtr		species = GetSpecies();
	Double				*molarMass = species->GetMolarMass()->vec;
	int					i, j, k, l;
	int					nOfSpecies = species->GetNOfSpecies();
	int					nSpeciesIn = species->GetNSpeciesInSystem();
	int					nSootMoments;
	time_t				theDate;
	char				buffer[80];
	char				**names = species->GetNames();
	Flag				fpOpen = FALSE;
	Double				*x = fSolGrid->vec;
	Double				**Y = fMassFracsWork->mat;
	Double				**theMom;
	Double				*temp = fTempWork->vec;
	double * prog = fProgWork->vec;
	double * enthun = fEnthWork->vec;
	Double				theTime = GetCurrentTime();
	char				tmpName[127];
	VectorPtr 			ZBilgerVec = NewVector( fNGridPoints + 2 );
	Double				*ZBilger = &ZBilgerVec->vec[kNext];
	VectorPtr ZBilgerSrcVec = NewVector(fNGridPoints + 2);
	double * ZBilgerSrc = &ZBilgerSrcVec->vec[kNext];
	VectorPtr HCRatVec = NewVector(fNGridPoints + 2);
	double * HCRat = &HCRatVec->vec[kNext];
	VectorPtr ProgSrcVec = NewVector(fNGridPoints + 2);
	double * ProgSrc = &ProgSrcVec->vec[kNext];
	VectorPtr DimerRateVec = NewVector(fNGridPoints+2);
	double * DimerRate = &DimerRateVec->vec[kNext];
	VectorPtr DimerNbrCVec = NewVector(fNGridPoints+2);
	double * DimerNbrC = &DimerNbrCVec->vec[kNext];
	VectorPtr DimerNbrHVec = NewVector(fNGridPoints+2);
	double * DimerNbrH = &DimerNbrHVec->vec[kNext];
	VectorPtr SootSurfCoeffVec = NewVector(fNGridPoints+2);
	double * SootSurfCoeff = &SootSurfCoeffVec->vec[kNext];
	VectorPtr SootOxCoeffVec = NewVector(fNGridPoints+2);
	double * SootOxCoeff = &SootOxCoeffVec->vec[kNext];
	VectorPtr SootOxCoeff_O2Vec = NewVector(fNGridPoints+2);
	double * SootOxCoeff_O2 = &SootOxCoeff_O2Vec->vec[kNext];
	VectorPtr sqrtTVec = NewVector(fNGridPoints+2);
	double * sqrtT = &sqrtTVec->vec[kNext];
	VectorPtr T_MUVec = NewVector(fNGridPoints+2);
	double * T_MU = &T_MUVec->vec[kNext];
	VectorPtr MUsqrtW_RHOsqrtTVec = NewVector(fNGridPoints+2);
	double * MUsqrtW_RHOsqrtT = &MUsqrtW_RHOsqrtTVec->vec[kNext];
// 	VectorPtr EnthFluxVec = NewVector(fNGridPoints + 2);
// 	double * EnthFlux = &EnthFluxVec->vec[kNext];
// 	double enthprev, enthcurr, enthnext;

	static const Double	molarMassC = 12.01, 
						molarMassO = 16.0,
						molarMassH = 1.008;
	static Double		elemMassLastCInit = 10000.0;
	static Double		elemMassLastHInit = 10000.0;
	static Double		elemMassLastOInit = 10000.0;

	double * nrhodot = rhodot->vec;
	double * nenthdot = enthdot->vec;

	if ( fSoot ) {
		nSootMoments = fSoot->GetNSootMoments();
		theMom = fSootMomentsWork->mat;
	}

	SetOutSolution();

	if ( fZREnd >= 0.9999 ) {
		elemMassLastCInit = GetElementMassFraction( Y[fNGridPoints], "C", molarMassC );
		elemMassLastHInit = GetElementMassFraction( Y[fNGridPoints], "H", molarMassH );
		elemMassLastOInit = GetElementMassFraction( Y[fNGridPoints], "O", molarMassO );
	}

/*	fprintf( stderr, "ZC1 = %g\n", elemMassLastCInit );*/
/*	fprintf( stderr, "ZH1 = %g\n", elemMassLastHInit );*/
/*	fprintf( stderr, "ZO1 = %g\n", elemMassLastOInit );*/

	if ( !fp ) {
		fpOpen = TRUE;
		fp = GetOutputFile( theTime, head, tail, TFlame::kText );
	}

	// write header
	fprintf( fp, "header\n\n" );

	fprintf( fp, "title = \"transient flamelet\"\n" );
	fprintf( fp, "mechanism = \"%s\"\n", fInputData->fReactionFile );
	fprintf( fp, "author = \"%s\"\n", GetAuthor() );
	time( &theDate );
	strcpy( buffer, ctime( &theDate ) );
	if ( buffer[strlen(buffer)-1] == '\n' )
		buffer[strlen(buffer)-1] = '\0';
	fprintf( fp, "date = \"%s\"\n\n", buffer );
	for ( i = 0; i < GetNFuels(); ++i ) {
		fprintf( fp, "fuel = \"%s\"\n", fVariableNames[fFirstSpecies+GetFuelIndex( i )] );
	}
	fprintf( fp, "time = %g [ms]\n", theTime * 1.0e3 );
	fprintf( fp, "pressure = %g [bar]\n", GetPressure( theTime ) / 1.0e5 );
	fprintf( fp, "chi_ref = %g [1/s]\n", GetRefDissRate( theTime ) );


	fprintf( fp, "Tmax = %g [K]\n", temp[LocationOfMax( fNGridPoints+2
												, &temp[kPrev] )-1] );
	fprintf( fp, "ZofTmax = %g [K]\n", x[LocationOfMax( fNGridPoints+2
												, &temp[kPrev] )-1] );
	Double	zR = Interpol( theTime, fZRStart, fTStart, fZREnd, fTEnd );
	fprintf( fp, "ZR = %g\n", zR );

	int gp = fNGridPoints;
	while ( gp >= 0 && temp[gp-1] > temp[gp] ) --gp;
	fprintf( fp, "LocTmax = %g [K]\n", temp[gp] );

	Double	EIFuel = 0.0;
	for ( i = 0; i < GetNFuels(); ++i ) {
		EIFuel += ComputeEmissionIndex( theTime, GetFuelIndex( i ), NULL, NULL );
	}
	fprintf( fp, "EmissionIndexFuel = %g [kg/m^3s]\n", EIFuel );
	int	index = species->FindSpecies( "NO" );
	if ( index > -1 ) {
		fprintf( fp, "EmissionIndexNO = %g [kg/sm^3]\n"
				, ComputeEmissionIndex( theTime, index, NULL, NULL ) / EIFuel );
	}
	
	index = species->FindSpecies( "NO2" );
	if ( index > -1 && !species->IsSteadyState(index)) {
		fprintf( fp, "EmissionIndexNO2 = %g [kg/sm^3]\n"
				, ComputeEmissionIndex( theTime, index, NULL, NULL ) / EIFuel );
	}
	
	index = species->FindSpecies( "N2O" );
	if ( index > -1 && !species->IsSteadyState(index)) {
		fprintf( fp, "EmissionIndexN2O = %g [kg/sm^3]\n"
				, ComputeEmissionIndex( theTime, index, NULL, NULL ) );
	}

	if ( fSoot ) {
		fprintf( fp, "EmissionIndexSoot = %g [kg/sm^3]\n"
				, ComputeEmissionIndexSoot( theTime, 1, NULL, NULL ) * 24 );
	}

	fprintf( fp, "CurrTimeStep = %g [s]\n", *fActualTimeStepSize );

	fprintf( fp, "\nFuelSide\n" );
	fprintf( fp, "begin\n" );
	fprintf( fp, "\tTemperature = %g [K]\n", fTempFuelEnd );
	for ( i = 0; i < nSpeciesIn; ++i ) {
		if ( fabs( Y[fNGridPoints][i] ) > 1.0e-10 ) {
			fprintf( fp, "\tMassfraction-%s = %g\n", names[i], Y[fNGridPoints][i] );
		}
	}
	fprintf( fp, "end\n\n" );

	fprintf( fp, "OxidizerSide\n" );
	fprintf( fp, "begin\n" );
	fprintf( fp, "\tTemperature = %g [K]\n", fTempOxEnd );
	for ( i = 0; i < nSpeciesIn; ++i ) {
		if ( fabs( Y[kPrev][i] ) > 1.0e-10 ) {
			fprintf( fp, "\tMassfraction-%s = %g\n", names[i], Y[kPrev][i] );
		}
	}
	fprintf( fp, "end\n\n" );

	fprintf( fp, "numOfSpecies = %d\n", nOfSpecies );
	fprintf( fp, "gridPoints = %d\n\n", fNGridPoints+2 );

	fprintf( fp, "body\n" );

	// write solution
	Double	*mixfrac = New1DArray( fNGridPoints+2 );

	for ( k = -1; k <= fNGridPoints; ++k ) {
		mixfrac[k+1] = x[k] * zR;
		if ( mixfrac[k+1] < 1.0e-10 ) {
			mixfrac[k+1] = 0.0;
		}
	}

	PrintFlameletVector( fNGridPoints+2, mixfrac, "Z", fp );
	PrintFlameletVector( fNGridPoints+2, &x[kPrev], "zeta", fp );

	Free1DArray( mixfrac );

	PrintFlameletVector( fNGridPoints+2, &temp[kPrev], "temperature [K]", fp );

	// write massfractions of species
	for ( i = 0; i < nOfSpecies; ++i ) {
		sprintf( tmpName, "massfraction-%s", names[i] );
		PrintFlameletVector( fNGridPoints+2, &Y[kPrev][i], tmpName
				, fp, fMassFracsWork->phys_rows );
	}

	if ( fPrintMolarFractions ) {

		MatrixPtr	moleMat = NewMatrix( nOfSpecies, fNGridPoints+2, kColumnPointers );
		Double		**molefracs = &moleMat->mat[kNext];
		for ( k = -1; k <= fNGridPoints; ++k )
		{
			fProperties->ComputeMixtureMolarMass( props->GetMixMolarMassRef()
					, Y[k], molarMass, nSpeciesIn );
			for ( i = 0; i < nOfSpecies; ++i ) {
				molefracs[k][i] = Y[k][i] * props->GetMixMolarMass() / molarMass[i];
			}
		}
		for ( i = 0; i < nOfSpecies; ++i ) {
			sprintf( tmpName, "molarfraction-%s", names[i] );
			PrintFlameletVector( fNGridPoints+2, &molefracs[kPrev][i], tmpName
					, fp, moleMat->phys_rows );
		}
		DisposeMatrix( moleMat );
	}
	
	// write progvar 
	int	CO2 = fSpecies->FindSpecies( "CO2" );
	int	CO = fSpecies->FindSpecies( "CO" );
	int	H2O = fSpecies->FindSpecies( "H2O" );
	int	H2 = fSpecies->FindSpecies( "H2" );
	Double	*prog_old = New1DArray( fNGridPoints+2 );
	sprintf( tmpName, "ProgVar" );
	for ( k = 0; k <= fNGridPoints+1; ++k ) {
		prog_old[k] = Y[k-1][CO2] + Y[k-1][CO] + Y[k-1][H2O] + Y[k-1][H2];
	}
	PrintFlameletVector( fNGridPoints+2, prog_old, "ProgVar", fp );
	Free1DArray( prog_old );

// element mixture fraction
	Double	*ZElem = New1DArray( fNGridPoints+2 );
	Double	elemMassFirst;
//	ZElem[0] = 0.0;
//	ZElem[fNGridPoints+1] = 1.0;
	
	
// C
	elemMassFirst = GetElementMassFraction( Y[kPrev], "C", molarMassC );
	for ( k = -1; k <= fNGridPoints; ++k ) {
		ZElem[k+1] = ( GetElementMassFraction( Y[k], "C", molarMassC ) - elemMassFirst ) / MAX( elemMassLastCInit - elemMassFirst, 1.0e-30 );
	}
	PrintFlameletVector( fNGridPoints+2, ZElem, "CMassFrac", fp );
	
// H
	elemMassFirst = GetElementMassFraction( Y[kPrev], "H", molarMassH );
	for ( k = -1; k <= fNGridPoints; ++k ) {
		ZElem[k+1] = ( GetElementMassFraction( Y[k], "H", molarMassH ) - elemMassFirst ) / MAX( elemMassLastHInit - elemMassFirst, 1.0e-30 );
	}
	PrintFlameletVector( fNGridPoints+2, ZElem, "HMassFrac", fp );

// O
	elemMassFirst = GetElementMassFraction( Y[kPrev], "O", molarMassO );
	for ( k = -1; k <= fNGridPoints; ++k ) {
		ZElem[k+1] = ( GetElementMassFraction( Y[k], "O", molarMassO ) - elemMassFirst ) / MAX( elemMassLastOInit - elemMassFirst, 1.0e-30 );
	}
	PrintFlameletVector( fNGridPoints+2, ZElem, "OMassFrac", fp );

	Free1DArray( ZElem );

// Mixture fraction following Barlows definition
	Double	*ZBarlow = New1DArray( fNGridPoints+2 );
	Double	elemMassFirstH = GetElementMassFraction( Y[kPrev], "H", molarMassH );
	Double	elemMassFirstC = GetElementMassFraction( Y[kPrev], "C", molarMassC );
	Double	denomBarlow = 0.5 * ( elemMassLastHInit - elemMassFirstH ) / molarMassH
				+ 2.0 * ( elemMassLastCInit - elemMassFirstC ) / molarMassC;
	for ( k = -1; k <= fNGridPoints; ++k ) {
		ZBarlow[k+1] = ( 0.5 * ( GetElementMassFraction( Y[k], "H", molarMassH ) - elemMassFirstH ) / molarMassH
						+ 2.0 * ( GetElementMassFraction( Y[k], "C", molarMassC ) - elemMassFirstC ) / molarMassC )
						/ denomBarlow;
	}

	PrintFlameletVector( fNGridPoints+2, ZBarlow, "ZBarlow", fp );

	Free1DArray( ZBarlow );


// density
	Double	*density = New1DArray( fNGridPoints+2 );
	Double	*locchi = New1DArray( fNGridPoints+2 );
//	Double	*locchi = New1DArray( fNGridPoints+2 );
	for ( k = -1; k <= fNGridPoints; ++k ) {
		fProperties->ComputeMixtureMolarMass( props->GetMixMolarMassRef()
				, Y[k], molarMass, nSpeciesIn );
		density[k+1] = GetPressure( theTime ) * props->GetMixMolarMass() / ( RGAS * temp[k] );
		locchi[k+1] = GetDissRate( theTime, x[k] );
	}

	PrintFlameletVector( fNGridPoints+2, density, "density [kg/m^3]", fp );
	PrintFlameletVector( fNGridPoints+2, locchi, "chi [1/s]", fp );

	Double	*dummy = New1DArray( fNGridPoints+2 );

	for ( k = -1; k <= fNGridPoints; ++k ) {
		dummy[k+1] = k+1;
	}

//	PrintFlameletVector( fNGridPoints+2, dummy, "GridNum", fp );

	Free1DArray( dummy );

		Double	*sootrad;
		Double	*gasrad;
		double *sootradcoeff;
		if ( fProperties->GetRadiation() ) {
			sootrad = New1DArray( fNGridPoints+2 );
			gasrad = New1DArray( fNGridPoints+2 );
			sootradcoeff = New1DArray( fNGridPoints+2 );
		}
		Double	*summh = New1DArray( fNGridPoints+2 );
		Double	*prodRate = fSpecies->GetProductionRate()->vec;
		Double	**prodRateGrid = fProdRate->mat;
		Double	*enth = fSpecies->GetEnthalpy()->vec;
		Double	press = GetPressure( theTime );
		Double	*cpMix = fHeatCpMix->vec;
		Double	*mu = fViscosity->vec;
		Double	*mixMolarMass = fMolarMassMix->vec;
		Double	*lambdaOverCp = fLambdaOverCpMix->vec;
		for ( k = -1; k <= fNGridPoints; ++k )
		{
#ifdef DELTAINEW
			UpdateThermoProps( k, Y[k], temp[k], press, fProperties->GetDensityRef()
					   , kDensFromPress, (fSoot) ? theMom[k] : NULL, &nrhodot[k], &nenthdot[k] );
#else 
			T0DFlame::UpdateThermoProps( Y[k], temp[k], press, fProperties->GetDensityRef()
						     , kDensFromPress, (fSoot) ? theMom[k] : NULL, &nrhodot[k], &nenthdot[k] );
#endif
		
			summh[k+1] = 0.0;
			for ( i = 0; i < nSpeciesIn; ++i ) {
				summh[k+1] += enth[i] * prodRate[i];
				prodRateGrid[k][i] = prodRate[i];
			}
			cpMix[k] = fProperties->GetMixHeatCapacity();
			mu[k] = fProperties->GetMixViscosity();
			lambdaOverCp[k] = fProperties->GetMixConductivity() / fProperties->GetMixHeatCapacity();
			mixMolarMass[k] = fProperties->GetMixMolarMass();

			if ( fProperties->GetRadiation() ) {
				gasrad[k+1] = fProperties->GetRadiation()->GetRadiation();
				if ( fSoot ) {
#ifdef ROSSELAND
				  if (k==-1||k==fNGridPoints) {
					sootrad[k+1] = 0.0;
				  }
				  else {
					Double rfM = density[k] * GetDissRate( theTime, x[k-1] ) * lambdaOverCp[k-1];
					Double rfC = density[k+1] * GetDissRate( theTime, x[k] ) * lambdaOverCp[k];
					Double rfP = density[k+2] * GetDissRate( theTime, x[k+1] ) * lambdaOverCp[k+1];

					sootrad[k+1] = -GetRosseRadiation( k, &temp[k], &theMom[k], rfM, rfC, rfP );
				  }
#else

					sootrad[k+1] = -fSoot->GetSootRadiation( temp[k], theMom[k] );
					sootradcoeff[k+1] = fSoot->GetSootRadiationCoeff( temp[k] );
#endif
				}
			}
		}
	
		if ( fProperties->GetRadiation() ) {
#ifdef NOSOOTRAD
#else
			if ( fSoot ) {
				PrintFlameletVector( fNGridPoints+2, sootrad, "SootRadiation [J/m^3s]", fp );
			}
#endif
			PrintFlameletVector( fNGridPoints+2, gasrad, "GasRadiation [J/m^3s]", fp );
			PrintFlameletVector( fNGridPoints+2, sootradcoeff, "SootRadCoeff [J/m^3s]", fp );
		}
		PrintFlameletVector( fNGridPoints+2, summh, "SumMH [J/m^3s]", fp );
		PrintFlameletVector( fNGridPoints+2, &cpMix[kPrev], "Cp", fp );
		PrintFlameletVector( fNGridPoints+2, &mu[kPrev], "mu", fp );
		PrintFlameletVector( fNGridPoints+2, &lambdaOverCp[kPrev], "lambdaOverCp", fp );
		PrintFlameletVector( fNGridPoints+2, &mixMolarMass[kPrev], "W", fp );
	
	
//  write source term ProdVar
	fprintf( fp, "ProdRateProgVar [kg/m^3s]\n" );
	fprintf( fp, "\t%-.6e", 0.0 );
	for ( k = 0; k < fNGridPoints; ++k ) {
		fprintf( fp, "\t%-.6e", prodRateGrid[k][CO2] + prodRateGrid[k][CO] + prodRateGrid[k][H2O] + prodRateGrid[k][H2] );
		if ( (k+2) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	fprintf( fp, "\t%-.6e\n", 0.0 );
	
//  write source term CO2
	fprintf( fp, "ProdRateCO2 [kg/m^3s]\n" );
	fprintf( fp, "\t%-.6e", 0.0 );
	for ( k = 0; k < fNGridPoints; ++k ) {
		fprintf( fp, "\t%-.6e", prodRateGrid[k][CO2] );
		if ( (k+2) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	fprintf( fp, "\t%-.6e\n", 0.0 );


//  write source term H2O
	fprintf( fp, "ProdRateH2O [kg/m^3s]\n" );
	fprintf( fp, "\t%-.6e", 0.0 );
	for ( k = 0; k < fNGridPoints; ++k ) {
		fprintf( fp, "\t%-.6e", prodRateGrid[k][H2O] );
		if ( (k+2) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	fprintf( fp, "\t%-.6e\n", 0.0 );


//  write source term H2
	fprintf( fp, "ProdRateH2 [kg/m^3s]\n" );
	fprintf( fp, "\t%-.6e", 0.0 );
	for ( k = 0; k < fNGridPoints; ++k ) {
		fprintf( fp, "\t%-.6e", prodRateGrid[k][H2] );
		if ( (k+2) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	fprintf( fp, "\t%-.6e\n", 0.0 );


//  write source term CO
	fprintf( fp, "ProdRateCO [kg/m^3s]\n" );
	fprintf( fp, "\t%-.6e", 0.0 );
	for ( k = 0; k < fNGridPoints; ++k ) {
		fprintf( fp, "\t%-.6e", prodRateGrid[k][CO] );
		if ( (k+2) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	fprintf( fp, "\t%-.6e\n", 0.0 );

	ZBilger[kPrev] = x[kPrev];
	ZBilger[fNGridPoints] = x[fNGridPoints];
	for ( k = 0; k < fNGridPoints; ++k ) {
		ZBilger[k] = ComputeZBilger( Y[k], Y[fNGridPoints], Y[kPrev] );
	}

	ZBilgerSrc[kPrev] = 0.0;
	ZBilgerSrc[fNGridPoints] = 0.0;
	for (k = 0; k < fNGridPoints; ++k)
	{
	  if (fSoot)
	    ZBilgerSrc[k] = ComputeZBilgerSource(prodRateGrid[k], Y[fNGridPoints], Y[-1], nrhodot[k]);
	  else
	    ZBilgerSrc[k] = ComputeZBilgerSource(prodRateGrid[k], Y[fNGridPoints], Y[-1], 0.0);
	}

	HCRat[kPrev] = 0.0;
	for (k = 0; k <= fNGridPoints; ++k)
	{
	  HCRat[k] = ComputeHC(Y[k], Y[fNGridPoints], Y[-1]);
	}

	double C_MAX;
	ProgSrc[kPrev] = 0.0;
	ProgSrc[fNGridPoints] = 0.0;
	for (k = 0; k < fNGridPoints; ++k)
	{
	  C_MAX = ComputeCMAX(Y[k], Y[fNGridPoints], Y[-1]);
	  ProgSrc[k] = (prodRateGrid[k][CO2] + prodRateGrid[k][CO] + prodRateGrid[k][H2O] + prodRateGrid[k][H2]) / C_MAX;
	}

// 	EnthFlux[kPrev] = 0.0;
// 	EnthFlux[fNGridPoints] = 0.0;
// 	for (k = 0; k < fNGridPoints; ++k)
// 	{
// 	  EnthFlux[k] = 0.0;
// 	  for (i = 0; i < nOfSpecies; ++i)
// 	  {
// 	    fSpecies->ComputeSpeciesProperties(temp[k-1]);
// 	    enthprev = fSpecies->GetEnthalpy()->vec[i];
// 	    fSpecies->ComputeSpeciesProperties(temp[k  ]);
// 	    enthcurr = fSpecies->GetEnthalpy()->vec[i];
// 	    fSpecies->ComputeSpeciesProperties(temp[k+1]);
// 	    enthnext = fSpecies->GetEnthalpy()->vec[i];
// 	    EnthFlux[k] += density[k+1]*locchi[k+1] / 2.0 * enthcurr * 
// 	      (fWMinus[k]*Y[k-1][i] + fWCurr[k]*Y[k][i] + fWPlus[k]*Y[k+1][i]);
// 	    EnthFlux[k] += density[k+1]*locchi[k+1] / 2.0 * (1.0/lambdaOverCp[k]) * 
// 	      (fFDWMinus[k]*lambdaOverCp[k-1]*enthprev + fFDWCurr[k]*lambdaOverCp[k]*enthcurr + fFDWPlus[k]*lambdaOverCp[k+1]*enthnext) *
// 	      (fFDWMinus[k]*Y[k-1][i] + fFDWCurr[k]*Y[k][i] + fFDWPlus[k]*Y[k+1][i]);
// 	  }
// 	}

	PrintFlameletVector( fNGridPoints+2, &prog[kPrev], "ProgRat", fp);
	PrintFlameletVector( fNGridPoints+2, &ProgSrc[kPrev], "ProgSrc", fp);
// 	for ( int k = -1; k <= fNGridPoints; ++k)
// 	  enthun[k] = enthun[k] * 1.0e20;
	PrintFlameletVector( fNGridPoints+2, &enthun[kPrev], "EnthLoss", fp);
// 	for ( int k = -1; k <= fNGridPoints; ++k)
// 	  enthun[k] = enthun[k] / 1.0e20;
	PrintFlameletVector( fNGridPoints+2, &ZBilger[kPrev], "ZBilger", fp );
	PrintFlameletVector( fNGridPoints+2, &ZBilgerSrc[kPrev], "ZBilgerSrc", fp );
	PrintFlameletVector( fNGridPoints+2, &HCRat[kPrev], "H_C", fp);
// 	PrintFlameletVector( fNGridPoints+2, &EnthFlux[kPrev], "EnthFlux [J/m^3 s]", fp );

		Free1DArray( summh );
		if ( fProperties->GetRadiation() ) {
			Free1DArray( gasrad );
			Free1DArray( sootrad );
			Free1DArray( sootradcoeff );
		}

	if (fSoot)
	{
	  for (i = 0; i < nSootMoments; ++i)
	  {
	    j = fSoot->Geti(i);
#ifdef V
	    sprintf(tmpName, "conc-SOOT%d", j);
#endif
#ifdef VS
	    l = fSoot->Getj(i);
	    sprintf(tmpName, "conc-SOOT%d-%d", j, l);
#endif

	    PrintFlameletVector(fNGridPoints+2, &theMom[kPrev][i], tmpName, fp, fSootMomentsWork->phys_rows);
	  }
		
	  //WriteSootInfo( theTime, temp, Y, GetPressure( theTime ), theMom, fp );		

	  PrintFlameletVector(fNGridPoints+2, &nrhodot[kPrev], "RhoDot [kg/m3-s]", fp);
	  PrintFlameletVector(fNGridPoints+2, &nenthdot[kPrev], "EnthDot [J/m3-s]", fp);

	  //DimerRate[kPrev] = 0.0;
	  //DimerRate[fNGridPoints] = 0.0;
	  for (k = -1; k <= fNGridPoints; ++k)
	  {
	    fSoot->CalcDimerProductionRate(species, GetReaction(), &DimerRate[k], &DimerNbrC[k], &DimerNbrH[k], density[k+1], Y[k], temp[k], theMom[k]);
	    SootSurfCoeff[k] = fSoot->CalcSurfaceGrowthCoeff(GetReaction(), Y[k], temp[k], density[k+1], molarMass, mixMolarMass[k]);
	    SootOxCoeff[k] = fSoot->CalcOxidationCoeff(GetReaction(), Y[k], temp[k], density[k+1], molarMass, mixMolarMass[k]);
	    SootOxCoeff_O2[k] = fSoot->CalcOxidationCoeff_O2(GetReaction(), Y[k], temp[k], density[k+1], molarMass, mixMolarMass[k]);
	    sqrtT[k] = sqrt(temp[k]);
	    T_MU[k] = temp[k] / mu[k];
	    MUsqrtW_RHOsqrtT[k] = mu[k] / density[k+1] * sqrt(mixMolarMass[k]/1000.0 / temp[k]);
	  }
	  PrintFlameletVector(fNGridPoints+2, &DimerRate[kPrev], "Dimer_ProdRate [kmol/m3-s]", fp);
	  PrintFlameletVector(fNGridPoints+2, &DimerNbrC[kPrev], "Dimer_nbrC", fp);
	  PrintFlameletVector(fNGridPoints+2, &DimerNbrH[kPrev], "Dimer_nbrH", fp);
	  PrintFlameletVector(fNGridPoints+2, &SootSurfCoeff[kPrev], "SootSurfGrowthCoeff", fp);
	  PrintFlameletVector(fNGridPoints+2, &SootOxCoeff[kPrev], "SootOxCoeff", fp);
	  PrintFlameletVector(fNGridPoints+2, &SootOxCoeff_O2[kPrev], "SootOxCoeff_O2", fp);
	  PrintFlameletVector(fNGridPoints+2, &sqrtT[kPrev], "sqrtT", fp);
	  PrintFlameletVector(fNGridPoints+2, &T_MU[kPrev], "T_MU", fp);
	  PrintFlameletVector(fNGridPoints+2, &MUsqrtW_RHOsqrtT[kPrev], "MUsqrtW_RHOsqrtT", fp);

	  PrintProdRatePAH(theTime, fp);
	}

	Free1DArray( density );
	Free1DArray( locchi );

	Double	*totEnt = New1DArray( fNGridPoints+2 );
	Double	totentstart, totentend;
	for ( k = -1; k <= fNGridPoints; ++k )
	{
//		totentstart = Interpol( x[k] , fTotEntStart->vec[kPrev], x[kPrev], fTotEntStart->vec[fNGridPoints], x[fNGridPoints] );
//		totentend = Interpol( x[k] , fTotEntEnd->vec[kPrev], x[kPrev] , fTotEntEnd->vec[fNGridPoints], x[fNGridPoints] );
//		totEnt[k+1] = Interpol( theTime, totentstart, fTStart, totentend, fTEnd );
		totEnt[k+1] = GetTotEnt( k, x, theTime );
	}
	PrintFlameletVector( fNGridPoints+2, totEnt, "TotalEnthalpy [J/kg]", fp );

	for ( k = -1; k <= fNGridPoints; ++k )
	{
		totEnt[k+1] = GetTotEnt( temp[k], Y[k] );
	}
	PrintFlameletVector( fNGridPoints+2, totEnt, "TotalEnthalpy2 [J/kg]", fp );

	Free1DArray( totEnt );

	PrintFlameletVector( fNGridPoints+2, &fMonFct[kPrev], "MonitorFunc", fp );

	// Production Rate of NO
	index = species->FindSpecies( "NO" );
	if ( index > -1 ) {
		PrintProdRate( theTime, index, fp );
	}
	// Production Rate of NO
	index = species->FindSpecies( "CH4" );
	if ( index > -1 ) {
		PrintProdRate( theTime, index, fp );
	}
	
#ifdef WRITEGLOBALPRODRATE
	PrintProdRateGlobalReac( theTime );	
#endif
	
	fprintf( fp, "trailer\n" );
	
#ifdef LEWISCHANGE
	Double					*LeOrig = fSpecies->GetLewisNumber()->vec;
	Double					*Le = new Double[nSpeciesIn];
	Double					LeFunc;
	Double					switchTime = LEWISSWITCHTIME;
	Double					switchPeriod = LEWISSWITCHPERIOD;
	const Double			Pi = 4.0 * atan( 1.0 );
	Double					omega = 2.0 * Pi / switchPeriod;
	Double					tStar = theTime - ( switchTime - 0.5 * switchPeriod );

	if ( tStar >= 0.0 ) {
		if ( tStar < switchPeriod ) {
			LeFunc = 0.5 * ( cos( omega * tStar ) + 1.0 );
			for ( int iLe = 0; iLe < nSpeciesIn; ++iLe ) {
				Le[iLe] = LeFunc * ( LeOrig[iLe] - 1.0 ) + 1.0;
			}
		}
		else {
			for ( int iLe = 0; iLe < nSpeciesIn; ++iLe ) {
				Le[iLe] = 1.0;
			}
		}
	}
	else {
		for ( int iLe = 0; iLe < nSpeciesIn; ++iLe ) {
			Le[iLe] = LeOrig[iLe];
		}
	}
#else
	Double					*Le = fSpecies->GetLewisNumber()->vec;
#endif
	for ( i = 0; i < nSpeciesIn; ++i ) {
		fprintf( fp, "%s\t%g\n", names[i], Le[i] );
	}
	if ( fpOpen ) {
		fclose( fp );
	}
#ifdef LEWISCHANGE
	delete Le;
#endif
	DisposeVector( DimerRateVec );
	DisposeVector( DimerNbrHVec );
	DisposeVector( DimerNbrCVec );
	DisposeVector( SootSurfCoeffVec );
	DisposeVector( SootOxCoeffVec );
	DisposeVector( SootOxCoeff_O2Vec );
	DisposeVector( sqrtTVec );
	DisposeVector( T_MUVec );
	DisposeVector( MUsqrtW_RHOsqrtTVec );
	DisposeVector( ZBilgerVec );
	DisposeVector( ZBilgerSrcVec );
	DisposeVector( ProgSrcVec );
	DisposeVector( HCRatVec );
// 	DisposeVector( EnthFluxVec );
}

void TTransFlameSolver::PrintProdRate( Double t, int speciesIndex, FILE *fp )
{
	int			k;
	Double		**Y = fSolMassFracs->mat;
	Double		*T = fSolTemp->vec;
	Double		**M = ( fSoot ) ? fSolSootMoments->mat : NULL;
	double * nrhodot = rhodot->vec;
	double * nenthdot = enthdot->vec;
	Double		*prodRate = fSpecies->GetProductionRate()->vec;
	Double		*reactionRate = fReaction->GetReactionRate()->vec;
	Double		pressure = GetPressure( t );
	Double		*prodOfZ = New1DArray( fNGridPoints+2 );
	Double		*prodPlusOfZ = New1DArray( fNGridPoints+2 );
	Double		*prodMinusOfZ = New1DArray( fNGridPoints+2 );
	Double		*prodMinusOverYOfZ = New1DArray( fNGridPoints+2 );
	char		buffer[128];
	sprintf( buffer, "ProdRate-%s", fSpecies->GetNames()[speciesIndex] );

	prodOfZ[0] = prodOfZ[fNGridPoints+1] = 0.0;
	prodPlusOfZ[0] = prodPlusOfZ[fNGridPoints+1] = 0.0;
	prodMinusOfZ[0] = prodMinusOfZ[fNGridPoints+1] = 0.0;
	prodMinusOverYOfZ[0] = prodMinusOverYOfZ[fNGridPoints+1] = 0.0;
	for ( k = 0; k < fNGridPoints; ++k )
	{
#ifdef DELTAINEW
		UpdateThermoProps( k, Y[k], T[k], pressure, fProperties->GetDensityRef()
				   , kDensFromPress, ( !fSoot ) ? NULL : M[k], &nrhodot[k], &nenthdot[k] );
#else 
		T0DFlame::UpdateThermoProps( Y[k], T[k], pressure, fProperties->GetDensityRef()
										, kDensFromPress, ( !fSoot ) ? NULL : M[k] );
#endif
		prodOfZ[k+1] = prodRate[speciesIndex];
		prodPlusOfZ[k+1] = fSpecies->GetPosNegProductionRate( speciesIndex, reactionRate, TRUE );
		prodMinusOfZ[k+1] = fSpecies->GetPosNegProductionRate( speciesIndex, reactionRate, FALSE );
		prodMinusOverYOfZ[k+1] = prodMinusOfZ[k+1] / MAX(1.0e-10, Y[k][speciesIndex] );
	}
	PrintFlameletVector( fNGridPoints+2, prodOfZ, buffer, fp );
	sprintf( buffer, "ProdRatePos-%s", fSpecies->GetNames()[speciesIndex] );
	PrintFlameletVector( fNGridPoints+2, prodPlusOfZ, buffer, fp );
	sprintf( buffer, "ProdRateNeg-%s", fSpecies->GetNames()[speciesIndex] );
	PrintFlameletVector( fNGridPoints+2, prodMinusOfZ, buffer, fp );
	sprintf( buffer, "ProdRateNegOverYNO-%s", fSpecies->GetNames()[speciesIndex] );
	PrintFlameletVector( fNGridPoints+2, prodMinusOverYOfZ, buffer, fp );
	
	Free1DArray( prodMinusOverYOfZ );
	Free1DArray( prodPlusOfZ );
	Free1DArray( prodMinusOfZ );
	Free1DArray( prodOfZ );
}

void TTransFlameSolver::PrintProdRatePAH(double t, FILE * fp)
{
  int k, s;
  int nPAH = 8;
  int * aromIndex = New1DIntArray(nPAH);
  double * Y_PAH = New1DArray(fNGridPoints+2);
  double * prodOfZ = New1DArray(fNGridPoints+2);
  double * prodPlusOfZ = New1DArray(fNGridPoints+2);
  double * prodMinusOfZ = New1DArray(fNGridPoints+2);
  double * prodMinusOverYOfZ = New1DArray(fNGridPoints+2);
  double ** Y = fSolMassFracs->mat;
  double * T = fSolTemp->vec;
  double ** M = (fSoot) ? fSolSootMoments->mat : NULL;
  double * reactionRate = fReaction->GetReactionRate()->vec;
  double * prodRate = fSpecies->GetProductionRate()->vec;
  double * nrhodot = rhodot->vec;
  double * nenthdot = enthdot->vec;
  double pressure = GetPressure(t);
  char buffer[128];

  Y_PAH[0] = Y_PAH[fNGridPoints+1] = 0.0;
  prodOfZ[0] = prodOfZ[fNGridPoints+1] = 0.0;
  prodPlusOfZ[0] = prodPlusOfZ[fNGridPoints+1] = 0.0;
  prodMinusOfZ[0] = prodMinusOfZ[fNGridPoints+1] = 0.0;
  prodMinusOverYOfZ[0] = prodMinusOverYOfZ[fNGridPoints+1] = 0.0;

  aromIndex[0] = fSpecies->FindSpecies("A2-C10H8");
  aromIndex[1] = fSpecies->FindSpecies("A2R5-C12H8");
  aromIndex[2] = fSpecies->FindSpecies("P2-C12H10");
  aromIndex[3] = fSpecies->FindSpecies("A3-C14H10");
  aromIndex[4] = fSpecies->FindSpecies("A3R5-C16H10");
  aromIndex[5] = fSpecies->FindSpecies("A4-C16H10");
  aromIndex[6] = fSpecies->FindSpecies("FLTN-C16H10");
  aromIndex[7] = fSpecies->FindSpecies("A4R5-C18H10");

  for (s = 0; s < nPAH; ++s)
  {
    Y_PAH[0] += MAX(Y[-1][aromIndex[s]],1.0e-60); //Mueller
    Y_PAH[fNGridPoints+1] += MAX(Y[fNGridPoints][aromIndex[s]],1.0e-60); //Mueller
  }

  for (k = 0; k < fNGridPoints; ++k)
  {
#ifdef DELTAINEW
    UpdateThermoProps(k, Y[k], T[k], pressure, fProperties->GetDensityRef(), kDensFromPress, (!fSoot) ? NULL : M[k], &nrhodot[k], &nenthdot[k]);
#else
    T0DFlame::UpdateThermoProps(Y[k], T[k], pressure, fProperties->GetDensityRef(), kDensFromPress, (!fSoot) ? NULL : M[k]);
#endif
    Y_PAH[k+1] = 0.0;
    for (s = 0; s < nPAH; ++s)
    {
      if (aromIndex[s] > -1)
      {
	prodOfZ[k+1] += prodRate[aromIndex[s]];
	prodPlusOfZ[k+1] += fSpecies->GetPosNegProductionRate(aromIndex[s], reactionRate, TRUE);
	prodMinusOfZ[k+1] += fSpecies->GetPosNegProductionRate(aromIndex[s], reactionRate, FALSE);
	Y_PAH[k+1] += MAX(Y[k][aromIndex[s]],1.0e-60); //Mueller
      }
    }
    prodMinusOverYOfZ[k+1] = prodMinusOfZ[k+1] / MAX(1.0e-60, Y_PAH[k+1]);
  }

  sprintf(buffer, "Y-PAH");
  PrintFlameletVector(fNGridPoints+2, Y_PAH, buffer, fp);
  sprintf(buffer, "ProdRate-PAH [kg/m^3s]");
  PrintFlameletVector(fNGridPoints+2, prodOfZ, buffer, fp);
  sprintf(buffer, "ProdRatePos-PAH [kg/m^3s]");
  PrintFlameletVector(fNGridPoints+2, prodPlusOfZ, buffer, fp);
  sprintf(buffer, "ProdRateNeg-PAH [kg/m^3s]");
  PrintFlameletVector(fNGridPoints+2, prodMinusOfZ, buffer, fp);
  sprintf(buffer, "ProdRateNegOverY-PAH [kg/m^3s]");
  PrintFlameletVector(fNGridPoints+2, prodMinusOverYOfZ, buffer, fp);

  Free1DArray(prodMinusOverYOfZ);
  Free1DArray(prodPlusOfZ);
  Free1DArray(prodMinusOfZ);
  Free1DArray(prodOfZ);
  Free1DArray(Y_PAH);
}

/*void TTransFlameSolver::PrintProdRate( Double t, int speciesIndex )*/
/*{*/
/*	char		buffer[128];*/
/*	sprintf( buffer, "PR_%s", fSpecies->GetNames()[speciesIndex] );*/
/*	FILE		*fp = GetOutputFile( t, buffer, NULL, TFlame::kText );*/
/*	int			k;*/
/*	Double		**Y = fSolMassFracs->mat;*/
/*	Double		*T = fSolTemp->vec;*/
/*	Double		**M = ( fSoot ) ? fSolSootMoments->mat : NULL;*/
/*	Double		*Z = fSolGrid->vec;*/
/*	Double		*prodRate = fSpecies->GetProductionRate()->vec;*/
/*	Double		pressure = GetPressure( t );*/
/**/
/*	fprintf( fp, "*\nZ\t%s\n%g\t%g\n", fSpecies->GetNames()[speciesIndex], Z[-1], 0.0 );*/
/*	for ( k = 0; k < fNGridPoints; ++k ) {*/
/*		UpdateThermoProps( k, Y[k], T[k], pressure, fProperties->GetDensityRef()*/
/*										, kDensFromPress, ( !fSoot ) ? NULL : M[k] );*/
/*		fprintf( fp, "%g\t%g\n", Z[k], prodRate[speciesIndex] );*/
/*	}*/
/*	fclose( fp );*/
/*}*/

Double TTransFlameSolver::ComputeEmissionIndex( Double t, int speciesIndex, Double **PDF, Double *timePDF )
{
	int			k;
	Double		**Y = fSolMassFracs->mat;
	Double		*T = fSolTemp->vec;
	Double		**M = ( fSoot ) ? fSolSootMoments->mat : NULL;
	Double		*Z = fSolGrid->vec;
	Double		*prodRate = fSpecies->GetProductionRate()->vec;
	Double		pressure = GetPressure( t );
	char		**names = fSpecies->GetNames();
	Double		sum;
	int 		indPDF = 0;

	Double	zR = Interpol( t, fZRStart, fTStart, fZREnd, fTEnd );

	if ( PDF ) {
		while ( timePDF[indPDF] < t ) ++indPDF;
	}

	sum = 0.0;  // m[0] = 0, m[L] = 0
	for ( k = 1; k < fNGridPoints; ++k )
	{
#ifdef DELTAINEW
		UpdateThermoProps( k, Y[k], T[k], pressure, fProperties->GetDensityRef()
										, kDensFromPress, ( !fSoot ) ? NULL : M[k] );
#else 
		T0DFlame::UpdateThermoProps( Y[k], T[k], pressure, fProperties->GetDensityRef()
										, kDensFromPress, ( !fSoot ) ? NULL : M[k] );
#endif
		sum += prodRate[speciesIndex] * ( Z[k+1] - Z[k-1] )*zR
		 	* ( ( PDF ) ? Interpol( t, PDF[indPDF-1][k], timePDF[indPDF-1]
			, PDF[indPDF][k], timePDF[indPDF] ) : 1.0 );
	}
			
	return 0.5 * sum;
}

Double TTransFlameSolver::ComputeEmissionIndexSoot( Double t, int which, Double **PDF, Double *timePDF )
{
	int			k;
	Double		**Y = fSolMassFracs->mat;
	Double		*T = fSolTemp->vec;
	Double		**M = ( fSoot ) ? fSolSootMoments->mat : NULL;
	Double		*Z = fSolGrid->vec;
	Double		pressure = GetPressure( t );
	char		**names = fSpecies->GetNames();
	Double		sum;
	int 		indPDF = 0;

	Double	zR = Interpol( t, fZRStart, fTStart, fZREnd, fTEnd );
	
	if ( PDF ) {
		while ( timePDF[indPDF] < t ) ++indPDF;
	}

	sum = 0.0;  // m[0] = 0, m[L] = 0
	for ( k = 1; k < fNGridPoints; ++k )
	{
#ifdef DELTAINEW
		UpdateThermoProps( k, Y[k], T[k], pressure, fProperties->GetDensityRef()
										, kDensFromPress, ( !fSoot ) ? NULL : M[k] );
#else 
		T0DFlame::UpdateThermoProps( Y[k], T[k], pressure, fProperties->GetDensityRef()
										, kDensFromPress, ( !fSoot ) ? NULL : M[k] );
#endif
		sum += fSoot->NucleationNew( which, T[k], fSoot->GetPAHMoments()->vec ) 
				* ( Z[k+1] - Z[k-1] )*zR
				* ( ( PDF ) ? Interpol( t, PDF[indPDF-1][k], timePDF[indPDF-1], PDF[indPDF][k], timePDF[indPDF] ) : 1.0 );
	}
	
	return 0.5 * sum;
}

Double TTransFlameSolver::TurbMeanZBarlow( Double t, Double ZMean, Double ZVar )
{
	int			k;
	Double		zR = Interpol( t, fZRStart, fTStart, fZREnd, fTEnd );
	int			nj;
	Double		*Z; // Z[0] = 0, Z[nj-2] = ZR, Z[nj-1] = 1.0
	Double		*zeta;
	Flag		newZ = FALSE;
	
	if ( zR < 1.0 ) {
		nj = fNGridPoints + 2 + 1;
		Z = New1DArray( nj ); // Z[0] = 0, Z[nj-2] = ZR, Z[nj-1] = 1.0
		newZ = TRUE;
		zeta = fSolGrid->vec;
		for ( k = -1; k <= fNGridPoints; ++k ) {
			Z[k+1] = MAX(0.0,zeta[k] * zR);
		}
		Z[nj-1] = 1.0;
	}
	else {
		nj = fNGridPoints + 2;
		Z = &fSolGrid->vec[kPrev]; // Z[0] = 0, Z[nj-1] = 1.0
		zeta = NULL;
	}


	Double		*pdf = New1DArray( nj );
	Double		sum;
	Double		**Y = fSolMassFracs->mat;

// Mixture fraction following Barlows definition
	static const Double	molarMassC = 12.01, 
//						molarMassO = 16.0,
						molarMassH = 1.008;
	static Double		elemMassLastCInit = 10000.0;
	static Double		elemMassLastHInit = 10000.0;
//	static Double		elemMassLastOInit = 10000.0;
	if ( fZREnd >= 0.9999 ) {
		elemMassLastCInit = GetElementMassFraction( Y[fNGridPoints], "C", molarMassC );
		elemMassLastHInit = GetElementMassFraction( Y[fNGridPoints], "H", molarMassH );
//		elemMassLastOInit = GetElementMassFraction( Y[fNGridPoints], "O", molarMassO );
	}
	Double	*ZBarlow = New1DArray( fNGridPoints+2 );
	ZBarlow = &ZBarlow[kNext];
	Double	elemMassFirstH = GetElementMassFraction( Y[kPrev], "H", molarMassH );
	Double	elemMassFirstC = GetElementMassFraction( Y[kPrev], "C", molarMassC );
	Double	denomBarlow = 0.5 * ( elemMassLastHInit - elemMassFirstH ) / molarMassH
				+ 2.0 * ( elemMassLastCInit - elemMassFirstC ) / molarMassC;
	for ( k = -1; k <= fNGridPoints; ++k ) {
		ZBarlow[k] = ( 0.5 * ( GetElementMassFraction( Y[k], "H", molarMassH ) - elemMassFirstH ) / molarMassH
						+ 2.0 * ( GetElementMassFraction( Y[k], "C", molarMassC ) - elemMassFirstC ) / molarMassC )
						/ denomBarlow;
	}


	ZVar = MAX( 1e-10, ZVar );
	
	int		i, j;
	int		newPoints = ( fNGridPoints + 1 ) * 5 + 1;
	Double	*newGrid = New1DArray( newPoints );
	Double	*newZBarlow = New1DArray( newPoints );
	Double	*newPdf = New1DArray( newPoints );
	Double	*ZBarlowBase0 = &ZBarlow[kPrev];
	
	for ( i = 0; i < fNGridPoints + 1; ++i ) {
		for ( j = 0; j < 5; ++j ) {
			newGrid[i*5+j] = Z[i] + j * 0.2 * ( Z[i+1] - Z[i] );
			newZBarlow[i*5+j] = ZBarlowBase0[i] + j * 0.2 * ( ZBarlowBase0[i+1] - ZBarlowBase0[i] );
		}
	}
    newGrid[0]=0.0;
	newGrid[i*5] = Z[fNGridPoints + 1];
	newZBarlow[i*5] = ZBarlowBase0[fNGridPoints + 1];

#ifdef NEWBETAPDF
	F77BETAPDF( newPdf, newGrid, &ZMean, &ZVar, &newPoints, &newPoints );
#else
	PDFZBETA( &ZMean, &ZVar, newGrid, newPdf, &newPoints );
#endif

	sum = 0.0;
	for ( k = 0; k < newPoints; ++k ) {
		sum += newPdf[k] * newZBarlow[k];
	}
		
	if ( newZ ) {
		Free1DArray( Z );
	}
	Free1DArray( newPdf );
	Free1DArray( newZBarlow );
	Free1DArray( newGrid );
	Free1DArray( pdf );

	ZBarlow = &ZBarlow[kPrev];
	Free1DArray( ZBarlow );
	return sum;
}

Double TTransFlameSolver::TurbMeanTemp( Double t, Double ZMean, Double ZVar, Double *tempVar )
{
	int			k;
	Double		zR = Interpol( t, fZRStart, fTStart, fZREnd, fTEnd );
	int			nj;
	Double		*Z; // Z[0] = 0, Z[nj-2] = ZR, Z[nj-1] = 1.0
	Double		*zeta;
	Flag		newZ = FALSE;
	
	if ( zR < 1.0 ) {
		nj = fNGridPoints + 2 + 1;
		Z = New1DArray( nj ); // Z[0] = 0, Z[nj-2] = ZR, Z[nj-1] = 1.0
		newZ = TRUE;
		zeta = fSolGrid->vec;
		for ( k = -1; k <= fNGridPoints; ++k ) {
			Z[k+1] = MAX(0.0,zeta[k] * zR);
		}
        Z[0]=0.0;
		Z[nj-1] = 1.0;
	}
	else {
		nj = fNGridPoints + 2;
		Z = &fSolGrid->vec[kPrev]; // Z[0] = 0, Z[nj-1] = 1.0
		zeta = NULL;
	}
	Double		*T = fSolTemp->vec;
	Double		*pdf = New1DArray( nj );
	Double		sum;


	ZVar = MAX( 1e-10, ZVar );
	
#ifdef FINEPDF
	int		i, j;
	int		newPoints = ( fNGridPoints + 1 ) * 5 + 1;
	Double	*newGrid = New1DArray( newPoints );
	Double	*newTemp = New1DArray( newPoints );
	Double	*newPdf = New1DArray( newPoints );
	Double	*TBase0 = &T[kPrev];
	
	for ( i = 0; i < fNGridPoints + 1; ++i ) {
		for ( j = 0; j < 5; ++j ) {
			newGrid[i*5+j] = Z[i] + j * 0.2 * ( Z[i+1] - Z[i] );
			newTemp[i*5+j] = TBase0[i] + j * 0.2 * ( TBase0[i+1] - TBase0[i] );
		}
	}
    newGrid[0]=0.0;
	newGrid[i*5] = Z[fNGridPoints + 1];
	newTemp[i*5] = TBase0[fNGridPoints + 1];

#ifdef NEWBETAPDF
	F77BETAPDF( newPdf, newGrid, &ZMean, &ZVar, &newPoints, &newPoints );
#else
	PDFZBETA( &ZMean, &ZVar, newGrid, newPdf, &newPoints );
#endif

#undef DEBUGNEWTURB
#ifdef DEBUGNEWTURB
	FILE *fp = GetOutputFile( t, "NewTurbGrid", NULL, kData );
	
	fprintf( fp, "*\nZ\tT\tPDFZ\n" );
	for ( i = 0; i < newPoints; ++i ) {
		fprintf( fp, "%g\t%g\t%g\n", newGrid[i], newTemp[i], newPdf[i] );		
	}
	
	fclose( fp );
#endif

	sum = 0.0;
	for ( k = 0; k < newPoints; ++k ) {
		sum += newPdf[k] * newTemp[k];
	}
	
	if ( tempVar ) {
		Double	sumVar = 0.0;
		for ( k = 0; k < newPoints; ++k ) {
			sumVar += newPdf[k] * ( newTemp[k] - sum ) * ( newTemp[k] - sum );
		}	
		*tempVar = sumVar;
	}
	
	if ( newZ ) {
		Free1DArray( Z );
	}
	Free1DArray( newPdf );
	Free1DArray( newTemp );
	Free1DArray( newGrid );
	Free1DArray( pdf );

	return sum;
#else

#ifdef NEWBETAPDF
	F77BETAPDF( pdf, Z, &ZMean, &ZVar, &nj, &nj );
#else
	PDFZBETA( &ZMean, &ZVar, Z, pdf, &nj );
#endif
	
	sum = 0.0;
	for ( k = -1; k <= fNGridPoints; ++k ) {
		sum += pdf[k+1] * T[k];
	}
	
	if ( tempVar ) {
		Double	sumVar = 0.0;
		for ( k = -1; k <= fNGridPoints; ++k ) {
			sumVar += pdf[k+1] * ( T[k] - sum ) * ( T[k] - sum );
		}	
		*tempVar = sumVar;
	}
	
	if ( newZ ) {
		Free1DArray( Z );
	}
	Free1DArray( pdf );

	return sum;
#endif
}

Double TTransFlameSolver::TurbMeanTotEnergy( Double t, Double ZMean, Double ZVar )
{
	int			k;
	Double		zR = Interpol( t, fZRStart, fTStart, fZREnd, fTEnd );
	int			nj;
	Double		*Z; // Z[0] = 0, Z[nj-2] = ZR, Z[nj-1] = 1.0
	Double		*zeta;
	Flag		newZ = FALSE;

	if ( zR < 1.0 ) {
		nj = fNGridPoints + 2 + 1;
		Z = New1DArray( nj ); // Z[0] = 0, Z[nj-2] = ZR, Z[nj-1] = 1.0
		newZ = TRUE;
		zeta = fSolGrid->vec;
		for ( k = -1; k <= fNGridPoints; ++k ) {
			Z[k+1] = MAX(0.0,zeta[k] * zR);
		}
		Z[nj-1] = 1.0;
	}
	else {
		nj = fNGridPoints + 2;
		Z = &fSolGrid->vec[kPrev]; // Z[0] = 0, Z[nj-1] = 1.0
		zeta = NULL;
	}
	Double		*T = fSolTemp->vec;
	Double		**Y = fSolMassFracs->mat;
	Double		*pdf = New1DArray( nj );
	Double		sum;
	int 		nSpeciesIn = fSpecies->GetNSpeciesInSystem();
	Double		*molarMass = fSpecies->GetMolarMass()->vec;
	Double		dens;
	
	Double		*totEnt = New1DArray( fNGridPoints+2 );
	totEnt = &totEnt[kNext];

	for ( k = -1; k <= fNGridPoints; ++k )
	{
		fProperties->ComputeMixtureMolarMass( fProperties->GetMixMolarMassRef()
				, Y[k], molarMass, nSpeciesIn );
		dens = GetPressure( t ) * fProperties->GetMixMolarMass() / ( RGAS * T[k] );
		totEnt[k] = dens * ( GetTotEnt( k, fSolGrid->vec, t ) - GetTotEnt( T[k], Y[k] ) );
	}

	ZVar = MAX( 1e-10, ZVar );

#ifdef FINEPDF
	int		i, j;
	int		newPoints = ( fNGridPoints + 1 ) * 5 + 1;
	Double	*newGrid = New1DArray( newPoints );
	Double	*newTotEnt = New1DArray( newPoints );
	Double	*newPdf = New1DArray( newPoints );
	Double	*HBase0 = &totEnt[kPrev];
	
	for ( i = 0; i < fNGridPoints + 1; ++i ) {
		for ( j = 0; j < 5; ++j ) {
			newGrid[i*5+j] = Z[i] + j * 0.2 * ( Z[i+1] - Z[i] );
			newTotEnt[i*5+j] = HBase0[i] + j * 0.2 * ( HBase0[i+1] - HBase0[i] );
		}
	}
    newGrid[0]=0.0;
	newGrid[i*5] = Z[fNGridPoints + 1];
	newTotEnt[i*5] = HBase0[fNGridPoints + 1];

#ifdef NEWBETAPDF
	F77BETAPDF( newPdf, newGrid, &ZMean, &ZVar, &newPoints, &newPoints );
#else
	PDFZBETA( &ZMean, &ZVar, newGrid, newPdf, &newPoints );
#endif

#undef DEBUGNEWTURB
#ifdef DEBUGNEWTURB
	FILE *fp = GetOutputFile( t, "NewTurbGrid", NULL, kData );
	
	fprintf( fp, "*\nZ\tTotEnt\tPDFZ\n" );
	for ( i = 0; i < newPoints; ++i ) {
		fprintf( fp, "%g\t%g\t%g\n", newGrid[i], newTotEnt[i], newPdf[i] );		
	}
	
	fclose( fp );
#endif

	sum = 0.0;
	for ( k = 0; k < newPoints; ++k ) {
		sum += newPdf[k] * newTotEnt[k];
	}
	
	if ( newZ ) {
		Free1DArray( Z );
	}
	
	Free1DArray( newPdf );
	Free1DArray( newTotEnt );
	Free1DArray( newGrid );
	Free1DArray( &totEnt[kPrev] );
	Free1DArray( pdf );

	return sum;
#else

#ifdef NEWBETAPDF
	F77BETAPDF( pdf, Z, &ZMean, &ZVar, &nj, &nj );
#else
	PDFZBETA( &ZMean, &ZVar, Z, pdf, &nj );
#endif
	
	sum = 0.0;
	for ( k = -1; k <= fNGridPoints; ++k ) {
		sum += pdf[k+1] * totEnt[k];
	}

	if ( newZ ) {
		Free1DArray( Z );
	}
	Free1DArray( &totEnt[kPrev] );
	Free1DArray( pdf );

	return sum;
#endif
}

Double TTransFlameSolver::TurbMeanX( Double t, int speciesIndex, Double ZMean, Double ZVar, Double *specVar )
{
	int			k;
	Double		zR = Interpol( t, fZRStart, fTStart, fZREnd, fTEnd );
	int			nj;
	Double		*Z; // Z[0] = 0, Z[nj-2] = ZR, Z[nj-1] = 1.0
	Double		*zeta;
	Flag		newZ = FALSE;

	if ( zR < 1.0 ) {
		nj = fNGridPoints + 2 + 1;
		Z = New1DArray( nj ); // Z[0] = 0, Z[nj-2] = ZR, Z[nj-1] = 1.0
		newZ = TRUE;
		zeta = fSolGrid->vec;
		for ( k = -1; k <= fNGridPoints; ++k ) {
			Z[k+1] = zeta[k] * zR;
		}
		Z[nj-1] = 1.0;
	}
	else {
		nj = fNGridPoints + 2;
		Z = &fSolGrid->vec[kPrev]; // Z[0] = 0, Z[nj-1] = 1.0
		zeta = NULL;
	}
	int 		nSpeciesIn = fSpecies->GetNSpeciesInSystem();
	Double		*molarMass = fSpecies->GetMolarMass()->vec;
	Double		**Y = fSolMassFracs->mat;
	Double		*pdf = New1DArray( nj );
	Double		sum;

	ZVar = MAX( 1e-10, ZVar );

#ifdef NEWBETAPDF
	F77BETAPDF( pdf, Z, &ZMean, &ZVar, &nj, &nj );
#else
	PDFZBETA( &ZMean, &ZVar, Z, pdf, &nj );
#endif
	
	sum = 0.0;
	for ( k = -1; k <= fNGridPoints; ++k ) {
		fProperties->ComputeMixtureMolarMass( fProperties->GetMixMolarMassRef()
				, Y[k], molarMass, nSpeciesIn );
		sum += pdf[k+1] * Y[k][speciesIndex] * fProperties->GetMixMolarMass() / molarMass[speciesIndex];
	}

	if ( specVar ) {
		Double X;
		Double sumVar = 0.0;
		for ( k = -1; k <= fNGridPoints; ++k ) {
			X = Y[k][speciesIndex] * fProperties->GetMixMolarMass() / molarMass[speciesIndex];
			sumVar += pdf[k+1] * ( X - sum ) * ( X - sum );
		}	
		*specVar = sumVar;
	}

	if ( newZ ) {
		Free1DArray( Z );
	}
	Free1DArray( pdf );

	return sum;
}

Double TTransFlameSolver::TurbMeanY( Double t, int speciesIndex, Double ZMean, Double ZVar, Double *specVar )
{
	int			k;
	Double		zR = Interpol( t, fZRStart, fTStart, fZREnd, fTEnd );
	int			nj;
	Double		*Z; // Z[0] = 0, Z[nj-2] = ZR, Z[nj-1] = 1.0
	Double		*zeta;
	Flag		newZ = FALSE;

	if ( zR < 1.0 ) {
		nj = fNGridPoints + 2 + 1;
		Z = New1DArray( nj ); // Z[0] = 0, Z[nj-2] = ZR, Z[nj-1] = 1.0
		newZ = TRUE;
		zeta = fSolGrid->vec;
		for ( k = -1; k <= fNGridPoints; ++k ) {
			Z[k+1] = MAX(0.0,zeta[k] * zR);
		}
        Z[0]=0.0;
		Z[nj-1] = 1.0;
	}
	else {
		nj = fNGridPoints + 2;
		Z = &fSolGrid->vec[kPrev]; // Z[0] = 0, Z[nj-1] = 1.0
		zeta = NULL;
	}
	int 		nSpeciesIn = fSpecies->GetNSpeciesInSystem();
	Double		*molarMass = fSpecies->GetMolarMass()->vec;
	Double		**Y = fSolMassFracs->mat;
	Double		*pdf = New1DArray( nj );
	Double		sum;

	ZVar = MAX( 1e-10, ZVar );

#ifdef FINEPDF
	int		i, j;
	int		newPoints = ( fNGridPoints + 1 ) * 5 + 1;
	Double	*newGrid = New1DArray( newPoints );
	Double	*newY = New1DArray( newPoints );
	Double	*newPdf = New1DArray( newPoints );
	Double	**YBase0 = &Y[kPrev];
	
	for ( i = 0; i < fNGridPoints + 1; ++i ) {
		for ( j = 0; j < 5; ++j ) {
			newGrid[i*5+j] = Z[i] + j * 0.2 * ( Z[i+1] - Z[i] );
			newY[i*5+j] = YBase0[i][speciesIndex] + j * 0.2 * ( YBase0[i+1][speciesIndex] - YBase0[i][speciesIndex] );
		}
	}
    newGrid[0]=0.0;
	newGrid[i*5] = Z[fNGridPoints + 1];
	newY[i*5] = YBase0[fNGridPoints + 1][speciesIndex];

#ifdef NEWBETAPDF
	F77BETAPDF( newPdf, newGrid, &ZMean, &ZVar, &newPoints, &newPoints );
#else
	PDFZBETA( &ZMean, &ZVar, newGrid, newPdf, &newPoints );
#endif

/*#define DEBUGNEWTURB*/
/*#ifdef DEBUGNEWTURB*/
/*	FILE *fp = GetOutputFile( t, "NewTurbGrid", NULL, kData );*/
/*	*/
/*	fprintf( fp, "*\nZ\tT\tPDFZ\n" );*/
/*	for ( i = 0; i < newPoints; ++i ) {*/
/*		fprintf( fp, "%g\t%g\t%g\n", newGrid[i], newY[i], newPdf[i] );		*/
/*	}*/
/*	*/
/*	fclose( fp );*/
/*#endif*/

	sum = 0.0;
	for ( k = 0; k < newPoints; ++k ) {
		sum += newPdf[k] * newY[k];
	}
	
	if ( specVar ) {
		Double	sumVar = 0.0;
		for ( k = 0; k < newPoints; ++k ) {
			sumVar += newPdf[k] * ( newY[k] - sum ) * ( newY[k] - sum );
		}	
		*specVar = sumVar;
	}
	
	if ( newZ ) {
		Free1DArray( Z );
	}
	Free1DArray( newPdf );
	Free1DArray( newY );
	Free1DArray( newGrid );
	Free1DArray( pdf );

	return sum;
#else

#ifdef NEWBETAPDF
	F77BETAPDF( pdf, Z, &ZMean, &ZVar, &nj, &nj );
#else
	PDFZBETA( &ZMean, &ZVar, Z, pdf, &nj );
#endif
	
	sum = 0.0;
	for ( k = -1; k <= fNGridPoints; ++k ) {
		sum += pdf[k+1] * Y[k][speciesIndex];
	}

	if ( specVar ) {
		Double	sumVar = 0.0;
		for ( k = -1; k <= fNGridPoints; ++k ) {
			sumVar += pdf[k+1] * ( Y[k][speciesIndex] - sum ) * ( Y[k][speciesIndex] - sum );
		}	
		*specVar = sumVar;
	}

	if ( newZ ) {
		Free1DArray( Z );
	}
	Free1DArray( pdf );

	return sum;
#endif
}

Double TTransFlameSolver::TurbMeanSoot( Double t, int sootIndex, Double ZMean, Double ZVar, Double *sootVar )
{
	int			k;
	Double		zR = Interpol( t, fZRStart, fTStart, fZREnd, fTEnd );
	int			nj;
	Double		*Z; // Z[0] = 0, Z[nj-2] = ZR, Z[nj-1] = 1.0
	Double		*zeta;
	Flag		newZ = FALSE;

	if ( zR < 1.0 ) {
		nj = fNGridPoints + 2 + 1;
		Z = New1DArray( nj ); // Z[0] = 0, Z[nj-2] = ZR, Z[nj-1] = 1.0
		newZ = TRUE;
		zeta = fSolGrid->vec;
		for ( k = -1; k <= fNGridPoints; ++k ) {
		  Z[k+1] = MAX(0.0,zeta[k] * zR);
		}
		Z[0] = 0.0;
		Z[nj-1] = 1.0;
	}
	else {
		nj = fNGridPoints + 2;
		Z = &fSolGrid->vec[kPrev]; // Z[0] = 0, Z[nj-1] = 1.0
		zeta = NULL;
	}
	int 		nSpeciesIn = fSpecies->GetNSpeciesInSystem();
	Double		*molarMass = fSpecies->GetMolarMass()->vec;
	Double		*pdf = New1DArray( nj );
	Double		sum;
	Double		**moments;

	if ( fSoot ) {
		moments = fSolSootMoments->mat;
//		sootOff = fSoot->GetOffsetSootMoments();
	}
	else {
		return 0.0;
	}

	ZVar = MAX( 1e-10, ZVar );

#ifdef FINEPDF
	int		i, j;
	int		newPoints = ( fNGridPoints + 1 ) * 5 + 1;
	Double	*newGrid = New1DArray( newPoints );
	Double	*newMom = New1DArray( newPoints );
	Double	*newPdf = New1DArray( newPoints );
	Double	**MomBase0 = &moments[kPrev];
	
	for ( i = 0; i < fNGridPoints + 1; ++i ) {
		for ( j = 0; j < 5; ++j ) {
			newGrid[i*5+j] = Z[i] + j * 0.2 * ( Z[i+1] - Z[i] );
			newMom[i*5+j] = MomBase0[i][sootIndex] + j * 0.2 * ( MomBase0[i+1][sootIndex] - MomBase0[i][sootIndex] );
		}
	}
	newGrid[0]=0.0;
	newGrid[i*5] = Z[fNGridPoints + 1];
	newMom[i*5] = MomBase0[fNGridPoints + 1][sootIndex];

#ifdef NEWBETAPDF
	F77BETAPDF( newPdf, newGrid, &ZMean, &ZVar, &newPoints, &newPoints );
#else
	PDFZBETA( &ZMean, &ZVar, newGrid, newPdf, &newPoints );
#endif
	sum = 0.0;
	for ( k = 0; k < newPoints; ++k ) {
		sum += newPdf[k] * newMom[k];
	}
	
	if ( sootVar ) {
		Double	sumVar = 0.0;
		for ( k = 0; k < newPoints; ++k ) {
			sumVar += newPdf[k] * ( newMom[k] - sum ) * ( newMom[k] - sum );
		}	
		*sootVar = sumVar;
	}
	
	if ( newZ ) {
		Free1DArray( Z );
	}
	Free1DArray( newPdf );
	Free1DArray( newMom );
	Free1DArray( newGrid );
	Free1DArray( pdf );

	return sum;
#else

#ifdef NEWBETAPDF
	F77BETAPDF( pdf, Z, &ZMean, &ZVar, &nj, &nj );
#else
	PDFZBETA( &ZMean, &ZVar, Z, pdf, &nj );
#endif
	
	sum = 0.0;
	for ( k = -1; k <= fNGridPoints; ++k ) {
		sum += pdf[k+1] * moments[k][sootIndex];
	}

	if ( sootVar ) {
		Double	sumVar = 0.0;
		for ( k = -1; k <= fNGridPoints; ++k ) {
			sumVar += pdf[k+1] * ( moments[k][sootIndex] - sum ) * ( moments[k][sootIndex] - sum );
		}	
		*sootVar = sumVar;
	}

	if (isnan(sum)) 
	{
	  cerr << "Mean,Variance" << ZMean << " " << ZVar << endl;
	  for (k=-1; k<=fNGridPoints; ++k)
	    cerr << Z[k+1] << " " << pdf[k+1] << endl;

	  double wait;
	  cin >> wait;
	}

	if ( newZ ) {
		Free1DArray( Z );
	}
	Free1DArray( pdf );

	return sum;
#endif
}

FILE *TTransFlameSolver::GetOutputFile( Double theTime, const char *head, const char *tail, FileType type )
{
	int 		nOfSpeciesIn = fSpecies->GetNSpeciesInSystem();
	char		*name = new char[64];
	FILE 		*fp;
	
	sprintf( name, "%s%s%.8s_p%.2dt%9.3ems%s"
					, ( head ) ? head : "", ( head ) ? "_" : ""
					, fSpecies->GetNames()[GetFuelIndex()]
					, ( int ) floor( GetPressure( theTime ) * 1.0e-5 + 0.5 )	// in [bar]
					, theTime * 1000.0
					, ( tail ) ? tail : "" );
	
	fp = GetOutfile( name, type );
	delete name;
	
	return fp;
}

FILE *TTransFlameSolver::GetOutputFile( char */*head*/, char */*tail*/, FileType /*type*/ )
{ 
	fprintf( fOutFilePtr, "wrong instance of function 'TTransFlameSolver::GetOutputFile'\n" );
	exit( 2 );

	return NULL;
}


void TTransFlameSolver::DoExit( void )
{
	FILE	*fp = GetOutfile( "interrupt", kText );
	WriteFlameletFile( fp, NULL, NULL );
	fprintf( fOutFilePtr, "\nprogram stopped by user\noutput written to %s\n"
				, GetOutfileName( "interrupt", kText ) );
/*	cerr << NEWL << "program stopped by user\noutput written to " 
				<< GetOutfileName( "interrupt", kText ) << NEWL;*/
	fclose( fp );
	exit( 2 );
}

Double TTransFlameSolver::GetCurrentTime( void )
{
	int	theTime = GetActualPoint( fTEnd );

	return ( theTime >= 0.0 ) ? fSolTime->vec[theTime] : fTEnd;
}

#ifdef NOCCLINK
extern "C" int _main();
#endif

void GETFLAMELETSOLVER( int **object, int *objectSize, char *fInputName, int *nChars )
{
	getflameletsolver( object, objectSize, fInputName, nChars );
}

void getflameletsolver( int **object, int *objectSize, char *fInputName, int *nChars )
{
#ifdef NOCCLINK
	static Flag	first = TRUE;
	
	if ( first ) {
		_main();
		first = FALSE;
	}
#endif

	fprintf( stderr, "\n***  RifMan Version %-7s ***\n", (char *)VERSION );
	fprintf( stderr, "*** written by Heinz Pitsch ***\n" );
	fprintf( stderr, "***     ITM RWTH Aachen     ***\n" );

	if ( *objectSize != sizeof( int * ) ) {
		fprintf( stderr, "#error: size of object is currently %d, but has to be %d\n"
		"        set variable 'objectSize' to %d and ensure that this is the physical size of variable 'object'\n"
						, *objectSize, sizeof( int * ), sizeof( int * ) );
		exit( 2 );
	}

	char *inputFile = FortranToCChar( fInputName, *nChars );

	FirstInputPtr input = new FirstInput( inputFile );
   	if ( !input ) FatalError( "memory allocation of FirstInputPtr failed" );
	TTransFlameSolverPtr flame = new TTransFlameSolver( input );
	*object = (int *)flame;
	delete input;
}

void INITFLAMELET( int **object, Double *timeStart
					, char *names, int *nChars
					, Double *startSolution, int *vars, int *offsetVars
					, Double *grid, int *gridPoints, int *offsetGridPoints
					, Double *pressureStart, Double *scalarDissRateStart
					, Double *firstTimeStep )
{
	initflamelet( object, timeStart
					, names, nChars
					, startSolution, vars, offsetVars
					, grid, gridPoints, offsetGridPoints
					, pressureStart, scalarDissRateStart
					, firstTimeStep );
}

void initflamelet( int **object, Double *timeStart
					, char *names, int *nChars
					, Double *startSolution, int *vars, int *offsetVars
					, Double *grid, int *gridPoints, int *offsetGridPoints
					, Double *pressureStart, Double *scalarDissRateStart
					, Double *firstTimeStep )
{
	TTransFlameSolverPtr flame = *( TTransFlameSolverPtr *)object;

	MatrixPtr	sol = FortranToCMat( startSolution, *vars, *offsetVars
													, *gridPoints, *offsetGridPoints );
	char		**newNames = FortranToCCharArray( names, *nChars, *vars );

	flame->Initialize( *timeStart, newNames, sol->mat, *vars, grid, *gridPoints
						, *pressureStart, *scalarDissRateStart, *firstTimeStep, -1, 1 );
					
	DisposeFToCCharArray( newNames, *vars );
	DisposeFToCMat( sol );
}

void INITFLAMELETZREF( int **object, Double *timeStart
					, char *names, int *nChars
					, Double *startSolution, int *vars, int *offsetVars
					, Double *grid, int *gridPoints, int *offsetGridPoints
					, Double *pressureStart, Double *scalarDissRateStart
					, Double *firstTimeStep, Double *zRef )
{
	initflameletzref( object, timeStart
					, names, nChars
					, startSolution, vars, offsetVars
					, grid, gridPoints, offsetGridPoints
					, pressureStart, scalarDissRateStart
					, firstTimeStep, zRef );
}

void initflameletzref( int **object, Double *timeStart
					, char *names, int *nChars
					, Double *startSolution, int *vars, int *offsetVars
					, Double *grid, int *gridPoints, int *offsetGridPoints
					, Double *pressureStart, Double *scalarDissRateStart
					, Double *firstTimeStep, Double *zRef )
{
	TTransFlameSolverPtr flame = *( TTransFlameSolverPtr *)object;

	MatrixPtr	sol = FortranToCMat( startSolution, *vars, *offsetVars
													, *gridPoints, *offsetGridPoints );
	char		**newNames = FortranToCCharArray( names, *nChars, *vars );

	flame->Initialize( *timeStart, newNames, sol->mat, *vars, grid, *gridPoints
					, *pressureStart, *scalarDissRateStart, *firstTimeStep, *zRef, 1 );
					
	DisposeFToCCharArray( newNames, *vars );
	DisposeFToCMat( sol );
}

void GETSPECIESDATA( int **object, char *names, int *nChars
					, int *vars, int *offsetVars
					, Double *high, Double *low, Double *molarMass )
{
	getspeciesdata( object, names, nChars
					, vars, offsetVars
					, high, low, molarMass );
}

void getspeciesdata( int **object, char *names, int *nChars
					, int *vars, int *offsetVars
					, Double *high, Double *low, Double *molarMass )
{
	TTransFlameSolverPtr flame = *( TTransFlameSolverPtr *)object;

	MatrixPtr	highMat = FortranToCMat( high, NOFCOEFFS, NOFCOEFFS, *vars, *offsetVars );
	MatrixPtr	lowMat = FortranToCMat( low, NOFCOEFFS, NOFCOEFFS, *vars, *offsetVars );
	char		**newNames = FortranToCCharArray( names, *nChars, *vars );

	flame->GetSpeciesData( newNames, highMat, lowMat, molarMass );

	DisposeFToCCharArray( newNames, *vars );
	DisposeFToCMat( lowMat );
	DisposeFToCMat( highMat );
}

int SOLVEFLAMELETZR( int **object, Double *timeEnd, Double *pressureEnd
					, Double *scalarDissRateEnd
					, Double *temOxEnd, Double *tempFuelEnd, Double *zR )
{
	return solveflameletzr( object, timeEnd, pressureEnd
					, scalarDissRateEnd
					, temOxEnd, tempFuelEnd, zR );
}

int SOLVEFLAMELET( int **object, Double *timeEnd, Double *pressureEnd
					, Double *scalarDissRateEnd
					, Double *temOxEnd, Double *tempFuelEnd )
{
	return solveflamelet( object, timeEnd, pressureEnd
					, scalarDissRateEnd
					, temOxEnd, tempFuelEnd );
}

int solveflameletzr( int **object, Double *timeEnd, Double *pressureEnd
					, Double *scalarDissRateEnd
					, Double *temOxEnd, Double *tempFuelEnd, Double *zR )
{
	TTransFlameSolverPtr flame = *( TTransFlameSolverPtr *)object;
		
	return ( flame->Solve( *timeEnd, *pressureEnd, *scalarDissRateEnd
			, *temOxEnd, *tempFuelEnd, *zR, flame->GetInputData()->fNOutputs ) ) 
		? 1 : 0;
}

int solveflamelet( int **object, Double *timeEnd, Double *pressureEnd
					, Double *scalarDissRateEnd
					, Double *temOxEnd, Double *tempFuelEnd )
{
	TTransFlameSolverPtr flame = *( TTransFlameSolverPtr *)object;
		
	return ( flame->Solve( *timeEnd, *pressureEnd, *scalarDissRateEnd
			, *temOxEnd, *tempFuelEnd, flame->GetInputData()->fNOutputs ) ) 
		? 1 : 0;
}

void GETFLAMELETSOLUTIONEXT( int **object, char *names, int *nChars, Double *outSol, Double *grid
					, int *gridPoints, int *offsetGridPoints, int *vars, int *offsetVars, Double *density )
{
	getflameletsolutionext( object, names, nChars, outSol, grid
					, gridPoints, offsetGridPoints, vars, offsetVars, density );
}

void GETFLAMELETSOLUTION( int **object, char *names, int *nChars, Double *outSol, Double *grid
					, int *gridPoints, int *offsetGridPoints, int *vars, int *offsetVars )
{
	getflameletsolution( object, names, nChars, outSol, grid
					, gridPoints, offsetGridPoints, vars, offsetVars );
}

void getflameletsolution( int **object, char *names, int *nChars, Double *outSol, Double *grid
					, int *gridPoints, int *offsetGridPoints, int *vars, int *offsetVars )
{
	TTransFlameSolverPtr flame = *( TTransFlameSolverPtr *)object;

	MatrixPtr	sol = FortranToCMat( outSol, *vars, *offsetVars
									, *gridPoints, *offsetGridPoints );
	char	**newNames = FortranToCCharArray( names, *nChars, *vars );

	flame->GetSolution( newNames, sol->mat, grid, *gridPoints, *vars );

	DisposeFToCCharArray( newNames, *vars );
	DisposeFToCMat( sol );
}

void getflameletsolutionext( int **object, char *names, int *nChars, Double *outSol, Double *grid
					, int *gridPoints, int *offsetGridPoints, int *vars, int *offsetVars, Double *density )
{
	TTransFlameSolverPtr flame = *( TTransFlameSolverPtr *)object;

	MatrixPtr	sol = FortranToCMat( outSol, *vars, *offsetVars
									, *gridPoints, *offsetGridPoints );
	char	**newNames = FortranToCCharArray( names, *nChars, *vars );

	flame->GetSolution( newNames, sol->mat, grid, *gridPoints, *vars, density );

	DisposeFToCCharArray( newNames, *vars );
	DisposeFToCMat( sol );
}

void GETSOLUTIONINFO( int **object, Double *timeStep )
{
	getsolutioninfo( object, timeStep );
}

void getsolutioninfo( int **object, Double *timeStep )
{
	TTransFlameSolverPtr flame = *( TTransFlameSolverPtr *)object;

	flame->GetSolutionInfo( timeStep );
}

void DISPOSEFLAMELET( int **object )
{
	disposeflamelet( object );
}

void disposeflamelet( int **object )
{
	TTransFlameSolverPtr flame = *( TTransFlameSolverPtr *)object;
	delete flame;
}


void PRINTSOLUTION( void **object, Double *timeCurr, char *names, int *nChars
			, Double *outSol, Double *grid
			, int *gridPoints, int *offsetGridPoints, int *vars, int *offsetVars )
{
	printsolution( object, timeCurr, names, nChars
			, outSol, grid
			, gridPoints, offsetGridPoints, vars, offsetVars );
}

void printsolution( void **object, Double *timeCurr, char *names, int *nChars
			, Double *outSol, Double *grid
			, int *gridPoints, int *offsetGridPoints, int *vars, int *offsetVars )
{
	TTransFlameSolverPtr flame = *( TTransFlameSolverPtr *)object;

	FILE *fp = flame->GetOutputFile( *timeCurr, NULL, NULL, TFlame::kData );
	
	char		**newNames = FortranToCCharArray( names, *nChars, *vars );
	MatrixPtr	sol = FortranToCMat( outSol, *vars, *offsetVars
								, *gridPoints, *offsetGridPoints
													 );
	flame->PrintSolution( fp, *gridPoints-2, *vars, &grid[1]
									, &sol->mat[1], newNames );	
	DisposeFToCCharArray( newNames, *vars );
	DisposeFToCMat( sol );
	fclose( fp );
}

void GETDISSRATE( void **object, Double *dissRateReq, Double *ZReq, Double *dissRate, Double *z )
{
	getdissrate( object, dissRateReq, ZReq, dissRate, z );
}

void getdissrate( void **object, Double *dissRateReq, Double *ZReq, Double *dissRate, Double *z )
{
	TTransFlameSolverPtr flame = *( TTransFlameSolverPtr *)object;

	*dissRateReq = flame->GetDissRateReq( *ZReq, *dissRate, *z );
}

Double TTransFlameSolver::GetDissRateReq( Double ZReq, Double dissRate, Double z )
{
	Double	zReqStar = ( ZReq - fZl ) / ( fZr - fZl );
	Double	zStar = ( z - fZl ) / ( fZr - fZl );
	Double	zReqNew = zReqStar - 0.5;
	Double	zReqNew2 = zReqNew * zReqNew;
	Double	zReqNew4 = zReqNew2 * zReqNew2;
	Double	zReqNew6 = zReqNew4 * zReqNew2;
	Double	zReqNew8 = zReqNew6 * zReqNew2;
	Double	zNew = zStar - 0.5;
	Double	zNew2 = zNew * zNew;
	Double	zNew4 = zNew2 * zNew2;
	Double	zNew6 = zNew4 * zNew2;
	Double	zNew8 = zNew6 * zNew2;
	Double	fitZReq = 1.00539 * ( 12.9041 - 
					82.123 * zReqNew2 +
                    115.29 * zReqNew4 -
                    201.898 * zReqNew6 +
                    912.136 * zReqNew8 );
	Double	fitZ = 1.00539 * ( 12.9041 - 
                    82.123  * zNew2 +
                    115.29  * zNew4 -
                    201.898 * zNew6 +
                    912.136 * zNew8 );
	
	
	return dissRate * fitZReq / fitZ;
}

void GETDISSFUNC( void **object, Double *a, Double *z )
{
	getdissfunc( object, a, z );
}

void getdissfunc( void **object, Double *a, Double *z )
{
	TTransFlameSolverPtr flame = *( TTransFlameSolverPtr *)object;

	*a = flame->GetDissFunc( *z );
}

#ifdef MOVEZRIGHT
Double TTransFlameSolver::GetDissFunc( Double z )
{
	Double	ZR = Interpol( GetCurrentTime(), fZRStart, fTStart, fZREnd, fTEnd );
	
	
	if ( z < 1.0e-10 || z > ZR ) {
		return 0.0;
	}

	return - ( z * z * log( z / ZR ) );
}
#else
Double TTransFlameSolver::GetDissFunc( Double z )
{
	Double	zStar = ( z - fZl ) / ( fZr - fZl );
	Double	zNew = zStar - 0.5;
	Double	zNew2 = zNew * zNew;
	Double	zNew4 = zNew2 * zNew2;
	Double	zNew6 = zNew4 * zNew2;
	Double	zNew8 = zNew6 * zNew2;

	return 1.00539 * ( 12.9041 - 
                    82.123  * zNew2 +
                    115.29  * zNew4 -
                    201.898 * zNew6 +
                    912.136 * zNew8 );
}
#endif

#ifdef WRITEGLOBALPRODRATE
void TTransFlameSolver::PrintProdRateGlobalReac( Double t )
{
#include "nHeptaneNO.red4.h"
fprintf( stderr, "#error code -13: old kinetic data\n" );
exit(2);
	int			k;
	Double		**Y = fSolMassFracs->mat;
	Double		*T = fSolTemp->vec;
	Double		*Z = fSolGrid->vec;
	Double		**M = ( fSoot ) ? fSolSootMoments->mat : NULL;
	Double		pressure = GetPressure( t );
	Double		*w = GetReaction()->GetReactionRate()->vec;
	Double		reacRate;
	FILE		*fp = GetOutputFile( t, "GlRa", NULL, TFlame::kData );
	
//	header
	fprintf( fp, "*\n" );
	fprintf( fp, "%s", "y" );
	fprintf( fp, "\t%s", "I.    N\\d2\\n + O\\d2\\n = 2 NO" );
	fprintf( fp, "\t%s", "IIa.  N\\d2\\n + C\\d2\\nH\\d4\\n + O\\d2\\n + OH = NO + HCN + H\\d2\\nO + CO H\\d2\\n" );
	fprintf( fp, "\t%s", "IIb.  HCN + O\\d2\\n + H\\d2\\nO = NO + CO + H\\d2\\n + OH" );
	fprintf( fp, "\t%s", "IIIa. N\\d2\\n + 2 OH + M = N\\d2\\nO + H\\d2\\nO + M" );
	fprintf( fp, "\t%s", "IIIb. N\\d2\\nO + 2OH = 2 NO + H\\d2\\nO" );
	fprintf( fp, "\t%s", "IV.   NO + H\\d2\\nO = NO\\d2\\n + H\\d2\\n" );
	fprintf( fp, "\t%s", "Va.   NO + C\\d2\\nH\\d4\\n + OH = HNCO + CO + 2H\\d2\\n" );
	fprintf( fp, "\t%s", "Vb.   HNCO + NO = N\\d2\\n + CO + OH" );
	fprintf( fp, "\t%s", "VI.   NH\\d3\\n + NO = N\\d2\\n + H\\d2\\n + OH" );
	fprintf( fp, "\n%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g"
		, Z[kPrev], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );

	for ( k = 0; k < fNGridPoints; ++k ) {
#ifdef DELTAINEW
		UpdateThermoProps( k, Y[k], T[k], pressure, fProperties->GetDensityRef()
										, kDensFromPress, ( !fSoot ) ? NULL : M[k] );
#else 
		T0DFlame::UpdateThermoProps( Y[k], T[k], pressure, fProperties->GetDensityRef()
										, kDensFromPress, ( !fSoot ) ? NULL : M[k] );
#endif
		fprintf( fp, "\n%-12E", Z[k] );
		// I
		reacRate = w[rN7f] - w[rN7b] - w[rN101f] + w[rN101b] - w[rN104f] + w[rN104b] 
				 - w[rN37f] + w[rN37b];
		fprintf( fp, "\t%-12E", reacRate );
		// IIa
		reacRate = w[rN110f] - w[rN110b];
		fprintf( fp, "\t%-12E", reacRate );
		
		// IIb
		reacRate = - w[rN111] - w[rN113f] + w[rN113b] + w[rN120f] - w[rN120b]
			 	   - w[rN130f] + w[rN130b] - w[rN131f] + w[rN131b];
		fprintf( fp, "\t%-12E", reacRate );
		
		// IIIa
		reacRate = w[rN6f] - w[rN6b] + w[rN9f] - w[rN9b] - w[rN25f] + w[rN25b]
				 + w[rN40f] - w[rN40b] - w[rN45f] + w[rN45b] - w[rN46yf] + w[rN46yb]
				 - w[rN49f] + w[rN49b] + w[rN60f] - w[rN60b] + w[rN63f] - w[rN63b];
		fprintf( fp, "\t%-12E", reacRate );
		
		// IIIb
		reacRate = w[rN17f] - w[rN17b] + w[rN18f] - w[rN18b] + w[rN19f] - w[rN19b]
				 + w[rN20f] - w[rN20b] + w[rN47f] - w[rN47b] - w[rN71f] + w[rN71b];
		fprintf( fp, "\t%-12E", reacRate );
		
		// IV
		reacRate = - w[rN84f] + w[rN84b] + w[rN85f] - w[rN85b] - w[rN87f] + w[rN87b]
				   - w[rN88f] + w[rN88b];
		fprintf( fp, "\t%-12E", reacRate );
		
		// Va
		reacRate = - w[rN62f] + w[rN62b] - w[rN65f] + w[rN65b] + w[rN74f] - w[rN74b]
				   + w[rN112f] - w[rN112b];
		fprintf( fp, "\t%-12E", reacRate );
		
		// Vb
		reacRate = w[rN60f] - w[rN60b] + w[rN61f] - w[rN61b] + w[rN63f] - w[rN63b];
		fprintf( fp, "\t%-12E", reacRate );
		
		// VI
		reacRate = w[rN1f] - w[rN1b] + w[rN2f] - w[rN2b] + w[rN3f] - w[rN3b] + 
				   w[rN4f] - w[rN4b];
		fprintf( fp, "\t%-12E", reacRate );
		
	}
	fprintf( fp, "\n%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g"
		, Z[fNGridPoints], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );
	
	fclose( fp );
}
#endif

void TTransFlameSolver::WriteSootInfo( Double theTime, Double *temp, Double **Y, Double pressure, Double **moments, FILE *fp )
{
	if ( fSoot ) {
		int		i, k;
		int		nSootMoments = fSoot->GetNSootMoments();
		const int	nPAHMolecules = 5;
		char	tmpName[127];
		
		// write soot volume fraction
		Double	fvFact = 24.0/1800.0; // fMolarMassSoot / fSootDensity;
		fprintf( fp, "fv\n" );
		for ( k = 0; k < fNGridPoints+2; ++k ) {
			fprintf( fp, "\t%-.6e", moments[k-1][1] * fvFact );
			if ( (k+1) % 5 == 0 ) {
				fprintf( fp, "\n" );
			}
		}
		if ( k % 5 ) {
			fprintf( fp, "\n" );
		}
				
		// write pah moments
		int			nPAHMoments = fSoot->GetNPAHMoments();
		int 		nSpeciesIn = fSpecies->GetNSpeciesInSystem();
		VectorPtr	conv1Vec = NewVector( fNGridPoints+2 );
		Double		*conv1 = &conv1Vec->vec[kNext];
		VectorPtr	conv2Vec = NewVector( fNGridPoints+2 );
		Double		*conv2 = &conv2Vec->vec[kNext];
		VectorPtr	theCSootVec = NewVector( fNGridPoints+2 );
		Double		*theCSoot = &theCSootVec->vec[kNext];
		VectorPtr	theOxFactVec = NewVector( fNGridPoints+2 );
		Double		*theOxFact = &theOxFactVec->vec[kNext];
		MatrixPtr	thePAHMomMat = NewMatrix( nPAHMoments, fNGridPoints+2, kColumnPointers );
		Double		**thePAHMom = &thePAHMomMat->mat[kNext];
		TensorPtr	sourcesMat = NewTensor( nSootMoments, kNSootSources+2, fNGridPoints+2, kColumnPointers );
		Double		***sources = sourcesMat->tensor;
		MatrixPtr	pahConcsMat1 = NewMatrix( nPAHMolecules, fNGridPoints+2, kColumnPointers );
		Double		**pahConcs1 = &pahConcsMat1->mat[kNext];
		MatrixPtr	pahConcsMat2 = NewMatrix( nPAHMolecules, fNGridPoints+2, kColumnPointers );
		Double		**pahConcs2 = &pahConcsMat2->mat[kNext];
		Double		**si = NULL;
		Double		*sik = NULL;
		Double		*pahMoments = fSoot->GetPAHMoments()->vec;
		Double		**pahConcs = fSoot->GetPij()->mat;
		Double		*molarMass = fSpecies->GetMolarMass()->vec;
		Double		*Z = fSolGrid->vec;
		Double		density;

		for ( k = -1; k <= fNGridPoints; ++k )
		{
#ifdef DELTAINEW
			UpdateThermoProps( k, Y[k], temp[k], pressure, density, kDensFromPress, moments[k] ); 
#else
			T0DFlame::UpdateThermoProps( Y[k], temp[k], pressure, density, kDensFromPress, moments[k] );
#endif
			// save pah moments
			copy( nPAHMoments, pahMoments, 1, thePAHMom[k], 1 );
			// save pah concentrations
			copy( nPAHMolecules, &pahConcs[0][1], 1, pahConcs1[k], 1 );
			copy( nPAHMolecules, &pahConcs[1][1], 1, pahConcs2[k], 1 );
			// save oxidation factor
			theOxFact[k] = fSoot->GetSootOxCoeff( Y[k], density, molarMass );
			theCSoot[k] = fSoot->GetCCSootStar();
			// save souce terms
			for ( i = 0; i < nSootMoments; ++i ) {
				sik = sources[i][k+1];
				if ( fSoot->WithNucleation() ) {
					sik[kNuc] = fSoot->NucleationNew( i, temp[k], pahMoments ) / density;
				}
				if ( fSoot->WithCoagulation() ) {
					sik[kCoag] = fSoot->SourceCoagulationNew( i, temp[k], moments[k] ) / density;
				}
				if ( fSoot->WithCondensation() ) {
//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//
					sik[kCond] = fSoot->SourceCondensationNew( i, temp[k], fSoot->GetPAHMoments()->vec
									, moments[k], Y[k], density, molarMass ) / density;
//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//
				}
				if ( fSoot->WithSurfaceGrowth() ) {
					sik[kSG] = fSoot->SourceSurfGrowthNew( i, moments[k], Y[k], density, molarMass ) / density;
				}
				if ( fSoot->WithSurfaceOxidation() ) {
					sik[kOx] = fSoot->SourceSootOxidationNew( i, moments[k], Y[k], density, molarMass ) / density;
				}
				if ( k > -1 && k < fNGridPoints )
				{
					Double	lambdaOverCpCurr = fProperties->GetMixConductivity() / fProperties->GetMixHeatCapacity();
					Double	fracIndex = i - 2.0 / 3.0;
					Double	rhoPrev, rhoNext, mm;
					Double	lambdaOverCpPrev, lambdaOverCpNext;
					fProperties->ComputeMixtureMolarMass( mm, Y[k-1], molarMass, nSpeciesIn );
					rhoPrev = pressure * mm / ( RGAS * temp[k-1] );
					fSpecies->ComputeSpeciesProperties( temp[k-1] );
					ComputeDeltaI( k-1, Y[k-1], temp[k-1] );
					fProperties->CompMixtureProps( fSpecies->GetHeatCapacity()->vec
					, fSpecies->GetConductivity()->vec, fSpecies->GetViscosity()->vec, Y[k-1], temp[k-1]
					, pressure, fProperties->GetDensityRef()
					, kDensFromPress, nSpeciesIn, fSpecies );
					lambdaOverCpPrev = fProperties->GetMixConductivity() / fProperties->GetMixHeatCapacity();

					fProperties->ComputeMixtureMolarMass( mm, Y[k+1], molarMass, nSpeciesIn );
					rhoNext = pressure * mm / ( RGAS * temp[k+1] );
					fSpecies->ComputeSpeciesProperties( temp[k+1] );
					ComputeDeltaI( k+1, Y[k+1], temp[k+1] );
					fProperties->CompMixtureProps( fSpecies->GetHeatCapacity()->vec
					, fSpecies->GetConductivity()->vec, fSpecies->GetViscosity()->vec, Y[k+1], temp[k+1]
					, pressure, fProperties->GetDensityRef()
					, kDensFromPress, nSpeciesIn, fSpecies );
					lambdaOverCpNext = fProperties->GetMixConductivity() / fProperties->GetMixHeatCapacity();

					fProperties->ComputeMixtureMolarMass( mm, Y[k], molarMass, nSpeciesIn );
					rhoNext = pressure * mm / ( RGAS * temp[k+1] );
					fSpecies->ComputeSpeciesProperties( temp[k+1] );
					ComputeDeltaI( k+1, Y[k+1], temp[k+1] );
					fProperties->CompMixtureProps( fSpecies->GetHeatCapacity()->vec
					, fSpecies->GetConductivity()->vec, fSpecies->GetViscosity()->vec, Y[k+1], temp[k+1]
					, pressure, fProperties->GetDensityRef()
					, kDensFromPress, nSpeciesIn, fSpecies );

#ifdef SOOTCONVECTION
					conv1[k] = 0.25 / density * ( 
										fFDWMinus[k] * rhoPrev * GetDissRate( theTime, Z[k-1] ) 
										+ fFDWCurr[k] * density * GetDissRate( theTime, Z[k] )
										+ fFDWPlus[k] * rhoNext * GetDissRate( theTime, Z[k+1] ) );
					conv2[k] = 0.25 / density * ( 
										density * GetDissRate( theTime, Z[k] ) / lambdaOverCpCurr
											* ( fFDWMinus[k] * lambdaOverCpPrev
												+ fFDWCurr[k] * lambdaOverCpCurr
												+ fFDWPlus[k] * lambdaOverCpNext ) );
					Double	convVelo = 0.25 / density * ( 
										fFDWMinus[k] * rhoPrev * GetDissRate( theTime, Z[k-1] ) 
										+ fFDWCurr[k] * density * GetDissRate( theTime, Z[k] )
										+ fFDWPlus[k] * rhoNext * GetDissRate( theTime, Z[k+1] )
										+ density * GetDissRate( theTime, Z[k] ) / lambdaOverCpCurr
											* ( fFDWMinus[k] * lambdaOverCpPrev
												+ fFDWCurr[k] * lambdaOverCpCurr
												+ fFDWPlus[k] * lambdaOverCpNext ) );
#	ifdef UPWINDSOOT
//					fprintf( stderr, "y[gpOff+sootOff]-stuff has to be changed to moments\n" );
				if ( convVelo > 0.0 ) {
					sik[kNSootSources+1] = -convVelo * ( ( moments[k][i] / density - moments[k-1][i] / rhoPrev ) 
							/ ( Z[k] - Z[k-1] )
							- 1.0 / fSoot->GetLewis1() 
							* FRACMOMFACT * ( fSoot->FracMom2( fracIndex, moments[k+1] ) / rhoNext 
								- fSoot->FracMom2( fracIndex, moments[k] ) / density ) 
								/ ( Z[k+1] - Z[k] ) );
				}
				else {
					sik[kNSootSources+1] = -convVelo * ( ( moments[k+1][i] / rhoNext - moments[k][i] / density ) 
							/ ( Z[k+1] - Z[k] )
							- 1.0 / fSoot->GetLewis1() 
							* FRACMOMFACT * ( fSoot->FracMom2( fracIndex, moments[k] ) / density
							- fSoot->FracMom2( fracIndex, moments[k-1] ) / rhoPrev ) 
								/ ( Z[k] - Z[k-1] ) );
				}
#	else
				sik[kNSootSources+1] = -convVelo * ( 
					1.0 / fSoot->GetLewis1() 
					* FRACMOMFACT * ( fFDWMinus[k] * fSoot->FracMom2( fracIndex, moments[k-1] ) / rhoPrev
					+ fFDWCurr[k] * fSoot->FracMom2( fracIndex, moments[k] ) / density
					+ fFDWPlus[k] * fSoot->FracMom2( fracIndex, moments[k+1] ) / rhoNext )
										
					- ( fFDWMinus[k] * moments[k-1][i] / rhoPrev
					+ fFDWCurr[k] * moments[k][i] / density
					+ fFDWPlus[k] * moments[k+1][i] / rhoNext )
										);
#	endif
#endif
				}
				else {
					sik[kNSootSources+1] = 0.0;
				}

//				if ( fSoot->WithThermoPhoresis() ) {
					if ( k > -1 && k < fNGridPoints ) {
#	ifdef PRONE
						Double	PrCurr = 1.0;
#	else
						Double	PrCurr = fProperties->GetMixViscosity() / ( fProperties->GetMixConductivity() / fProperties->GetMixHeatCapacity() );
#	endif
						sik[kTherPhor] = 0.275 * moments[k][i] / density / ( temp[k] ) * PrCurr * GetDissRate( theTime, Z[k] )
										 * ( fWMinus[k] * temp[k-1]
											+fWCurr[k] * temp[k]
											+fWPlus[k] * temp[k+1] );
					}
//				}
				if ( k > -1 && k < fNGridPoints ) {
					Double	rhoPrev, rhoNext, mm;
					fProperties->ComputeMixtureMolarMass( mm, Y[k-1], molarMass, nSpeciesIn );
					rhoPrev = pressure * mm / ( RGAS * temp[k-1] );
					fProperties->ComputeMixtureMolarMass( mm, Y[k+1], molarMass, nSpeciesIn );
					rhoNext = pressure * mm / ( RGAS * temp[k+1] );
#ifdef SIZEDEPDIFFUSION
					sik[kNSootSources] = GetDissRate( theTime, Z[k] ) / ( 2.0 * fSoot->GetLewis1() ) 
						* FRACMOMFACT * ( fWMinus[k] * fSoot->FracMom2( i - 2.0 / 3.0, moments[k-1] ) / rhoPrev
						+ fWCurr[k] * fSoot->FracMom2( i - 2.0 / 3.0, moments[k] ) / density
						+ fWPlus[k] * fSoot->FracMom2( i - 2.0 / 3.0, moments[k+1] ) / rhoNext );
#else
					sik[kNSootSources] = GetDissRate( theTime, Z[k] ) / ( 2.0 * fSoot->GetLewis1() ) 
						* ( fWMinus[k] * moments[k-1][i] / rhoPrev
						+ fWCurr[k] * moments[k][i] / density
						+ fWPlus[k] * moments[k+1][i] / rhoNext );
#endif
				} else {
					sik[kNSootSources] = 0.0;
				}
			}
		}

		PrintFlameletVector( fNGridPoints+2, theOxFact, "OxidationFactor"
				, fp, 1 );
		PrintFlameletVector( fNGridPoints+2, theCSoot, "CSootStar"
				, fp, 1 );

		PrintFlameletVector( fNGridPoints+2, conv1, "ConvVelo1", fp, 1 );
		PrintFlameletVector( fNGridPoints+2, conv2, "ConvVelo2", fp, 1 );
		
		for ( i = 0; i < nPAHMoments; ++i ) {
			sprintf( tmpName, "MPAH%d", i );
			PrintFlameletVector( fNGridPoints+2, &thePAHMom[kPrev][i], tmpName
					, fp, thePAHMomMat->phys_rows );
		}

		for ( i = 0; i < nSootMoments; ++i ) {
			si = sources[i];
			sprintf( tmpName, "So_Nuc_%d", i );
			PrintFlameletVector( fNGridPoints+2, &si[kCurr][kNuc], tmpName, fp, sourcesMat->phys_rows );
			if ( fSoot->WithCoagulation() ) {
				sprintf( tmpName, "So_Coag_%d", i );
				PrintFlameletVector( fNGridPoints+2, &si[kCurr][kCoag], tmpName, fp, sourcesMat->phys_rows );
			}
			if ( fSoot->WithCondensation() ) {
				sprintf( tmpName, "So_Cond_%d", i );
				PrintFlameletVector( fNGridPoints+2, &si[kCurr][kCond], tmpName, fp, sourcesMat->phys_rows );
			}
			if ( fSoot->WithSurfaceGrowth() ) {
				sprintf( tmpName, "So_SG_%d", i );
				PrintFlameletVector( fNGridPoints+2, &si[kCurr][kSG], tmpName, fp, sourcesMat->phys_rows );
			}
			if ( fSoot->WithSurfaceOxidation() ) {
				sprintf( tmpName, "So_Ox_%d", i );
				PrintFlameletVector( fNGridPoints+2, &si[kCurr][kOx], tmpName, fp, sourcesMat->phys_rows );
			}
//			if ( fSoot->WithThermoPhoresis() ) {
				sprintf( tmpName, "So_TherPhor_%d", i );
				PrintFlameletVector( fNGridPoints+2, &si[kCurr][kTherPhor], tmpName, fp, sourcesMat->phys_rows );
//			}
#ifdef SOOTCONVECTION
			sprintf( tmpName, "So_Conv_%d", i );
			PrintFlameletVector( fNGridPoints+2, &si[kCurr][kNSootSources+1], tmpName, fp, sourcesMat->phys_rows );
#endif
			sprintf( tmpName, "So_Diff_%d", i );
			PrintFlameletVector( fNGridPoints+2, &si[kCurr][kNSootSources], tmpName, fp, sourcesMat->phys_rows );
		}

		for ( i = 0; i <= nPAHMolecules; ++i ) {
			sprintf( tmpName, "P0%d", i+1 );
			PrintFlameletVector( fNGridPoints+2, &pahConcs1[kPrev][i], tmpName
					, fp, pahConcsMat1->phys_rows );
		}

		for ( i = 0; i <= nPAHMolecules; ++i ) {
			sprintf( tmpName, "P1%d", i+1 );
			PrintFlameletVector( fNGridPoints+2, &pahConcs2[kPrev][i], tmpName
					, fp, pahConcsMat2->phys_rows );
		}

		DisposeMatrix( pahConcsMat2 );
		DisposeMatrix( pahConcsMat1 );
		DisposeTensor( sourcesMat );
		DisposeMatrix( thePAHMomMat );
		DisposeVector( theOxFactVec );
		DisposeVector( theCSootVec );
		DisposeVector( conv1Vec );
		DisposeVector( conv2Vec );
	}
}

#ifdef LOGNORMALCHI
//#include "RandomNums.h"
void TTransFlameSolver::SetRandomNumber( void )
{
	static int	idum = 0;

	fRandomNumber = exp( gasdev( &idum ) );
}
#endif

int GETSOOTSOURCES( int **object, int *whichMoment, char *names, int *nChars, Double *outSol, Double *grid
					, int *gridPoints, int *offsetGridPoints, int *vars, int *offsetVars )
{
	return getsootsources( object, whichMoment, names, nChars, outSol, grid
					, gridPoints, offsetGridPoints, vars, offsetVars );
}

int getsootsources( int **object, int *whichMoment, char *names, int *nChars, Double *sootSources, Double *grid
					, int *gridPoints, int *offsetGridPoints, int *nSources, int *offsetSources )
{
	TTransFlameSolverPtr flame = *( TTransFlameSolverPtr *)object;
	int error;

	MatrixPtr	sources = FortranToCMat( sootSources, *nSources, *offsetSources
									, *gridPoints, *offsetGridPoints );
									
	char	**newNames = FortranToCCharArray( names, *nChars, *nSources );

	error = flame->GetSootSources( *whichMoment, newNames, sources->mat, grid, *gridPoints, *nSources );

	DisposeFToCCharArray( newNames, *nSources );
	DisposeFToCMat( sources );
	
	return error;
}

int TTransFlameSolver::GetSootSources( int whichMoment, ConstStringArray names, Double **sources, Double *grid, int gridPointsA, int vars )
{
	if ( !fSoot ) {
		return 1;
	}
	if ( whichMoment >= fSoot->GetNSootMoments() || whichMoment < 0 ) {
		fprintf( stderr, "###error in function 'GetSootSources': soot moment no. %d not allowed\n", whichMoment );
	}
	
	int			i, k, ind;
	int			nOfSpecies = fSpecies->GetNOfSpecies();
	int			nSootMoments = fSoot->GetNSootMoments();
	MMDataBag	bag( kNSootSources );
	Double		**Y = fMassFracsWork->mat;
	Double		*temp = fTempWork->vec;
	Double		**moments = fSootMomentsWork->mat;
	Double		theTime = GetCurrentTime();
	MatrixPtr	sourceWork = NewMatrix( kNSootSources, fNGridPoints+2, kColumnPointers );
	Double		**theSources = sourceWork->mat;
	VectorPtr	fOutSol = NewVector( gridPointsA );
	Flag		*outSet = new Flag[vars];
	if ( !outSet ) fprintf( stderr, "new failed\n" );
	char		theNames[kNSootSources][20];
	for ( i = 0; i < vars; ++i ) {
		outSet[i] = FALSE;
	}
	// init names
	strcpy( theNames[kNuc], "Nucleation" );
	strcpy( theNames[kCoag], "Coagulation" );
	strcpy( theNames[kCond], "Condensation" );
	strcpy( theNames[kSG], "SurfaceGrowth" );
	strcpy( theNames[kOx], "SurfaceOxidation" );
	strcpy( theNames[kTherPhor], "Thermophoresis" );

	if ( fSoot ) {
//		momOff = fSoot->GetOffsetSootMoments();
//		theMomOff = fSootMomentsWork->rows;
	}
	
// copy solution to bag
	bag.Initialize();
	bag.SetOldInpedVar( &fSolGrid->vec[kPrev], fSolGrid->len, 1, "xIn" );

	SetOutSolution();
	
	SetSootSources( whichMoment, sourceWork
				, theTime, temp, Y
				, GetPressure( theTime ), moments );	
	
	for ( i = 0; i < kNSootSources; ++i ) {
		bag.Insert( &theSources[0][i], fNGridPoints+2, kNSootSources, theNames[i] );
	}
	
	bag.SetNewInpedVar( grid, gridPointsA, 1, "xNew" );

#ifdef DEBUGBAG
	cout << bag;
#endif

	Double		**newSource = &sources[kNext];
	Double		*outWork = &fOutSol->vec[kNext];
	for ( i = 0; i < bag.NumElems(); ++i ) {
		ind = GetVariableIndex( bag[i].Name(), names, vars );
		if ( ind >= 0 ) {
			bag[i].Map( fOutSol );
			for ( k = -1; k < gridPointsA-1; ++k ) {
				newSource[k][ind] = outWork[k];
			}
			outSet[ind] = TRUE;
		}
		else {
//			cerr << "#warning: no match for source term " << bag[i].Name() << NEWL;
		}
	}
	
	for ( i = 0; i < vars; ++i ) {
		if ( outSet[i] == FALSE ) {
			fprintf( fOutFilePtr, "#warning from function 'GetSolution': no match for source term %s\n", names[i] );
//			cerr << "#warning from function 'GetSolution': no match for variable " 
//				<< names[i] << NEWL;
		}
	}
	
#ifdef DEBUGSOURCES
	FILE	*fp = GetOutfile( "SootSources", kData );
	PrintSolution( fp, gridPointsA-2, vars, &grid[kNext], newSource, names );
	fclose( fp );
#endif

	// clean up
	delete outSet;
	DisposeVector( fOutSol );
	DisposeMatrix( sourceWork );
	
	return 0;
}

void TTransFlameSolver::SetSootSources( int whichMoment, MatrixPtr sourcesMat
				, Double theTime, Double *temp, Double **Y
				, Double pressure, Double **moments )
{
	int			k;
	int			nPAHMoments = fSoot->GetNPAHMoments();
	int 		nSpeciesIn = fSpecies->GetNSpeciesInSystem();
	MatrixPtr	thePAHMomMat = NewMatrix( nPAHMoments, fNGridPoints+2, kColumnPointers );
	Double		**sources = &sourcesMat->mat[kNext];
	Double		*sik = NULL;
	Double		*pahMoments = fSoot->GetPAHMoments()->vec;
	Double		*molarMass = fSpecies->GetMolarMass()->vec;
	Double		*Z = fSolGrid->vec;
	Double		density;

	for ( k = -1; k <= fNGridPoints; ++k )
	{
#ifdef DELTAINEW
		UpdateThermoProps( k, Y[k], temp[k], pressure
						, density, kDensFromPress, moments[k] );
#else 
		T0DFlame::UpdateThermoProps( Y[k], temp[k], pressure
						, density, kDensFromPress, moments[k] );
#endif

		// save souce terms
		sik = sources[k];
		sik[kNuc] = fSoot->NucleationNew( whichMoment, temp[k], pahMoments );
		
		if ( fSoot->WithCoagulation() ) {
			sik[kCoag] = fSoot->SourceCoagulationNew( whichMoment, temp[k], moments[k] );
		}
		else {
			sik[kCoag] = 0.0;
		}
		
		if ( fSoot->WithCondensation() ) {
//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//
			sik[kCond] = fSoot->SourceCondensationNew( whichMoment, temp[k], fSoot->GetPAHMoments()->vec
							, moments[k], Y[k], density, molarMass ); 
//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//
		}
		else {
			sik[kCond] = 0.0;
		}
		
		if ( fSoot->WithSurfaceGrowth() ) {
			sik[kSG] = fSoot->SourceSurfGrowthNew( whichMoment, moments[k], Y[k], density, molarMass );
		}
		else {
			sik[kSG] = 0.0;
		}
		
		if ( fSoot->WithSurfaceOxidation() ) {
			sik[kOx] = fSoot->SourceSootOxidationNew( whichMoment, moments[k], Y[k], density, molarMass );
		}
		else {
			sik[kOx] = 0.0;
		}
		
		if ( fSoot->WithThermoPhoresis() ) {
			if ( k > -1 && k < fNGridPoints )
			{
				Double	fracIndex = whichMoment - 2.0 / 3.0;
				Double	rhoPrev, rhoNext, mm;
				fProperties->ComputeMixtureMolarMass( mm, Y[k-1], molarMass, nSpeciesIn );
				rhoPrev = pressure * mm / ( RGAS * temp[k-1] );
				fProperties->ComputeMixtureMolarMass( mm, Y[k+1], molarMass, nSpeciesIn );
				rhoNext = pressure * mm / ( RGAS * temp[k+1] );
				Double	convVelo = 0.25 * ( fFDWMinus[k] * rhoPrev * GetDissRate( theTime, Z[k-1] ) 
									+ fFDWCurr[k] * density * GetDissRate( theTime, Z[k] )
									+ fFDWPlus[k] * rhoNext * GetDissRate( theTime, Z[k+1] ) );
				sik[kTherPhor] = convVelo * ( fFDWMinus[k] * moments[k-1][whichMoment] / rhoPrev
									+ fFDWCurr[k] * moments[k][whichMoment] / density
									+ fFDWPlus[k] * moments[k+1][whichMoment] / rhoNext
									- 1.0 / fSoot->GetLewis1() 
									* FRACMOMFACT * ( fFDWMinus[k] * fSoot->FracMom2( fracIndex, moments[k-1] ) / rhoPrev
									+ fFDWCurr[k] * fSoot->FracMom2( fracIndex, moments[k] ) / density
									+ fFDWPlus[k] * fSoot->FracMom2( fracIndex,moments[k+1] ) / rhoNext ) );
			}
			else {
				sik[kTherPhor] = 0.0;
			}
		}
	}
}

void READSTARTPROF( int **object
					, int *gridPoints, int *offsetGridPoints, int *vars, int *offsetVars
					, Double *grid, Double *startSolution, Double *pressure, Double *chi
					, Double *theTime, Double *currTimeStep
		    , int *tempOff, int *progOff, int *enthOff, int *speciesOff, int *sootOff )
{
	readstartprof( object
					, gridPoints, offsetGridPoints, vars, offsetVars
					, grid, startSolution, pressure, chi
					, theTime, currTimeStep
		       , tempOff, progOff, enthOff, speciesOff, sootOff );
}

void readstartprof( int **object
					, int *gridPoints, int *offsetGridPoints, int *vars, int *offsetVars
					, Double *grid, Double *startSolution, Double *pressure, Double *chi
					, Double *theTime, Double *currTimeStep
		    , int *tempOff, int *progOff, int *enthOff, int *speciesOff, int *sootOff )
{
	TTransFlameSolverPtr flame = *( TTransFlameSolverPtr *)object;

	MatrixPtr	sol = FortranToCMat( startSolution, *vars, *offsetVars
													, *gridPoints, *offsetGridPoints );

	flame->ReadStartProfiles( flame->GetInputData()
				, *gridPoints, *vars, grid, sol->mat, pressure
				, chi, theTime, currTimeStep
				  , *tempOff, *progOff, *enthOff, *speciesOff, *sootOff );
					
	DisposeFToCMat( sol );
}

void TTransFlameSolver::ReadStartProfiles( TInputDataPtr inp
				, int nGridPoints, int nVars, Double *x, Double **y, Double *pressure
				, Double *chi, Double *theTime, Double *currTimeStep
					   , int tempOff, int progOff, int enthOff, int speciesOff, int sootOff )
{
	StartProfilePtr	sp = NULL;
	char 	*insp = inp->fStartProfileFile;
	FILE	*fpS = NULL;
	char	*fileName = GetFullPath( insp, kFileName );

	sp = new StartProfile;	
	if ( !sp ) FatalError( "new failed" );
	if ( !insp  || ( fpS = fopen( fileName, "r" ) ) == NULL ) {
		fprintf( stderr, "startprofiles file '%s' not found\n", fileName );
		exit( 2 );
	} 
	else {
		fprintf( stderr, "read initial solution from file '%s'\n", fileName );
		::ReadStartProfiles( sp, fpS );
		SetInitialValues( fInputData, sp, nGridPoints, nVars, x, y
						, pressure, chi, theTime, currTimeStep
				  , tempOff, progOff, enthOff, speciesOff, sootOff );
		CleanReadStartProfile();
		delete sp;
		fclose( fpS );
	}
	delete fileName;
}

void TTransFlameSolver::SetInitialValues( TInputDataPtr inp, StartProfilePtr sp
				, int nGridPoints, int nVars, Double *x, Double **y, Double *pressure, Double *chi
				, Double *theTime, Double *currTimeStep
					  , int tempOff, int progOff, int enthOff, int speciesOff, int sootOff )
{
	x = &x[kNext];
	y = &y[kNext];

	int 				i, j, k;
	int					variable, speciesIndex;
	Flag				ZSet = FALSE;
	int					gridPointsIn = sp->gridPoints;	// including boundaries
	int					nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();
	Double				*yInFloat = sp->data;
	char				*string = sp->labels;
	struct _parameter	*param;
	
	SetInitialValues( nGridPoints-2, nVars, x, y ); // default values
	
	if ( nGridPoints < gridPointsIn ) {
		fprintf( stderr, "nGridPointsIn = %d > nGridPoints = %d\n", gridPointsIn, nGridPoints );
		exit(2);
	}

// find independent coordinate
	for ( i = 0; i < sp->variables; ++i ) {
		if ( strncmp( string, "z", 1 ) == 0 && ZSet == FALSE ) {
			cerr << "choose inputGrid" << NEWL;
			for ( j = -1; j <= gridPointsIn-2; ++j ) {
				x[j] = yInFloat[i*gridPointsIn + j+1];		// implicit cast from float to Double
			}
			ZSet = TRUE;
		}
		string += strlen(string) + 1;
	}

// error checking
	if ( !ZSet ) {
		cerr << "error: can't find coordinate 'Z'" << NEWL;
		exit(2);
	}

// reset string
	string = sp->labels;
	
	for ( i = 0; i < sp->variables; ++i ) {
		if ( strncmp( string, "temperature", 11 ) == 0 ) {
			variable = tempOff;
		}
		else if ( strncmp( string, "prograt", 7 ) == 0 ) {
		  variable = progOff;
		}
		else if ( strncmp( string, "enthloss", 8 ) == 0 ) {
		  variable = enthOff;
		}
		else if ( strncmp( string, "massfraction-", 13 ) == 0 ){
			string += 13;
			UpperString( string );
			if ( ( speciesIndex = inp->FindSpecies( string ) ) >= 0 ) {
				if ( speciesIndex < nSpeciesInSystem ) {
					variable = speciesOff + speciesIndex;
				}
				else {
					string += strlen(string) + 1;
					continue;
				}
			}
			else {
				cerr << "warning: no match for species " << string << NEWL;
				string += strlen(string) + 1;
				continue;
			}
		}
		else if ( sootOff >= 0 && GetSoot() && strncmp( string, "m", 1 ) == 0 ) {
			string += 1;
			int	num = atoi( string );
			if ( isdigit( string[0] ) && num < GetSoot()->GetNSootMoments() ) {
				variable = num + sootOff;
			}
			else {
				string += strlen(string) + 1;
				continue;
			}
		}
		else {
			string += strlen(string) + 1;
			continue;
		}

		string += strlen(string) + 1;
		for ( k = -1; k <= gridPointsIn-2; ++k ) {
			y[k][variable] = yInFloat[i*gridPointsIn + k + 1];	// copy workspace to vector of solution
		}
	}
		
	param = GetParameter( "pressure" );
	if ( param ) {
		*pressure = (Double)param->what.quantity.value;
		if ( strcmp( param->what.quantity.unit, "bar" ) == 0 ) {
			*pressure *= 1.0e5;
		}
	}

	param = GetParameter( "chi" );
	if ( param ) {
		*chi = (Double)param->what.quantity.value;
	}
	else {
		param = GetParameter( "chi_ref" );
		if ( param ) {
			*chi = (Double)param->what.quantity.value;
		}
	}

	param = GetParameter( "time" );
	if ( param ) {
		*theTime = (Double)param->what.quantity.value;
		if ( strcmp( param->what.quantity.unit, "ms" ) == 0 ) {
			*theTime *= 1.0e-3;
		}
	}

	param = GetParameter( "currtimestep" );
	if ( param ) {
		*currTimeStep = (Double)param->what.quantity.value;
		if ( strcmp( param->what.quantity.unit, "ms" ) == 0 ) {
			*currTimeStep *= 1.0e-3;
		}
	}
	
/*	if ( fSoot ) {*/
/*		// solver works with M/rho*/
/*		int		nSootMoments = fSoot->GetNSootMoments();*/
/*		int		momOff = fSoot->GetOffsetSootMoments();*/
/*		Double	mixMolarMass;*/
/*		Double	*solk;*/
/*		Double	density;*/
/*		for ( k = -1; k <= fNGridPoints; ++k ) {*/
/*			solk = y[k];*/
/*			fProperties->ComputeMixtureMolarMass( mixMolarMass*/
/*					, &solk[speciesOff], fSpecies->GetMolarMass()->vec, nSpeciesInSystem );*/
/*			density = ( GetPressure( fTime->vec[0] ) * mixMolarMass ) */
/*						/ ( RGAS * solk[tempOff] );*/
/*			for ( int i = 0; i < nSootMoments; ++i ) {*/
/*				solk[momOff+i] /= density;*/
/*			}*/
/*		}*/
/*	}*/

/*	FILE	*fp = GetOutfile( "initialguess", TFlame::kData );
	PrintSolution( fp, fNGridPoints, fNOfEquations, fSolGrid->vec, fSolution->mat, fVariableNames );
	fclose(fp);*/

/*// save initial solution to fSol*/
/*	for ( i = -1; i <= fNGridPoints; ++i ) {*/
/*		SaveSolution( i, fTime->vec[0], fSolution->mat[i] );*/
/*		// second call puts solution down to fSolOld*/
/*		SaveSolution( i, fTime->vec[0], fSolution->mat[i] ); */
/*	}*/
}

void TTransFlameSolver::SetInitialValues( int nGridPoints, int nVars, Double *x, Double **y )
{
	Double	*yLeft = y[kPrev];
	Double	*yRight = y[nGridPoints];
	Double	left = x[kPrev];
	Double	right = x[nGridPoints];

	for ( int k = 0; k < nGridPoints; ++k ) {
		for ( int j = 0; j < nVars; ++j ) {
			y[k][j] = Interpol( x[k], yLeft[j], left, yRight[j], right );
		}
	}
}

Double TTransFlameSolver::GetTotEnt( int k, Double *Z, Double t )
{
	Double	totentstart = Interpol( Z[k] , fTotEntStart->vec[kPrev], Z[kPrev], fTotEntStart->vec[fNGridPoints], Z[fNGridPoints] );
	Double	totentend = Interpol( Z[k] , fTotEntEnd->vec[kPrev], Z[kPrev] , fTotEntEnd->vec[fNGridPoints], Z[fNGridPoints] );

	return Interpol( t, totentstart, fTStart, totentend, fTEnd );
//	return Interpol( t, fTotEntStart->vec[k], fTStart, fTotEntEnd->vec[k], fTEnd );
}

#ifdef INGOOPT
void TTransFlameSolver::UpdateThermoProps( int k, Double *Y, Double temp, Double &pressure
								, Double &density, EqOfState what, Double *sootMoments )
{
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
	int		i;
	int		nSpeciesIn = fSpecies->GetNSpeciesInSystem();
	Flag	specChanged = FALSE;
		
	for ( i = 0; i < nSpeciesIn; ++i ) {
		if ( Y[i] != fYGSave[k][i] ) {
			specChanged = TRUE;
			copy( nSpeciesIn, Y, 1, fYGSave[k], 1 );
			break;
		}
	}

	if ( specChanged || temp != fTempGSave[k] ) {
		ComputeDeltaI( k, Y, temp );
		T0DFlame::UpdateThermoProps( Y, temp, pressure, density, what, sootMoments );
		return;
	}

	copy( nSpeciesIn, fDeltaI[k], 1, fSpecies->GetDeltaI()->vec, 1 );
	T0DFlame::UpdateThermoProps( Y, temp, pressure, density, what, sootMoments );
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo   
}

void TTransFlameSolver::UpdateThermoProps(int k, double *Y, double temp, double &pressure, double &density, EqOfState what, double *sootMoments, double * rhodot, double * enthdot)
{
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
	int		i;
	int		nSpeciesIn = fSpecies->GetNSpeciesInSystem();
	Flag	specChanged = FALSE;
		
	for ( i = 0; i < nSpeciesIn; ++i ) {
		if ( Y[i] != fYGSave[k][i] ) {
			specChanged = TRUE;
			copy( nSpeciesIn, Y, 1, fYGSave[k], 1 );
			break;
		}
	}

	if ( specChanged || temp != fTempGSave[k] ) {
		ComputeDeltaI( k, Y, temp );
		T0DFlame::UpdateThermoProps( Y, temp, pressure, density, what, sootMoments, rhodot, enthdot);
		return;
	}

	copy( nSpeciesIn, fDeltaI[k], 1, fSpecies->GetDeltaI()->vec, 1 );
	T0DFlame::UpdateThermoProps( Y, temp, pressure, density, what, sootMoments, rhodot, enthdot);
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo   
}

#else
void TTransFlameSolver::UpdateThermoProps( int k, Double *Y, Double temp, Double &pressure
								, Double &density, EqOfState what, Double *sootMoments )
{
	ComputeDeltaI( k, Y, temp );
	T0DFlame::UpdateThermoProps( Y, temp, pressure, density, what, sootMoments ); 
}
#endif

void TTransFlameSolver::ComputeDeltaI( int k, Double *Y, Double temp )
{
	int 	i;
	int		nSpeciesIn = fSpecies->GetNSpeciesInSystem();
	Double	*deltaI = fSpecies->GetDeltaI()->vec;
	
	CompDeltaIG( k, temp );

	for ( i = 0; i < nSpeciesIn; ++i ) {
		deltaI[i] = CompOneDeltaI( i, k, Y );
	}
	copy( nSpeciesIn, deltaI, 1, fDeltaI[k], 1 );
}

Double TTransFlameSolver::CompOneDeltaI( int i, int k, Double *Y )
{
	int		nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();
	Double	*M = fSpecies->GetMolarMass()->vec;
	Double	Delta_i = 0;

	for ( int j = 0; j < nSpeciesInSystem; ++j ){
		Delta_i += fG_ij[k][i][j] / M[j] * Y[j];
	}
	
	return Delta_i * M[i];
}

void TTransFlameSolver::CompDeltaIG( int k, Double temp )
{
	int		nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();
	Double	*mu = fSpecies->GetViscosity()->vec;
	
	if ( temp == fTempGSave[k] ) {
		return;
	}
	
	fTempGSave[k] = temp;

	fSpecies->ComputeSpeciesProperties( temp );

	Double	**sqrtMjOverMi = fSpecies->GetSqrtMjOverMi();
	Double	**sqrtMiOverMjPlusOne = fSpecies->GetSqrtMiOverMjPlusOne();
	for ( int i = 0; i < nSpeciesInSystem; ++i ) {
		for ( int j = 0; j < nSpeciesInSystem; ++j ) {
//#ifdef SPEED
		  fG_ij[k][i][j] = sqrt( sqrtMjOverMi[j][i] * mu[i] / mu[j] ) + 1.0; //Bug fixed here!!!!!
			fG_ij[k][i][j] *= fG_ij[k][i][j] * sqrtMiOverMjPlusOne[i][j];
//#else		
//			fprintf( stderr, "CompDeltaIG not yet implemented\n" );
//			exit( 2 );
//#endif
		}
	}
}

Double TTransFlameSolver::GetU( Double t, Double Z )
{
	int	i = 1;
	int	j = 1;
	int	zlen = fZIn->len;
	int	tlen = fTimeIn->len;
	Double	*zIn = fZIn->vec;
	Double	*tIn = fTimeIn->vec;
	Double	**UIn = fUstOverU->mat;
	Double	UAtZt1, UAtZt2;

	while ( j < zlen-1 && Z > zIn[j] ) ++j;
	if ( j == zlen-1 && Z > zIn[j]+1.0e-8 ) {
		fprintf( stderr, "Z = %g out of range[%g,%g]\n", Z, zIn[0]
								, zIn[zlen-1] );
	}
	
	while ( i < tlen-1 && t > tIn[i] ) ++i;
	if ( i == tlen-1 && t > tIn[i] ) {
		fprintf( stderr, "t = %g out of range[%g,%g]: linear extrapolation\n", Z, tIn[0]
								, tIn[tlen-1] );
	}
	
	if ( zIn[j] == zIn[j-1] ) {
		fprintf( stderr, "something's wrong\n" );
	}
	else {
		UAtZt1 = Interpol( Z, UIn[i-1][j-1], zIn[j-1], UIn[i-1][j], zIn[j] );
		UAtZt2 = Interpol( Z, UIn[i][j-1], zIn[j-1], UIn[i][j], zIn[j] );
	}
	
	if ( tIn[i] == tIn[i-1] ) {
		fprintf( stderr, "something else's wrong\n" );
		return 0.0;
	}
	else {
		return Interpol( t, UAtZt1, tIn[i-1], UAtZt2, tIn[i] );
	}
}

#define M1 259200
#define IA1 7141
#define IC1 54773
#define RM1 (1.0/M1)
#define M2 134456
#define IA2 8121
#define IC2 28411
#define RM2 (1.0/M2)
#define M3 243000
#define IA3 4561
#define IC3 51349

Double ran1( int *idum )
{
	static int		iff = 0;
	static long		ix1, ix2, ix3;
	static Double	r[98];
	Double			temp;
	int				j;

	if (*idum < 0 || iff == 0) {
		iff=1;
		ix1=(IC1-(*idum)) % M1;
		ix1=(IA1*ix1+IC1) % M1;
		ix2=ix1 % M2;
		ix1=(IA1*ix1+IC1) % M1;
		ix3=ix1 % M3;
		for (j=1;j<=97;j++) {
			ix1=(IA1*ix1+IC1) % M1;
			ix2=(IA2*ix2+IC2) % M2;
			r[j]=(ix1+ix2*RM2)*RM1;
		}
		*idum=1;
	}
	ix1=(IA1*ix1+IC1) % M1;
	ix2=(IA2*ix2+IC2) % M2;
	ix3=(IA3*ix3+IC3) % M3;
	j=1 + ((97*ix3)/M3);
	if (j > 97 || j < 1) fprintf(stderr,"RAN1: This cannot happen.");
	temp=r[j];
	r[j]=(ix1+ix2*RM2)*RM1;
	return temp;
}

#undef M1
#undef IA1
#undef IC1
#undef RM1
#undef M2
#undef IA2
#undef IC2
#undef RM2
#undef M3
#undef IA3
#undef IC3


Double ran1( int *idum );
Double expdev( int *idum )
{
	return -log(ran1(idum));
}


Double gasdev( int *idum )
{
	static int iset=0;
	static Double gset;
	Double fac,r,v1,v2;

	if  (iset == 0) {
		do {
			v1=2.0*ran1(idum)-1.0;
			v2=2.0*ran1(idum)-1.0;
			r=v1*v1+v2*v2;
		} while (r >= 1.0);
		fac=sqrt(-2.0*log(r)/r);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}

Double TTransFlameSolver::GetRosseRadiation( int k, Double *nTemp, Double **moments, Double rfPrev
											, Double rfCurr, Double rfNext )
{
// rf is rho * chi * cp / lambda
  Double rossecutoff = 5.0e-7;
  if ( moments[kCurr][1]*24.0/1800.0 > rossecutoff ) { 
					// coeffs negative
	Double	coeff = fSoot->GetSootRadRossCoeff( nTemp[kCurr], moments[kCurr], rossecutoff );
	Double	coeffPlus = fSoot->GetSootRadRossCoeff( nTemp[kNext], moments[kNext], rossecutoff );
	Double	coeffMinus = fSoot->GetSootRadRossCoeff( nTemp[kPrev], moments[kPrev], rossecutoff );
	Double	delq_R = 0.5 * coeff * rfCurr
	  * ( fWMinus[k] * nTemp[kPrev]
		 + fWCurr[k] * nTemp[kCurr]
		 + fWPlus[k] * nTemp[kNext] )
		+ 0.25 
		  * ( fFDWMinus[k] * rfPrev * coeffMinus
			 + fFDWCurr[k] * rfCurr * coeff 
			 + fFDWPlus[k] * rfNext * coeffPlus )
			+ 0.25 * rfCurr
			  * ( fFDWMinus[k] * coeffMinus 
				 + fFDWCurr[k] * coeff 
				 + fFDWPlus[k] * coeffPlus );
	return delq_R;
  }
  else {
	return fSoot->GetSootRadiation( nTemp[kCurr], moments[kCurr] );
  }
}


