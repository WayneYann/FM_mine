#ifdef SUN
#define UNDERSCORE
#elif defined SILGRAPH
#define UNDERSCORE
#elif defined DIGITAL
#define UNDERSCORE
#elif defined SUNGXX
#define UNDERSCORE
#elif defined LINUXGXX
#define UNDERSCORE
#endif

#ifdef UNDERSCORE
#define ADIABFLAMETEMP adiabflametemp_
#endif

extern "C" {
	void ADIABFLAMETEMP(Double *z, int *komp, Double *c0, Double *temp0, Double *p, char *symbol, 
                                 Double *xmol_i, Double *Enth0,Double *ponalo, Double *ponahi,
                                 int *nndd,
                                 Double *adiaT,Double *adiaY);
}

//*************************************************************************************************************
void Res0DIsoChor( Double *T, Double *y, Double *yPrime, Double *delta
			, int *iRes, Double *rPar, int *iPar );
void Res0DIsoBar( Double *T, Double *y, Double *yPrime, Double *delta
			, int *iRes, Double *rPar, int *iPar );

class T0DIsoChor : public T0DFlame, public TPremixed {
friend void Res0DIsoChor( Double *T, Double *y, Double *yPrime, Double *delta
			, int *iRes, Double *rPar, int *iPar );
friend void Res0DIsoBar( Double *T, Double *y, Double *yPrime, Double *delta
			, int *iRes, Double *rPar, int *iPar );
public:
	T0DIsoChor( FirstInputPtr firstInp ) : 
			fFirstSpecies(0), fTemperature(fFirstSpecies + fInputData->GetCounter()->species - fInputData->GetCounter()->steadyStates),
			fVariablesWithoutSpecies(fInputData->fVariablesWithoutSpecies),
			fSootMoments(fTemperature+1),
			T0DFlame( firstInp ), TPremixed( fInputData ) { InitT0DIsoChor();};	
	virtual ~T0DIsoChor( void );

	int			GetSolution( void );
	int			GetOffsetFirstSpecies( void ) { return fFirstSpecies; };
	int			GetOffsetTemperature( void ) { return fTemperature; };
	int			GetOffsetSootMoments( void ) { return fSootMoments; };
	int			GetVariablesWithoutSpecies( void ) { return fVariablesWithoutSpecies; };
	ConstStringArray	GetVariableNames( void ) { return fVariableNames; };

	void		Solve( void );

private:
	enum ReInitType { kTemp, kPhi, kPress, kSensAnal };
	
  	void	InitT0DIsoChor( void );
	void	SetInitialConditions( Double *y, TInputDataPtr inp );
	Flag	FirstStep( void );
	Flag	OneStep( void );
	void	SaveSolution( int count, Double time, Double *y, Flag doIt = FALSE );
	void	RestoreSolution( int count, Double &time, Double *y );
	void	SaveAdditional( int count );
	void	Reduce( void );
	void	ReduceAdditional( int toLength );
	void	SaveMoleFracDiff( int jTemp );
	void	WriteMoleFracDiff( char *tail = NULL );
	Double	GetTIgnition( void );
	void	WriteSolution( char *tail = NULL );
	void	PrintAdditional( char *tail );
	void	PrintHeatRelease( char *tail );
	void	PrintEpsilonM1( char *tail );
	void	PrintConsPro( char *tail );
	void	PrintDetProdRate( char *tail );
	void	PrintIntProdRate( char *tail );
	void	PrintReacRateSave( char *tail );
	void	PrintEquilibrium( char *tail );
	void    WriteFlameletFile( char *tail );

	//PP
	void	PrintHeatRelSave( char *tail );	
	//PP

	Double	GetEquilibrium( Double temp, Double *massFracs );
	void	WriteIgnDelTimes( void );
	void	ReInit( ReInitType what );
	void	InitDasslWorkSpace( void );
	void	ShowSolution( void );
	void	SetInfo( void );
	void	SetDensityPres( Double *Y, Double temp );
	void	SetInitialPressure( void ) { fInitialPressure = fPressure->vec[fPressure->len]; };
	int		GetLenOfDetProdRate( char **names, int quant );
	FILE 	*GetOutputFile( char *head, char *tail, FileType type );
	Double	GetPressure( void ) { cerr << "#error: function GetPressure not allowed" << NEWL; return -1; };
	Flag	IsLastFreqFac( void );
	void	NextRate( void );
	void	SaveSensCoeff( void );
	void	WriteSensCoeffs( void );
	void    WriteSortedSensCoeffs( void );
	void	(*HomoResFunc)( Double *T, Double *y, Double *yPrime, Double *delta
			, int *iRes, Double *rPar, int *iPar );
	Double	GetDPDt( Double time );


	const int	fFirstSpecies;
	const int	fTemperature;
	const int	fVariablesWithoutSpecies;
	const int	fSootMoments;

	// input
	FlameType	fType;
	Flag		fKeepMassFracs;
	Flag		fTempDiff;
	Flag		fWriteSolution;
	Flag		fAdditionalOutput;
	Flag		fEquidistant;
	Flag		fPrintMolarFractions;
	int		fNOutput;
	Double		fTStart;
	Double		fTEnd;
	Flag            fImposeTEnd;
	Double		fRTol;
	Double		fATol;
	
	Double		fEGR;
	Double		fCHMax;
	int		fCH;
	Flag            fError;
	Flag		fInvInitialTemp;
	int		fMaxOrd;
	char		**fVariableNames;
	int		fDeltaStepSave;
	Double		fMaxStepSize;
	Double		fTotEnt;
	int		*fNActualOrd;
	int		*fNActualStep;
	int		fNOfEquations;
	int		fActLength;
	Double		fInitialPressure;
	Double		fInitialTemp;
	Double		fArtificialSource;
	Double		fDensity;
	Double		fDeltaT;
	Double		fTout;
	Double		fTime;
	VectorPtr	fSolution;	// workspace for vector of solution
	VectorPtr	fSolPrime;	// workspace for partial derivatives

	VectorPtr	fSolTime;
	MatrixPtr	fSolMassFracs;
	VectorPtr	fSolTemp;
	MatrixPtr	fSolSootMoments;

	MatrixPtr	fConsPro;	// last element is fConsPro[timesteps][3 * nOfSpeciesIn]
	MatrixPtr	fEpsilonM1;	// means eps^(-1) = |cons|/|cons-prod|; 
							// last element is fEpsilonM1[timesteps][nOfSpeciesIn]
	VectorPtr	fHeatRelease;
	MatrixPtr	fDetProdRate;// lenght sum_{i over nObj}( nUsedReactions[i] )
	MatrixPtr	fReacRateSave;// lenght is nOfReactions

	//PP
	MatrixPtr	fHeatRelSave;// lenght is nOfReactions
	//PP
	
	MatrixPtr	fSaveMoleFracDiff;
	VectorPtr	fSaveTempDiff;

	TensorPtr	fIgnDelTimes;

//	variables used for dassl
	int		*fInfo;
	int		fIdid;
	int		fLRW;
	Double		*fRWork;	// workspace for dassl
	int		fLIW;
	int		*fIWork;	// workspace for dassl

//	variables used for udassl
	int		*fIsen;
	Double		fRTolS;
	Double		fRPar;
	Double		fATolS;
	int		fLSW;
	Double		*fSWork;	// workspace for udassl
	int		fLISW;
	int		*fISWork;	// workspace for udassl

// variables used for sensitivity analysis
	int		fDistFreqFac;
	int		fFirstSensRate;
	Double		fSavedFreqFac;
	Double		fSavedFreqFacBack;
	Double		fIgnDelOrig;
	VectorPtr	fSaveAllRates;
	VectorPtr       fDataRef;            // Vector storing results from first computation
	MatrixPtr       fSensCoeffs;         // Matric containing sensibility coefficients
	VectorPtr       fStoreTIg;           // Store Tig for each analysis


};

typedef T0DIsoChor *T0DIsoChorPtr;
