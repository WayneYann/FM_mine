void CountPremSimJacFirst( void *object, NodeInfoPtr nodeInfo );
void CountPremSimJacRest( void *object, NodeInfoPtr nodeInfo );
void CountPremSimRHSRest( void *object, NodeInfoPtr nodeInfo, RHSMode rhsMode );
void CountPremSimOutput( void *object, FILE *fp, char* tail );
int CountPremSimPostIter( void *object );
void CountPremSimUpdateLeftBoundary( void  *object );
void CountPremSimUpdateRightBoundary( void *object );
void SetCountPremSimNodeInfo( int k, void *object );
void CountPremSimPostConv( void *object );
ConstStringArray GetCountPremSimVarNames( void *object );
// the following function is used as utility function for dfdyUpwind
// the void * is an object of the class T1DFlame
// Double Numerical_rhoDk( int equation, NodeInfoPtr nodeInfo, void *object, Flag CalcNewProperties = FALSE );
Double EntFluxFuncPremSim( int equation, NodeInfoPtr nodeInfo, void *object, Flag theFlag = FALSE );


//*************************************************************************************************************
class TCountPremFlameSim : public T1DFlame, public TPremixed {
friend void CountPremSimJacFirst( void *object, NodeInfoPtr nodeInfo );
friend void CountPremSimJacRest( void *object, NodeInfoPtr nodeInfo );
friend void CountPremSimRHSRest( void *object, NodeInfoPtr nodeInfo, RHSMode rhsMode );
friend void CountPremSimOutput( void *object, FILE *fp, char* tail );
friend int CountPremSimPostIter( void *object );
friend void CountPremSimUpdateLeftBoundary( void  *object );
friend void CountPremSimUpdateRightBoundary( void *object );
friend void CountPremSimPostConv( void *object );
friend Double EntFluxFuncPremSim( int equation, NodeInfoPtr nodeInfo, void *object, Flag theFlag  );

public:
	TCountPremFlameSim( FirstInputPtr firstInp ) : fVVelocity(0), fUVelocity(1), 
								fFirstSpecies(2), fTemperature(fFirstSpecies + fInputData->GetCounter()->species - fInputData->GetCounter()->steadyStates),
								fSootMoments(fTemperature+1), fVariablesWithoutSpecies(fInputData->fVariablesWithoutSpecies),
								T1DFlame( firstInp ), TPremixed( fInputData ) { InitCountPremFlameSim();};	
	~TCountPremFlameSim( void );

	int		GetOffsetVVelocity( void );
	int		GetOffsetUVelocity( void );
	int		GetOffsetTemperature( void );
	int		GetOffsetFirstSpecies( void );
	int		GetVariablesWithoutSpecies( void );
	ConstStringArray	GetVariableNames( void );
// not used by TCountPremFlame but required by TFlame as pure virtual functions
	int		GetOffsetMixFrac( void );
	
	VectorPtr	GetV( void ) { return fSolV; };
	VectorPtr	GetU( void ) { return fSolU; };

	void	SaveSolution( void );
	void	RestoreSolution( void );
	void	PrintRHSSpecies( TNewtonPtr bt );
	Double	GetBurningVelocity( void );

protected:
	void	PrintRHSTemp( TNewtonPtr bt );

	const int			fVVelocity;
	const int			fUVelocity;
	const int			fFirstSpecies;
	const int			fTemperature;
	const int			fVariablesWithoutSpecies;
	char				**fVariableNames;
	Flag				fPrintMolarFractions;
	TMassFractionPtr	fMassFraction;
	const int			fSootMoments;
	
//	vectors of solution
	VectorPtr			fSolV;			// length of these vectors is nGridPoints + 2
	VectorPtr			fSolU;			// len is set to nGridPoints
										// first element is fSol->vec[-1]

	VectorPtr			fSavedV;			// length of these vectors is nGridPoints + 2
	VectorPtr			fSavedU;			// len is set to nGridPoints
	VectorPtr			fSavedMixFrac;		// first element is fSol->vec[-1]

	void	UpdateSolution( MatrixPtr yMat, VectorPtr yLeftVec, VectorPtr yRightVec );
	void	UpdateSolutionOnePoint( Double *y, int gridPoint );

private:
  	void	InitCountPremFlameSim( void );
	void	SetInitialBC( TGridPtr grid, TInputDataPtr inp );
	void	SetInitialValues( TInputDataPtr inp, StartProfilePtr sp );
	void	FilldTdxdYdxOverTT( Double coeff, NodeInfoPtr nodeInfo );
	void	FilldTdxdYdxOverTTUpwind( Double coeff, NodeInfoPtr nodeInfo );
	void	FillImplicitSpeciesDiffusion( int speciesIndexk, NodeInfoPtr nodeInfo );
	void	FillJacSpeciesDiffusion( int nVariable, Double constCoeff, NodeInfoPtr nodeInfo, Flag sign = kPositive );
	void	FillJacFullDiffusion( int lInd, int kInd, NodeInfoPtr nodeInfo );
	Double	NewDiffCorr( int nVariable, NodeInfoPtr nodeInfo );
	Double	FullSpeciesDiffusion( int lInd, int kInd, NodeInfoPtr nodeInfo );
	Double	ImplicitSpeciesDiffusion( int nVariable, NodeInfoPtr nodeInfo );
	Double	SpeciesDiffusion( int nVariable, NodeInfoPtr nodeInfo );
	Double	SecondDerivSpeciesDiffusion( int nVariable, NodeInfoPtr nodeInfo );
	void	PrintRHSSpecies( TNewtonPtr bt, NodeInfoPtr nodeInfo, Double physX, FILE *fp );
	void	PrintRHSTemp( TNewtonPtr bt, NodeInfoPtr nodeInfo, Double physX, FILE *fp );
	Double	ComputeAbsStrainRate( int k1, Double f1, Double f2, Double fp1, Double fp2, Double h );
	void	FillJacDiffCorr( int nVariable, Double constCoeff, NodeInfoPtr nodeInfo, Flag sign = kPositive );
	Double	Numerical_dRhoDkdeta( int equation, NodeInfoPtr nodeInfo, Flag CalcNewProperties );
	void	EtaToX( TNewtonPtr bt, VectorPtr xPhysVec );
	FILE	*GetOutputFile( char *head, char *tail, FileType type );
	void	UpdateDimensions( int len );
	void	SolutionToSolver( void );
	int		CheckComputationalDomain( void );
};

typedef TCountPremFlameSim *TCountPremFlameSimPtr;
