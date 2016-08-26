//#define FOOL_SOFTBENCH( x ) 

#include "FlameMaster.h"
#include "ListTool.h"
#include "Spline.h"
#include "TCountPremFlameSim.h"

#undef IMPLICITDIFFUSION

#define UPWINDCONVECTION

#undef TESTCONSISTENCE

#undef PRINTRHSTERMS

#define DIFFUSIVITYCORRECTION

#define ENTHALPYFLUX

#undef ENTFLUXINJAC

void TCountPremFlameSim::InitCountPremFlameSim( void )
{
	int i;
	TBVPSolverPtr	solver = GetSolver();
	TNewtonPtr		bt = solver->bt;
	TGridPtr		fine = bt->GetGrid()->GetFine();
	TGridPtr		coarse = bt->GetGrid()->GetCoarse();
	int				nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();
	int				maxGridPoints = bt->GetMaxGridPoints();

	if ( fSoot ) {
		fSoot->SetMomentsOffset( fSootMoments );
	}

	fVariableNames = new String[fVariablesWithoutSpecies + nSpeciesInSystem];

	fVariableNames[fVVelocity] = new char[2];
	strcpy( fVariableNames[fVVelocity], "f" );
	fVariableNames[fUVelocity] = new char[3];
	strcpy( fVariableNames[fUVelocity], "f'" );
	fVariableNames[fTemperature] = new char[2];
	strcpy( fVariableNames[fTemperature], "T" );
	for ( i = 0; i < nSpeciesInSystem; ++i ) {
		fVariableNames[fFirstSpecies + i] = new char[strlen( fSpecies->GetNames()[i] ) + 1];
		strcpy( fVariableNames[fFirstSpecies + i], fSpecies->GetNames()[i] );
	}
	if ( fSoot ) {
		int	offset = fSoot->GetOffsetSootMoments();
		for ( i = 0; i < fSoot->GetNSootMoments(); ++i ) {
			fVariableNames[offset + i] = new char[8];
			sprintf( fVariableNames[offset + i], "M%d/rho", i );
		}
	}
	fPrintMolarFractions = fInputData->fPrintMolarFractions;
	
//	vectors of solution
	fSolV = NewVector( maxGridPoints + 2 );
	fSolU = NewVector( maxGridPoints + 2 );

	fSolV->vec = &fSolV->vec[kNext];
	fSolU->vec = &fSolU->vec[kNext];

	fSolV->len -= 2;
	fSolU->len -= 2;

//	saved solution
	fSavedV = NewVector( maxGridPoints + 2 );
	fSavedU = NewVector( maxGridPoints + 2 );

	fSavedV->vec = &fSavedV->vec[kNext];
	fSavedU->vec = &fSavedU->vec[kNext];

	fSavedV->len -= 2;
	fSavedU->len -= 2;

	bt->SetUtFuncs( CountPremSimJacFirst, CountPremSimJacRest, CountPremSimJacRest
					, CountPremSimRHSRest, CountPremSimRHSRest, CountPremSimRHSRest
					, CountPremSimOutput, CountPremSimPostIter
					, SetCountPremSimNodeInfo, CountPremSimPostConv
					, GetCountPremSimVarNames
					, CountPremSimUpdateLeftBoundary, CountPremSimUpdateRightBoundary );
	SetInitialBC( fine, fInputData );
	SetInitialBC( coarse, fInputData );
	cerr << "initial equivalence ratio is " << GetPhi() << NEWL;
/*	fMassFraction = new TMassFraction( this );*/
/*	if ( !fMassFraction ) FatalError( "memory allocation of TMassFraction failed" );*/
	//	fMassFraction->PrintMassFraction();

	ReadStartProfiles( fInputData );
	CheckBC();
	CheckInitialGuess();
	UpdateSolution( fine->GetY(), fine->GetYLeft(), fine->GetYRight() );
}

TCountPremFlameSim::~TCountPremFlameSim( void )
{
	int	nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();

	delete fMassFraction;

	fSavedU->vec = &fSavedU->vec[kPrev];
	fSavedV->vec = &fSavedV->vec[kPrev];

	DisposeVector( fSavedU );
	DisposeVector( fSavedV );

	fSolV->vec = &fSolV->vec[kPrev];
	fSolU->vec = &fSolU->vec[kPrev];

	DisposeVector( fSolU );
	DisposeVector( fSolV );

	for ( int i = 0; i < nSpeciesInSystem+fVariablesWithoutSpecies; ++i ) {
		delete fVariableNames[i];
	}
	delete fVariableNames;
}

void CountPremSimJacFirst( void *object, NodeInfoPtr nodeInfo )
{
	TCountPremFlameSimPtr	flame = ( TCountPremFlameSimPtr )object;
	TFlameNodePtr	flameNode = flame->fFlameNode;
	int 	fFirstSpecies = flame->GetOffsetFirstSpecies();
	int 	fTemperature = flame->GetOffsetTemperature();
	int 	fUVelocity = flame->GetOffsetUVelocity();
	int 	fVVelocity = flame->GetOffsetVVelocity();
	int		M = nodeInfo->nOfEquations;
	int		nSpeciesInSystem = flame->GetSpecies()->GetNSpeciesInSystem();
	int		speciesEq, speciesVar, speciesIndexEq, speciesIndexVar;
	int		lastSpeciesEq = nSpeciesInSystem + fFirstSpecies;
    Double  h = nodeInfo->h;
    Double  hnenn = nodeInfo->hnenn;
	Double	**a = nodeInfo->a;
	Double	*temp = flameNode->temp;
	Double	*y = nodeInfo->y;
	Double	strainRate = flame->GetStrainRate();
	Double	*mixDensity = flameNode->mixDensity;
	Double	*mixHeatCapacity = flameNode->mixHeatCapacity;
	Double	*mixViscosity = flameNode->mixViscosity;
	Double	*mixConductivity = flameNode->mixConductivity;
	Double	*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
	Double	*enthalpy = flameNode->enthalpy;
	Double	*heatCapacity = flameNode->heatCapacity;
	Double	*diffusivity = flameNode->diffusivity;
	Double	*diffusivityPrev = flameNode->diffusivityPrev;
	Double	*productionRate = flameNode->productionRate;
	Double	mixDensityInf = flameNode->rhoInf;
	Double	mixViscosityInf = flameNode->viscosityInf;
	Double	constMassDiffCoeff = 1.0 / ( mixDensityInf * mixViscosityInf );
	Double	constThermDiffCoeff = constMassDiffCoeff / mixHeatCapacity[kCurr];
	Double	oneOverATimesRho = 1.0 / ( strainRate * mixDensity[kCurr] );
	Double	mixMolarMass = *flameNode->mixMolarMass;
	Double	idealGasCoeff = flame->GetPressure() * mixMolarMass / RGAS; // rho = idealGasCoeff / T
//	Double	ROverAP = RGAS / ( strainRate * flame->GetPressure() );
//	Double	ROverAPM = ROverAP / mixMolarMass;
//	Double	energySourceCoeff = ROverAPM / mixHeatCapacity[kCurr];
	Double	energySourceCoeff = oneOverATimesRho / mixHeatCapacity[kCurr];
	Double	**dMdY = flameNode->dMdY;
	Double	*dMdT = flameNode->dMdY[nSpeciesInSystem];

	flame->FilldMdYOnePoint( flameNode );
	flame->FilldMdTOnePoint( flameNode );

	if ( flame->GetSoot() ) {
		flame->GetSoot()->UpdateJacobian( flame );
		flame->GetSoot()->FillJacobi( flame, nodeInfo, kSimilarity );
	}

// first fill all convection terms
	// first equation ( mass )
	FillJacFirstDerivUp( fVVelocity, fVVelocity, nodeInfo );

#ifdef UPWINDCONVECTION
	// second equation ( momentum )
	FillJacNonlinearConvectUpwind( fVVelocity, fUVelocity, nodeInfo, 1.0, FALSE );
	// third to three + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		FillJacNonlinearConvectUpwind( fVVelocity, speciesEq, nodeInfo, 1.0, FALSE );
	}
	if ( fTemperature < M ) {
		FillJacNonlinearConvectUpwind( fVVelocity, fTemperature, nodeInfo, 1.0, FALSE );
	}
#else
	// second equation ( momentum )
	FillJacNonlinearConvectCentral( fVVelocity, fUVelocity, nodeInfo );
	a[fUVelocity][fUVelocity] -= h * h * y[fVVelocity];
	// third to three + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		FillJacNonlinearConvectCentral( fVVelocity, speciesEq, nodeInfo );
		a[speciesEq][speciesEq] -= h * h * y[fVVelocity];
	}
	if ( fTemperature < M ) {
		FillJacNonlinearConvectCentral( fVVelocity, fTemperature, nodeInfo );
		a[fTemperature][fTemperature] -= h * h * y[fVVelocity];
	}
#endif

	
// first equation ( mass )
	a[fUVelocity][fVVelocity] -= hnenn;
	
// second equation ( momentum )
	flame->FillJacDiffusion( fUVelocity, fUVelocity, constMassDiffCoeff, flameNode->mixViscosity, nodeInfo );
	a[fUVelocity][fUVelocity] += constMassDiffCoeff * ( mixViscosity[kPrev] * mixDensity[kPrev]
								 + mixViscosity[kCurr] * mixDensity[kCurr] ) * h;
	a[fUVelocity][fUVelocity] -= 2.0 * y[fUVelocity] * hnenn;
	a[fTemperature][fUVelocity] += mixDensityInf / ( mixDensity[kCurr] * temp[kCurr] ) * hnenn;
	
	Double	theConst = mixDensityInf * mixMolarMass / mixDensity[kCurr] * hnenn;
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq; ++speciesEq ) {
		a[speciesEq][fUVelocity] += theConst / molarMass[speciesEq-fFirstSpecies];
	}
	
// third to three + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		speciesIndexEq = speciesEq - fFirstSpecies;
		flame->FillJacSpeciesDiffusion( speciesEq, constMassDiffCoeff, nodeInfo );
		a[speciesEq][speciesEq] += constMassDiffCoeff * ( diffusivityPrev[speciesIndexEq] * mixDensity[kPrev] * mixDensity[kPrev]
					+ diffusivity[speciesIndexEq] * mixDensity[kCurr] * mixDensity[kCurr] ) * h;

#ifdef DIFFUSIVITYCORRECTION
		if ( flame->UseDiffCorr() ) {
			flame->FillJacDiffCorr( speciesEq, constMassDiffCoeff, nodeInfo, kNegative );
		}
#endif

		//	implicit source term
		for ( speciesVar = fFirstSpecies; speciesVar < lastSpeciesEq && speciesVar < M; ++ speciesVar ) {
			speciesIndexVar = speciesVar - fFirstSpecies;
			a[speciesVar][speciesEq] += ( dMdY[speciesIndexVar][speciesIndexEq]
					+ productionRate[speciesIndexEq] * mixMolarMass / molarMass[speciesIndexVar] )
					* oneOverATimesRho * hnenn;
		}
		a[fTemperature][speciesEq] += productionRate[speciesIndexEq] / y[fTemperature]
										* oneOverATimesRho * hnenn;
		a[fTemperature][speciesEq] += dMdT[speciesIndexEq]
										* oneOverATimesRho * hnenn;
	}

// last equation ( temperature )
	if ( fTemperature < M ) {
		int 	i, l;
		flame->FillJacDiffusion( fTemperature, fTemperature, constThermDiffCoeff, flameNode->mixConductivity, nodeInfo );
		a[fTemperature][fTemperature] += constThermDiffCoeff * ( mixConductivity[kPrev] * mixDensity[kPrev] 
										 + mixConductivity[kCurr] * mixDensity[kCurr] ) * h;

		for ( i = 0; i < nSpeciesInSystem; ++i ) {
			//	d/dT sum (m_i * h_i ) 
			a[fTemperature][fTemperature] -= ( dMdT[i] * enthalpy[i] + 
						productionRate[i] * ( enthalpy[i] / temp[kCurr] + heatCapacity[i] ) ) 
						* energySourceCoeff * hnenn;
			for ( l = 0; l < nSpeciesInSystem; ++l ) {
				//	d/dY_l sum (m_i * h_i ) 
				a[l + fFirstSpecies][fTemperature] -= ( dMdY[l][i]
							+ productionRate[i] * mixMolarMass / molarMass[l] ) 
							* enthalpy[i] * energySourceCoeff * hnenn;
			}
		}
		
#ifdef ENTFLUXINJAC
#ifdef ENTHALPYFLUX
		Double	sumCpY = 0.0;
		Double	sumCpDdYdx = 0.0;
		Double	entFluxCoeff = constThermDiffCoeff * mixDensity[kCurr] * mixDensity[kCurr]
					* FirstDeriv( yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h );
		Double	entTerm;
		Double	**Y = flameNode->Y;
		for ( i = 0; i < nSpeciesInSystem; ++i ) {
			sumCpDdYdx += heatCapacity[i] * diffusivity[i] 
						* FirstDeriv( Y[kPrev][i], Y[kCurr][i], Y[kNext][i], hm, h );
#	ifdef DIFFUSIVITYCORRECTION
			if ( flame->UseDiffCorr() ) {
				sumCpY += heatCapacity[i] * Y[kCurr][i];
			}
#	endif
		}
		entTerm = sumCpDdYdx;
#	ifdef DIFFUSIVITYCORRECTION
		if ( flame->UseDiffCorr() ) {
			entTerm -= flameNode->diffCorr[kCurr] * sumCpY;
		}
#	endif
		entTerm *= entFluxCoeff;
		entTerm *= 2.0 * hnenn;
		for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq; ++speciesEq ) {
			i = speciesEq - fFirstSpecies;
			a[speciesEq][fTemperature] += dfdyUpwind( speciesEq, fTemperature, EntFluxFuncPremSim, nodeInfo, flame ) * hnenn;
			a[speciesEq][fTemperature] -= entTerm 
					* mixMolarMass / molarMass[speciesEq-fFirstSpecies];
		}
#endif
#endif
		
		if ( flame->fProperties->GetRadiation() ) {
			flame->fProperties->GetRadiation()->FillJacRadiation( energySourceCoeff, flame, nodeInfo );
		}
	}
}

void CountPremSimJacRest( void *object, NodeInfoPtr nodeInfo )
{
	TCountPremFlameSimPtr	flame = ( TCountPremFlameSimPtr )object;
	TFlameNodePtr	flameNode = flame->fFlameNode;
	int 	fFirstSpecies = flame->GetOffsetFirstSpecies();
	int 	fTemperature = flame->GetOffsetTemperature();
	int 	fUVelocity = flame->GetOffsetUVelocity();
	int 	fVVelocity = flame->GetOffsetVVelocity();
	int		M = nodeInfo->nOfEquations;
	int		nSpeciesInSystem = flame->GetSpecies()->GetNSpeciesInSystem();
	int		speciesEq, speciesVar, speciesIndexEq, speciesIndexVar;
	int		lastSpeciesEq = nSpeciesInSystem + fFirstSpecies;
    Double  hnenn = nodeInfo->hnenn;
	Double	**a = nodeInfo->a;
	Double	*temp = flameNode->temp;
	Double	*y = nodeInfo->y;
	Double	strainRate = flame->GetStrainRate();
	Double	mixDensity = *flameNode->mixDensity;
	Double	mixHeatCapacity = *flameNode->mixHeatCapacity;
	Double	*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
	Double	*enthalpy = flameNode->enthalpy;
	Double	*heatCapacity = flameNode->heatCapacity;
	Double	*productionRate = flameNode->productionRate;
	Double	mixDensityInf = flameNode->rhoInf;
	Double	mixViscosityInf = flameNode->viscosityInf;
	Double	constMassDiffCoeff = 1.0 / ( mixDensityInf * mixViscosityInf );
	Double	constThermDiffCoeff = constMassDiffCoeff / mixHeatCapacity;
	Double	oneOverATimesRho = 1.0 / ( strainRate * mixDensity );
	Double	mixMolarMass = *flameNode->mixMolarMass;
	Double	idealGasCoeff = flame->GetPressure() * mixMolarMass / RGAS; // rho = idealGasCoeff / T
//	Double	ROverAP = RGAS / ( strainRate * flame->GetPressure() );
//	Double	ROverAPM = ROverAP / mixMolarMass;
//	Double	energySourceCoeff = ROverAPM / mixHeatCapacity;
	Double	energySourceCoeff = oneOverATimesRho / mixHeatCapacity;
	Double	**dMdY = flameNode->dMdY;
	Double	*dMdT = flameNode->dMdY[nSpeciesInSystem];

	flame->FilldMdYOnePoint( flameNode );
	flame->FilldMdTOnePoint( flameNode );

	if ( flame->GetSoot() ) {
		flame->GetSoot()->UpdateJacobian( flame );
		flame->GetSoot()->FillJacobi( flame, nodeInfo, kSimilarity );
	}

// first fill all convection terms
	// first equation ( mass )
	FillJacFirstDerivUp( fVVelocity, fVVelocity, nodeInfo );

#ifdef UPWINDCONVECTION
	// second equation ( momentum )
	FillJacNonlinearConvectUpwind( fVVelocity, fUVelocity, nodeInfo, 1.0, FALSE );
	// third to three + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		FillJacNonlinearConvectUpwind( fVVelocity, speciesEq, nodeInfo, 1.0, FALSE );
	}
	if ( fTemperature < M ) {
		FillJacNonlinearConvectUpwind( fVVelocity, fTemperature, nodeInfo, 1.0, FALSE );
	}
#else
	// second equation ( momentum )
	FillJacNonlinearConvectCentral( fVVelocity, fUVelocity, nodeInfo );
	// third to three + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		FillJacNonlinearConvectCentral( fVVelocity, speciesEq, nodeInfo );
	}
	if ( fTemperature < M ) {
		FillJacNonlinearConvectCentral( fVVelocity, fTemperature, nodeInfo );
	}
#endif

	
// first equation ( mass )
	a[fUVelocity][fVVelocity] -= hnenn;
	
// second equation ( momentum )
	flame->FillJacDiffusion( fUVelocity, fUVelocity, constMassDiffCoeff, flameNode->mixViscosity, nodeInfo );
	a[fUVelocity][fUVelocity] -= 2.0 * y[fUVelocity] * hnenn;
	a[fTemperature][fUVelocity] += mixDensityInf / ( mixDensity * temp[kCurr] ) * hnenn;
	
	Double	theConst = mixDensityInf * mixMolarMass / mixDensity * hnenn;
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq; ++speciesEq ) {
		a[speciesEq][fUVelocity] += theConst / molarMass[speciesEq-fFirstSpecies];
	}
	
// third to three + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		speciesIndexEq = speciesEq - fFirstSpecies;
		flame->FillJacSpeciesDiffusion( speciesEq, constMassDiffCoeff, nodeInfo );
#ifdef DIFFUSIVITYCORRECTION
		if ( flame->UseDiffCorr() ) {
			flame->FillJacDiffCorr( speciesEq, constMassDiffCoeff, nodeInfo, kNegative );
		}
#endif
		//	implicit source term
		for ( speciesVar = fFirstSpecies; speciesVar < lastSpeciesEq && speciesVar < M; ++ speciesVar ) {
			speciesIndexVar = speciesVar - fFirstSpecies;
			a[speciesVar][speciesEq] += ( dMdY[speciesIndexVar][speciesIndexEq]
					+ productionRate[speciesIndexEq] * mixMolarMass / molarMass[speciesIndexVar] )
					* oneOverATimesRho * hnenn;
		}
		a[fTemperature][speciesEq] += productionRate[speciesIndexEq] / y[fTemperature]
										* oneOverATimesRho * hnenn;
		a[fTemperature][speciesEq] += dMdT[speciesIndexEq]
										* oneOverATimesRho * hnenn;
	}

// last equation ( temperature )
	if ( fTemperature < M ) {
		int 	i, l;
		flame->FillJacDiffusion( fTemperature, fTemperature, constThermDiffCoeff, flameNode->mixConductivity, nodeInfo );

		for ( i = 0; i < nSpeciesInSystem; ++i ) {
			//	d/dT sum (m_i * h_i ) 
			a[fTemperature][fTemperature] -= ( dMdT[i] * enthalpy[i] + 
						productionRate[i] * ( enthalpy[i] / temp[kCurr] + heatCapacity[i] ) ) 
						* energySourceCoeff * hnenn;
			for ( l = 0; l < nSpeciesInSystem; ++l ) {
				//	d/dY_l sum (m_i * h_i ) 
				a[l + fFirstSpecies][fTemperature] -= ( dMdY[l][i]
							+ productionRate[i] * mixMolarMass / molarMass[l] ) 
							* enthalpy[i] * energySourceCoeff * hnenn;
			}
		}
		
#ifdef ENTFLUXINJAC
#ifdef ENTHALPYFLUX
		Double	sumCpY = 0.0;
		Double	sumCpDdYdx = 0.0;
		Double	entFluxCoeff = constThermDiffCoeff * mixDensity * mixDensity
					* FirstDeriv( yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h );
		Double	entTerm;
		Double	**Y = flameNode->Y;
		for ( i = 0; i < nSpeciesInSystem; ++i ) {
			sumCpDdYdx += heatCapacity[i] * diffusivity[i] 
						* FirstDeriv( Y[kPrev][i], Y[kCurr][i], Y[kNext][i], hm, h );
#	ifdef DIFFUSIVITYCORRECTION
			if ( flame->UseDiffCorr() ) {
				sumCpY += heatCapacity[i] * Y[kCurr][i];
			}
#	endif
		}
		entTerm = sumCpDdYdx;
#	ifdef DIFFUSIVITYCORRECTION
		if ( flame->UseDiffCorr() ) {
			entTerm -= flameNode->diffCorr[kCurr] * sumCpY;
		}
#	endif
		entTerm *= entFluxCoeff;
		entTerm *= 2.0 * hnenn;
		for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq; ++speciesEq ) {
			i = speciesEq - fFirstSpecies;
			a[speciesEq][fTemperature] += dfdyUpwind( speciesEq, fTemperature, EntFluxFuncPremSim, nodeInfo, flame ) * hnenn;
			a[speciesEq][fTemperature] -= entTerm 
					* mixMolarMass / molarMass[speciesEq-fFirstSpecies];
		}
#endif
#endif
		
		if ( flame->fProperties->GetRadiation() ) {
			flame->fProperties->GetRadiation()->FillJacRadiation( energySourceCoeff, flame, nodeInfo );
		}
	}
}
/*
void CountPremSimRHSRest( void *object, NodeInfoPtr nodeInfo, RHSMode rhsMode )
{
	TCountPremFlameSimPtr	flame = ( TCountPremFlameSimPtr )object;
	if ( !flame->RHSAction( nodeInfo, rhsMode ) ) {
		return;
	}
	
	TFlameNodePtr	flameNode = flame->fFlameNode;
	TFlameNodePtr	flameNodeSaved = flame->fFlameNodeSaved;
	int 	fFirstSpecies = flame->GetOffsetFirstSpecies();
	int 	fTemperature = flame->GetOffsetTemperature();
	int 	fUVelocity = flame->GetOffsetUVelocity();
	int 	fVVelocity = flame->GetOffsetVVelocity();
	int		eqLoop, speciesEq;
	int		M = nodeInfo->nOfEquations;
	int		nSpeciesInSystem = flame->GetSpecies()->GetNSpeciesInSystem();
	int		lastSpeciesEq = nSpeciesInSystem + fFirstSpecies;
    Double  h = nodeInfo->h;
    Double  hm = nodeInfo->hm;
    Double  hnenn = nodeInfo->hnenn;
	Double	*rhs = nodeInfo->rhs;
	Double	*YPrev = flameNode->Y[kPrev];
	Double	*Y = flameNode->Y[kCurr];
	Double	*YNext = flameNode->Y[kNext];
	Double	*temp = flameNode->temp;
	Double	*yPrev = nodeInfo->yPrev;
	Double	*y = nodeInfo->y;
	Double	*yNext = nodeInfo->yNext;
	Double	yMax;
	Double	strainRate = flame->GetStrainRate();
	Double	*enthalpy = flameNode->enthalpy;
	Double	*mixDensity = flameNode->mixDensity;
	Double	mixDensityInf = flameNode->rhoInf;
	Double	mixViscosityInf = flameNode->viscosityInf;
	Double	*mixViscosity = flameNode->mixViscosity;
	Double	*productionRate = flameNode->productionRate;
	Double	*diffusivity = flameNode->diffusivity;
	Double	*heatCapacity = flameNode->heatCapacity;
	Double	oneOverATimesRho = 1.0 / ( strainRate * mixDensity[kCurr] );
	Double	idealGasCoeff = flame->GetPressure() * *flameNode->mixMolarMass / RGAS; // rho = idealGasCoeff / T
	Double	ROverAPM = 1.0 / ( strainRate * idealGasCoeff );
	Double	constMassDiffCoeff = 1.0 / ( mixDensityInf * mixViscosityInf );
	Double	mixHeatCapacity = *flameNode->mixHeatCapacity;
	Double	constThermDiffCoeff = constMassDiffCoeff / flameNode->mixHeatCapacity[kCurr];
	Double	*diffCorr = flameNode->diffCorr;
	Double	sumCpDdYdx;
	Double	sumMH;

	if ( flame->GetSoot() ) {
		flame->GetSoot()->FillRHS( flame, nodeInfo, kSimilarity );
	}
	
// first fill all convection terms
	// first equation ( mass )
	rhs[fVVelocity] += FirstDerivUpwind( y[fVVelocity], yPrev[fVVelocity], hm );

#ifdef UPWINDCONVECTION
	// second equation ( momentum )
	rhs[fUVelocity] += NonlinearConvectUpwind( y[fVVelocity], yPrev[fUVelocity], y[fUVelocity], yNext[fUVelocity], hm, h, FALSE );
	// third to three + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		yMax = nodeInfo->yMax[speciesEq-fFirstSpecies];
		rhs[speciesEq] += yMax * NonlinearConvectUpwind( y[fVVelocity], yPrev[speciesEq], y[speciesEq], yNext[speciesEq], hm, h, FALSE );
	}
	if ( fTemperature < M ) {
		rhs[fTemperature] += NonlinearConvectUpwind( y[fVVelocity], yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h, FALSE );
	}
#else
	// second equation ( momentum )
	rhs[fUVelocity] += NonlinearConvectCentral( y[fVVelocity], yPrev[fUVelocity], y[fUVelocity], yNext[fUVelocity], hm, h );
	// third to three + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		rhs[speciesEq] += NonlinearConvectCentral( y[fVVelocity], yPrev[speciesEq], y[speciesEq], yNext[speciesEq], hm, h );
	}
	if ( fTemperature < M ) {
		rhs[fTemperature] += NonlinearConvectCentral( y[fVVelocity], yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h );
	}
#endif

// mass equation
	rhs[fVVelocity] -= y[fUVelocity];

// momentum equation
	rhs[fUVelocity] += mixDensityInf / idealGasCoeff * y[fTemperature];
	rhs[fUVelocity] -= y[fUVelocity] * y[fUVelocity];
	rhs[fUVelocity] += constMassDiffCoeff * flame->StandardDiffusion( fUVelocity, mixViscosity, nodeInfo );
	
// fFirstSpecies to fFirstSpecies + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		eqLoop = speciesEq - fFirstSpecies;
		yMax = nodeInfo->yMax[eqLoop];
		rhs[speciesEq] += yMax * constMassDiffCoeff * flame->SecondDerivSpeciesDiffusion( speciesEq, nodeInfo );
#ifdef DIFFUSIVITYCORRECTION
		if ( flame->UseDiffCorr() ) {
			rhs[speciesEq] -= yMax * constMassDiffCoeff * flame->NewDiffCorr( speciesEq, nodeInfo );
		}
#endif
		if ( flame->fThermoDiffusion ) {
			rhs[speciesEq] += yMax * constMassDiffCoeff * flame->ThermoDiffusion( eqLoop, kSimilarity, nodeInfo );
		}
		rhs[speciesEq] += productionRate[eqLoop] * oneOverATimesRho;
	}
	
// energy equation
	if ( fTemperature < M ) {
		Double	sumCpY = 0.0;
		Double	corrCoeff = 0.0;
		Double	oneOverARhoCp = 1.0 / ( strainRate * mixDensity[kCurr] * mixHeatCapacity );
		sumCpDdYdx = 0.0;
		sumMH = 0.0;
		rhs[fTemperature] += constThermDiffCoeff * flame->StandardDiffusion( fTemperature, flameNode->mixConductivity, nodeInfo );
#ifdef DIFFUSIVITYCORRECTION
		if ( flame->UseDiffCorr() ) {
			corrCoeff = constMassDiffCoeff * mixDensity[kCurr] * mixDensity[kCurr]
									* diffCorr[kCurr] / mixHeatCapacity;
		}
#endif

		for ( eqLoop = 0; eqLoop < nSpeciesInSystem; ++eqLoop ) {
			sumCpDdYdx += heatCapacity[eqLoop] * diffusivity[eqLoop] 
						* FirstDeriv( YPrev[eqLoop], Y[eqLoop], YNext[eqLoop], hm, h );
#ifdef DIFFUSIVITYCORRECTION
			if ( flame->UseDiffCorr() ) {
				sumCpY += heatCapacity[eqLoop] * Y[eqLoop];
			}
#endif
			sumMH += productionRate[eqLoop] * enthalpy[eqLoop];
		}
#ifdef ENTHALPYFLUX
		rhs[fTemperature] +=  constThermDiffCoeff * sumCpDdYdx * mixDensity[kCurr] * mixDensity[kCurr]
						* FirstDeriv( yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h );
#	ifdef DIFFUSIVITYCORRECTION
		if ( flame->UseDiffCorr() ) {
			rhs[fTemperature] -= corrCoeff * sumCpY
						* FirstDeriv( yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h );
		}
#	endif
#endif
		rhs[fTemperature] -= sumMH * oneOverARhoCp;
		
		if ( flame->fProperties->GetRadiation() ) {
			rhs[fTemperature] += flameNode->radiation[kCurr] * oneOverATimesRho / mixHeatCapacity;
			if ( flame->GetSoot() ) {
				rhs[fTemperature] -= oneOverATimesRho / mixHeatCapacity * flame->GetSoot()->GetSootRadiation( y[fTemperature], flameNode->moments );
			}
		}
	}

//	for ( eqLoop = 0; eqLoop < M; ++eqLoop ){
//		rhs[eqLoop] *= - hnenn;
//	}
}
*/
void CountPremSimRHSRest( void *object, NodeInfoPtr nodeInfo, RHSMode rhsMode )
{
	TCountPremFlameSimPtr	flame = ( TCountPremFlameSimPtr )object;
	if ( !flame->RHSAction( nodeInfo, rhsMode ) ) {
		return;
	}
	
	TFlameNodePtr	flameNode = flame->fFlameNode;
	int 	fFirstSpecies = flame->GetOffsetFirstSpecies();
	int 	fTemperature = flame->GetOffsetTemperature();
	int 	fUVelocity = flame->GetOffsetUVelocity();
	int 	fVVelocity = flame->GetOffsetVVelocity();
	int		eqLoop, speciesEq;
	int		M = nodeInfo->nOfEquations;
	int		nSpeciesInSystem = flame->GetSpecies()->GetNSpeciesInSystem();
	int		lastSpeciesEq = nSpeciesInSystem + fFirstSpecies;
    Double  h = nodeInfo->h;
    Double  hm = nodeInfo->hm;
    Double  hnenn = nodeInfo->hnenn;
	Double	*rhs = nodeInfo->rhs;
	Double	*YPrev = flameNode->Y[kPrev];
	Double	*Y = flameNode->Y[kCurr];
	Double	*YNext = flameNode->Y[kNext];
	Double	*yPrev = nodeInfo->yPrev;
	Double	*y = nodeInfo->y;
	Double	*yNext = nodeInfo->yNext;
	Double	strainRate = flame->GetStrainRate();
	Double	*enthalpy = flameNode->enthalpy;
	Double	*mixDensity = flameNode->mixDensity;
	Double	mixDensityInf = flameNode->rhoInf;
	Double	mixViscosityInf = flameNode->viscosityInf;
	Double	*mixViscosity = flameNode->mixViscosity;
	Double	*productionRate = flameNode->productionRate;
	Double	*diffusivity = flameNode->diffusivity;
	Double	*heatCapacity = flameNode->heatCapacity;
	Double	oneOverATimesRho = 1.0 / ( strainRate * mixDensity[kCurr] );
	Double	idealGasCoeff = flame->GetPressure() * *flameNode->mixMolarMass / RGAS; // rho = idealGasCoeff / T
	Double	constMassDiffCoeff = 1.0 / ( mixDensityInf * mixViscosityInf );
	Double	mixHeatCapacity = *flameNode->mixHeatCapacity;
	Double	constThermDiffCoeff = constMassDiffCoeff / flameNode->mixHeatCapacity[kCurr];
	Double	*diffCorr = flameNode->diffCorr;
	Double	sumCpDdYdx;
	Double	sumMH;

	if ( flame->GetSoot() ) {
		flame->GetSoot()->FillRHS( flame, nodeInfo, kSimilarity );
	}
	
// first fill all convection terms
	// first equation ( mass )
	rhs[fVVelocity] += FirstDerivUpwind( y[fVVelocity], yPrev[fVVelocity], hm );

#ifdef UPWINDCONVECTION
	// second equation ( momentum )
	rhs[fUVelocity] += NonlinearConvectUpwind( y[fVVelocity], yPrev[fUVelocity], y[fUVelocity], yNext[fUVelocity], hm, h, FALSE );
	// third to three + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		rhs[speciesEq] += NonlinearConvectUpwind( y[fVVelocity], yPrev[speciesEq], y[speciesEq], yNext[speciesEq], hm, h, FALSE );
	}
	if ( fTemperature < M ) {
		rhs[fTemperature] += NonlinearConvectUpwind( y[fVVelocity], yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h, FALSE );
	}
#else
	// second equation ( momentum )
	rhs[fUVelocity] += NonlinearConvectCentral( y[fVVelocity], yPrev[fUVelocity], y[fUVelocity], yNext[fUVelocity], hm, h );
	// third to three + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		rhs[speciesEq] += NonlinearConvectCentral( y[fVVelocity], yPrev[speciesEq], y[speciesEq], yNext[speciesEq], hm, h );
	}
	if ( fTemperature < M ) {
		rhs[fTemperature] += NonlinearConvectCentral( y[fVVelocity], yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h );
	}
#endif

// mass equation
	rhs[fVVelocity] -= y[fUVelocity];

// momentum equation
	rhs[fUVelocity] += mixDensityInf / idealGasCoeff * y[fTemperature];
	rhs[fUVelocity] -= y[fUVelocity] * y[fUVelocity];
	rhs[fUVelocity] += constMassDiffCoeff * flame->StandardDiffusion( fUVelocity, mixViscosity, nodeInfo );
	
// fFirstSpecies to fFirstSpecies + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		eqLoop = speciesEq - fFirstSpecies;
		rhs[speciesEq] += constMassDiffCoeff * flame->SecondDerivSpeciesDiffusion( speciesEq, nodeInfo );
#ifdef DIFFUSIVITYCORRECTION
		if ( flame->UseDiffCorr() ) {
			rhs[speciesEq] -= constMassDiffCoeff * flame->NewDiffCorr( speciesEq, nodeInfo );
		}
#endif
		if ( flame->fThermoDiffusion ) {
			rhs[speciesEq] += constMassDiffCoeff * flame->ThermoDiffusion( eqLoop, kSimilarity, nodeInfo );
		}
		rhs[speciesEq] += productionRate[eqLoop] * oneOverATimesRho;
	}
	
// energy equation
	if ( fTemperature < M ) {
		Double	sumCpY = 0.0;
		Double	corrCoeff = 0.0;
		Double	oneOverARhoCp = 1.0 / ( strainRate * mixDensity[kCurr] * mixHeatCapacity );
		sumCpDdYdx = 0.0;
		sumMH = 0.0;
		rhs[fTemperature] += constThermDiffCoeff * flame->StandardDiffusion( fTemperature, flameNode->mixConductivity, nodeInfo );
#ifdef DIFFUSIVITYCORRECTION
		if ( flame->UseDiffCorr() ) {
			corrCoeff = constMassDiffCoeff * mixDensity[kCurr] * mixDensity[kCurr]
									* diffCorr[kCurr] / mixHeatCapacity;
		}
#endif

		for ( eqLoop = 0; eqLoop < nSpeciesInSystem; ++eqLoop ) {
			sumCpDdYdx += heatCapacity[eqLoop] * diffusivity[eqLoop] 
						* FirstDeriv( YPrev[eqLoop], Y[eqLoop], YNext[eqLoop], hm, h );
#ifdef DIFFUSIVITYCORRECTION
			if ( flame->UseDiffCorr() ) {
				sumCpY += heatCapacity[eqLoop] * Y[eqLoop];
			}
#endif
			sumMH += productionRate[eqLoop] * enthalpy[eqLoop];
		}
#ifdef ENTHALPYFLUX
		rhs[fTemperature] +=  constThermDiffCoeff * sumCpDdYdx * mixDensity[kCurr] * mixDensity[kCurr]
						* FirstDeriv( yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h );
#	ifdef DIFFUSIVITYCORRECTION
		if ( flame->UseDiffCorr() ) {
			rhs[fTemperature] -= corrCoeff * sumCpY
						* FirstDeriv( yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h );
		}
#	endif
#endif
		rhs[fTemperature] -= sumMH * oneOverARhoCp;
		
		if ( flame->fProperties->GetRadiation() ) {
			rhs[fTemperature] += flameNode->radiation[kCurr] * oneOverATimesRho / mixHeatCapacity;
			if ( flame->GetSoot() ) {
				rhs[fTemperature] -= oneOverATimesRho / mixHeatCapacity * flame->GetSoot()->GetSootRadiation( y[fTemperature], flameNode->moments );
			}
		}
	}

	for ( eqLoop = 0; eqLoop < M; ++eqLoop ){
		rhs[eqLoop] *= - hnenn;
	}
}

int CountPremSimPostIter( void *object )
{
	TCountPremFlameSimPtr	flame = ( TCountPremFlameSimPtr )object;
	int			fUVelocity = flame->GetOffsetUVelocity();
	int			fVVelocity = flame->GetOffsetVVelocity();
	int			fTemperature = flame->GetOffsetTemperature();
	int			fFirstSpecies = flame->GetOffsetFirstSpecies();
	int			i;
	TNewtonPtr	bt = flame->GetSolver()->bt;
	NodeInfoPtr nodeInfo = bt->GetNodeInfo();
	TGridPtr 	currGrid = bt->GetGrid()->GetCurrentGrid();
	int			nGridPoints = currGrid->GetNGridPoints();
	int			*bcFlagLeft = currGrid->GetBcFlagLeft();
	int			mixtureSpecificationLeft = flame->GetMixtureSpecificationLeft();
	MatrixPtr	yMat = currGrid->GetY();
	Double		**y = yMat->mat;
	VectorPtr	xVec = currGrid->GetX();
	Double		*x = xVec->vec;
	VectorPtr	yLeftVec = currGrid->GetYLeft();
	VectorPtr	yRightVec = currGrid->GetYRight();
	Double		hLast = bt->GetRight() - x[nGridPoints-1];
	Double		*bcLeft = currGrid->GetBcLeft()->vec;
	Double		hFirst = x[0] - bt->GetLeft();
	T1DSpeciesPtr	species = flame->GetSpecies();
	T1DPropertiesPtr	properties = flame->GetProperties();
	Double		*diffusivity = species->GetDiffusivity()->mat[0];
	Double		*density =  properties->GetDensity()->vec;
	Double		viscosityInf = properties->GetViscosity()->vec[nGridPoints];
	Double		pressure = flame->GetPressure();
	Double		parameter = bt->GetParameter();
	Double		*temp = flame->GetTemperature()->vec;
	Double		**Y = flame->GetMassFracs()->mat;
	int			variables = bt->GetNVariables();
	int			nSpeciesInSystem = species->GetNSpeciesInSystem();
	
	for ( i = 0; i < nGridPoints; ++i ) {
		if ( flame->CheckSolution( y[i][fTemperature], &y[i][fFirstSpecies], nSpeciesInSystem ) ) {
			return 1;
		}
	}

	if ( flame->GetSoot() ) {
		flame->GetSoot()->PostIter( flame );
	}

	CountPremSimUpdateLeftBoundary( flame );
	CountPremSimUpdateRightBoundary( flame );
/*
// left boundary
	for ( int j = 0; j < variables; ++j ) {
		yLeft[j] = y[0][j];
	}
	yLeft[fVVelocity] = 0.0;

//	update properties
	flame->UpdateSolution( yMat, yLeftVec, yRightVec );
	flame->SetFlameNode( kPrev );
	flame->ComputeProperties( flame->fFlameNode, temp[kPrev], Y[kPrev], pressure );
	flame->SetFlameNode( nGridPoints );
	flame->ComputeProperties( flame->fFlameNode, temp[nGridPoints], Y[nGridPoints], pressure );
	
// right boundary 
	yRight[fUVelocity] = 1.0;
	yRight[fVVelocity] = yLast[fVVelocity] + hLast * yRight[fUVelocity];
*/
	flame->UpdateSolution( yMat, yLeftVec, yRightVec );
	
	flame->UpdateThermoProps();

	return 0;
}

#include "TofZ.h"

void TCountPremFlameSim::UpdateDimensions( int len )
{
	T1DFlame::UpdateDimensions( len );
	fSolV->len = len;
	fSolU->len = len;
}

void TCountPremFlameSim::UpdateSolution( MatrixPtr yMat, VectorPtr yLeftVec, VectorPtr yRightVec )
{
	int		nGridPoints = yMat->cols;
	Double	*V = fSolV->vec;
	Double	*U = fSolU->vec;
	Double	**y = yMat->mat;
	Double	*yLeft = yLeftVec->vec;
	Double	*yRight = yRightVec->vec;

	UpdateDimensions( nGridPoints );

	T1DFlame::UpdateSolution( yMat, yLeftVec, yRightVec );
	
	V[kPrev] = yLeft[fVVelocity];
	U[kPrev] = yLeft[fUVelocity];
	for ( int k = 0; k < nGridPoints; ++k ) {
		V[k] = y[k][fVVelocity];
		U[k] = y[k][fUVelocity];
	}
	V[nGridPoints] = yRight[fVVelocity];
	U[nGridPoints] = yRight[fUVelocity];
}

void TCountPremFlameSim::UpdateSolutionOnePoint( Double *y, int gridPoint /*, Double *yMax */ )
{
	int		i;
	int		nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();

	fSolU->vec[gridPoint] = y[fUVelocity];
	fSolV->vec[gridPoint] = y[fVVelocity];
	fSolTemp->vec[gridPoint] = y[fTemperature];
	for ( i = 0; i < nSpeciesInSystem; ++i ) {
		fSolMassFracs->mat[gridPoint][i] = y[i+fFirstSpecies];
//		fSolMassFracs->mat[gridPoint][i] = y[i+fFirstSpecies] * yMax[i];
	}
}

void TCountPremFlameSim::SolutionToSolver( void )
{
	TNewtonPtr	bt = fSolver->bt;
	TGridPtr	grid = bt->GetGrid()->GetFine();
	int		nGridPoints = fSolV->len;
	Double	*V = fSolV->vec;
	Double	*U = fSolU->vec;
	Double	**y = grid->GetY()->mat;

	T1DFlame::SolutionToSolver();
	
	for ( int k = 0; k < nGridPoints; ++k ) {
		y[k][fVVelocity] = V[k];
		y[k][fUVelocity] = U[k];
	}
	
	CountPremSimPostIter( this );
}

void TCountPremFlameSim::SaveSolution( void )
{
	int		k;
	int		len = fSolV->len;
	Double	*v = fSolV->vec;
	Double	*saveV = fSavedV->vec;
	Double	*u = fSolU->vec;
	Double	*saveU = fSavedU->vec;

	T1DFlame::SaveSolution();
	fSavedV->len = fSolV->len;
	fSavedU->len = fSolU->len;

	for ( k = -1; k <= len; ++k ) {
		saveV[k] = v[k];
		saveU[k] = u[k];
	}
}

void TCountPremFlameSim::RestoreSolution( void )
{
	int		k;
	int		len = fSolV->len;
	Double	*v = fSolV->vec;
	Double	*saveV = fSavedV->vec;
	Double	*u = fSolU->vec;
	Double	*saveU = fSavedU->vec;

	UpdateDimensions( len );

	T1DFlame::RestoreSolution();

	for ( k = -1; k <= len; ++k ) {
		v[k] = saveV[k];
		u[k] = saveU[k];
	}
	
	SolutionToSolver();
}

int	TCountPremFlameSim::GetOffsetVVelocity( void )
{
	return fVVelocity; 
}

int	TCountPremFlameSim::GetOffsetUVelocity( void )
{
	return fUVelocity; 
}

int	TCountPremFlameSim::GetOffsetTemperature( void )
{
	return fTemperature; 
}

int	TCountPremFlameSim::GetOffsetFirstSpecies( void ) 
{
	return fFirstSpecies;
}

ConstStringArray TCountPremFlameSim::GetVariableNames( void )
{
	return ( ConstStringArray ) fVariableNames;
}

int TCountPremFlameSim::GetVariablesWithoutSpecies( void )
{
	return fVariablesWithoutSpecies;
}

void TCountPremFlameSim::SetInitialValues( TInputDataPtr inp, StartProfilePtr sp )
{
	int 				i, j, k;
	TBVPSolverPtr		solver = GetSolver();
	TNewtonPtr			bt = solver->bt;
	NodeInfoPtr			nodeInfo = bt->GetNodeInfo();
	TAdaptiveGridPtr	adapGrid = bt->GetGrid();
	TGridPtr			grid = adapGrid->GetFine();
	int					variables = bt->GetNVariables();
	int					nGridPoints = grid->GetNGridPoints();
	int					maxGridPoints = bt->GetMaxGridPoints();
	int					initialGridPoints = bt->GetInitialGridpoints();
	int					nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();
	int					speciesIndex;
	MatrixPtr			yMat = grid->GetY();
	VectorPtr			yLeftVec = grid->GetYLeft();
	VectorPtr			yRightVec = grid->GetYRight();
	Double			 	*yLeft = yLeftVec->vec;
	Double 				*yRight = yRightVec->vec;
	Double				left = inp->fLeft;
	Double				right = inp->fRight;
	Double				*locX = grid->GetX()->vec;
	int					gridPointsIn = sp->gridPoints;	// including boundaries
	Double				*yWork = adapGrid->GetWorkVector()->vec;
	Flag				etaSet = FALSE;
	Flag				oxidizerFound = FALSE;
	Flag				chooseInputGrid = FALSE;
	Double				*xIn = new Double[gridPointsIn];
	if ( !xIn ) FatalError( "memory allocation of TCountPremFlameSim failed" );
	Double				*yIn =  new Double[gridPointsIn];
	if ( !yIn ) FatalError( "memory allocation of TCountPremFlameSim failed" );
	Double				*yInFloat = sp->data;
	Double				**y = grid->GetY()->mat;
	int					variable;
	char				*string = sp->labels;
	SplinePtr			theSpline = NULL;
	Double				leftSlope;
	Double				rightSlope;
	int					oxidizerSide; // this program assumes that oxidizerSide = kRight
	FILE				*fp;
	Double				*temp = GetTemperature()->vec;
	Double				**Y = GetMassFracs()->mat;
	Double				strainRateIn;
	struct _parameter	*param = GetParameter( "strainrate" );
	
	if ( !fStrainRate ) {
	// get strainrate from input file
		if ( param ) {
			strainRateIn = (Double)param->what.quantity.value;
		}
		else { // choose default
			cerr << "#warning: no value for 'strainrate' in inputfile" << NEWL;
			strainRateIn = 10.0;
		}
		
		fStrainRate = NewVector( 1 );
		fStrainRate->vec[0] = strainRateIn;
		fStrainRate->len = 0;
		cerr << "initial strainrate is " << GetStrainRate() << NEWL;
	}
	
//	choose grid from input or equidistant
	if ( gridPointsIn <= maxGridPoints+2 && gridPointsIn > initialGridPoints+2 && ( gridPointsIn % 2 ) != 0 ) {
		grid->AdjustNGridPoints( gridPointsIn-2 );
		solver->UpdateAllDimensions( gridPointsIn-2 );	
		nGridPoints = grid->GetNGridPoints();
		chooseInputGrid = TRUE;
	}
	else {
		nGridPoints = initialGridPoints;
		grid->AdjustNGridPoints( nGridPoints );
		solver->UpdateAllDimensions( nGridPoints );	
		chooseInputGrid = FALSE;
	}
	
// find independent coordinate and oxidizerSide
	for ( i = 0; i < sp->variables; ++i ) {
		if ( strcmp( string, "eta" ) == 0 ) {
			if ( chooseInputGrid ) {
				cerr << "choose inputGrid" << NEWL;
				for ( j = 0; j < gridPointsIn-2; ++j ) {
					locX[j] = yInFloat[i*gridPointsIn + j+1];		// implicit cast from float to Double
				}
				bt->SetLeft( yInFloat[i*gridPointsIn] );
				bt->SetRight( yInFloat[(i+1)*gridPointsIn-1] );
				left = bt->GetLeft();
				right = bt->GetRight();
				etaSet = TRUE;
			}
			else { // choose own Grid, but read x for interpolation
				cerr << "choose own Grid" << NEWL;
				for ( j = 0; j < gridPointsIn; ++j ) {
					xIn[j] = yInFloat[i*gridPointsIn + j];		// implicit cast from float to Double
				}
				bt->SetLeft( left );
				bt->SetRight( right );
				if ( xIn[0] > left ) {
					xIn[0] = left;
				}
				if ( xIn[gridPointsIn-1] < right ) {
					xIn[gridPointsIn-1] = right;
				}
				grid->Make_equi_Grid();
				etaSet = TRUE;
			}
		}
		else if ( strcmp( string, "massfraction-o2" ) == 0 ) {
			if ( yInFloat[i * gridPointsIn] < yInFloat[(i+1) * gridPointsIn - 1] ) {
				oxidizerSide = kRight;
			}
			else {
				oxidizerSide = kLeft;
			}
			oxidizerFound = TRUE;
		}
		string += strlen(string) + 1;
	}

// set default values
	for ( i = 0; i < nGridPoints; ++i ) {
		for ( j = 0; j < variables; ++j ) {
			y[i][j] = yLeft[j] + ( yRight[j] - yLeft[j] ) / ( right - left ) * locX[i];
		}
	}

// error checking
	if ( !etaSet ) {
		cerr << "error: can't find coordinate 'eta'" << NEWL;
		exit(2);
	}
	if ( !oxidizerFound ) {
		cerr << "error: can't find massfraction of oxidizer 'massfraction-o2'" << NEWL;
		exit(2);
	}
	
// reset string
	string = sp->labels;
	
	for ( i = 0; i < sp->variables; ++i ) {
		if ( ( strcmp( string, "df/deta" ) == 0 ) || ( strcmp( string, "u" ) == 0 )
					|| ( strcmp( string, "f'" ) == 0 ) ) {
			variable = fUVelocity;
		}
		else if ( strcmp( string, "f" ) == 0 ) {
			variable = fVVelocity;
		}
		else if ( strncmp( string, "temperature", 11 ) == 0 ) {
			variable = fTemperature;
		}
		else if ( GetSoot() && strncmp( string, "conc-soot", 9 ) == 0 ) {
			string += 9;
//			cerr << "string = " << string << NEWL;
			int	num = atoi( string );
			if ( num < GetSoot()->GetNSootMoments() ) {
				variable = num + GetSoot()->GetOffsetSootMoments();
			}
			else {
				string += strlen(string) + 1;
				continue;
			}
		}
		else if ( strncmp( string, "massfraction-", 13 ) == 0 ){
			string += 13;
			UpperString( string );
			if ( ( speciesIndex = inp->FindSpecies( string ) ) >= 0 ) {
				if ( speciesIndex < nSpeciesInSystem ) {
					variable = fFirstSpecies + speciesIndex;
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
		else {
			string += strlen(string) + 1;
			continue;
		}

		string += strlen(string) + 1;
		if ( chooseInputGrid ) {
			if ( oxidizerSide == kRight ) {
				for ( k = 0; k < gridPointsIn-2; ++k ) {	// transform to current coordinate system
					if ( variable != fVVelocity || yInFloat[i*gridPointsIn] < yInFloat[i*gridPointsIn + 3] ) {
						y[k][variable] = yInFloat[i*gridPointsIn + k+1];	// copy workspace to vector of solution
					}
					else {
						y[k][variable] = -yInFloat[i*gridPointsIn + k+1];	// copy workspace to vector of solution and transform f
					}
				}
			}
			else {
				for ( k = 0; k < gridPointsIn-2; ++k ) { // turn vector and transform to current coordinate system
					if ( variable != fVVelocity || yInFloat[i*gridPointsIn] > yInFloat[i*gridPointsIn + 3] ) {
						y[k][variable] = yInFloat[(i+1)*gridPointsIn - k-2];	// copy workspace to vector of solution
					}
					else {
						y[k][variable] = -yInFloat[(i+1)*gridPointsIn - k-2];	// copy workspace to vector of solution
					}
				}
			}
		}
		else {
			for ( k = 0; k < gridPointsIn; ++k ) {	// store vector in workspace
				yIn[k] = yInFloat[i * gridPointsIn + k];	// implicit cast from float to Double
			}
		
			if ( variable == fVVelocity ) {
				yIn[0] = yIn[1] - ( yIn[2] - yIn[1] ) / ( xIn[2] - xIn[1] ) * ( xIn[1] - xIn[0] );
				yIn[gridPointsIn-1] = yIn[gridPointsIn-2] + ( yIn[gridPointsIn-2] - yIn[gridPointsIn-3] ) 
										/ ( xIn[gridPointsIn-2] - xIn[gridPointsIn-3] ) 
										* ( xIn[gridPointsIn-1] - xIn[gridPointsIn-2] );
			}
					
			leftSlope = ( yIn[1] - yIn[0] ) / ( xIn[1] - xIn[0] );
			rightSlope = ( yIn[gridPointsIn-1] - yIn[gridPointsIn-2] ) / ( xIn[gridPointsIn-1] - xIn[gridPointsIn-2] );
			theSpline = ComputeSimpleSpline( xIn, yIn, gridPointsIn, FALSE, leftSlope, FALSE, rightSlope, NULL, TRUE );
			SplineInterpolate( theSpline, locX, yWork, nGridPoints );
			if ( oxidizerSide == kRight ) {
				for ( k = 0; k < nGridPoints; ++k ) {	// transform to current coordinate system
					if ( variable != fVVelocity || yInFloat[i*gridPointsIn] < yInFloat[i*gridPointsIn + 3] ) {
						y[k][variable] = yWork[k];	// copy workspace to vector of solution
					}
					else {
						y[k][variable] = -yWork[k];	// copy workspace to vector of solution and transform f
					}
				}
			}
			else {
				for ( k = 0; k < nGridPoints; ++k ) { // turn vector and transform to current coordinate system
					if ( variable != fVVelocity || yInFloat[i*gridPointsIn] > yInFloat[i*gridPointsIn + 3] ) {
						y[k][variable] = yWork[nGridPoints-k-1];	// copy workspace to vector of solution
					}
					else {
						y[k][variable] = -yWork[nGridPoints-k-1];	// copy workspace to vector of solution
					}
				}
			}
		}
		//	set bc
		if ( oxidizerSide == kRight ) {
			if ( variable == fTemperature ) {
				if ( yLeft[fTemperature] <= 0.0 ) {
					yLeft[fTemperature] = yInFloat[i*gridPointsIn];
				}
				if ( yRight[fTemperature] <= 0.0 ) {
					yRight[fTemperature] = yInFloat[(i+1)*gridPointsIn-1];
				}
			}
		}
		else {
			if ( variable == fTemperature ) {
				if ( yLeft[fTemperature] <= 0.0 ) {
					yLeft[fTemperature] = yInFloat[(i+1)*gridPointsIn-1];
				}
				if ( yRight[fTemperature] <= 0.0 ) {
					yRight[fTemperature] = yInFloat[i*gridPointsIn];
				}
			}
		}
	}

// set left hand side of the grid to zero
	if (  oxidizerSide == kRight ) {
		for ( k = 0; k < nGridPoints; ++k ) {
			locX[k] -= bt->GetLeft();
		}
		bt->SetRight( bt->GetRight() - bt->GetLeft() );
		bt->SetLeft( 0.0 );
	}
	else {
		for ( k = 0; k < nGridPoints; ++k ) {
			yWork[k] = bt->GetRight() - locX[nGridPoints-k-1]; 	// turn sign and sort for ascending values
		}
		for ( k = 0; k < nGridPoints; ++k ) {
			locX[k] = yWork[k];					// copy workspace to x
		}
		bt->SetRight( bt->GetRight() - bt->GetLeft() );
		bt->SetLeft( 0.0 );
	}
	
// set initial Boundary values and pressure

/*	for ( j = 0; j < variables; ++j ) {*/
/*		yLeft[j] = y[0][j];*/
/*	}*/
	if ( GetPressure() <= 0.0 ) {
		param = GetParameter( "pressure" );
		Double	thePressure;
		if ( param ) {
			thePressure = (Double)param->what.quantity.value;
			if ( strcmp( param->what.quantity.unit, "bar" ) == 0 ) {
				thePressure *= 1.0e5;
			}
			SetPressure( thePressure );
		}
		else { // exit
			cerr << "#error: no value for 'pressure' in inputfile" << NEWL;
			exit(2);
		}
	}
/*	for ( j = 0; j < variables; ++j ) {*/
/*		yLeft[j] = y[0][j];*/
/*	}*/
/*	yLeft[fVVelocity] = 0.0;*/
/*	
//	scale mass fractions in solver
	Double	*yMax = bt->fYMax->vec;
	for ( i = 0; i < nSpeciesInSystem; ++i ) {
		yMax[i] = y[0][i+fFirstSpecies];
		for ( k = 1; k < nGridPoints; ++k ) {
			if ( y[k][i+fFirstSpecies] > yMax[i] ) yMax[i] = y[k][i+fFirstSpecies];
		}
		for ( k = 0; k < nGridPoints; ++k ) {
			y[k][i+fFirstSpecies] /= yMax[i];
		}
		yLeft[i+fFirstSpecies] /= yMax[i];
		yRight[i+fFirstSpecies] /= yMax[i];
	}

*/
//	update properties
	UpdateSolution( yMat, yLeftVec, yRightVec );
	SetFlameNode( kPrev );
	ComputeProperties( fFlameNode, temp[kPrev], Y[kPrev], GetPressure() );
	SetFlameNode( nGridPoints );
	ComputeProperties( fFlameNode, temp[nGridPoints], Y[nGridPoints], GetPressure() );
	
	yRight[fUVelocity] = 1.0;
	yRight[fVVelocity] = y[nGridPoints-1][fVVelocity] + ( bt->GetRight() - locX[nGridPoints-1] ) 
					   * yRight[fUVelocity];

	if ( GetSpecies()->IsConstantLewisNumber() ) {
		CompLewisNumbers( GetSpecies()->GetLewisNumberFile() );
	}
	CountPremSimPostIter( this );

	FreeSpline( theSpline );
	delete yIn;
	delete xIn;

	adapGrid->SetSolutionScaler();
	
	fp = GetOutfile( "initialguess", TFlame::kData );
	bt->PrintSolution( locX, y, GetVariableNames(), fp );
	fclose(fp);
}

void TCountPremFlameSim::FillJacSpeciesDiffusion( int nVariable, Double constCoeff, NodeInfoPtr nodeInfo, Flag sign )
{
// fills the jacobian with     constCoeff * d/dy ( rho * rho * diffusivity * df/dy )

	int		speciesIndex = nVariable - fFirstSpecies;
	Double	diff = fFlameNode->diffusivity[speciesIndex];
	Double	diffNext = fFlameNode->diffusivityNext[speciesIndex];
	Double	diffPrev = fFlameNode->diffusivityPrev[speciesIndex];
	Double	*rho = fFlameNode->mixDensity;
//	Double	yPrev = nodeInfo->yPrev[nVariable];
//	Double	y = nodeInfo->y[nVariable];
//	Double	yNext = nodeInfo->yNext[nVariable];

	if ( sign == kNegative ) {
		constCoeff *= -1.0;
	}
	Double	diffPlusHm = nodeInfo->hm * constCoeff * ( diff * rho[kCurr] * rho[kCurr]
										+ diffNext * rho[kNext] * rho[kNext] );
	Double	diffMinusH = nodeInfo->h * constCoeff * ( diffPrev * rho[kPrev] * rho[kPrev]
					+ diff * rho[kCurr] * rho[kCurr] );

	nodeInfo->a[nVariable][nVariable] -= ( diffPlusHm + diffMinusH );
//	nodeInfo->a[fTemperature][nVariable] -= 2.0 / fFlameNode->temp[kCurr] 
//			* constCoeff * diff * rho[kCurr] * rho[kCurr] 
//			* ( hm * ( yNext - y ) + h * ( yPrev - y ) );
	if ( !nodeInfo->lastPoint ) {
		nodeInfo->b[nVariable][nVariable] += diffPlusHm;
//		nodeInfo->b[fTemperature][nVariable] -= 2.0 / fFlameNode->temp[kNext] 
//				* constCoeff * diffNext * rho[kNext] * rho[kNext]
//				* hm * ( yNext - y );
	}
	if ( !nodeInfo->firstPoint ) {
		nodeInfo->c[nVariable][nVariable] += diffMinusH;
//		nodeInfo->c[fTemperature][nVariable] -= 2.0 / fFlameNode->temp[kPrev] 
//				* constCoeff * diffPrev * rho[kPrev] * rho[kPrev]
//				* h * ( yPrev - y );
	}
}

Double TCountPremFlameSim::SecondDerivSpeciesDiffusion( int nVariable, NodeInfoPtr nodeInfo )
{
	int		speciesIndex = nVariable - fFirstSpecies;
	Double	yPrev = nodeInfo->yPrev[nVariable];
	Double	y = nodeInfo->y[nVariable];
	Double	yNext = nodeInfo->yNext[nVariable];
	Double	*diffusivity = fFlameNode->diffusivity;
	Double	*diffusivityNext = fFlameNode->diffusivityNext;
	Double	*diffusivityPrev = fFlameNode->diffusivityPrev;
	Double	*mixDensity = fFlameNode->mixDensity;
	Double	diffPlus = diffusivity[speciesIndex] * mixDensity[kCurr] * mixDensity[kCurr]
					+ diffusivityNext[speciesIndex] * mixDensity[kNext] * mixDensity[kNext];
	Double	diffMinus = diffusivityPrev[speciesIndex] * mixDensity[kPrev] * mixDensity[kPrev]
					+ diffusivity[speciesIndex] * mixDensity[kCurr] * mixDensity[kCurr];
	Double	hm = nodeInfo->hm;
	Double	h = nodeInfo->h;
	
	return ( diffPlus * hm * ( yNext - y ) + diffMinus * h * ( yPrev - y ) ) 
				/ nodeInfo->hnenn;
}

void TCountPremFlameSim::PrintRHSTemp( TNewtonPtr bt )
{
	TAdaptiveGridPtr	grid = bt->GetGrid();
    TGridPtr			currentGrid = grid->GetCurrentGrid();
	NodeInfoPtr			nodeInfo = bt->GetNodeInfo();
    int         		k;
    int         		N = currentGrid->GetNGridPoints();
	VectorPtr 			physXVec = NewVector( N + 2 );
	Double				*physX = physXVec->vec;
	FILE				*fp = GetOutfile( "tempRHS", TFlame::kData );
	
#if defined (applec) || defined (powerc)
    RotateCursor( 32 * bt->GetNIter() );
#endif
	
	UpdateThermoProps();
		
	EtaToX( bt, physXVec );
		
	fprintf( fp, "*\n%-12s\t%-12s\t%-12s\t%-12s", "eta", "y", "convection", "diffusion" );
#ifdef ENTHALPYFLUX
	fprintf( fp, "\t%-12s", "entflux" );
#	ifdef DIFFUSIVITYCORRECTION
	if ( UseDiffCorr() ) {
		fprintf( fp, "\t%-12s", "diffCorr" );
	}
#	endif
#endif
	fprintf( fp, "\t%-12s", "production" );
	if ( fProperties->GetRadiation() ) {
		fprintf( fp, "\t%-12s", "radiation" );
		if ( GetSoot() ) {
			fprintf( fp, "\t%-12s", "sootRad" );
		}
	}
	fprintf( fp, "\n" );
	
	for ( k = 0; k < N; ++k ){
#if defined (applec) || defined (powerc)
	RotateCursor( 32 );
#endif
		bt->SetNodeInfo( this, k );
		PrintRHSTemp( bt, nodeInfo, physX[k+1], fp );
	}
    fclose( fp );
}

void TCountPremFlameSim::PrintRHSTemp( TNewtonPtr bt, NodeInfoPtr nodeInfo, Double physX, FILE *fp )
{
	int		eqLoop, speciesEq;
    int		M = bt->GetNEquations();
	int		nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();
    Double  h = nodeInfo->h;
    Double  hm = nodeInfo->hm;
	Double	*yPrev = nodeInfo->yPrev;
	Double	*y = nodeInfo->y;
	Double	*yNext = nodeInfo->yNext;
	Double	strainRate = GetStrainRate();
	Double	*enthalpy = fFlameNode->enthalpy;
	Double	mixDensity = *fFlameNode->mixDensity;
	Double	mixDensityInf = fFlameNode->rhoInf;
	Double	mixViscosityInf = fFlameNode->viscosityInf;
	Double	*productionRate = fFlameNode->productionRate;
	Double	*diffusivity = fFlameNode->diffusivity;
	Double	*heatCapacity = fFlameNode->heatCapacity;
	Double	oneOverATimesRho = 1.0 / ( strainRate * mixDensity );
	Double	idealGasCoeff = GetPressure() * *fFlameNode->mixMolarMass / RGAS; // rho = idealGasCoeff / T
	Double	constMassDiffCoeff = 1.0 / ( mixDensityInf * mixViscosityInf );
	Double	mixHeatCapacity = *fFlameNode->mixHeatCapacity;
	Double	constThermDiffCoeff = constMassDiffCoeff / fFlameNode->mixHeatCapacity[kCurr];
	Double	sumCpDdYdx;
	Double	sumMH;
	Double	*diffCorr = fFlameNode->diffCorr;

	fprintf( fp, "%-.6e\t%-.6e", *nodeInfo->x, physX );
#ifdef UPWINDCONVECTION
	fprintf( fp, "\t%-.6e", NonlinearConvectUpwind( y[fVVelocity], yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h ) );
#else
	fprintf( fp, "\t%-.6e", NonlinearConvectCentral( y[fVVelocity], yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h ) );
#endif

// energy equation
	Double	sumCpY = 0.0;
	Double	corrCoeff = 0.0;
	Double	oneOverARhoCp = 1.0 / ( strainRate * mixDensity * mixHeatCapacity );
	sumCpDdYdx = 0.0;
	sumMH = 0.0;

	fprintf( fp, "\t%-.6e", constThermDiffCoeff * StandardDiffusion( fTemperature, fFlameNode->mixConductivity, nodeInfo ) );
#ifdef DIFFUSIVITYCORRECTION
	if ( UseDiffCorr() ) {
		corrCoeff = constMassDiffCoeff * mixDensity * mixDensity
								* diffCorr[kCurr] / mixHeatCapacity;
	}
#endif
	for ( eqLoop = 0; eqLoop < nSpeciesInSystem; ++eqLoop ) {
		speciesEq = fFirstSpecies+eqLoop;
		sumCpDdYdx += heatCapacity[eqLoop] * diffusivity[eqLoop] 
					* FirstDeriv( yPrev[speciesEq], y[speciesEq], yNext[speciesEq], hm, h );
#ifdef DIFFUSIVITYCORRECTION
		if ( UseDiffCorr() ) {
			sumCpY += heatCapacity[eqLoop] * y[speciesEq];
		}
#endif
		sumMH += productionRate[eqLoop] * enthalpy[eqLoop];
	}
#ifdef ENTHALPYFLUX
	fprintf( fp, "\t%-.6e", constThermDiffCoeff * sumCpDdYdx * mixDensity * mixDensity
						* FirstDeriv( yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h ) );
#	ifdef DIFFUSIVITYCORRECTION
	if ( UseDiffCorr() ) {
		fprintf( fp, "\t%-.6e", -corrCoeff * sumCpY
					* FirstDeriv( yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h ) );
	}
#	endif
#endif
	fprintf( fp, "\t%-.6e", -sumMH * oneOverARhoCp );
	if ( fProperties->GetRadiation() ) {
		fprintf( fp, "\t%-.6e", fFlameNode->radiation[kCurr] * oneOverATimesRho / mixHeatCapacity );
		if ( GetSoot() ) {
			fprintf( fp, "\t%-.6e", -oneOverATimesRho / mixHeatCapacity * GetSoot()->GetSootRadiation( y[fTemperature], fFlameNode->moments ) );
		}
	}
	fprintf( fp, "\n" );
}

void TCountPremFlameSim::PrintRHSSpecies( TNewtonPtr bt )
{
	TAdaptiveGridPtr	grid = bt->GetGrid();
    TGridPtr			currentGrid = grid->GetCurrentGrid();
	NodeInfoPtr			nodeInfo = bt->GetNodeInfo();
    int         		i, k;
    int         		N = currentGrid->GetNGridPoints();
	VectorPtr 			physXVec = NewVector( N + 2 );
	Double				*physX = physXVec->vec;
	char				**names = GetSpecies()->GetNames();
	int					nSpecIn = fSpecies->GetNSpeciesInSystem();
	int					num = 7;
	int					start = 0, end;
	int					fileCount = 1;
	FILE				*fp;
	
	do {
		end = ( int )MIN( start + 255 / num, nSpecIn );
		sprintf( GetOutFileBuff(), "%sspeciesRHS%d.dout", GetOutputPath(), fileCount++ );
		if ( !( fp = fopen( GetOutFileBuff(), "w") ) ) { 
			cerr << "#warning: unable to open file " << GetOutFileBuff() << NEWL;
			exit(2);
		}
#if defined (applec) || defined (powerc)
  	  RotateCursor( 32 * bt->GetNIter() );
#endif
	
		UpdateThermoProps();
		
		EtaToX( bt, physXVec );
			
		fprintf( fp, "*\n%-12s\t%-12s", "eta", "y" );
		for ( i = start; i < end; ++i ) {
			fprintf( fp, "\tConv_%-7s\tDiff_%-7s", names[i], names[i] );
#ifdef DIFFUSIVITYCORRECTION
			if ( UseDiffCorr() ) {
				fprintf( fp, "\tDiCo_%-7s", names[i] );
			}
#endif
			if ( fThermoDiffusion ) {
				fprintf( fp, "\tThDi_%-7s", names[i] );
			}
			fprintf( fp, "\tProd_%-7s\tCons_%-7s\tSource_%-7s", names[i], names[i], names[i] );
		}
		fprintf( fp, "\n" );
			
		for ( k = 0; k < N; ++k ){
#if defined (applec) || defined (powerc)
		RotateCursor( 32 );
#endif
			bt->SetNodeInfo( this, k );
			PrintRHSSpecies( bt, nodeInfo, physX[k+1], fp );
		}
		fclose( fp );
		start = end;
	} while ( end < nSpecIn );
}

void TCountPremFlameSim::PrintRHSSpecies( TNewtonPtr bt, NodeInfoPtr nodeInfo, Double physX, FILE *fp )
{
	int		i, j;
	int		speciesEq;
    int		M = bt->GetNEquations();
	int		nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();
    Double  h = nodeInfo->h;
    Double  hm = nodeInfo->hm;
	Double	*yPrev = nodeInfo->yPrev;
	Double	*y = nodeInfo->y;
	Double	*yNext = nodeInfo->yNext;
	Double	strainRate = GetStrainRate();
	Double	mixDensityInf = fFlameNode->rhoInf;
	Double	mixViscosityInf = fFlameNode->viscosityInf;
	Double	idealGasCoeff = GetPressure() * *fFlameNode->mixMolarMass / RGAS; // rho = idealGasCoeff / T
	Double	ROverAPM = 1.0 / ( strainRate * idealGasCoeff );
	Double	constMassDiffCoeff = 1.0 / ( mixDensityInf * mixViscosityInf );
	Double	diffTerm;
	Double	*reactionRate = fFlameNode->reactionRate;
	Double	productionRate;
	Double	source;
	Double	sink;
	int		*nOfUsedReactions = fSpecies->GetNOfUsedReactions()->vec;
	VectorPtr	*nu = fSpecies->GetNu();
	IntVectorPtr	*usedReactions = fSpecies->GetUsedReactions();
	Double	*molarMass =  fSpecies->GetMolarMass()->vec;

	fprintf( fp, "%-.6e\t%-.6e", *nodeInfo->x, physX );

	for ( i = 0; i < nSpeciesInSystem; ++i ) {
		speciesEq = fFirstSpecies+i;

#ifdef UPWINDCONVECTION
		fprintf( fp, "\t%-.6e", NonlinearConvectUpwind( y[fVVelocity], yPrev[speciesEq], y[speciesEq], yNext[speciesEq], hm, h ) );
#else
		fprintf( fp, "\t%-.6e", NonlinearConvectCentral( y[fVVelocity], yPrev[speciesEq], y[speciesEq], yNext[speciesEq], hm, h ) );
#endif

		diffTerm = constMassDiffCoeff * SecondDerivSpeciesDiffusion( speciesEq, nodeInfo );
		fprintf( fp, "\t%-.6e", diffTerm );
#ifdef DIFFUSIVITYCORRECTION
		if ( UseDiffCorr() ) {
			diffTerm = -constMassDiffCoeff * NewDiffCorr( speciesEq, nodeInfo );
			fprintf( fp, "\t%-.6e", diffTerm );
		}
#endif
		if ( fThermoDiffusion ) {
			diffTerm = constMassDiffCoeff * ThermoDiffusion( i, kSimilarity, nodeInfo );
			fprintf( fp, "\t%-.6e", diffTerm );
		}

		sink = source = 0.0;
		for ( j = 0; j < nOfUsedReactions[i]; ++j ) {
			productionRate = nu[i]->vec[j] * reactionRate[usedReactions[i]->vec[j]];
			if ( productionRate > 0.0 ) {
				sink += productionRate;
			}
			else {
				source += productionRate;
			}
		}
		sink *= - molarMass[i] * y[fTemperature] * ROverAPM;
		source *= - molarMass[i] * y[fTemperature] * ROverAPM;
		fprintf( fp, "\t%-.6e\t%-.6e\t%-.6e", source, -sink, source+sink );
	}
	fprintf( fp, "\n" );
}

void CountPremSimOutput( void *object, FILE *fp, char* tail )
{
	TCountPremFlameSimPtr	flame = ( TCountPremFlameSimPtr )object;
	TNewtonPtr		bt = flame->GetSolver()->bt;
	T1DPropertiesPtr	props = flame->GetProperties();
	TSpeciesPtr		species = flame->GetSpecies();
	NodeInfoPtr		nodeInfo = bt->GetNodeInfo();
	Double			*rho = props->GetDensity()->vec;
	TGridPtr		currentGrid = bt->GetGrid()->GetCurrentGrid();
	Double			*x = currentGrid->GetX()->vec;
	Double			*mixMolarMass = props->GetMolarMass()->vec;
	Double			*molarMass = species->GetMolarMass()->vec;
	Double			**massFracs = flame->GetMassFracs()->mat;
	Double			*temp = flame->GetTemperature()->vec;
	Double			*V = flame->GetV()->vec;
	Double			*U = flame->GetU()->vec;
	Double			**y = currentGrid->GetY()->mat;
	Double			*yLeft = currentGrid->GetYLeft()->vec,
					*yRight = currentGrid->GetYRight()->vec;
	int				i, k;
	int				gridPoints = currentGrid->GetNGridPoints();
	int				nOfSpecies = flame->GetSpecies()->GetNOfSpecies();
	int				nOfVariables = bt->GetNVariables();
	int				nOfEquations = bt->GetNEquations();
	int				firstSpecies = flame->GetOffsetFirstSpecies();
	int				tempOffset = flame->GetOffsetTemperature();
	int				fUVelocity = flame->GetOffsetUVelocity();
	int				fVVelocity = flame->GetOffsetVVelocity();
	time_t			theDate;
	char			buffer[80];
	ConstStringArray	varNames = flame->GetVariableNames();
	char			**names = species->GetNames();
	VectorPtr 		physXVec = NewVector( gridPoints + 2 );
	Double			*physX = physXVec->vec;
	Flag			fpOpen = FALSE;
	
	if ( !fp ) {
		fpOpen = TRUE;
		fp = flame->GetOutputFile( NULL, tail, TFlame::kNone );
	}

	flame->EtaToX( bt, physXVec );
	
// write header
	fprintf( fp, "header\n\n" );

	fprintf( fp, "title = \"planar counterflow premixed flame\"\n" );
	fprintf( fp, "author = \"%s\"\n", flame->GetAuthor() );
	time( &theDate );
	strcpy( buffer, ctime( &theDate ) );
	if ( buffer[strlen(buffer)-1] == '\n' )
		buffer[strlen(buffer)-1] = '\0';
	fprintf( fp, "date = \"%s\"\n\n", buffer );
	for ( i = 0; i < flame->GetNFuels(); ++i ) {
		fprintf( fp, "\tfuel = \"%s\"\n", varNames[firstSpecies+flame->GetFuelIndex( i )] );
	}
	fprintf( fp, "pressure = %g [bar]\n", flame->GetPressure() / 1.0e5 );
	fprintf( fp, "strainrate = %g [1/s]\n", flame->GetStrainRate() );
	fprintf( fp, "fuel-air-equivalence-ratio = %g\n", flame->GetPhi() );

	if ( species->IsConstantLewisNumber() ) {
		fprintf( fp, "ConstantLewisNumbers = \"True\"\n" );
	}
	
	fprintf( fp, "burningVelocity = %g [cm/s]\n", flame->GetBurningVelocity() * 100.0 );
	fprintf( fp, "Tmax = %g [K]\n", temp[LocationOfMax( gridPoints+2, &temp[kPrev] ) - 1] );
	
	Double	EIFuel = 0.0;
	for ( i = 0; i < flame->GetNFuels(); ++i ) {
		EIFuel += flame->ComputeEmissionIndex( flame->GetFuelIndex( i ), &physX[1] );
	}

	int	indNO = species->FindSpecies( "NO" );
	if ( indNO > -1 ) {
		fprintf( fp, "EmissionIndexNO = %g [g/kg]\n"
				, -1000.0 * flame->ComputeEmissionIndex( indNO, &physX[1] )
				/ EIFuel );
	}
	
	int	indNO2 = species->FindSpecies( "NO2" );
	if ( indNO2 > -1 ) {
		fprintf( fp, "EmissionIndexNO2 = %g [g/kg]\n"
				, -1000.0 * flame->ComputeEmissionIndex( indNO2, &physX[1] ) 
				/ EIFuel );
	}
	
//  write bc
	Double	locMoleMass;
	fprintf( fp, "\nunburnt\n" );
	fprintf( fp, "begin\n" );
	fprintf( fp, "\tTemperature = %g [K]\n", temp[gridPoints] );
	for ( i = 0; i < nOfSpecies; ++i ) { // write Y_i
		if ( fabs( massFracs[gridPoints][i] ) > 1.0e-5 ) {
			fprintf( fp, "\tMassfraction-%s = %g\n", names[i], massFracs[gridPoints][i] );
		}
	}
	for ( i = 0; i < nOfSpecies; ++i ) { // write X_i
		locMoleMass = massFracs[gridPoints][i] * mixMolarMass[gridPoints] / molarMass[i];
		if ( fabs( locMoleMass ) > 1.0e-5 ) {
			fprintf( fp, "\tMolefraction-%s = %g\n", names[i], locMoleMass );
		}
	}
	fprintf( fp, "end\n\n" );

	fprintf( fp, "numOfSpecies = %d\n", nOfSpecies );
	fprintf( fp, "gridPoints = %d\n\n", gridPoints+2 );

	fprintf( fp, "body\n" );

// write independent coordinate
	fprintf( fp, "eta\n" );
	fprintf( fp, "\t%-.6e", bt->GetLeft() );
	for ( k = 0; k < gridPoints; ++k ) {
		fprintf( fp, "\t%-.6e", x[k] );
		if ( (k+2) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	fprintf( fp, "\t%-.6e\n", bt->GetRight() );
	
// write physical coordinate
	
	fprintf( fp, "y [m]\n" );
	for ( k = 0; k < gridPoints+2; ++k ) {
		fprintf( fp, "\t%-.6e", physX[k] );
		if ( (k+1) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	if ( k % 5 ) {
		fprintf( fp, "\n" );
	}
		
// write solution
	// write f-Velocity, f'-Velocity, mixture fraction and temperature
	flame->PrintFlameletVector( gridPoints+2, &V[kPrev], "f", fp );
	flame->PrintFlameletVector( gridPoints+2, &U[kPrev], "df/deta", fp );
	flame->PrintFlameletVector( gridPoints+2, &temp[kPrev], "temperature [K]", fp );

	// write massfractions of species
	for ( i = 0; i < nOfSpecies; ++i ) {
		fprintf( fp, "massfraction-%s\n", names[i] );
		for ( k = 0; k < gridPoints+2; ++k ) {
			fprintf( fp, "\t%-.6e", massFracs[k-1][i] );
			if ( (k+1) % 5 == 0 ) {
				fprintf( fp, "\n" );
			}
		}
		if ( k % 5 ) {
			fprintf( fp, "\n" );
		}
	}

	if ( flame->fPrintMolarFractions ) {
		Double	locMolarMass;
		for ( i = 0; i < nOfSpecies; ++i ) {
//			varIndex = i + firstSpecies;
			locMolarMass = molarMass[i];
			fprintf( fp, "molefraction-%s\n", names[i] );
			fprintf( fp, "\t%-.6e", massFracs[kPrev][i] * mixMolarMass[-1] / locMolarMass );
			for ( k = 0; k < gridPoints; ++k ) {
				fprintf( fp, "\t%-.6e", massFracs[k][i] * mixMolarMass[k] / locMolarMass );
				if ( (k+2) % 5 == 0 ) {
					fprintf( fp, "\n" );
				}
			}
			fprintf( fp, "\t%-.6e\n", massFracs[gridPoints][i] * mixMolarMass[gridPoints] / locMolarMass );
		}
	}

	if ( flame->fSoot ) {
		flame->GetSoot()->PrintFlameletFile( gridPoints, flame, fp );
	}
	
//	write V
	Double	coeff = -sqrt( flame->GetStrainRate() * flame->fFlameNode->rhoInf * flame->fFlameNode->viscosityInf );
	fprintf( fp, "V [kg/(s*m^2)]\n" );
	fprintf( fp, "\t%-.6e", coeff * V[kPrev] );
	for ( k = 0; k < gridPoints; ++k ) {
		fprintf( fp, "\t%-.6e", coeff * V[k] );
		if ( (k+2) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	fprintf( fp, "\t%-.6e\n", coeff * V[gridPoints] );
	
		
/*//	write local strainrate
	Double			absStrainRate;
	fprintf( fp, "localStrainRate\n" );
	absStrainRate = flame->ComputeAbsStrainRate( -1, yLeft[fVVelocity], y[0][fVVelocity], 
										 yLeft[fUVelocity], y[0][fUVelocity], x[0] - bt->GetLeft() );
	fprintf( fp, "\t%-.6e", absStrainRate );
	for ( k = 0; k < gridPoints-1; ++k ) {
		absStrainRate = flame->ComputeAbsStrainRate( k, y[k][fVVelocity], y[k+1][fVVelocity], 
											 y[k][fUVelocity], y[k+1][fUVelocity], x[k+1] - x[k] );
		fprintf( fp, "\t%-.6e", absStrainRate );
		if ( (k+2) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	absStrainRate = flame->ComputeAbsStrainRate( gridPoints-1, y[gridPoints-1][fVVelocity], yRight[fVVelocity], 
										 y[gridPoints-1][fUVelocity], yRight[fUVelocity], bt->GetRight() - x[gridPoints-1] );
	fprintf( fp, "\t%-.6e", absStrainRate );
	if ( (k+2) % 5 ) {
		fprintf( fp, "\n" );
	}
	fprintf( fp, "\t%-.6e\n", absStrainRate );
*/
	
	fprintf( fp, "trailer\n" );
	if ( nOfEquations < nOfVariables) {
		fprintf( fp, "number of converged equations is %d\n", nOfEquations );
	}
	
	DisposeVector( physXVec );
	
	if ( fpOpen ) {
		fclose( fp );
	}
}

Double TCountPremFlameSim::ComputeAbsStrainRate( int k1, Double f1, Double f2, Double fp1, Double fp2, Double h )
{
	Double	*rho = GetProperties()->GetDensity()->vec;
	Double	rho1 = rho[k1];
	Double	rho2 = rho[k1+1];
	Double	strainRate = GetStrainRate();
	Double	fpm = 0.5 * ( fp1 + fp2 );
	
	Double fact = 0.5 * ( rho1 + rho2 ) 
					* ( f2/rho2 - f1/rho1 ) / h;
	return strainRate * sqrt( fpm * fpm + fact * fact ); 
}

void TCountPremFlameSim::FillJacDiffCorr( int nVariable, Double constCoeff, NodeInfoPtr nodeInfo, Flag sign )
{
// fills the jacobian with     constCoeff * sum_j ( d/dy(rho^2 Y_k D_j dY_j/dy) )
// this function implies that dYi/dx = 0 at both boundaries

	if ( sign == kNegative ) {
		constCoeff *= -1.0;
	}

	int		i, lVar;
	int		nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();
	Double	hm = nodeInfo->hm;
	Double	h = nodeInfo->h;
	Double	**a = nodeInfo->a;
	Double	**b = nodeInfo->b;
	Double	**c = nodeInfo->c;
	Double	*Y = fFlameNode->Y[kCurr];
	Double	*YPrev = fFlameNode->Y[kPrev];
	Double	*YNext = fFlameNode->Y[kNext];
//	Double	*temp = fFlameNode->temp;
	Double	*density = fFlameNode->mixDensity;
	Double	*diffusivity = fFlameNode->diffusivity;
	Double	*diffusivityPrev = fFlameNode->diffusivityPrev;
	Double	*diffusivityNext = fFlameNode->diffusivityNext;
	lVar = nVariable-fFirstSpecies;
	Double	coeffCurr = constCoeff * density[kCurr] * density[kCurr] * Y[lVar];
	Double	coeffPrev = constCoeff * density[kPrev] * density[kPrev] * YPrev[lVar];
	Double	coeffNext = constCoeff * density[kNext] * density[kNext] * YNext[lVar];
	Double	diffPlusHm, diffMinusH;
	Double	sumYCurr = 0.0;
	Double	sumYNext = 0.0;
	Double	sumYPrev = 0.0;
	
	for ( i = 0; i < nSpeciesInSystem; ++i ) {
		sumYCurr += Y[i];
		sumYPrev += YPrev[i];
		sumYNext += YNext[i];
	}
	coeffCurr /= sumYCurr;
	coeffPrev /= sumYPrev;
	coeffNext /= sumYNext;
	
	for ( i = 0; i < nSpeciesInSystem; ++i ) {
		lVar = fFirstSpecies + i;
		// d/dY_l
		diffPlusHm = hm * ( coeffCurr * diffusivity[i] + coeffNext * diffusivityNext[i] );
		diffMinusH = h * ( coeffCurr * diffusivity[i] + coeffPrev * diffusivityPrev[i] );
		if ( !nodeInfo->lastPoint ) {
			a[lVar][nVariable] -= diffPlusHm;
			b[lVar][nVariable] += diffPlusHm;
		}
		if ( !nodeInfo->firstPoint ) {
			a[lVar][nVariable] -= diffMinusH;
			c[lVar][nVariable] += diffMinusH;
		}
	}
}

Double TCountPremFlameSim::NewDiffCorr( int nVariable, NodeInfoPtr nodeInfo )
{
// returns     sum_j ( d/dy(rho^2 Y_k D_j dY_j/dy) )

	int		i;
	int		nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();
	Double	hm = nodeInfo->hm;
	Double	h = nodeInfo->h;
	Double	*Y = fFlameNode->Y[kCurr];
	Double	*YPrev = fFlameNode->Y[kPrev];
	Double	*YNext = fFlameNode->Y[kNext];
	Double	*density = fFlameNode->mixDensity;
	Double	*diffusivity = fFlameNode->diffusivity;
	Double	*diffusivityPrev = fFlameNode->diffusivityPrev;
	Double	*diffusivityNext = fFlameNode->diffusivityNext;
	i = nVariable-fFirstSpecies;
	Double	coeffCurr = density[kCurr] * density[kCurr] * Y[i];
	Double	coeffPrev = density[kPrev] * density[kPrev] * YPrev[i];
	Double	coeffNext = density[kNext] * density[kNext] * YNext[i];
	Double	diffPlus, diffMinus;
	Double	value = 0.0;
	Double	sumYCurr = 0.0;
	Double	sumYNext = 0.0;
	Double	sumYPrev = 0.0;
	
	for ( i = 0; i < nSpeciesInSystem; ++i ) {
		sumYCurr += Y[i];
		sumYPrev += YPrev[i];
		sumYNext += YNext[i];
	}
	coeffCurr /= sumYCurr;
	coeffPrev /= sumYPrev;
	coeffNext /= sumYNext;
	
	for ( i = 0; i < nSpeciesInSystem; ++i ) {
		diffPlus = coeffCurr * diffusivity[i] + coeffNext * diffusivityNext[i];
		diffMinus = coeffCurr * diffusivity[i] + coeffPrev * diffusivityPrev[i];
		value += ( diffPlus * hm * ( YNext[i] - Y[i] ) 
					+ diffMinus * h * ( YPrev[i] - Y[i] ) );
	}

	return value / nodeInfo->hnenn;
}

void SetCountPremSimNodeInfo( int k, void *object )
{
	TCountPremFlameSimPtr	flame = ( TCountPremFlameSimPtr )object;
	
	flame->SetFlameNode( k );
}

void CountPremSimPostConv( void *object )
{
	TCountPremFlameSimPtr	flame = ( TCountPremFlameSimPtr )object;
	TNewtonPtr 				bt = flame->GetSolver()->bt;
	TAdaptiveGridPtr		grid = bt->GetGrid();
	TGridPtr				fine = grid->GetFine();
	int						fFirstSpecies = flame->GetOffsetFirstSpecies();
	int						nSpeciesInSystem = flame->fSpecies->GetNSpeciesInSystem();
	int						isConverged = bt->GetConvergeNewton();
	Double					coeffSpecies = 1.0 / flame->GetStrainRate();
	Double					coeffTemp = 1.0 / flame->GetStrainRate();
	Flag					leaveContinPrem;
	Double					*bcRightFine = fine->GetBcRight()->vec;
	
	
	// check length of the computational domain
	if ( isConverged && flame->CheckComputationalDomain() ) {
		flame->fSolver->ReInit();
		bt->WriteOutput( object, NULL, "" );
		return;
	}	
	
	flame->PostConvergence( object );
	if ( bt->GetLeaveContin() ) {
		if ( isConverged ) {
//				flame->GetSpecies()->PrintProdRateTerms( "CH3", flame );
				bt->WriteOutput( object, NULL, "" );
				if ( flame->fReactionFluxes ) {
					flame->ReactionFluxes( kSimilarity );
				}
				if ( flame->fSensAnal ) {
					flame->SensitivityAnalysis( coeffSpecies, coeffTemp, kSimilarity );
				}
				if ( flame->fPrintRHSSpecies ) {
					flame->PrintRHSSpecies( flame->GetSolver()->bt );
				}
				if ( flame->fPrintRHSTemp ) {
					flame->PrintRHSTemp( flame->GetSolver()->bt );
				}
		}
		leaveContinPrem = flame->PostConvTPremixed( isConverged );
		if ( !leaveContinPrem ) {
			flame->SetMassFracsOfPhi( flame, flame->GetPhi()
				, &fine->GetBcRight()->vec[fFirstSpecies], nSpeciesInSystem
				, flame->fSpecies->GetMolarMass()->vec, flame->fSpecies->GetNames() );
			flame->fSolver->ReInit();
		}
	}
	CountPremSimPostIter( flame );
}

ConstStringArray GetCountPremSimVarNames( void *object )
{
	TCountPremFlameSimPtr	flame = ( TCountPremFlameSimPtr )object;
	
	return flame->GetVariableNames();
}

FILE *TCountPremFlameSim::GetOutputFile( char *head, char *tail, FileType type )
{
	int				fuelIndex = GetFuelIndex();
	char			*name = new char[64];
	FILE			*fp;
	char			**speciesNames = fSpecies->GetNames();
	int				tOxidizer = ( int ) fSolTemp->vec[fSolTemp->len];
	Double			phi = GetPhi();
		
	sprintf( name, "%s%s%.8s_p%.2dphi%.1d-%.4da%.5dt%.4d%s"
					, ( head ) ? head : "", ( head ) ? "_" : ""
					, speciesNames[fuelIndex]
					, ( int ) floor( GetPressure() * 1.0e-5 + 0.5 )	// in [bar]
					, ( int ) floor( phi ) , ( int ) ( ( phi - floor( phi ) ) * 1.e4 )
					, ( int ) floor( GetStrainRate() + 0.5 )			// in [1/s]
					, ( int )( tOxidizer ) 						// in [K]
					, ( tail ) ? tail : "" );

	fp = GetOutfile( name, type );
	delete name;

	return fp;
}

void TCountPremFlameSim::EtaToX( TNewtonPtr bt, VectorPtr xPhysVec )
{
	int			k;
	int			gridPoints = bt->GetCurrentGridPoints();
	NodeInfoPtr nodeInfo = bt->GetNodeInfo();
	TGridPtr	grid = bt->GetGrid()->GetCurrentGrid();
#undef OLD

#ifdef OLD
	Double		factor = 2.0 * sqrt( fFlameNode->rhoInf * fFlameNode->viscosityInf / GetStrainRate() );
#else
	Double		factor = 0.5 * sqrt( fFlameNode->rhoInf * fFlameNode->viscosityInf / GetStrainRate() );
#endif
	Double		left = bt->GetLeft();
	Double		right = bt->GetRight();
	Double		*rho = GetProperties()->GetDensity()->vec;
	Double		*x = grid->GetX()->vec;
	Double		**y = grid->GetY()->mat;
	Double		*xPhys = xPhysVec->vec;
	
	if ( xPhysVec->len < gridPoints + 2 ) {
		cerr << "#warning: Vector xPhys too short, values for physical grid are not computed" << NEWL;
		return;
	}
	
	xPhys[0] = 0.0;
#ifdef OLD
	xPhys[1] = xPhys[0] + factor * ( x[0] - left ) / ( rho[-1] + rho[0] );
#else
	xPhys[1] = xPhys[0] + factor * ( x[0] - left ) * ( 1.0 / rho[-1] + 1.0 / rho[0] );
#endif
	for ( k = 1; k < gridPoints; ++k ) {
#ifdef OLD
		xPhys[k+1] = xPhys[k] + factor * ( x[k] - x[k-1] ) / ( rho[k-1] + rho[k] );
#else
		xPhys[k+1] = xPhys[k] + factor * ( x[k] - x[k-1] ) * ( 1.0 / rho[k-1] + 1.0 / rho[k] );
#endif
	}
#ifdef OLD
	xPhys[gridPoints+1] = xPhys[gridPoints] + factor * ( right - x[gridPoints-1] ) / ( rho[gridPoints-1] + rho[gridPoints] );
#else
	xPhys[gridPoints+1] = xPhys[gridPoints] + factor * ( right - x[gridPoints-1] ) * ( 1.0 / rho[gridPoints-1] + 1.0 / rho[gridPoints] );
#endif
	
}


Double EntFluxFuncPremSim( int /*equation*/, NodeInfoPtr nodeInfo, void *object, Flag /*theFlag*/  )
{
	Double	entFlux = 0.0;
#ifdef ENTHALPYFLUX
	TCountPremFlameSimPtr	flame = ( TCountPremFlameSimPtr )object;
	TFlameNodePtr			flameNode = flame->fFlameNode;
	int 	fFirstSpecies = flame->GetOffsetFirstSpecies();
	int 	fTemperature = flame->GetOffsetTemperature();
	int		eqLoop, speciesEq;
	int		nSpeciesInSystem = flame->GetSpecies()->GetNSpeciesInSystem();
	Double	*yPrev = nodeInfo->yPrev;
	Double	*y = nodeInfo->y;
	Double	*yNext = nodeInfo->yNext;
	Double	*mixDensity = flameNode->mixDensity;
	Double	mixDensityInf = flameNode->rhoInf;
	Double	mixViscosityInf = flameNode->viscosityInf;
	Double	*diffusivity = flameNode->diffusivity;
	Double	*heatCapacity = flameNode->heatCapacity;
	Double	constThermDiffCoeff = 1.0 / ( *flameNode->mixHeatCapacity * mixDensityInf * mixViscosityInf );
	Double	*diffCorr = flameNode->diffCorr;
	Double	h = nodeInfo->h;
	Double	hm = nodeInfo->hm;
	Double	sumCpY = 0.0;
	Double	entCoeff = 0.0;
	Double	sumCpDdYdx = 0.0;
			
	entCoeff = constThermDiffCoeff * mixDensity[kCurr] * mixDensity[kCurr]
					* FirstDeriv( yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h );

	for ( eqLoop = 0; eqLoop < nSpeciesInSystem; ++eqLoop ) {
		speciesEq = eqLoop + fFirstSpecies;
		sumCpDdYdx += heatCapacity[eqLoop] * diffusivity[eqLoop] 
					* FirstDeriv( yPrev[speciesEq], y[speciesEq], y[speciesEq], hm, h );
#	ifdef DIFFUSIVITYCORRECTION
		if ( flame->UseDiffCorr() ) {
			sumCpY += heatCapacity[eqLoop] * y[speciesEq];
		}
#	endif
	}
	entFlux = sumCpDdYdx;
#	ifdef DIFFUSIVITYCORRECTION
	if ( flame->UseDiffCorr() ) {
		entFlux -= diffCorr[kCurr] * sumCpY;
	}
#	endif
	entFlux *= entCoeff;
#endif
	
	return entFlux;
}

int	TCountPremFlameSim::GetOffsetMixFrac( void )
{
	cerr << "#error: class has no member fMixtureFraction" << NEWL;
	exit( 2 );
	return 0; 
}

void CountPremSimUpdateLeftBoundary( void  *object )
{
	TCountPremFlameSimPtr	flame = ( TCountPremFlameSimPtr )object;
	int				fVVelocity = flame->GetOffsetVVelocity();
	TNewtonPtr		bt = flame->GetSolver()->bt;
	TGridPtr 		currGrid = bt->GetGrid()->GetCurrentGrid();
	MatrixPtr		yMat = currGrid->GetY();
	Double			**y = yMat->mat;
	Double			*yLeft = currGrid->GetYLeft()->vec;
	Double			*yMax = bt->GetNodeInfo()->yMax;
	Double			pressure = flame->GetPressure();
	int				variables = bt->GetNVariables();

	for ( int j = 0; j < variables; ++j ) {
		yLeft[j] = y[0][j];
	}
	yLeft[fVVelocity] = 0.0;
	
	flame->UpdateSolutionOnePoint( yLeft, kPrev );
	flame->SetFlameNode( kPrev );
	flame->ComputeProperties( flame->fFlameNode, flame->fFlameNode->temp[kCurr]
							, flame->fFlameNode->Y[kCurr], pressure );
}

void CountPremSimUpdateRightBoundary( void *object )
{
	TCountPremFlameSimPtr	flame = ( TCountPremFlameSimPtr )object;
	int				fUVelocity = flame->GetOffsetUVelocity();
	int				fVVelocity = flame->GetOffsetVVelocity();
	int				fSpecies = flame->GetOffsetFirstSpecies();
	TNewtonPtr		bt = flame->GetSolver()->bt;
	TGridPtr 		currGrid = bt->GetGrid()->GetCurrentGrid();
	int				nGridPoints = currGrid->GetNGridPoints();
	MatrixPtr		yMat = currGrid->GetY();
	Double			**y = yMat->mat;
	Double			*bcRight = currGrid->GetBcRight()->vec;
	Double			*yRight = currGrid->GetYRight()->vec;
	Double			*yLast = y[nGridPoints-1];
	Double			*yMax = bt->GetNodeInfo()->yMax;
	Double			*x = currGrid->GetX()->vec;
	Double			hLast = bt->GetRight() - x[nGridPoints-1];
	Double			pressure = flame->GetPressure();
	int				variables = bt->GetNVariables();
	int				nSpecies = flame->fSpecies->GetNOfSpecies();

	yRight[fUVelocity] = 1.0;
	yRight[fVVelocity] = yLast[fVVelocity] + hLast * yRight[fUVelocity];
	for ( int i = 0; i < nSpecies; ++i ) {
		yRight[fSpecies+i] = bcRight[fSpecies+i];
	}

	flame->UpdateSolutionOnePoint( yRight, nGridPoints );
	flame->SetFlameNode( nGridPoints );
	flame->ComputeProperties( flame->fFlameNode, flame->fFlameNode->temp[kCurr]
							, flame->fFlameNode->Y[kCurr], pressure );
}

void TCountPremFlameSim::SetInitialBC( TGridPtr grid, TInputDataPtr inp )
{
	int					i;
	Double				mixMolarMass;
	int					nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();
	SpeciesPtr			species = inp->GetSpecies();
	BoundaryInputPtr	left = inp->leftBoundary;
	int					inpVOffset = inp->fVVelocityOffset;
	int					inpUOffset = inp->fUVelocityOffset;
	int					inpTOffset = inp->fTemperatureOffset;
	int					*speciesIndexRight = NULL;
	int					leftSpecifiedSpecies = left->fSpecifiedSpeciesBCs;
	int					*bcFlagRight = grid->GetBcFlagRight();
	Double				*yright = grid->GetYRight()->vec;
	Double				*bcRight = grid->GetBcRight()->vec;
	Double				*yleft = grid->GetYLeft()->vec;
	Double				*bcleft = grid->GetBcLeft()->vec;
	Double				*phi = GetPhiVector()->vec;
	
	// map boundary conditions for the unburnt from left to right

	//	allocate memory for speciesIndex
	speciesIndexRight = new int[left->fSpecifiedSpeciesBCs];
	if ( !speciesIndexRight ) FatalError( "memory allocation of TCountPremFlameSim failed" );

	//	set speciesIndex
	for ( i = 0; i < leftSpecifiedSpecies; ++i ) {
		speciesIndexRight[i] = inp->FindSpecies( left->speciesName[i] );
	}
	
	// set fMixtureSpecification
	SetMixtureSpecificationRight( left->fMixtureSpecification );
	
	// set BCFlags
	bcFlagRight[fVVelocity] = left->fBcFlag[inpVOffset];
	bcFlagRight[fUVelocity] = left->fBcFlag[inpUOffset];
	bcFlagRight[fTemperature] = left->fBcFlag[inpTOffset];
	for ( i = fFirstSpecies; i < nSpeciesInSystem+fFirstSpecies; ++i ) {
		bcFlagRight[i] = left->fBcFlagSpecies;
	}

	// set value
	yright[fVVelocity] = left->fValue[inpVOffset];
	yright[fUVelocity] = left->fValue[inpUOffset];
	yright[fTemperature] = left->fValue[inpTOffset];

	bcRight[fVVelocity] = left->fValue[inpVOffset];
	bcRight[fUVelocity] = left->fValue[inpUOffset];
	bcRight[fTemperature] = left->fValue[inpTOffset];

	for ( i = 0; i < leftSpecifiedSpecies; ++i ) {
		if ( speciesIndexRight[i] >= nSpeciesInSystem ) {
			fprintf( stderr, "%s%s%s\n", "#warning: value at left boundary for steady state species '"
				, fSpecies->GetNames()[i], "' specified, make no use of it"  );
		}
		else {
			yright[speciesIndexRight[i]+fFirstSpecies] = left->fValueSpecies[i];
			bcRight[speciesIndexRight[i]+fFirstSpecies] = left->fValueSpecies[i];
		}
	}

	if ( left->fMixtureSpecification == kMolarFraction ) {
		// first compute molar mass of mixture
		for ( i = 0, mixMolarMass = 0; i < nSpeciesInSystem; ++i ) {
			mixMolarMass += species[i].molarMass * yright[i+fFirstSpecies];
		}
		// compute massfractions
		for ( i = 0; i < nSpeciesInSystem; ++i ) {
			yright[i+fFirstSpecies] *= species[i].molarMass / mixMolarMass;
			bcRight[i+fFirstSpecies] = yright[i+fFirstSpecies];
		}
		for ( i = fFirstSpecies; i < nSpeciesInSystem+fFirstSpecies; ++i ) {
			bcFlagRight[i] = kMassFraction;
		}
	}
	
	if ( phi[0] != 0 ) {
		SetMassFracsOfPhi( this, phi[0], &bcRight[fFirstSpecies]
			, nSpeciesInSystem, fSpecies->GetMolarMass()->vec, fSpecies->GetNames() );
		SetMassFracsOfPhi( this, phi[0], &yright[fFirstSpecies]
			, nSpeciesInSystem, fSpecies->GetMolarMass()->vec, fSpecies->GetNames() );
		SetMixtureSpecificationRight( kMassFlux );

		SetMassFracsOfPhi( this, phi[0], &bcleft[fFirstSpecies]
			, nSpeciesInSystem, fSpecies->GetMolarMass()->vec, fSpecies->GetNames() );
		SetMassFracsOfPhi( this, phi[0], &yleft[fFirstSpecies]
			, nSpeciesInSystem, fSpecies->GetMolarMass()->vec, fSpecies->GetNames() );
		SetMixtureSpecificationLeft( kMassFlux );
	}
	else if ( leftSpecifiedSpecies > 0 ) {
		SetPhiOfMassFracs( this, phi, &bcRight[fFirstSpecies], fSpecies->GetMolarMass()->vec, fSpecies->GetNames() );
	}
	else {
		cerr << "#error: boundary conditions for the unburnt missing!" << NEWL;
		exit( 2 );
	}

	delete speciesIndexRight;
}

int TCountPremFlameSim::CheckComputationalDomain( void )
{
	int				last;
	TBVPSolverPtr	solver = GetSolver();
	TNewtonPtr	bt = solver->bt;
	TGridPtr 	grid = bt->GetGrid()->GetFine();
	VectorPtr	xVec = grid->GetX();
	Double		*x = xVec->vec;
	Double		**y = grid->GetY()->mat;
	int			nGridPoints = grid->GetNGridPoints();
	Double		*temp = fSolTemp->vec;
	Double		right =  bt->GetRight();
	Double		*yRight = grid->GetYRight()->vec;
	Double		slope;
	Double		slope_min = 1.0e-6;	
	Double		slope_max = 1.0e-3;	
	Double		alpha = 0.1;
	Double		deltaX;
		
	// check right boundary
	slope = fabs( FirstDerivUpwind( temp[nGridPoints-1], temp[nGridPoints-2], x[nGridPoints-1] - x[nGridPoints-2]) / temp[nGridPoints-1] );
	fprintf( stderr, "slope = %g\n", slope );
	if ( slope > slope_max ) {
		if ( !bt->GetGrid()->GetOneSolOneGrid() ) {
			cerr << "#warning: cannot enlarge computational domain for equidistant grid" << NEWL;
			return 0;
		}
		else {
			deltaX = right - x[nGridPoints-1];
			right *= 1.0 + alpha;
			bt->SetRight( right );
			yRight[fVVelocity] = y[nGridPoints-1][fVVelocity] + 
								 ( yRight[fVVelocity] - y[nGridPoints-1][fVVelocity] ) / deltaX *
								 ( right - x[nGridPoints-1] );
			cerr << "enlarge the computational domain by " << alpha * 100.0 << " percent" << NEWL << NEWL;
			return 1;
		}
	}
	else if ( slope < slope_min ) {
		if ( !bt->GetGrid()->GetOneSolOneGrid() ) {
			cerr << "#warning: cannot cut computational domain for equidistant grid" << NEWL;
			return 0;
		}
		else {
			right *= 1.0 - alpha;
			bt->SetRight( right );
			for ( last = nGridPoints - 1; x[last] > right; --last );
			grid->AdjustNGridPoints( last );
			solver->UpdateAllDimensions( last );
			CountPremSimPostIter( this );
			cerr << "cut the computational domain by " << alpha * 100.0 << " percent" << NEWL << NEWL;
			return 1;
		}
	}
	return 0;
}

Double TCountPremFlameSim::GetBurningVelocity( void )
{
	Double		*V = GetV()->vec;
	Double		*temp = GetTemperature()->vec;
	TGridPtr	currentGrid = GetSolver()->bt->GetGrid()->GetCurrentGrid();
	Double		*x = currentGrid->GetX()->vec;
	int			gridPoints = currentGrid->GetNGridPoints();
	int			inflectionPoint = LocationOfMaxSlope( temp, x, gridPoints );
	Double		rho_u = GetProperties()->GetDensity()->vec[gridPoints-1];
	Double		coeff = sqrt( GetStrainRate() * fFlameNode->rhoInf * fFlameNode->viscosityInf );

	Double		velocity = coeff * V[inflectionPoint] / rho_u;
	
	fprintf( stderr, "T0 = %g\tsL = %g\n", temp[inflectionPoint], velocity );

	return velocity;
}
