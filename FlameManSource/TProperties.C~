//#define FOOL_SOFTBENCH( x ) 

#include "FlameMaster.h"

#undef STEP4

//#define FABIANSSOOTCORRECTION 0.33 //Radiation in only one direction: 1/3
#define FABIANSSOOTCORRECTION 1


#define RADCAL
#define NEWRADCAL
#define RADCALCONST 1.0 //Radiation in only one direction: 1/3

#undef DEBUG

#define DELTAINEW

#define NGACOND

#undef SMOOKECONDUCTIVITY

/*#define SPEED*/

void TProperties::InitProperties( TInputDataPtr /*input*/, TSpeciesPtr species )
{
	int nSpecies = species->GetNOfSpecies();

//	if ( inher == kZeroDimension ) {
//		if ( input->fWithRadiation ) {
//			if ( !( fRadiation = new TRadiation( input ) ) ) FatalError( "memory allocation of TRadiation failed" );
//		}
//	}
//	else {
//		fRadiation = NULL;
//	}
	
/*#ifdef SPEED*/
/*	sqrtMjOverMi = New2DArray( nSpecies, nSpecies );*/
/*	sqrtMiOverMjPlusOne = New2DArray( nSpecies, nSpecies );*/
/**/
/*	FillProperties( species );*/
/*#else*/
/*	sqrtMjOverMi = NULL;*/
/*	sqrtMiOverMjPlusOne = NULL;*/
/*#endif*/
}

TProperties::~TProperties( void )
{
/*	Free2DArray( sqrtMiOverMjPlusOne );*/
/*	Free2DArray( sqrtMjOverMi );*/
}

/*void TProperties::FillProperties( TSpeciesPtr species )*/
/*{*/
/*	int 	i, j;*/
/*	int		nSpecies = species->GetNOfSpecies();*/
/*	Double	*M = species->GetMolarMass()->vec;*/
/**/
/*	for ( i = 0; i < nSpecies; ++i ) {*/
/*		for ( j = 0; j < nSpecies; ++j ) {*/
/*			sqrtMjOverMi[j][i] = sqrt( M[j] / M[i] );*/
/*			sqrtMiOverMjPlusOne[i][j] = 0.3535534 / sqrt( M[i] / M[j] + 1.0 );*/
/*		}*/
/*	}*/
/*}*/

void TProperties::ComputeDCpDT( Double &dCpdT, Double *Y, Double temp, TSpeciesPtr species )
{
	int 	nSpeciesInSystem = species->GetNSpeciesInSystem(); 
	
	dCpdT = 0;
	for ( int i = 0; i < nSpeciesInSystem; ++i ) {
		dCpdT += Y[i] * species->ComputeDCpiDT( i, temp );
	}
}

void TProperties::ComputeMixtureMolarMass( Double &mixMolarMass, Double *Y, Double *molarMass, int nSpeciesInSystem )
{
	int i;	
   
   	for ( i = 0, mixMolarMass = 0.0; i < nSpeciesInSystem; ++i ) {
		mixMolarMass += Y[i] / molarMass[i];
	}

	mixMolarMass = 1.0 / mixMolarMass;
#ifdef DEBUG
	cerr << "molar mass of mixture = " << mixMolarMass << NEWL;
#endif
}

void T0DProperties::InitT0DProperties( TInputDataPtr input )
{
	if ( input->fPressure ) {
		fPressure = input->fPressure->vec[0];
/*		if ( input->fPressure->len > 1 ) {
			cerr << "#warning: more than one pressure specified, use" 
					<< fPressure << NEWL;
		}*/
		cerr << "initial pressure is " << fPressure / 1.0e5 << " bar" << NEWL;
	}
	else {
		fPressure = -1.0;
//		cerr << "#error: no pressure specified" << NEWL;
//		exit( 2 );
	}
	fRadiation = NULL;
	if ( input->fWithRadiation ) {
 	  	if ( !( fRadiation = new T0DRadiation( input ) ) ) {
			FatalError( "memory allocation of TRadiation failed" );
		}
	}
}

T0DProperties::~T0DProperties( void )
{
	if ( fRadiation ) {
		delete fRadiation;
	}
}

void T0DProperties::CompMixtureProps( Double *heatCapacity, Double *conductivity, Double *mu, Double *Y, Double temp
			, Double &pressure, Double &density, EqOfState what, int nSpeciesIn, T0DSpeciesPtr species )
{
//  calculate c_p and either pressure or density, depending on 'EqOfState what'

	int 	i;
	Double	*deltaI = species->GetDeltaI()->vec;
	
// ATTENTION hp
// the following is a very dirty fix of the problem that
// in homogeneous combustion transport data is not needed, hence deltaI
// not computed. In 1D unsteady, however deltaI is computed and transport
// data is needed
	if ( deltaI[0] ) {
		fMixHeatCapacity = Y[0] * heatCapacity[0];
		fMixConductivity = Y[0] / deltaI[0] * conductivity[0];
		fMixViscosity = Y[0] / deltaI[0] * mu[0];
		for ( i = 1; i < nSpeciesIn; ++i ){
			fMixHeatCapacity += Y[i] * heatCapacity[i];
			fMixConductivity += Y[i] / deltaI[i] * conductivity[i];
			fMixViscosity += Y[i] / deltaI[i] * mu[i];
		}
	}
	else {
		fMixHeatCapacity = Y[0] * heatCapacity[0];
		for ( i = 1; i < nSpeciesIn; ++i ){
			fMixHeatCapacity += Y[i] * heatCapacity[i];
		}
	}

	switch ( what ) {
		case kDensFromPress:
			density = pressure * fMixMolarMass / ( RGAS * temp );
			break;
		case kPressFromDens:
			pressure = density * RGAS * temp / fMixMolarMass;
			break;
		default:
			cerr << "#error in function T0DProperties::CompMixtureProps" << NEWL;
			exit(2);
	}
//	fprintf( stderr, "cp = %g\trho = %g\n", fMixHeatCapacity, density );
}

void TRadiation::InitRadiation( TInputDataPtr input )
{
	fH2OIndex = input->fH2OIndex;
	fCO2Index = input->fCO2Index;
	fCH4Index = input->FindSpecies( "CH4" );
	fCOIndex = input->FindSpecies( "CO" );
}

TRadiation::~TRadiation( void )
{
}

#ifdef RADCAL
void TRadiation::ComputeRadiationOnePoint( Double *radiation, Double temp, Double *Y, Double *molarMass, Double density )
{
	static Double	constant = RADCALCONST * 4.0 * STEFBOLTZ * RGAS;

	Double			alphaCO2, alphaH2O, alphaCH4, alphaCO;
	Double			C_CO2 = ( fCO2Index == -1 ) ? 0.0 : Y[fCO2Index] / molarMass[fCO2Index];
	Double			C_H2O = ( fH2OIndex == -1 ) ? 0.0 : Y[fH2OIndex] / molarMass[fH2OIndex];
	Double			C_CH4 = ( fCH4Index == -1 ) ? 0.0 : Y[fCH4Index] / molarMass[fCH4Index];
	Double			C_CO = ( fCOIndex == -1 ) ? 0.0 : Y[fCOIndex] / molarMass[fCOIndex];
	Double			Tm1 = 1000.0/temp;
	Double			T2 = temp * temp;
	Double			c0, c1, c2, c3, c4, c5;
	Double			atmm1ToPam1 = 1.0 / 101330.0; // 1/atm to 1/Pa

#	ifdef NEWRADCAL
/* corrected version 06/15/00*/	
// first H2O
	c0 = -0.23093;
	c1 = -1.12390;
	c2 =  9.41530;
	c3 = -2.99880;
	c4 =  0.51382;
	c5 = -1.86840E-05;
	alphaH2O = c0+Tm1*(c1+Tm1*(c2+Tm1*(c3+Tm1*(c4+Tm1*c5)))); // (atm m)^-1

// then CO2
	c0 =  18.741;
	c1 = -121.310;
	c2 =  273.500;
	c3 = -194.050;
	c4 =  56.310;
	c5 = -5.8169;
	alphaCO2 = c0+Tm1*(c1+Tm1*(c2+Tm1*(c3+Tm1*(c4+Tm1*c5)))); // (atm m)^-1

#	else
	Double			logT = log( temp );
	Double			logT2 = logT * logT;
// first H2O
	c0 =  0.278713E+03;
	c1 = -0.153240E+03;
	c2 =  0.321971E+02;
	c3 = -0.300870E+01;
	c4 =  0.104055E+00;
	alphaH2O = exp( c0 + c1 * logT + c2 * logT2 + c3 * logT * logT2 + c4 * logT2 * logT2 ); // (atm m)^-1

// then CO2
	c0 =  0.96986E+03;
	c1 = -0.58838E+03;
	c2 =  0.13289E+03;
	c3 = -0.13182E+02;
	c4 =  0.48396E+00;
	alphaCO2 = exp( c0 + c1 * logT + c2 * logT2 + c3 * logT * logT2 + c4 * logT2 * logT2 ); // (atm m)^-1
#	endif

// then CH4
	alphaCH4 = 6.6334 - 0.0035686*temp + 1.6682e-08*T2 + 2.5611e-10*T2*temp - 2.6558e-14*T2*T2; // (atm m)^-1

// and CO
	if ( temp < 750.0 ) {
		c0 = 4.7869;
		c1 = -0.06953;
		c2 = 2.95775e-4;
		c3 = -4.25732e-7;
		c4 = 2.02894e-10;
	}
	else{
		c0 = 10.09;
		c1 = -0.01183;
		c2 = 4.7753e-6;
		c3 = -5.87209e-10;
		c4 = -2.5334e-14;
	}
	alphaCO =  c0 + temp * ( c1 + temp * ( c2 + temp * ( c3 + temp * c4) ) ); // (atm m)^-1

	//MM-Add fake absorption term to keep temperature above 298K
	//*radiation = -constant * density * pow( temp, 5.0 ) 
	//			* atmm1ToPam1 * ( alphaCO2 * C_CO2 + alphaH2O * C_H2O + alphaCH4 * C_CH4 + alphaCO * C_CO );
	 *radiation = -constant * density * temp * (pow( temp, 4.0 ) - pow(298.0,4.0))
	 			* atmm1ToPam1 * ( alphaCO2 * C_CO2 + alphaH2O * C_H2O + alphaCH4 * C_CH4 + alphaCO * C_CO );
// 	*radiation = -constant * density * temp * (MAX(pow( temp, 4.0 ) - pow(812.0,4.0),0.0))
// 				* atmm1ToPam1 * ( alphaCO2 * C_CO2 + alphaH2O * C_H2O + alphaCH4 * C_CH4 + alphaCO * C_CO );
 	//*radiation = -constant * density * temp * (MAX(pow( temp, 4.0 ) - pow(897.0,4.0),0.0))
 	//			* atmm1ToPam1 * ( alphaCO2 * C_CO2 + alphaH2O * C_H2O + alphaCH4 * C_CH4 + alphaCO * C_CO );
}
#else
void TRadiation::ComputeRadiationOnePoint( Double *radiation, Double temp, Double *Y, Double *molarMass, Double density )
{
	static Double	constant = FABIANSSOOTCORRECTION * 4.0 * STEFBOLTZ * RGAS;

	Double			alphaCO2, alphaH2O;
	Double			C_CO2 = ( fCO2Index == -1 ) ? 0.0 : Y[fCO2Index] / molarMass[fCO2Index];
	Double			C_H2O = ( fH2OIndex == -1 ) ? 0.0 : Y[fH2OIndex] / molarMass[fH2OIndex];
	
	alphaCO2 =  fAlphaCoeffCO2 * exp( fExpFuncCoeffCO2 * temp ); // (Pa m)^-1 eq. m^2/J
	alphaH2O =  fAlphaCoeffH2O * exp( fExpFuncCoeffH2O * temp ); // (Pa m)^-1 eq. m^2/J
	
	*radiation = -constant * density * pow( temp, 5.0 ) 
				* ( alphaCO2 * C_CO2 + alphaH2O * C_H2O );
}
#endif

#ifndef ZEROD
void T1DRadiation::InitT1DRadiation( TInputDataPtr input )
{
	fRadiation = NewVector( input->fMaxGridPoints );
}

T1DRadiation::~T1DRadiation( void )
{
	DisposeVector( fRadiation );
}

#ifdef RADCAL
void T1DRadiation::FillJacRadiation( Double coeff, T1DFlamePtr flame, NodeInfoPtr nodeInfo )
{
	static int init = 0;
	
	if ( !init ) {
		fprintf( stderr, "###warning: function T1DRadiation::FillJacRadiation not available for RADCAL radiation model\n" );
		fprintf( stderr, "###warning: might influence convergence, not results\n" );
		init = 1;
	}
}
#else
void T1DRadiation::FillJacRadiation( Double coeff, T1DFlamePtr flame, NodeInfoPtr nodeInfo )
{
	TFlameNodePtr	flameNode = flame->GetFlameNode();
	int				tempOff = flame->GetOffsetTemperature();
	int				specOff = flame->GetOffsetFirstSpecies();
	Double			temp = flameNode->temp[kCurr];
	Double			*Y = flameNode->Y[kCurr];
	Double			*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
	static Double	constant = FABIANSSOOTCORRECTION * 2.0 * STEFBOLTZ * RGAS;
	Double			**a = nodeInfo->a;
	Double			hnenn = nodeInfo->hnenn;

	Double			alphaCO2, alphaH2O;
	Double			mixDensity = *flameNode->mixDensity;
	Double			constRhoT5;
	
	constant *= coeff * hnenn;
	constRhoT5 = constant * mixDensity * pow( temp, 5.0 );

	alphaCO2 =  fAlphaCoeffCO2 * exp( fExpFuncCoeffCO2 * temp );
	alphaH2O =  fAlphaCoeffH2O * exp( fExpFuncCoeffH2O * temp );
	
	if ( fCO2Index != -1 ) {
		a[fCO2Index+specOff][tempOff] += -constRhoT5 * alphaCO2 / molarMass[fCO2Index];
	} 
	if ( fH2OIndex != -1 ) {
		a[fH2OIndex+specOff][tempOff] += -constRhoT5 * alphaH2O / molarMass[fH2OIndex]; 
	} 

	a[tempOff][tempOff] += 5.0 * coeff * flameNode->radiation[kCurr] / temp * hnenn;
	if ( fCO2Index != -1 ) {
		a[tempOff][tempOff] += -constRhoT5 
				* alphaCO2 * fExpFuncCoeffCO2 * Y[fCO2Index] / molarMass[fCO2Index];
	} 
	if ( fH2OIndex != -1 ) {
		a[tempOff][tempOff] += -constRhoT5 
				* alphaH2O * fExpFuncCoeffH2O * Y[fH2OIndex] / molarMass[fH2OIndex];
	} 
}
#endif

void T1DProperties::InitT1DProperties( TInputDataPtr input )
{
	int	maxGridPoints = input->fMaxGridPoints;

	fViscosity = NewVector( maxGridPoints+2 );
	fViscosity->vec = &fViscosity->vec[kNext];
	fConductivity = NewVector( maxGridPoints+2 );
	fConductivity->vec = &fConductivity->vec[kNext];
	fDensity = NewVector( maxGridPoints+2 );
	fDensity->vec = &fDensity->vec[kNext];
	fHeatCapacity = NewVector( maxGridPoints+2 );
	fHeatCapacity->vec = &fHeatCapacity->vec[kNext];
	fMolarMass = NewVector( maxGridPoints+2 );
	fMolarMass->vec = &fMolarMass->vec[kNext];
	fRadiation = NULL;
	if ( input->fWithRadiation ) {
 	  	if ( !( fRadiation = new T1DRadiation( input ) ) ) FatalError( "memory allocation of TRadiation failed" );
	}
	fdCpdT = NewVector( maxGridPoints );
}

T1DProperties::~T1DProperties( void )
{
	DisposeVector( fdCpdT );
	if ( fRadiation ) {
		delete fRadiation;
	}
	fMolarMass->vec = &fMolarMass->vec[kPrev];
	DisposeVector( fMolarMass );
	fHeatCapacity->vec = &fHeatCapacity->vec[kPrev];
	DisposeVector( fHeatCapacity );
	fDensity->vec = &fDensity->vec[kPrev];
	DisposeVector( fDensity );
	fConductivity->vec = &fConductivity->vec[kPrev];
	DisposeVector( fConductivity );
	fViscosity->vec = &fViscosity->vec[kPrev];
	DisposeVector( fViscosity );
}

void T1DProperties::CompMixtureProps( TFlameNodePtr flameNode, Double *Y, Double temp
					, Double pressure, TSpeciesPtr species )
{
//  calculate c_p, lambda and mu of the mixture

	int 	i;
	int		nSpeciesInSystem = species->GetNSpeciesInSystem();
	Double	y_over_delta = 0.0;
	Double	*c_p = flameNode->heatCapacity;
	Double	*lambda = flameNode->conductivity;
	Double	*mu = flameNode->viscosity;
	Double	*molarMass = species->GetMolarMass()->vec;
	Double	&mixHeatCapacity = *flameNode->mixHeatCapacity;
	Double	&mixConductivity = *flameNode->mixConductivity;
	Double	&mixViscosity = *flameNode->mixViscosity;
	Double	&mixDensity = *flameNode->mixDensity;
	Double	mixMolarMass = *flameNode->mixMolarMass;
	Double	*deltaI = flameNode->deltaI;
	

#ifdef STEP4
	mixDensity = pressure * mixMolarMass / ( RGAS * temp );
	mixHeatCapacity = 1429; /* [J / kg K]*/
	mixConductivity = 2.58e-5 * pow( ( *flameNode->temp / 298.0 ), 0.7 ) * mixHeatCapacity;
	mixViscosity = 0.75 * mixConductivity / mixHeatCapacity;
#else
	mixHeatCapacity = Y[0] * c_p[0];
	
//	y_over_delta = Y[0] / CompDeltaI( 0, nSpeciesInSystem, molarMass, mu, Y );
	y_over_delta = Y[0] / deltaI[0];

	mixConductivity = y_over_delta * lambda[0];
	mixViscosity = y_over_delta * mu[0];

	for ( i = 1; i < nSpeciesInSystem; ++i ){
		mixHeatCapacity += Y[i] * c_p[i];
	
//		y_over_delta = Y[i] / CompDeltaI( i, nSpeciesInSystem, molarMass, mu, Y );
		y_over_delta = Y[i] / deltaI[i];
		mixConductivity += y_over_delta * lambda[i];
		mixViscosity += y_over_delta * mu[i];
	}

#ifdef NGACOND
	//Recomp mixConductivity
	mixConductivity = 0.0;
	Double sum_1 = 0.0;
	Double sum_2 = 0.0;
	for ( i = 0; i < nSpeciesInSystem; ++i ){
	  sum_1=sum_1+Y[i]*mixMolarMass/(lambda[i]*molarMass[i]);
	  sum_2=sum_2+Y[i]*lambda[i]*mixMolarMass/molarMass[i];
	}
	mixConductivity = 0.5*(sum_2+(1.0/sum_1));
#endif

#	ifdef SMOOKECONDUCTIVITY
	mixConductivity = 2.58e-5 * pow( ( *flameNode->temp / 298.0 ), 0.7 ) * mixHeatCapacity;
#	endif
	mixDensity = pressure * mixMolarMass / ( RGAS * temp );
#endif
}

void T1DProperties::PrintProperties( int k )
{
	static Flag 	init = FALSE;
	FILE		*fp = NULL;
	
	if ( !init ) {
		if ( !( fp = fopen( "properties.tout", "w" ) ) ) { 
			cerr << "#warning: unable to open file 'properties.tout'" << NEWL;
			return;
		}
		init = TRUE;
	}
	else {
		if ( !( fp = fopen( "properties.tout", "a" ) ) ) { 
			cerr << "#warning: unable to open file 'properties.tout'" << NEWL;
			return;
		}
	}
	fprintf( fp, "Viscosity = %g\n", fViscosity->vec[k] );
	fprintf( fp, "Density = %g\n", fDensity->vec[k] );
	fprintf( fp, "HeatCapacity = %g\n", fHeatCapacity->vec[k] );
	fprintf( fp, "Conductivity = %g\n", fConductivity->vec[k] );
	fprintf( fp, "fMolarMass = %g\n", fMolarMass->vec[k] );
	
	fprintf( fp, "\n\n\n");
}

void T1DProperties::PrintProperties( TNewtonPtr bt )
{
	char		fName[32];
	FILE		*fp;
	TGridPtr	grid = bt->GetGrid()->GetCurrentGrid();
	Double		*x = grid->GetX()->vec;
	Double		*mu = fViscosity->vec;
	Double		*rho = fDensity->vec;
	Double		*cp = fHeatCapacity->vec;
	Double		*lambda = fConductivity->vec;
	Double		*molarMass = fMolarMass->vec;
	int			nGridPoints = bt->GetCurrentGridPoints();
	Double		left = bt->GetLeft();
	Double		right = bt->GetRight();
	static int	counter = 0;
	
	sprintf( fName, "props%d.dout", ++counter );
	if ( counter >= 10 ) counter = 0;
	if ( !( fp = fopen( fName, "w" ) ) ) { 
		cerr << "#warning: unable to open file " << fName << NEWL;
		return;
	}
	fprintf( fp, "*\n" );
	fprintf( fp, "%-12s\t%-12s\t%-12s\t%-12s\t%-12s\t%-12s\n", "eta", "mu", "rho", "cp", "lambda", "M" );
	fprintf( fp, "%-9E\t%-9E\t%-9E\t%-9E\t%-9E\t%-9E\n", left, mu[-1], rho[-1], cp[-1], lambda[-1], molarMass[-1] );
	for ( int i = 0; i < nGridPoints; ++i ) {
		fprintf( fp, "%-9E\t%-9E\t%-9E\t%-9E\t%-9E\t%-9E\n", x[i], mu[i], rho[i], cp[i], lambda[i], molarMass[i] );
	}
	fprintf( fp, "%-9E\t%-9E\t%-9E\t%-9E\t%-9E\t%-9E\n", right, mu[nGridPoints], rho[nGridPoints], cp[nGridPoints], lambda[nGridPoints], molarMass[nGridPoints] );
	fclose( fp );
}
#endif // ZEROD
