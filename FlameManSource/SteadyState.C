#include "FlameMaster.h"
#include "CH4.72GlobalFive.h"

void ComputeSteadyStates( Double *k, Double *c, Double *M )
{
	Double	w[rEnd];
	Double	K2 = k[r2f] / k[r2b];
	Double	K3 = k[r3f] / k[r3b];
	Double	K38 = k[r38f] / k[r38b];
	Double	OHFact = c[sH2O] / CatchZero( 
		 			c[sH2] * K3 );
	Double	OFact = c[sH2O] / CatchZero( 
		 			c[sH2] * c[sH2] * K2 * K3 );
	Double	CH3Fact = K38 * c[sCH4] / CatchZero( c[sH2] );

/*	c[sCH3] = CH3Fact * c[sH];*/
/*	c[sCH3] = SolveQuadratic( 2.0 * k[r36f]*/
/*					, CatchZero( k[r34] * c[sH] + k[r35] * c[sO] + k[r37f] * c[sO2] + k[r38b] * c[sH2] )*/
/*					, k[r38f] * c[sCH4] * c[sH] );*/
				
	
	c[sO] = OFact * c[sH] * c[sH];

	c[sOH] = OHFact * c[sH];

	w[r1f] = k[r1f] * c[sH] * c[sO2];
	w[r1b] = k[r1b] * c[sOH] * c[sO];
	w[r5f] = k[r5f] * c[sH] * c[sO2] * M[mM1];
	w[r18f] = k[r18f] * c[sCO] * c[sOH];
	w[r18b] = k[r18b] * c[sCO2] * c[sH];
	w[r34] = k[r34] * c[sCH3] * c[sH];
	w[r35] = k[r35] * c[sCH3] * c[sO];
	w[r36f] = k[r36f] * c[sCH3] * c[sCH3];
	w[r37f] = k[r37f] * c[sCH3] * c[sO2];
	w[r38f] = k[r38f] * c[sCH4] * c[sH];
	w[r38b] = k[r38b] * c[sCH3] * c[sH2];
	w[r40f] = k[r40f] * c[sCH4] * c[sOH];
	w[r40b] = k[r40b] * c[sCH3] * c[sH2O];

	k[rI] = -w[r34] + w[r38f] - w[r38b] /*+ w[r40f] - w[r40b]*/;	
	k[rII] = 0.5 * w[r35] + w[r36f] + 0.5 * w[r37f];
	k[rIII] = w[r18f]/* - w[r18b]*/;
	k[rIV] = w[r5f] + w[r34];
	k[rV] = w[r1f] - w[r1b] - 0.5 * w[r35] + 0.5 * w[r37f];

	for ( int j = r1f; j < rEnd; ++j ) {
		k[j] = 0.0;
	}
}

int SteadyStatesFunc( const VectorPtr /*x*/, VectorPtr /*fVec*/, void */*object*/ )
{
	return 0;
}

#ifndef MECHANISM
#define MECHANISM ""
#endif

void TReaction::CheckSteadyStatesMech( const char *mechName )
{
	char	*name = new char[strlen( MECHANISM ) + 5];
	sprintf( name, "/%s.pre", MECHANISM );
	if ( strstr( mechName, name ) == NULL ) {
		cerr << "#error: program linked with mechanism " << MECHANISM << ".mech"
					<< ", input mechanism is " << mechName << NEWL;
		exit(2);
	}
	else {
		cerr << "use reduced kinetics " << MECHANISM << NEWL;
	}
}
