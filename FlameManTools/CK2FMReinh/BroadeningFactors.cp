/*  The functions in this file compute the broadening factors.
 *  The names of the functions start with "Fc" followed by
 *  the label of the corresponding reaction.
 */

#include"BroadeningFactors.h"

/*	Mechanism file: "n_butanol_91.mech"	*/

BFFunction gBroadening[24] = { Fc18, Fc57, Fc76, Fc86, Fc97, Fc129, Fc149, Fc163, Fc200, Fc201, Fc208, Fc228, Fc229, Fc281, Fc288, Fc378, Fc447, Fc448, Fc453, Fc464, Fc564, Fc570, Fc591, Fc715 };

#ifndef MECHANISM
#define MECHANISM ""
#endif

void TReaction::CheckBroadeningFactors( const char *mechName )
{
	char	*name = new char[strlen( MECHANISM ) + 6];
	sprintf( name, "/%s.pre", MECHANISM );
	if ( strstr( mechName, name ) == NULL ) {
		for ( int i = 0; i < 24; ++i ) {
			gBroadening[i] = FcErr;
		}
	}
}

Double FcErr( Double /*T*/ )
{
	fprintf( stderr, "#error: wrong broadening factors (%s) linked to program\n", MECHANISM );
	exit( 2 );

	return 0;
}

Double Fc18( Double T )
{
#line 27 "n_butanol_91.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc57( Double T )
{
#line 70 "n_butanol_91.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc76( Double T )
{
#line 93 "n_butanol_91.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc86( Double T )
{
#line 107 "n_butanol_91.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc97( Double T )
{
#line 122 "n_butanol_91.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc129( Double T )
{
#line 158 "n_butanol_91.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc149( Double T )
{
#line 182 "n_butanol_91.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc163( Double T )
{
#line 200 "n_butanol_91.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc200( Double T )
{
#line 241 "n_butanol_91.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc201( Double T )
{
#line 246 "n_butanol_91.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc208( Double T )
{
#line 257 "n_butanol_91.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc228( Double T )
{
#line 280 "n_butanol_91.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc229( Double T )
{
#line 284 "n_butanol_91.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc281( Double T )
{
#line 340 "n_butanol_91.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc288( Double T )
{
#line 351 "n_butanol_91.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc378( Double T )
{
#line 445 "n_butanol_91.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc447( Double T )
{
#line 517 "n_butanol_91.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc448( Double T )
{
#line 521 "n_butanol_91.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc453( Double T )
{
#line 529 "n_butanol_91.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc464( Double T )
{
#line 543 "n_butanol_91.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc564( Double T )
{
#line 647 "n_butanol_91.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc570( Double T )
{
#line 657 "n_butanol_91.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc591( Double T )
{
#line 682 "n_butanol_91.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc715( Double T )
{
#line 810 "n_butanol_91.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

