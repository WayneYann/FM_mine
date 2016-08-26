#define MECHANISM "n_butanol_91"
#include "FlameMaster.h"

/*	Mechanism file: "n_butanol_91.mech"	*/

typedef Double (*BFFunction)(Double T);

/* prototypes */
Double Fc18( Double T );
Double Fc57( Double T );
Double Fc76( Double T );
Double Fc86( Double T );
Double Fc97( Double T );
Double Fc129( Double T );
Double Fc149( Double T );
Double Fc163( Double T );
Double Fc200( Double T );
Double Fc201( Double T );
Double Fc208( Double T );
Double Fc228( Double T );
Double Fc229( Double T );
Double Fc281( Double T );
Double Fc288( Double T );
Double Fc378( Double T );
Double Fc447( Double T );
Double Fc448( Double T );
Double Fc453( Double T );
Double Fc464( Double T );
Double Fc564( Double T );
Double Fc570( Double T );
Double Fc591( Double T );
Double Fc715( Double T );
Double FcErr( Double T );


extern BFFunction gBroadening[24];
