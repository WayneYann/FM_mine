#   File:       alligator.make
#   Target:     alligator
#   Sources:    alligator.c
#   Created:    Tuesday, January 14, 1997 06:10:41 PM


MAKEFILE     = alligator.make
¥MondoBuild¥ = {MAKEFILE}  # Make blank to avoid rebuilds when makefile is modified
Includes     =
Sym¥PPC      = 
ObjDir¥PPC   =

PPCCOptions  = {Includes} {Sym¥PPC} 

Objects¥PPC  = ¶
		"{ObjDir¥PPC}alligator.c.x"


alligator ÄÄ {¥MondoBuild¥} {Objects¥PPC}
	PPCLink ¶
		-o {Targ} {Sym¥PPC} ¶
		{Objects¥PPC} ¶
		-xm l
		move {Targ} {myLibs}alligator.lib


"{ObjDir¥PPC}alligator.c.x" Ä {¥MondoBuild¥} alligator.c
	{PPCC} alligator.c -o {Targ} {PPCCOptions}

