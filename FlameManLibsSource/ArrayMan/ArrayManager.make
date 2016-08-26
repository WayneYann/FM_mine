#   File:       ArrayManager.make
#   Target:     ArrayManager
#   Sources:    AM_alloc.c
#               AM_blas.c
#               AM_BlockTriDiagSolver.c
#               AM_Gauss.c
#               AM_Printing.c
#               ArrayManager.c
#   Created:    Wednesday, January 15, 1997 09:59:08 AM


MAKEFILE     = ArrayManager.make
¥MondoBuild¥ = {MAKEFILE}  # Make blank to avoid rebuilds when makefile is modified
Includes     =
Sym¥PPC      = 
ObjDir¥PPC   =

PPCCOptions  = {Includes} {Sym¥PPC} 

Objects¥PPC  = ¶
		"{ObjDir¥PPC}AM_alloc.c.x" ¶
		"{ObjDir¥PPC}AM_blas.c.x" ¶
		"{ObjDir¥PPC}AM_BlockTriDiagSolver.c.x" ¶
		"{ObjDir¥PPC}AM_Gauss.c.x" ¶
		"{ObjDir¥PPC}AM_Printing.c.x" ¶
		"{ObjDir¥PPC}ArrayManager.c.x"


ArrayManager ÄÄ {¥MondoBuild¥} {Objects¥PPC}
	PPCLink ¶
		-o {Targ} {Sym¥PPC} ¶
		{Objects¥PPC} ¶
		-xm l
		move {Targ} {myLibs}ArrayManager.lib
		


"{ObjDir¥PPC}AM_alloc.c.x" Ä {¥MondoBuild¥} AM_alloc.c
	{PPCC} AM_alloc.c -o {Targ} {PPCCOptions}

"{ObjDir¥PPC}AM_blas.c.x" Ä {¥MondoBuild¥} AM_blas.c
	{PPCC} AM_blas.c -o {Targ} {PPCCOptions}

"{ObjDir¥PPC}AM_BlockTriDiagSolver.c.x" Ä {¥MondoBuild¥} AM_BlockTriDiagSolver.c
	{PPCC} AM_BlockTriDiagSolver.c -o {Targ} {PPCCOptions}

"{ObjDir¥PPC}AM_Gauss.c.x" Ä {¥MondoBuild¥} AM_Gauss.c
	{PPCC} AM_Gauss.c -o {Targ} {PPCCOptions}

"{ObjDir¥PPC}AM_Printing.c.x" Ä {¥MondoBuild¥} AM_Printing.c
	{PPCC} AM_Printing.c -o {Targ} {PPCCOptions}

"{ObjDir¥PPC}ArrayManager.c.x" Ä {¥MondoBuild¥} ArrayManager.c
	{PPCC} ArrayManager.c -o {Targ} {PPCCOptions}

