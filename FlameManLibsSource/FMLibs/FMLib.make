#   File:       FMLib.make
#   Target:     FMLib
#   Sources:    Interrupt.c
#               LTBody.c
#               LTHeader.c
#               LTUtilities.c
#               MapMan.cp
#               NonlinearSystems.c
#               ReadStartProfile.c
#               SmallNewton.c
#               Spline.c
#   Created:    Wednesday, January 15, 1997 04:50:37 PM


MAKEFILE     = FMLib.make
¥MondoBuild¥ = {MAKEFILE}  # Make blank to avoid rebuilds when makefile is modified
Headers      = Interrupt.h ¶
				ListTool.h ¶
				MapMan.h ¶
				NonlinearSystems.h ¶
				paramtr.h ¶
				SmallNewton.h ¶
				Spline.h ¶
				TofZ.h
Includes     =
Sym¥PPC      = 
ObjDir¥PPC   =
FCompileOptions = -opt=3

PPCCOptions  = {Includes} {Sym¥PPC} -w 2
PPCCPlusOptions = {Includes} {Sym¥PPC} -w 2

Objects¥PPC  = ¶
		"{ObjDir¥PPC}Interrupt.c.x" ¶
		"{ObjDir¥PPC}LTBody.c.x" ¶
		"{ObjDir¥PPC}LTHeader.c.x" ¶
		"{ObjDir¥PPC}LTUtilities.c.x" ¶
		"{ObjDir¥PPC}MapMan.cp.x" ¶
		"{ObjDir¥PPC}NonlinearSystems.c.x" ¶
		"{ObjDir¥PPC}ReadStartProfile.c.x" ¶
		"{ObjDir¥PPC}SmallNewton.c.x" ¶
		"{ObjDir¥PPC}Spline.c.x" ¶
		"{ObjDir¥PPC}TofZ.c.x" ¶
		"{ObjDir¥PPC}Adiab_flam.f.o" ¶
		"{ObjDir¥PPC}PDFF.f.o"


FMLib ÄÄ {¥MondoBuild¥} {Objects¥PPC}
	PPCLink ¶
		-o {Targ} {Sym¥PPC} ¶
		{Objects¥PPC} ¶
		-xm l
	Move -y {Targ} "{myLibs}"FMLib.lib 
	Duplicate -y {Headers} "{myIncludes}"


"{ObjDir¥PPC}Interrupt.c.x" Ä {¥MondoBuild¥} Interrupt.c
	{PPCC} Interrupt.c -o {Targ} {PPCCOptions}

"{ObjDir¥PPC}LTBody.c.x" Ä {¥MondoBuild¥} LTBody.c
	{PPCC} LTBody.c -o {Targ} {PPCCOptions}

"{ObjDir¥PPC}LTHeader.c.x" Ä {¥MondoBuild¥} LTHeader.c
	{PPCC} LTHeader.c -o {Targ} {PPCCOptions}

"{ObjDir¥PPC}LTUtilities.c.x" Ä {¥MondoBuild¥} LTUtilities.c
	{PPCC} LTUtilities.c -o {Targ} {PPCCOptions}

"{ObjDir¥PPC}MapMan.cp.x" Ä {¥MondoBuild¥} MapMan.cp
	{PPCCPlus} MapMan.cp -o {Targ} {PPCCPlusOptions}

"{ObjDir¥PPC}NonlinearSystems.c.x" Ä {¥MondoBuild¥} NonlinearSystems.c
	{PPCC} NonlinearSystems.c -o {Targ} {PPCCOptions}

"{ObjDir¥PPC}ReadStartProfile.c.x" Ä {¥MondoBuild¥} ReadStartProfile.c
	{PPCC} ReadStartProfile.c -o {Targ} {PPCCOptions}

"{ObjDir¥PPC}SmallNewton.c.x" Ä {¥MondoBuild¥} SmallNewton.c
	{PPCC} SmallNewton.c -o {Targ} {PPCCOptions}

"{ObjDir¥PPC}Spline.c.x" Ä {¥MondoBuild¥} Spline.c
	{PPCC} Spline.c -o {Targ} {PPCCOptions}

"{ObjDir¥PPC}TofZ.c.x" Ä {¥MondoBuild¥} TofZ.c
	{PPCC} TofZ.c -o {Targ} {PPCCOptions}

"{ObjDir¥PPC}Adiab_flam.f.o" Ä {¥MondoBuild¥} Adiab_flam.f
	 FORTRAN.PPC Adiab_flam.f {FCompileOptions} 

"{ObjDir¥PPC}PDFF.f.o" Ä {¥MondoBuild¥} PDFF.f
	 FORTRAN.PPC PDFF.f {FCompileOptions} 


