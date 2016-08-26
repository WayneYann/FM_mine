#   File:       NewtonLib.make
#   Target:     NewtonLib
#   Sources:    AdaptiveGrid.cp
#               Newton.cp
#               NewtonUt.cp
#   Created:    Wednesday, January 15, 1997 04:46:55 PM


MAKEFILE     = NewtonLib.make
¥MondoBuild¥ = {MAKEFILE}  # Make blank to avoid rebuilds when makefile is modified
Includes     =
Sym¥PPC      = 
ObjDir¥PPC   =

PPCCPlusOptions = {Includes} {Sym¥PPC}  -w 2

Objects¥PPC  = ¶
		"{ObjDir¥PPC}AdaptiveGrid.cp.x" ¶
		"{ObjDir¥PPC}Newton.cp.x" ¶
		"{ObjDir¥PPC}NewtonUt.cp.x"


NewtonLib ÄÄ {¥MondoBuild¥} {Objects¥PPC}
	PPCLink ¶
		-o {Targ} {Sym¥PPC} ¶
		{Objects¥PPC} ¶
		-xm l
	move -y {Targ} "{myLibs}"NewtonLib.lib
	duplicate -y Newton.h "{myIncludes}"

"{ObjDir¥PPC}AdaptiveGrid.cp.x" Ä {¥MondoBuild¥} AdaptiveGrid.cp Newton.h
	{PPCCPlus} AdaptiveGrid.cp -o {Targ} {PPCCPlusOptions}

"{ObjDir¥PPC}Newton.cp.x" Ä {¥MondoBuild¥} Newton.cp Newton.h
	{PPCCPlus} Newton.cp -o {Targ} {PPCCPlusOptions}

"{ObjDir¥PPC}NewtonUt.cp.x" Ä {¥MondoBuild¥} NewtonUt.cp Newton.h
	{PPCCPlus} NewtonUt.cp -o {Targ} {PPCCPlusOptions}

