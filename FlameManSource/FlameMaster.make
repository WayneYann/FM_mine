#   File:       FlameMaster.make
#   Target:     FlameMaster
#   Sources:    FlameMasterMain.cp
#               FlameUtilities.cp
#               TFlame.cp
#               TInputData.cp
#               TProperties.cp
#               TReaction.cp
#               SteadyStateDefault.cp
#               TSoot.cp
#               TSpecies.cp
#               TPAH.cp
#               TCountDiffFlameCont.cp
#               TCountDiffFlameMix.cp
#               TCountDiffFlamePhys.cp
#               TCountDiffFlameSim.cp
#               TCountDiffPhysEigen.cp
#               TFlameSheet.cp
#               T0DIsoChor.cp
#               TTransFlamelet.cp
#               TTransFlameSolver.cp
#               TCountPremFlamePhys.cp
#               TCountPremFlameSim.cp
#               TUnstrPremFlamePhys.cp
#               TTrans1DIsoChor.cp
#               TTrans1DIsoChorSolver.cp
#   Created:    Wednesday, January 15, 1997 04:58:34 PM


MAKEFILE     = FlameMaster.make
¥MondoBuild¥ = {MAKEFILE} FlameMaster.h # Make blank to avoid rebuilds when makefile is modified
Includes     =
Sym¥PPC      = #-sym on
ObjDir¥PPC   =

PPCCPlusOptions = {Includes} {Sym¥PPC} -ansi strict -w 2

Objects¥PPC  = ¶
		"{ObjDir¥PPC}FlameMasterMain.cp.x" ¶
		"{ObjDir¥PPC}FlameUtilities.cp.x" ¶
		"{ObjDir¥PPC}TFlame.cp.x" ¶
		"{ObjDir¥PPC}TInputData.cp.x" ¶
		"{ObjDir¥PPC}TProperties.cp.x" ¶
		"{ObjDir¥PPC}TReaction.cp.x" ¶
		"{ObjDir¥PPC}SteadyStates.cp.x" ¶
		"{ObjDir¥PPC}TSoot.cp.x" ¶
		"{ObjDir¥PPC}TSpecies.cp.x" ¶
		"{ObjDir¥PPC}TPAH.cp.x" ¶
		"{ObjDir¥PPC}TCountDiffFlameCont.cp.x" ¶
		"{ObjDir¥PPC}TCountDiffFlameMix.cp.x" ¶
		"{ObjDir¥PPC}TCountDiffFlamePhys.cp.x" ¶
		"{ObjDir¥PPC}TCountDiffFlameSim.cp.x" ¶
		"{ObjDir¥PPC}TCountDiffPhysEigen.cp.x" ¶
		"{ObjDir¥PPC}TFlameSheet.cp.x" ¶
		"{ObjDir¥PPC}T0DIsoChor.cp.x" ¶
		"{ObjDir¥PPC}TTransFlamelet.cp.x" ¶
		"{ObjDir¥PPC}TTransFlameSolver.cp.x" ¶
		"{ObjDir¥PPC}TCountPremFlamePhys.cp.x" ¶
		"{ObjDir¥PPC}TCountPremFlameSim.cp.x" ¶
		"{ObjDir¥PPC}TUnstrPremFlamePhys.cp.x" ¶
		"{ObjDir¥PPC}TTrans1DIsoChor.cp.x" ¶
		"{ObjDir¥PPC}BroadeningFactors.cp.x" ¶
		"{ObjDir¥PPC}lex.yy.cp.x" ¶
		"{ObjDir¥PPC}TTrans1DIsoChorSolver.cp.x"


FlameMaster ÄÄ {¥MondoBuild¥} {Objects¥PPC}
	PPCLink ¶
		-d -o {Targ} {Sym¥PPC} ¶
		{Objects¥PPC} ¶
		-t 'MPST' ¶
		-c 'MPS ' ¶
		"{myLibs}"NewtonLib.lib ¶
		"{myLibs}"FMLib.lib ¶
		"{myLibs}"ArrayManager.lib ¶
		"{myLibs}"alligator.lib ¶
		"{myLibs}"blas.lib ¶
		"{myLibs}"linpack.lib ¶
		"{myLibs}"dassl.lib ¶
		"{myLibs}"UniStd.lib ¶
		"{SharedLibraries}InterfaceLib" ¶
		"{SharedLibraries}StdCLib" ¶
		"{SharedLibraries}MathLib" ¶
		"{PPCLibraries}StdCRuntime.o" ¶
		"{PPCLibraries}PPCCRuntime.o" ¶
		"{PPCLibraries}PPCToolLibs.o" ¶
		"{PPCLibraries}MrCPlusLib.o" ¶
		"{PPCLibraries}MrCIOStreams.o" ¶
		"{PPCFLibraries}FortranLibPPC.o"
	move -y FlameMaster "{myTools}"FlameMaster


"{ObjDir¥PPC}FlameMasterMain.cp.x" Ä {¥MondoBuild¥} FlameMasterMain.cp
	{PPCCPlus} FlameMasterMain.cp -o {Targ} {PPCCPlusOptions}

"{ObjDir¥PPC}FlameUtilities.cp.x" Ä {¥MondoBuild¥} FlameUtilities.cp
	{PPCCPlus} FlameUtilities.cp -o {Targ} {PPCCPlusOptions}

"{ObjDir¥PPC}TFlame.cp.x" Ä {¥MondoBuild¥} TFlame.cp
	{PPCCPlus} TFlame.cp -o {Targ} {PPCCPlusOptions}

"{ObjDir¥PPC}TInputData.cp.x" Ä {¥MondoBuild¥} TInputData.cp
	{PPCCPlus} TInputData.cp -o {Targ} {PPCCPlusOptions}

"{ObjDir¥PPC}TProperties.cp.x" Ä {¥MondoBuild¥} TProperties.cp
	{PPCCPlus} TProperties.cp -o {Targ} {PPCCPlusOptions}

"{ObjDir¥PPC}TReaction.cp.x" Ä {¥MondoBuild¥} TReaction.cp
	{PPCCPlus} TReaction.cp -o {Targ} {PPCCPlusOptions}

"{ObjDir¥PPC}SteadyStates.cp.x" Ä {¥MondoBuild¥} SteadyStates.cp
	{PPCCPlus} SteadyStates.cp -o {Targ} {PPCCPlusOptions}

"{ObjDir¥PPC}TSoot.cp.x" Ä {¥MondoBuild¥} TSoot.cp
	{PPCCPlus} TSoot.cp -o {Targ} {PPCCPlusOptions}

"{ObjDir¥PPC}TSpecies.cp.x" Ä {¥MondoBuild¥} TSpecies.cp
	{PPCCPlus} TSpecies.cp -o {Targ} {PPCCPlusOptions}

"{ObjDir¥PPC}TPAH.cp.x" Ä {¥MondoBuild¥} TPAH.cp
	{PPCCPlus} TPAH.cp -o {Targ} {PPCCPlusOptions}

"{ObjDir¥PPC}TCountDiffFlameCont.cp.x" Ä {¥MondoBuild¥} TCountDiffFlameCont.cp
	{PPCCPlus} TCountDiffFlameCont.cp -o {Targ} {PPCCPlusOptions}

"{ObjDir¥PPC}TCountDiffFlameMix.cp.x" Ä {¥MondoBuild¥} TCountDiffFlameMix.cp
	{PPCCPlus} TCountDiffFlameMix.cp -o {Targ} {PPCCPlusOptions}

"{ObjDir¥PPC}TCountDiffFlamePhys.cp.x" Ä {¥MondoBuild¥} TCountDiffFlamePhys.cp
	{PPCCPlus} TCountDiffFlamePhys.cp -o {Targ} {PPCCPlusOptions}

"{ObjDir¥PPC}TCountDiffFlameSim.cp.x" Ä {¥MondoBuild¥} TCountDiffFlameSim.cp
	{PPCCPlus} TCountDiffFlameSim.cp -o {Targ} {PPCCPlusOptions}

"{ObjDir¥PPC}TCountDiffPhysEigen.cp.x" Ä {¥MondoBuild¥} TCountDiffPhysEigen.cp
	{PPCCPlus} TCountDiffPhysEigen.cp -o {Targ} {PPCCPlusOptions}

"{ObjDir¥PPC}TFlameSheet.cp.x" Ä {¥MondoBuild¥} TFlameSheet.cp
	{PPCCPlus} TFlameSheet.cp -o {Targ} {PPCCPlusOptions}

"{ObjDir¥PPC}T0DIsoChor.cp.x" Ä {¥MondoBuild¥} T0DIsoChor.cp
	{PPCCPlus} T0DIsoChor.cp -o {Targ} {PPCCPlusOptions}

"{ObjDir¥PPC}TTransFlamelet.cp.x" Ä {¥MondoBuild¥} TTransFlamelet.cp
	{PPCCPlus} TTransFlamelet.cp -o {Targ} {PPCCPlusOptions}

"{ObjDir¥PPC}TTransFlameSolver.cp.x" Ä {¥MondoBuild¥} TTransFlameSolver.cp
	{PPCCPlus} TTransFlameSolver.cp -o {Targ} {PPCCPlusOptions}

"{ObjDir¥PPC}TCountPremFlamePhys.cp.x" Ä {¥MondoBuild¥} TCountPremFlamePhys.cp
	{PPCCPlus} TCountPremFlamePhys.cp -o {Targ} {PPCCPlusOptions}

"{ObjDir¥PPC}TCountPremFlameSim.cp.x" Ä {¥MondoBuild¥} TCountPremFlameSim.cp
	{PPCCPlus} TCountPremFlameSim.cp -o {Targ} {PPCCPlusOptions}

"{ObjDir¥PPC}TUnstrPremFlamePhys.cp.x" Ä {¥MondoBuild¥} TUnstrPremFlamePhys.cp
	{PPCCPlus} TUnstrPremFlamePhys.cp -o {Targ} {PPCCPlusOptions}

"{ObjDir¥PPC}TTrans1DIsoChor.cp.x" Ä {¥MondoBuild¥} TTrans1DIsoChor.cp
	{PPCCPlus} TTrans1DIsoChor.cp -o {Targ} {PPCCPlusOptions}

"{ObjDir¥PPC}BroadeningFactors.cp.x" Ä {¥MondoBuild¥} BroadeningFactors.cp
	{PPCCPlus} BroadeningFactors.cp -o {Targ} {PPCCPlusOptions}

"{ObjDir¥PPC}TTrans1DIsoChorSolver.cp.x" Ä {¥MondoBuild¥} TTrans1DIsoChorSolver.cp
	{PPCCPlus} TTrans1DIsoChorSolver.cp -o {Targ} {PPCCPlusOptions}

"{ObjDir¥PPC}dassl.cp.x" Ä {¥MondoBuild¥} dassl.cp
	{PPCCPlus} dassl.cp -o {Targ} {Includes} {Sym¥PPC} -w 2

lex.yy.cp Ä ReadAddData.flex {¥MondoBuild¥}
	Flex -t -s -i -8 ReadAddData.flex > lex.yy.cp

"{ObjDir¥PPC}lex.yy.cp.x" Ä {¥MondoBuild¥} lex.yy.cp
	{PPCCPlus} lex.yy.cp -o {Targ} {Includes} {Sym¥PPC} -w 2
