#include "FlameMaster.h"

#include "../FlameManLibsSource/DasslLib/dassl.h"

#undef DIFFFLAME

#undef V
#define VS
#undef VSH

#undef MOMIC
#define HMOM
#undef FRENKLACH

#undef AGG
#define BLEND

#undef LINDSTEDT

#undef FRAGTEST
#ifdef FRAGTEST
  #undef PHI_0_85
  #define PHI_0_95
  #undef PHI_1_05
  #undef PHI_1_10
  #undef PHI_1_15
#endif

#define SOLVELOG

#ifdef DIFFFLAME
#define NODIFF
#define MAGICSOURCE 1.0
#else
#define NODIFF
#define MAGICSOURCE 1.0
#endif

#undef FLUXBC
#define SMALLSOOT 1.0e-20
#define MAGIC 1.0 //Radiation in only one direction: 1/3

//Global rates of surface growth and oxidation
#undef GLOBAL_SOOT


//**************************************//
//************Initialization************//
//**************************************//

void TSoot::InitTSoot(TInputDataPtr input)
{
  CounterPtr	counter = input->GetCounter();
  ReactionPtr	sootReac = input->GetSootReactions();	
  int 		nSootReactions = counter->sootReactions;
  int		i;

  fSootDASSL = input->fSootDASSL; //Mueller-4/03/08
  fNucleation = input->fNucleation;
  fCondensation = input->fCondensation;
  fCoagulation = input->fCoagulation;
  fSurfaceGrowth = input->fSurfaceGrowth;
  fSurfaceOxidation = input->fSurfaceOxidation;
  fThermoPhoresis = input->fThermoPhoresis;
  fFragmentation = input->fFragmentation;
  fSootDiffusion = input->fSootDiffusion;
  fCoagFact = input->fCoagFact;
  fLewis1 = 1.0;
  fSootRadiation = input->fSootRadiation;
  fSootUpdateProdRates = input->fSootUpdateProdRates;
  fSizeDepDiff = input->fSizeDepDiff;
  fSurfDepCoag = input->fSurfDepCoag;

  InitPAH(input);

  //Initialize soot reactions
  if (!nSootReactions)
  {
    cerr << "error: number of soot reactions is zero" << NEWL;
    exit(2);
  }

  fASoot = NewVector(nSootReactions);
  fNSoot = NewVector(nSootReactions);
  fEOverRSoot = NewVector(nSootReactions);
  
  fSpecNumSoot = new IntVectorPtr[nSootReactions];
  if (!fSpecNumSoot)
    FatalError("memory allocation of TSoot failed");
  fNuSoot = new VectorPtr[nSootReactions];
  if (!fNuSoot)
    FatalError("memory allocation of TSoot failed");

  for (i = 0; i < nSootReactions; ++i)
  {
    fSpecNumSoot[i] = NewIntVector(sootReac[i].numberOfSpecies);
    fNuSoot[i] = NewVector(sootReac[i].numberOfSpecies);
  }

  fSootRateCoeffs = NewVector(nSootReactions);
  fSootReactionRate = NewVector(nSootReactions);

  if (!(fSootLabels = new String[nSootReactions]))
    FatalError("memory allocation of TSoot failed");
  if (!(fPAHSymbol = new char[strlen(input->fPAHSymbol)+1]))
    FatalError("memory allocation of TSoot failed");
  strcpy(fPAHSymbol, input->fPAHSymbol);
  if (!(fSootSymbol = new char[strlen(input->fSootSymbol)+1]))
    FatalError("memory allocation of TSoot failed");
  strcpy(fSootSymbol, input->fSootSymbol);

  FillTSoot(input);

//   CheckSootReactions();
//   CheckPAHReactions();

//   int *	ind = fPAHSpeciesIndex->vec;

//   if ((ind[1] = input->FindSpecies("A4-C18H10")) == -1)
//   {
//     cerr << "error: can't find the molecule A4-C18H10" << NEWL;
//     exit(2);
//   }
//   if ((ind[2] = input->FindSpecies("A4--C18H9")) == -1)
//   {
//     cerr << "error: can't find the molecule A4--C18H9" << NEWL;
//     exit(2);
//   }
//   if ((ind[3] = input->FindSpecies("A4C2H-C20H10")) == -1)
//   {
//     cerr << "error: can't find the molecule A4C2H-C20H10" << NEWL;
//     exit(2);
//   }
//   if ((ind[4] = input->FindSpecies("A4C2H*-C20H9")) == -1)
//   {
//     cerr << "error: can't find the molecule A4C2H*-C20H9" << NEWL;
//     exit(2);
//   }
//   if ((ind[5] = input->FindSpecies("A5--C22H11")) == -1)
//   {
//     cerr << "error: can't find the molecule A5--C22H11" << NEWL;
//     exit(2);
//   }
//   if ((ind[6] = input->FindSpecies("A5-C22H12")) == -1)
//   {
//     cerr << "error: can't find the molecule A5-C22H12" << NEWL;
//     exit(2);
//   }

  if ((f_H = input->FindSpecies("H")) == -1)
  {
    cerr << "error: can't find the molecule H" << NEWL;
    exit(2);
  }
  if ((f_H2 = input->FindSpecies("H2")) == -1)
  {
    cerr << "error: can't find the molecule H2" << NEWL;
    exit(2);
  }
  if ((f_O2 = input->FindSpecies("O2")) == -1)
  {
    cerr << "error: can't find the molecule O2" << NEWL;
    exit(2);
  }
  if ((f_OH = input->FindSpecies("OH")) == -1)
  {
    cerr << "error: can't find the molecule OH" << NEWL;
    exit(2);
  }
  if ((f_CO = input->FindSpecies("CO")) == -1)
  {
    cerr << "error: can't find the molecule f_CO" << NEWL;
    exit(2);
  }
//   if ((f_CH = input->FindSpecies("CH")) == -1)
//   {
//     cerr << "error: can't find the molecule f_CH" << NEWL;
//     exit(2);
//   }
//   if ((f_CHO = input->FindSpecies("CHO")) == -1)
//   {
//     if ((f_CHO = input->FindSpecies("HCO")) == -1)
//     {
//       cerr << "error: can't find the molecule f_CHO or f_HCO" << NEWL;
//       exit(2);
//     }
//   }
  if ((f_H2O = input->FindSpecies("H2O")) == -1)
  {
    cerr << "error: can't find the molecule H2O" << NEWL;
    exit(2);
  }
//   if ((f_C2H = input->FindSpecies("C2H")) == -1)
//   {
//     cerr << "error: can't find the molecule C2H" << NEWL;
//     exit(2);
//   }
  if ((f_C2H2 = input->FindSpecies("C2H2")) == -1)
  {
    cerr << "error: can't find the molecule C2H2" << NEWL;
    exit(2);
  }
//   if ((f_A3R5M = input->FindSpecies("A3R5--C16H9")) == -1)
//   {
//     cerr << "error: can't find the molecule A3R5--C16H9" << NEWL;
//     exit(2);
//   }
//   if ((f_A1 = input->FindSpecies("A1-C6H6")) == -1)
//   {
//     cerr << "error: can't find the molecule A1-C6H6" << NEWL;
//     exit(2);
//   }
//   if ((fFirstPAH = input->FindSpecies("A4-C18H10")) == -1)
//   {
//     cerr << "error: can't find the molecule A4-C18H10" << NEWL;
//     exit(2);
//   }

//   for (i = 1; i < fPAHSpeciesIndex->len ; ++i)
//   {
//     if ((ind[i] - fFirstPAH + 1) != i)
//     {
//       cerr << "error: PAH molecule no. " << i << TAB << (ind[i] - fFirstPAH + 1) << " is not in right order" << NEWL;
//     }
//   }
	
  fOffsetMoments = -1; //Value of fOffsetMoments is set later
}

TSoot::~TSoot(void)
{
  int i;
  double nSootReactions = GetNSootReactions();

  delete fSootSymbol;
  for (i = 0; i < nSootReactions; ++i)
  {
    delete fSootLabels[i];
  }
  delete fSootLabels;
  
  DisposeVector(fSootReactionRate);
  DisposeVector(fSootRateCoeffs);
  for (i = 0; i < GetNSootReactions(); ++i)
  {
    DisposeVector(fNuSoot[i]);
    DisposeIntVector(fSpecNumSoot[i]);
  }
  delete fNuSoot;
  delete fSpecNumSoot;
  
  DisposeVector(fEOverRSoot);
  DisposeVector(fNSoot);
  DisposeVector(fASoot);
  
//   DisposeVector(fPhiPAH);
//   DisposeVector(fFracMomPAH);

//   DisposeVector(fPhi);
//   DisposeVector(fFracMom);
}

void T0DSoot::InitT0DSoot(TInputDataPtr input)
{
  fMoments = NewVector(fNSootMoments);
}

T0DSoot::~T0DSoot(void)
{
//  DisposeVector(fReactionRate);

  DisposeVector(fMoments);
//  DisposeVector(fPAHMoments);
//  DisposeVector(fSumPi);
//  DisposeMatrix(fPij);
}

void T1DSoot::InitT1DSoot(TInputDataPtr input)
{
  int nOfGridPoints = input->fMaxGridPoints;
//   int nPAHReactions = input->GetCounter()->pahReactions;
  int nSpeciesIn = input->GetCounter()->species - input->GetCounter()->steadyStates;
	
//   fReactionRate = NewMatrix(nPAHReactions, nOfGridPoints, kColumnPointers);

//   fPij = NewTensor(nOfGridPoints+2, fNPAHMolecules+1, fNStages, kColumnPointers);
//   fPij->tensor = &fPij->tensor[kNext];
//   fSumPi = NewMatrix(fNPAHMolecules, nOfGridPoints+2, kColumnPointers);
//   fSumPi->mat = &fSumPi->mat[kNext];
//   fPAHMoments = NewMatrix(fNPAHMoments, nOfGridPoints+2, kColumnPointers);
//   fPAHMoments->mat = &fPAHMoments->mat[kNext];

  fMoments = NewMatrix(fNSootMoments, nOfGridPoints+2, kColumnPointers);
  fMoments->mat = &fMoments->mat[kNext];

  fSavedMoments = NewMatrix(fNSootMoments, nOfGridPoints+2, kColumnPointers);
  fSavedMoments->mat = &fSavedMoments->mat[kNext];

  fSootDiff = NewVector(nOfGridPoints+2);
  fSootDiff->vec = &fSootDiff->vec[kNext];

  //Density correction
  fRhoDot = NewVector(nOfGridPoints+2);
  fRhoDot->vec = &fRhoDot->vec[kNext];

  fEnthDot = NewVector(nOfGridPoints+2);
  fEnthDot->vec = &fEnthDot->vec[kNext];

  fDMdx = NewTensor(nOfGridPoints, fNSootMoments, nSpeciesIn+fNSootMoments+1, kColumnPointers);
  
  fSource = NewMatrix(fNSootMoments, nOfGridPoints, kColumnPointers);
	
  fSourceMod = NewVector(fNSootMoments);

  ySaved = NewMatrix(fNSootMoments, nOfGridPoints, kColumnPointers);
  ySavedPrev = NewMatrix(fNSootMoments, nOfGridPoints, kColumnPointers);
  ySavedNext = NewMatrix(fNSootMoments, nOfGridPoints, kColumnPointers);

  if (fSootDASSL)
  {
    //Allocate some well need arrays

    //Solution at previous iteration
    ySol = new double [fNSootMoments];
    yPrime = new double [fNSootMoments];
    src = new double [fNSootMoments];

    //Info array
    info = new int[15];

    info[0]  = 0; //New problem
    info[1]  = 0; //Accuracy as scalars
    info[2]  = 0;
    info[3]  = 0;
    info[4]  = 0;
    info[5]  = 0;
    info[6]  = 0;
    info[7]  = 0; //Do not set initial timestep
    info[8]  = 0;
    info[9]  = 0; //Positivity enforced
    info[10] = 0; //yPrime consistent
    info[11] = 0; //Not used
    info[12] = 0; //Not used
    info[13] = 0; //Not used
    info[14] = 0; //Not used

    //Precision
    //Force relative convergence and not absolute
//     if (input->fSootRelTol<0.0)
//     {
//       cerr << "error: SootRelTol not specified." << NEWL;
//       exit(2);
//     }
//     rtol = input->fSootRelTol;

    rtol = 1.0e-6; // 6
    atol = 0.0e-24;

    //Work arrays
    liw = 20+fNSootMoments;
    lrw = 40+(5+4)*fNSootMoments+fNSootMoments*fNSootMoments;
    iwork = new int[liw];
    rwork = new double[lrw];

    rwork[2] = 1.0e-24;
  }
}

T1DSoot::~T1DSoot(void)
{
  DisposeVector(fSourceMod);
  DisposeMatrix(fSource);
  DisposeTensor(fDMdx);

  fSootDiff->vec = &fSootDiff->vec[kPrev];
  DisposeVector(fSootDiff);

  //Density correction
  fRhoDot->vec = &fRhoDot->vec[kPrev];
  DisposeVector(fRhoDot);

  fEnthDot->vec = &fEnthDot->vec[kPrev];
  DisposeVector(fEnthDot);

  fSavedMoments->mat = &fSavedMoments->mat[kPrev];
  DisposeMatrix(fSavedMoments);
  fMoments->mat = &fMoments->mat[kPrev];
  DisposeMatrix(fMoments);

//   fPAHMoments->mat = &fPAHMoments->mat[kPrev];
//   DisposeMatrix(fPAHMoments);
//   fSumPi->mat = &fSumPi->mat[kPrev];
//   DisposeMatrix(fSumPi);
//   fPij->tensor = &fPij->tensor[kPrev];
//   DisposeTensor(fPij);

//   DisposeMatrix(fReactionRate);

  if (fSootDASSL)
  {
    //DASSL
    delete src;
    delete(ySol);
    delete(yPrime);
    delete(info);
    delete(iwork);
    delete(rwork);
  }
}


//***************************************************************************//
//**********************************Moments**********************************//
//***************************************************************************//

int TSoot::Geti(int n)
{
  switch (n)
  {
    case  0: return 0;
    case  1: return 6;
    case  2: return 0;
    case  3: return 50;
    case  4: return 6;
    case  5: return 0;
    case  6: return 50;
    case  7: return 12;
    case  8: return 6;
    case  9: return 0;
    case 10: return 50;
    case 11: return 18;
    case 12: return 12;
    case 13: return 6;
    case 14: return 0;
    case 15: return 50;
  }

  cerr << "More than sixteen moments not implemented.\n";
  return 0;
}

int TSoot::Getj(int n)
{
  switch (n)
  {
    case  0: return 0;
    case  1: return 0;
    case  2: return 6;
    case  3: return 50;
    case  4: return 6;
    case  5: return 12;
    case  6: return 50;
    case  7: return 6;
    case  8: return 12;
    case  9: return 18;
    case 10: return 50;
    case 11: return 6;
    case 12: return 12;
    case 13: return 18;
    case 14: return 24;
    case 15: return 50;
  }

  cerr << "More than sixteen moments not implemented.\n";
  return 0;
}

int TSoot::Getk(int n)
{
  switch (n)
  {
    case  0: return 0;
    case  1: return 0;
    case  2: return 0;
    case  3: return 6;
    case  4: return 50;
    case  5: return 0;
    case  6: return 0;
    case  7: return 0;
    case  8: return 0;
    case  9: return 0;
    case 10: return 0;
    case 11: return 0;
    case 12: return 0;
    case 13: return 0;
    case 14: return 0;
    case 15: return 0;
  }

  cerr << "More than sixteen moments not implemented.\n";
  return 0;
}


//***************************************************************************//
//*******************************Soot Reactions******************************//
//***************************************************************************//

void TSoot::CheckSootReactions(void)
{
  int	error = -1;
  int	offset = strlen(fSootSymbol);
	
  if (strcmp(fSootLabels[ks1f] + offset, "1f"))
    error = ks1f;
  if (strcmp(fSootLabels[ks1b] + offset, "1b"))
    error = ks1b;
  if (strcmp(fSootLabels[ks3f] + offset, "3f"))
    error = ks3f;
  if (strcmp(fSootLabels[ks3b] + offset, "3b"))
		error = ks3b;
// if (strcmp(fSootLabels[ks111] + offset, "111"))
//   error = ks111;
// if (strcmp(fSootLabels[ks112] + offset, "112"))
//   error = ks112;
// if (strcmp(fSootLabels[ks13f] + offset, "13f"))
//   error = ks13f;
// if (strcmp(fSootLabels[ks13b] + offset, "13b"))
//   error = ks13b;

  if (strcmp(fSootLabels[ks2f] + offset, "2f"))
    error = ks2f;
  if (strcmp(fSootLabels[ks2b] + offset, "2b"))
    error = ks2b;

  if (strcmp(fSootLabels[ks4] + offset, "4"))
    error = ks4;
//    if (strcmp(fSootLabels[ks10b] + offset, "10b"))
//      error = ks10b;

//    if (strcmp(fSootLabels[ks12] + offset, "12"))
//      error = ks12;

  if (strcmp(fSootLabels[ks5] + offset, "5"))
    error = ks5;
  
  if (strcmp(fSootLabels[ks6] + offset, "6"))
    error = ks6;

  if (error > -1)
  {
    cerr << "###error in soot reaction " << fSootLabels[error] << NEWL;
//     cerr << "ks12 =  " << (int)ks12 << NEWL;
//     cerr << "fSootLabels[ks12] + offset =  " << fSootLabels[ks12] + offset << NEWL;
    exit(2);
  }
}

void TSoot::PrintSootReactions(TFlamePtr flame, TSpeciesPtr species)
{
  char **	names = species->GetNames();
  int *	ind;
  FILE *	fp = flame->GetOutfile("sootreacs", flame->kText);

  for (int i = 0; i < GetNSootReactions(); ++i)
  {
    ind = fSpecNumSoot[i]->vec;
    fprintf(fp, "reaction %d:\n", i);
    for (int j = 0; j < fSpecNumSoot[i]->len; ++j)
    {
      fprintf(fp, "\t %g %s", fNuSoot[i]->vec[j], names[ind[j]]);
    }
    fprintf(fp, "\n\nA = %g\nn = %g\nE = %g\n\n\n", fASoot->vec[i], fNSoot->vec[i], fEOverRSoot->vec[i] * RGAS);
  }

  fclose(fp);
}

void TSoot::FillTSoot(TInputDataPtr input)
{
  int i, j;
  DimensionPtr dimension = input->GetDimension();
  int numberOfDimensions = input->GetCounter()->dimensions;
  double unitsConverterA = 0.0;
  double unitsConverterE = 0.0;
  double reacOrderConv = 0.0;
  double tempExpConv = 0.0;
  int orderOfReaction;
  const double	moleToKmole = 1000.0;

  //Set units converter
  //It is assumed that only A is a function of the order of reaction and the temperature exponent
  for (i = 0; i < numberOfDimensions; ++i )
  {
    if (dimension[i].name[0] == 'A')
    {
      if (unitsConverterA)
      {
	cerr << "#error: units converter for 'A' doubly defined or zero" << endl;
	exit(2);
      }
      else
      {
	unitsConverterA = dimension[i].value;
	if (dimension[i].orderOfReaction)
	{
	  reacOrderConv = dimension[i].orderOfReacValue;
	}
	if (dimension[i].tempExponent)
	{
	  tempExpConv = dimension[i].tempExpValue;
	}
      }
    }
    else if (dimension[i].name[0] == 'E')
    {
      if (unitsConverterE)
      {
	cerr << "#error: units converter for 'E' doubly defined or zero" << endl;
	exit(2);
      }
      else
      {
	unitsConverterE = dimension[i].value;
      }
    }
    else
    {
      cerr << "# warning: i don't need units for " << dimension->name << endl;
    }
  }

  //Set values for A, N, E, including the conversion to SI-units and mole to kmole for Soot reactions
  ReactionPtr reaction = input->GetSootReactions();

  Double *	a = fASoot->vec;
  Double *	n = fNSoot->vec;
  Double *	eOverR = fEOverRSoot->vec;

  for (i = 0; i < GetNSootReactions(); ++i)
  {
    //First get label
    fSootLabels[i] = new char[strlen(reaction[i].label)+1];
    if (!fSootLabels[i])
      FatalError("memory allocation of TSoot failed");
    strcpy(fSootLabels[i], reaction[i].label);
    
    //Then determine order of reaction
    orderOfReaction = (reaction[i].withThirdBody) ? 1: 0;
    for (j = 0; j < reaction[i].numberOfSpecies; ++j)
    {
      if (reaction[i].speciesCoeff[j] > 0.0)
      {
	orderOfReaction += (int)reaction[i].speciesCoeff[j];
      }
    }
    if (orderOfReaction > 3 && !reaction[i].partialEquilibrium)
    {
      cerr << "#warning: order of reaction no. " << i << " is " << orderOfReaction << NEWL;
    }
    
    //Then fill arrays
    //First convert A and E to SI-units, then change mole to kmole
    a[i] = unitsConverterA * myPow(reacOrderConv, orderOfReaction) * myPow(tempExpConv, reaction[i].n) * reaction[i].a * myPow(moleToKmole, orderOfReaction-1);
    n[i] = reaction[i].n;
    eOverR[i] = unitsConverterE * reaction[i].e / RGAS * 1000.0;
    
    int	speciesCounter = 0;
    for (j = 0; j < fNuSoot[i]->len; ++j)
    {
      if (reaction[i].speciesCoeff[j] > 0.0)
      {
	fSpecNumSoot[i]->vec[speciesCounter] = reaction[i].speciesNumber[j];
	fNuSoot[i]->vec[speciesCounter] = reaction[i].speciesCoeff[j];
	++speciesCounter;
      }
    }
    for (j = 0; j < fNuSoot[i]->len; ++ j)
    {
      if (reaction[i].speciesCoeff[j] < 0.0)
      {
	fSpecNumSoot[i]->vec[speciesCounter] = reaction[i].speciesNumber[j];
	fNuSoot[i]->vec[speciesCounter] = reaction[i].speciesCoeff[j];
	++speciesCounter;
      }
    }
  }
}

void TSoot::ComputeSootRateCoeffs(Double *k, Double temp, TReactionPtr reaction)
{
  double * a = fASoot->vec;
  double * n = fNSoot->vec;
  double * eOverR = fEOverRSoot->vec;
	
  for (int i = 0; i < GetNSootReactions(); ++i)
  {
    reaction->ComputeRateCoefficient(temp, k[i], a[i], n[i], eOverR[i]);
  }

  static const double GammaOH = 0.13; //0.13;
  k[ks6] = 8.94 * GammaOH * sqrt(temp) * AVOGADRO;
}


//**************************************//
//*************Source Terms*************//
//**************************************//

//----------Nucleation----------//
double TSoot::NucleationSource(double i, double temp)
{
  return 0.5 * GetBetaNucl(temp, dimer_nbrC2) * pow(2.0 * dimer_nbrC2, i) * dimer_conc * dimer_conc;
}

double TSoot::NucleationSource(double i, double j, double temp)
{
  return 0.5 * GetBetaNucl(temp, dimer_nbrC2) * pow(2.0 * dimer_nbrC2, i+(2.0/3.0)*j) * dimer_conc * dimer_conc;
}

double TSoot::NucleationSource(double i, double j, double k, double temp)
{
  return 0.5 * GetBetaNucl(temp, dimer_nbrC2) * pow(2.0 * dimer_nbrC2, i+(2.0/3.0)*j) * pow(2.0 * dimer_nbrH, k) * dimer_conc * dimer_conc;
}

//----------Coagulation----------//
double TSoot::CoagulationSource(double i, double temp, double * moments)
{
  //Only implemented for free-molecular regime

  //Coagulation does not change volume
  if (i < 1.01 && i > 0.99)
    return 0.0;

#ifdef MOMIC
  //Simplified forms for 0,1,2,3,4,5
  //Less error than general form below
  if (i == 0)
    return -0.5 * GetPsi(0.0, 0.0, moments, temp);
  if (i < 2.01 && i > 1.99)
    return GetPsi(1.0, 1.0, moments, temp);
  if (i < 3.01 && i > 2.99)
    return 3.0 * GetPsi(2.0, 1.0, moments, temp);
  if (i < 4.01 && i > 3.99)
    return 4.0 * GetPsi(3.0, 1.0, moments, temp) + 3.0 * GetPsi(2.0, 2.0, moments, temp);
  if (i < 5.01 && i > 4.99)
    return 5.0 * GetPsi(4.0, 1.0, moments, temp) + 10.0 * GetPsi(3.0, 2.0, moments, temp);

  //General form for all other moments
  return 0.5 * GetPhi(i, moments, temp) - GetPsi(i, 0.0, moments, temp);
#endif

#ifdef FRENKLACH
  //Simplified forms for 0,1,2,3,4,5
  //Less error than general form below
  if (i == 0)
    return -0.5 * GetPsi(0.0, 0.0, moments, temp);
  if (i < 2.01 && i > 1.99)
    return GetPsi(1.0, 1.0, moments, temp);
  if (i < 3.01 && i > 2.99)
    return 3.0 * GetPsi(2.0, 1.0, moments, temp);
  if (i < 4.01 && i > 3.99)
    return 4.0 * GetPsi(3.0, 1.0, moments, temp) + 3.0 * GetPsi(2.0, 2.0, moments, temp);
  if (i < 5.01 && i > 4.99)
    return 5.0 * GetPsi(4.0, 1.0, moments, temp) + 10.0 * GetPsi(3.0, 2.0, moments, temp);

  //General form for all other moments
  return 0.5 * GetPhi(i, moments, temp) - GetPsi(i, 0.0, moments, temp);
#endif

#ifdef HMOM
  double ss;
  double sl;
  double ll;

  double C = GetC(temp);
  double V0 = 2.0 * dimer_nbrC2;

  //Contribution of source term due to collision of small particles with small particles
  ss = C * pow(2.0, 2.5) * (pow(2.0, i-1.0) - 1.0) * pow(V0, i+(1.0/6.0)) * moments[fNSootMoments-1] * moments[fNSootMoments-1];

  //Contribution of source term due to collision of small particles with large particles
  if (i == 0)
    sl = - GetPsiSL(0.0, 0.0, moments, temp);
  else if (i < 2.01 && i > 1.99)
    sl = 2.0 * GetPsiSL(1.0, 1.0, moments, temp);
  else if (i < 3.01 && i > 2.99)
    sl = 3.0 * GetPsiSL(2.0, 1.0, moments, temp) + 3.0 * GetPsiSL(1.0, 2.0, moments, temp);
  else if (i < 4.01 && i > 3.99)
    sl = 4.0 * GetPsiSL(3.0, 1.0, moments, temp) + 6.0 * GetPsiSL(2.0, 2.0, moments, temp) + 4.0 * GetPsiSL(1.0, 3.0, moments, temp);
  else if (i < 5.01 && i > 4.99)
    sl = 5.0 * GetPsiSL(4.0, 1.0, moments, temp) + 10.0 * GetPsiSL(3.0, 2.0, moments, temp) + 10.0 * GetPsiSL(2.0, 3.0, moments, temp) + 5.0 * GetPsiSL(1.0, 4.0, moments, temp);
  //else
  //  sl = GetPhiSL(i, moments, temp) - GetPsiSL(i, 0.0, moments, temp) - GetPsiSL(0.0, i, moments, temp);

  //Contribution of source term due to collision of large particles with large particles
  if (i == 0)
    ll = -0.5 * GetPsiLL(0.0, 0.0, moments, temp);
  else if (i < 2.01 && i > 1.99)
    ll = GetPsiLL(1.0, 1.0, moments, temp);
  else if (i < 3.01 && i > 2.99)
    ll = 3.0 * GetPsiLL(2.0, 1.0, moments, temp);
  else if (i < 4.01 && i > 3.99)
    ll = 4.0 * GetPsiLL(3.0, 1.0, moments, temp) + 3.0 * GetPsiLL(2.0, 2.0, moments, temp);
  else if (i < 5.01 && i > 4.99)
    ll = 5.0 * GetPsiLL(4.0, 1.0, moments, temp) + 10.0 * GetPsiLL(3.0, 2.0, moments, temp);
  //else
  //  ll = 0.5 * GetPhiLL(i, moments, temp) - GetPsiLL(i, 0.0, moments, temp);

  return ss + sl + ll;
#endif
}

Double TSoot::CoagulationSource(Double i, Double j, Double temp, Double * moments)
{
  //Implemented for both free-molecular and continuum regimes

#ifdef MOMIC
  //Pure Aggregation
  return CoagCont(i, j, moments, temp);
#endif

 
#ifdef HMOM
  double ss;
  double ssfm;
  double sscont;

  double sl;
  double ll;

  double V0 = 2.0 * dimer_nbrC2;

  double C_fm = GetC(temp);

  double C_cont = 8.0 * RGAS/AVOGADRO * temp / (3.0*mu) * AVOGADRO;
  double lambda = 3.0*mu/rho*sqrt(PI*Wmix/(8.0*RGAS*temp)) / pow(6.0*fMolarMassSoot/(PI*fSootDensity*AVOGADRO),1.0/3.0);

  double FitC = 2.0/3.0;
  double FitExp = -0.2043;

 #ifdef AGG
  if (i < 1.01 && i > 0.99 && j == 0)
    return 0.0;
  if (i == 0 && j < 1.01 && j > 0.99)
    return 0.0;

  //Contribution of source term due to collision of small particles with small particles
  ssfm = C_fm * pow(2.0, 2.5) * (pow(2.0, i+j-1.0) - 1.0) * pow(V0, i+(2.0/3.0)*j+(1.0/6.0)) * moments[fNSootMoments-1] * moments[fNSootMoments-1];
  sscont = 4.0 * C_cont * (1 + 1.257*lambda*pow(V0,-1.0/3.0)) * (pow(2.0, i+j-1.0) - 1.0) * pow(V0, i+(2.0/3.0)*j) * moments[fNSootMoments-1] * moments[fNSootMoments-1];

  ss = (ssfm*sscont)/(ssfm+sscont);

  //Contribution of source term due to collision of small particles with large particles
  sl = CoagContSL(i,j,moments,temp);

  //Contribution of source term due to collision of large particles with large particles
  ll = CoagContLL(i,j,moments,temp);
 #endif

 #ifdef BLEND
  if (i < 1.01 && i > 0.99 && j == 0)
    return 0.0;

  //Contribution of source term due to collision of small particles with small particles
  //Small particles colliding with small particles is assumbed to be pure coalescence
  ssfm = C_fm * pow(2.0, 2.5) * (pow(2.0, i+(2.0/3.0)*j-1.0) - 1.0) * pow(V0, i+(2.0/3.0)*j+(1.0/6.0)) * moments[fNSootMoments-1] * moments[fNSootMoments-1];
  sscont = 4.0 * C_cont * (1 + 1.257*lambda*pow(V0,-1.0/3.0)) * (pow(2.0, i+(2.0/3.0)*j-1.0) - 1.0) * pow(V0, i+(2.0/3.0)*j) * moments[fNSootMoments-1] * moments[fNSootMoments-1];

  if (ssfm*sscont == 0.0)
    ss = 0.0;
  else
    ss = (ssfm*sscont)/(ssfm+sscont);

  //Contribution of source term due to collision of small particles with large particles
  //Small particles are assumed to 'splash' onto the large particles
  sl = CoagContSL(i,j,moments,temp);

  //Contribution of source term due to collision of large particles with large particles
  //Large particles colliding with large particles is assumed to be pure aggregation
  ll = CoagContLL(i,j,moments,temp);
 #endif
  
  return ss + sl + ll;
#endif
}

double TSoot::CoagCont(double i, double j, double * moments, double temp)
{
  if ((i < 1.01 && i > 0.99 && j == 0) || (i == 0 && j < 1.01 && j > 0.99))
    return 0.0;

  double fm;
  double cont;

  //Free Molecular
  if (i == 0 && j==0)
    fm = -0.5 * GetPsi(0.0, 0.0, 0.0, 0.0, moments, temp);
  else if (i < 2.01 && i > 1.99 && j == 0)
    fm =        GetPsi(1.0, 0.0, 1.0, 0.0, moments, temp);
  else if (i < 1.01 && i > 0.99 && j < 1.01 && j > 0.99)
    fm =        GetPsi(1.0, 0.0, 0.0, 1.0, moments, temp);
  else if (i == 0 && j < 2.01 && j > 1.99)
    fm =        GetPsi(0.0, 1.0, 0.0, 1.0, moments, temp);
  else if (i < 3.01 && i > 2.99 && j == 0)
    fm =  3.0 * GetPsi(2.0, 0.0, 1.0, 0.0, moments, temp);
  else if (i < 2.01 && i > 1.99 && j < 1.01 && j > 0.99)
    fm =        GetPsi(2.0, 0.0, 0.0, 1.0, moments, temp) +
          2.0 * GetPsi(1.0, 1.0, 1.0, 0.0, moments, temp);
  else if (i < 1.01 && i > 0.99 && j < 2.01 && j > 1.99)
    fm =        GetPsi(0.0, 2.0, 1.0, 0.0, moments, temp) +
          2.0 * GetPsi(1.0, 1.0, 0.0, 1.0, moments, temp);
  else if (i == 0 && j < 3.01 && j > 2.99)
    fm =  3.0 * GetPsi(0.0, 2.0, 0.0, 1.0, moments, temp);
  else if (i < 4.01 && i > 3.99 && j == 0)
    fm =  4.0 * GetPsi(3.0, 0.0, 1.0, 0.0, moments, temp) +
          3.0 * GetPsi(2.0, 0.0, 2.0, 0.0, moments, temp);
  else if (i < 3.01 && i > 2.99 && j < 1.01 && j > 0.99)
    fm =  3.0 * GetPsi(2.0, 1.0, 1.0, 0.0, moments, temp) +
          3.0 * GetPsi(1.0, 1.0, 2.0, 0.0, moments, temp) +
                GetPsi(0.0, 1.0, 3.0, 0.0, moments, temp);
  else if (i < 2.01 && i > 1.99 && j < 2.01 && j > 1.99)
    fm =  2.0 * GetPsi(2.0, 1.0, 0.0, 1.0, moments, temp) +
                GetPsi(2.0, 0.0, 0.0, 2.0, moments, temp) +
          2.0 * GetPsi(1.0, 2.0, 1.0, 0.0, moments, temp) +
          2.0 * GetPsi(1.0, 1.0, 1.0, 1.0, moments, temp);
  else if (i < 1.01 && i > 0.99 && j < 3.01 && j > 2.99)
    fm =  3.0 * GetPsi(1.0, 2.0, 0.0, 1.0, moments, temp) +
          3.0 * GetPsi(1.0, 1.0, 0.0, 2.0, moments, temp) +
                GetPsi(1.0, 0.0, 0.0, 3.0, moments, temp);
  else if (i == 0 && j < 4.01 && j > 3.99)
    fm =  4.0 * GetPsi(0.0, 3.0, 0.0, 1.0, moments, temp) +
          3.0 * GetPsi(0.0, 2.0, 0.0, 2.0, moments, temp);
  else
    fm = 0.5 * GetPhi(i, j, moments, temp) - GetPsi(i, j, 0.0, 0.0, moments, temp);

  //Continuum
  double Df = 1.8;
  double av = 1.0 - (2.0 / Df);
  double as = (3.0 / Df) - 1.0;

  double lambda = 3.0*mu/rho*sqrt(PI*Wmix/(8.0*RGAS*temp)) / pow(6.0*fMolarMassSoot/(PI*fSootDensity*AVOGADRO),1.0/3.0);
  double C = 8.0 * RGAS/AVOGADRO * temp / (3.0*mu) * AVOGADRO;

  if (i == 0 && j == 0)
    cont = -0.5 * C * (2.000*FracMom(0,0,moments)*FracMom(0,0,moments) + 
		     FracMom(av,as,moments)*FracMom(-av,-as,moments) +
		     FracMom(-av,-as,moments)*FracMom(av,as,moments) +
		     1.257*lambda*FracMom(-av,-as,moments)*FracMom(0,0,moments) +
		     1.257*lambda*FracMom(0,0,moments)*FracMom(-av,-as,moments) +
		     1.257*lambda*FracMom(-2.0*av,-2.0*as,moments)*FracMom(av,as,moments) +
		       1.257*lambda*FracMom(av,as,moments)*FracMom(-2.0*av,-2.0*as,moments));
  else if (i < 2.01 && i > 1.99 && j == 0)
    cont = C * (2.000*FracMom(1.0,0,moments)*FracMom(1.0,0,moments) + 
	      FracMom(av+1.0,as,moments)*FracMom(-av+1.0,-as,moments) +
	      FracMom(-av+1.0,-as,moments)*FracMom(av+1.0,as,moments) +
	      1.257*lambda*FracMom(-av+1.0,-as,moments)*FracMom(1.0,0,moments) +
	      1.257*lambda*FracMom(1.0,0,moments)*FracMom(-av+1.0,-as,moments) +
	      1.257*lambda*FracMom(-2.0*av+1.0,-2.0*as,moments)*FracMom(av+1.0,as,moments) +
		1.257*lambda*FracMom(av+1.0,as,moments)*FracMom(-2.0*av+1.0,-2.0*as,moments));
  else if (i < 1.01 && i > 0.99 && j < 1.01 && j > 0.99)
    cont = C * (2.000*FracMom(1.0,0,moments)*FracMom(0,1.0,moments) + 
	      FracMom(av+1.0,as,moments)*FracMom(-av,-as+1.0,moments) +
	      FracMom(-av+1.0,-as,moments)*FracMom(av,as+1.0,moments) +
	      1.257*lambda*FracMom(-av+1.0,-as,moments)*FracMom(0,1.0,moments) +
	      1.257*lambda*FracMom(1.0,0,moments)*FracMom(-av,-as+1.0,moments) +
	      1.257*lambda*FracMom(-2.0*av+1.0,-2.0*as,moments)*FracMom(av,as+1.0,moments) +
		1.257*lambda*FracMom(av+1.0,as,moments)*FracMom(-2.0*av,-2.0*as+1.0,moments));
  else if (i == 0 && j < 2.01 && j > 1.99)
    cont = C * (2.000*FracMom(0,1.0,moments)*FracMom(0,1.0,moments) + 
	      FracMom(av,as+1.0,moments)*FracMom(-av,-as+1.0,moments) +
	      FracMom(-av,-as+1.0,moments)*FracMom(av,as+1.0,moments) +
	      1.257*lambda*FracMom(-av,-as+1.0,moments)*FracMom(0,1.0,moments) +
	      1.257*lambda*FracMom(0,1.0,moments)*FracMom(-av,-as+1.0,moments) +
	      1.257*lambda*FracMom(-2.0*av,-2.0*as+1.0,moments)*FracMom(av,as+1.0,moments) +
		1.257*lambda*FracMom(av,as+1.0,moments)*FracMom(-2.0*av,-2.0*as+1.0,moments));
  if (i < 3.01 && i > 2.99 && j == 0)
    cont = 3.0 * C * (2.000*FracMom(2.0,0,moments)*FracMom(1.0,0,moments) + 
		     FracMom(av+2.0,as,moments)*FracMom(-av+1.0,-as,moments) +
		     FracMom(-av+2.0,-as,moments)*FracMom(av+1.0,as,moments) +
		     1.257*lambda*FracMom(-av+2.0,-as,moments)*FracMom(1.0,0,moments) +
		     1.257*lambda*FracMom(2.0,0,moments)*FracMom(-av+1.0,-as,moments) +
		     1.257*lambda*FracMom(-2.0*av+2.0,-2.0*as,moments)*FracMom(av+1.0,as,moments) +
		       1.257*lambda*FracMom(av+2.0,as,moments)*FracMom(-2.0*av+1.0,-2.0*as,moments));
  if (i < 2.01 && i > 1.99 && j < 1.01 && j > 0.99)
    cont = C * (
		(2.000*FracMom(2.0,0,moments)*FracMom(0,1.0,moments) + 
		     FracMom(av+2.0,as,moments)*FracMom(-av,-as+1.0,moments) +
		     FracMom(-av+2.0,-as,moments)*FracMom(av,as+1.0,moments) +
		     1.257*lambda*FracMom(-av+2.0,-as,moments)*FracMom(0,1.0,moments) +
		     1.257*lambda*FracMom(2.0,0,moments)*FracMom(-av,-as+1.0,moments) +
		     1.257*lambda*FracMom(-2.0*av+2.0,-2.0*as,moments)*FracMom(av,as+1.0,moments) +
		       1.257*lambda*FracMom(av+2.0,as,moments)*FracMom(-2.0*av,-2.0*as+1.0,moments)
		 ) + 2.0 * (
		 2.000*FracMom(1.0,1.0,moments)*FracMom(1.0,0,moments) + 
		     FracMom(av+1.0,as+1.0,moments)*FracMom(-av+1.0,-as,moments) +
		     FracMom(-av+1.0,-as+1.0,moments)*FracMom(av+1.0,as,moments) +
		     1.257*lambda*FracMom(-av+1.0,-as+1.0,moments)*FracMom(1.0,0,moments) +
		     1.257*lambda*FracMom(1.0,1.0,moments)*FracMom(-av+1.0,-as,moments) +
		     1.257*lambda*FracMom(-2.0*av+1.0,-2.0*as+1.0,moments)*FracMom(av+1.0,as,moments) +
		       1.257*lambda*FracMom(av+1.0,as+1.0,moments)*FracMom(-2.0*av+1.0,-2.0*as,moments)
		));
  if (i < 1.01 && i > 0.99 && j < 2.01 && j > 1.99)
    cont = C * (
		(2.000*FracMom(0,2.0,moments)*FracMom(1.0,0,moments) + 
		     FracMom(av,as+2.0,moments)*FracMom(-av+1.0,-as,moments) +
		     FracMom(-av,-as+2.0,moments)*FracMom(av+1.0,as,moments) +
		     1.257*lambda*FracMom(-av,-as+2.0,moments)*FracMom(1.0,0,moments) +
		     1.257*lambda*FracMom(0,2.0,moments)*FracMom(-av+1.0,-as,moments) +
		     1.257*lambda*FracMom(-2.0*av,-2.0*as+2.0,moments)*FracMom(av+1.0,as,moments) +
		       1.257*lambda*FracMom(av,as+2.0,moments)*FracMom(-2.0*av+1.0,-2.0*as,moments)
		 ) + 2.0 * (
		 2.000*FracMom(1.0,1.0,moments)*FracMom(0,1.0,moments) + 
		     FracMom(av+1.0,as+1.0,moments)*FracMom(-av,-as+1.0,moments) +
		     FracMom(-av+1.0,-as+1.0,moments)*FracMom(av,as+1.0,moments) +
		     1.257*lambda*FracMom(-av+1.0,-as+1.0,moments)*FracMom(0,1.0,moments) +
		     1.257*lambda*FracMom(1.0,1.0,moments)*FracMom(-av,-as+1.0,moments) +
		     1.257*lambda*FracMom(-2.0*av+1.0,-2.0*as+1.0,moments)*FracMom(av,as+1.0,moments) +
		       1.257*lambda*FracMom(av+1.0,as+1.0,moments)*FracMom(-2.0*av,-2.0*as+1.0,moments)
		));
  if (i == 0 && j < 3.01 && j > 2.99)
    cont = 3.0 * C * (2.000*FracMom(0,2.0,moments)*FracMom(0,1.0,moments) + 
		     FracMom(av,as+2.0,moments)*FracMom(-av,-as+1.0,moments) +
		     FracMom(-av,-as+2.0,moments)*FracMom(av,as+1.0,moments) +
		      1.257*lambda*FracMom(-av,-as+2.0,moments)*FracMom(0,1.0,moments) +
		     1.257*lambda*FracMom(0,2.0,moments)*FracMom(-av,-as+1.0,moments) +
		     1.257*lambda*FracMom(-2.0*av,-2.0*as+2.0,moments)*FracMom(av,as+1.0,moments) +
		       1.257*lambda*FracMom(av,as+2.0,moments)*FracMom(-2.0*av,-2.0*as+1.0,moments));

  return cont*fm/(cont+fm);
}

double TSoot::CoagContLL(double i, double j, double * moments, double temp)
{
  if ((i < 1.01 && i > 0.99 && j == 0) || (i == 0 && j < 1.01 && j > 0.99))
    return 0.0;

  double fm;
  double cont;

  //Free Molecular
  if (i == 0 && j == 0)
    fm = -0.5 * GetPsiLL(0.0, 0.0, 0.0, 0.0, moments, temp);
  else if (i < 2.01 && i > 1.99 && j == 0)
    fm =        GetPsiLL(1.0, 0.0, 1.0, 0.0, moments, temp);
  else if (i < 1.01 && i > 0.99 && j < 1.01 && j > 0.99)
    fm =        GetPsiLL(1.0, 0.0, 0.0, 1.0, moments, temp);
  else if (i == 0 && j < 2.01 && j > 1.99)
    fm =        GetPsiLL(0.0, 1.0, 0.0, 1.0, moments, temp);
  else if (i < 3.01 && i > 2.99 && j == 0)
    fm =  3.0 * GetPsiLL(2.0, 0.0, 1.0, 0.0, moments, temp);
  else if (i < 2.01 && i > 1.99 && j < 1.01 && j > 0.99)
    fm =        GetPsiLL(2.0, 0.0, 0.0, 1.0, moments, temp) +
          2.0 * GetPsiLL(1.0, 1.0, 1.0, 0.0, moments, temp);
  else if (i < 1.01 && i > 0.99 && j < 2.01 && j > 1.99)
    fm =        GetPsiLL(0.0, 2.0, 1.0, 0.0, moments, temp) +
          2.0 * GetPsiLL(1.0, 1.0, 0.0, 1.0, moments, temp);
  else if (i == 0 && j < 3.01 && j > 2.99)
    fm =  3.0 * GetPsiLL(0.0, 2.0, 0.0, 1.0, moments, temp);
  else if (i < 4.01 && i > 3.99 && j == 0)
    fm =  4.0 * GetPsiLL(3.0, 0.0, 1.0, 0.0, moments, temp) +
          3.0 * GetPsiLL(2.0, 0.0, 2.0, 0.0, moments, temp);
  else if (i < 3.01 && i > 2.99 && j < 1.01 && j > 0.99)
    fm =  3.0 * GetPsiLL(2.0, 1.0, 1.0, 0.0, moments, temp) +
          3.0 * GetPsiLL(1.0, 1.0, 2.0, 0.0, moments, temp) +
                GetPsiLL(0.0, 1.0, 3.0, 0.0, moments, temp);
  else if (i < 2.01 && i > 1.99 && j < 2.01 && j > 1.99)
    fm =  2.0 * GetPsiLL(2.0, 1.0, 0.0, 1.0, moments, temp) +
                GetPsiLL(2.0, 0.0, 0.0, 2.0, moments, temp) +
          2.0 * GetPsiLL(1.0, 2.0, 1.0, 0.0, moments, temp) +
          2.0 * GetPsiLL(1.0, 1.0, 1.0, 1.0, moments, temp);
  else if (i < 1.01 && i > 0.99 && j < 3.01 && j > 2.99)
    fm =  3.0 * GetPsiLL(1.0, 2.0, 0.0, 1.0, moments, temp) +
          3.0 * GetPsiLL(1.0, 1.0, 0.0, 2.0, moments, temp) +
                GetPsiLL(1.0, 0.0, 0.0, 3.0, moments, temp);
  else if (i == 0 && j < 4.01 && j > 3.99)
    fm =  4.0 * GetPsiLL(0.0, 3.0, 0.0, 1.0, moments, temp) +
          3.0 * GetPsiLL(0.0, 2.0, 0.0, 2.0, moments, temp);
  else
    fm = 0.5 * GetPhiLL(i, j, moments, temp) - GetPsiLL(i, j, 0.0, 0.0, moments, temp);

  //Continuum
  double Df = 1.8;
  double av = 1.0 - (2.0 / Df);
  double as = (3.0 / Df) - 1.0;

  double lambda = 3.0*mu/rho*sqrt(PI*Wmix/(8.0*RGAS*temp)) / pow(6.0*fMolarMassSoot/(PI*fSootDensity*AVOGADRO),1.0/3.0);
  double C = 8.0 * RGAS/AVOGADRO * temp / (3.0*mu) * AVOGADRO;

  if (i == 0 && j == 0)
    cont = -0.5 * C * (2.000*       FracMomLarge( 0     , 0     ,moments)*FracMomLarge( 0     , 0     ,moments) + 
		                    FracMomLarge(     av,     as,moments)*FracMomLarge(-    av,-    as,moments) +
		                    FracMomLarge(-    av,-    as,moments)*FracMomLarge(     av,     as,moments) +
		       1.257*lambda*FracMomLarge(-    av,-    as,moments)*FracMomLarge( 0     , 0     ,moments) +
		       1.257*lambda*FracMomLarge( 0     , 0     ,moments)*FracMomLarge(-    av,-    as,moments) +
		       1.257*lambda*FracMomLarge(-2.0*av,-2.0*as,moments)*FracMomLarge(     av,     as,moments) +
		       1.257*lambda*FracMomLarge(     av,     as,moments)*FracMomLarge(-2.0*av,-2.0*as,moments));
  else if (i < 2.01 && i > 1.99 && j == 0)
    cont = C * (2.000*       FracMomLarge(        1.0, 0     ,moments)*FracMomLarge(        1.0, 0     ,moments) + 
		             FracMomLarge(     av+1.0,     as,moments)*FracMomLarge(-    av+1.0,-    as,moments) +
	                     FracMomLarge(-    av+1.0,-    as,moments)*FracMomLarge(     av+1.0,     as,moments) +
	        1.257*lambda*FracMomLarge(-    av+1.0,-    as,moments)*FracMomLarge(        1.0, 0     ,moments) +
		1.257*lambda*FracMomLarge(        1.0, 0     ,moments)*FracMomLarge(-    av+1.0,-    as,moments) +
		1.257*lambda*FracMomLarge(-2.0*av+1.0,-2.0*as,moments)*FracMomLarge(     av+1.0,     as,moments) +
		1.257*lambda*FracMomLarge(     av+1.0,     as,moments)*FracMomLarge(-2.0*av+1.0,-2.0*as,moments));
  else if (i < 1.01 && i > 0.99 && j < 1.01 && j > 0.99)
    cont = C * (2.000*       FracMomLarge(        1.0, 0     ,moments)*FracMomLarge( 0     ,        1.0,moments) + 
		             FracMomLarge(     av+1.0,     as,moments)*FracMomLarge(-    av,-    as+1.0,moments) +
		             FracMomLarge(-    av+1.0,-    as,moments)*FracMomLarge(     av,     as+1.0,moments) +
		1.257*lambda*FracMomLarge(-    av+1.0,-    as,moments)*FracMomLarge( 0     ,        1.0,moments) +
		1.257*lambda*FracMomLarge(        1.0, 0     ,moments)*FracMomLarge(-    av,-    as+1.0,moments) +
		1.257*lambda*FracMomLarge(-2.0*av+1.0,-2.0*as,moments)*FracMomLarge(     av,     as+1.0,moments) +
		1.257*lambda*FracMomLarge(     av+1.0,     as,moments)*FracMomLarge(-2.0*av,-2.0*as+1.0,moments));
  else if (i == 0 && j < 2.01 && j > 1.99)
    cont = C * (2.000*       FracMomLarge( 0     ,        1.0,moments)*FracMomLarge( 0     ,        1.0,moments) + 
	                     FracMomLarge(     av,     as+1.0,moments)*FracMomLarge(-    av,-    as+1.0,moments) +
	                     FracMomLarge(-    av,-    as+1.0,moments)*FracMomLarge(     av,     as+1.0,moments) +
		1.257*lambda*FracMomLarge(-    av,-    as+1.0,moments)*FracMomLarge( 0     ,        1.0,moments) +
		1.257*lambda*FracMomLarge( 0     ,        1.0,moments)*FracMomLarge(-    av,-    as+1.0,moments) +
		1.257*lambda*FracMomLarge(-2.0*av,-2.0*as+1.0,moments)*FracMomLarge(     av,     as+1.0,moments) +
		1.257*lambda*FracMomLarge(     av,     as+1.0,moments)*FracMomLarge(-2.0*av,-2.0*as+1.0,moments));
  else if (i < 3.01 && i > 2.99 && j == 0)
    cont = 3.0 * C * (2.000*FracMomLarge(2.0,0,moments)*FracMomLarge(1.0,0,moments) + 
		     FracMomLarge(av+2.0,as,moments)*FracMomLarge(-av+1.0,-as,moments) +
		     FracMomLarge(-av+2.0,-as,moments)*FracMomLarge(av+1.0,as,moments) +
		     1.257*lambda*FracMomLarge(-av+2.0,-as,moments)*FracMomLarge(1.0,0,moments) +
		     1.257*lambda*FracMomLarge(2.0,0,moments)*FracMomLarge(-av+1.0,-as,moments) +
		     1.257*lambda*FracMomLarge(-2.0*av+2.0,-2.0*as,moments)*FracMomLarge(av+1.0,as,moments) +
		       1.257*lambda*FracMomLarge(av+2.0,as,moments)*FracMomLarge(-2.0*av+1.0,-2.0*as,moments));
  else if (i < 2.01 && i > 1.99 && j < 1.01 && j > 0.99)
    cont = C * (
		(2.000*FracMomLarge(2.0,0,moments)*FracMomLarge(0,1.0,moments) + 
		     FracMomLarge(av+2.0,as,moments)*FracMomLarge(-av,-as+1.0,moments) +
		     FracMomLarge(-av+2.0,-as,moments)*FracMomLarge(av,as+1.0,moments) +
		     1.257*lambda*FracMomLarge(-av+2.0,-as,moments)*FracMomLarge(0,1.0,moments) +
		     1.257*lambda*FracMomLarge(2.0,0,moments)*FracMomLarge(-av,-as+1.0,moments) +
		     1.257*lambda*FracMomLarge(-2.0*av+2.0,-2.0*as,moments)*FracMomLarge(av,as+1.0,moments) +
		       1.257*lambda*FracMomLarge(av+2.0,as,moments)*FracMomLarge(-2.0*av,-2.0*as+1.0,moments)
		 ) + 2.0 * (
		 2.000*FracMomLarge(1.0,1.0,moments)*FracMomLarge(1.0,0,moments) + 
		     FracMomLarge(av+1.0,as+1.0,moments)*FracMomLarge(-av+1.0,-as,moments) +
		     FracMomLarge(-av+1.0,-as+1.0,moments)*FracMomLarge(av+1.0,as,moments) +
		     1.257*lambda*FracMomLarge(-av+1.0,-as+1.0,moments)*FracMomLarge(1.0,0,moments) +
		     1.257*lambda*FracMomLarge(1.0,1.0,moments)*FracMomLarge(-av+1.0,-as,moments) +
		     1.257*lambda*FracMomLarge(-2.0*av+1.0,-2.0*as+1.0,moments)*FracMomLarge(av+1.0,as,moments) +
		       1.257*lambda*FracMomLarge(av+1.0,as+1.0,moments)*FracMomLarge(-2.0*av+1.0,-2.0*as,moments)
		));
  else if (i < 1.01 && i > 0.99 && j < 2.01 && j > 1.99)
    cont = C * (
		(2.000*FracMomLarge(0,2.0,moments)*FracMomLarge(1.0,0,moments) + 
		     FracMomLarge(av,as+2.0,moments)*FracMomLarge(-av+1.0,-as,moments) +
		     FracMomLarge(-av,-as+2.0,moments)*FracMomLarge(av+1.0,as,moments) +
		     1.257*lambda*FracMomLarge(-av,-as+2.0,moments)*FracMomLarge(1.0,0,moments) +
		     1.257*lambda*FracMomLarge(0,2.0,moments)*FracMomLarge(-av+1.0,-as,moments) +
		     1.257*lambda*FracMomLarge(-2.0*av,-2.0*as+2.0,moments)*FracMomLarge(av+1.0,as,moments) +
		       1.257*lambda*FracMomLarge(av,as+2.0,moments)*FracMomLarge(-2.0*av+1.0,-2.0*as,moments)
		 ) + 2.0 * (
		 2.000*FracMomLarge(1.0,1.0,moments)*FracMomLarge(0,1.0,moments) + 
		     FracMomLarge(av+1.0,as+1.0,moments)*FracMomLarge(-av,-as+1.0,moments) +
		     FracMomLarge(-av+1.0,-as+1.0,moments)*FracMomLarge(av,as+1.0,moments) +
		     1.257*lambda*FracMomLarge(-av+1.0,-as+1.0,moments)*FracMomLarge(0,1.0,moments) +
		     1.257*lambda*FracMomLarge(1.0,1.0,moments)*FracMomLarge(-av,-as+1.0,moments) +
		     1.257*lambda*FracMomLarge(-2.0*av+1.0,-2.0*as+1.0,moments)*FracMomLarge(av,as+1.0,moments) +
		       1.257*lambda*FracMomLarge(av+1.0,as+1.0,moments)*FracMomLarge(-2.0*av,-2.0*as+1.0,moments)
		));
  else if (i == 0 && j < 3.01 && j > 2.99)
    cont = 3.0 * C * (2.000*FracMomLarge(0,2.0,moments)*FracMomLarge(0,1.0,moments) + 
		     FracMomLarge(av,as+2.0,moments)*FracMomLarge(-av,-as+1.0,moments) +
		     FracMomLarge(-av,-as+2.0,moments)*FracMomLarge(av,as+1.0,moments) +
		     1.257*lambda*FracMomLarge(-av,-as+2.0,moments)*FracMomLarge(0,1.0,moments) +
		     1.257*lambda*FracMomLarge(0,2.0,moments)*FracMomLarge(-av,-as+1.0,moments) +
		     1.257*lambda*FracMomLarge(-2.0*av,-2.0*as+2.0,moments)*FracMomLarge(av,as+1.0,moments) +
		       1.257*lambda*FracMomLarge(av,as+2.0,moments)*FracMomLarge(-2.0*av,-2.0*as+1.0,moments));
  else
    cont = 4.0 * C * (1.0 + 1.257 * lambda * FracMomLarge(0.0, 0.0, moments) / FracMomLarge(av, as, moments)) * (pow(2.0, i+j-1.0)-1.0) *
      pow(FracMomLarge(1.0, 0.0, moments) / FracMomLarge(0.0, 0.0, moments), i) * pow(FracMomLarge(0.0, 1.0, moments) / FracMomLarge(0.0, 0.0, moments), j) *
      FracMomLarge(0.0, 0.0, moments) * FracMomLarge(0,0,moments);


  if (cont*fm == 0.0)
    return 0.0;
  else
    return cont*fm/(cont+fm);
}

double TSoot::CoagContSL(double i, double j, double * moments, double temp)
{
#ifdef AGG
  if ((i < 1.01 && i > 0.99 && j == 0) || (i == 0 && j < 1.01 && j > 0.99))
    return 0.0;
#endif
#ifdef BLEND
  if (i < 1.01 && i > 0.99 && j == 0)
    return 0.0;
#endif

  double fm;
  double cont;

  double V0 = 2.0 * dimer_nbrC2;

  double Df = 1.8;
  double av = 1.0 - (2.0 / Df);
  double as = (3.0 / Df) - 1.0;

  double FitC = 2.0/3.0;
  double FitExp = -0.2043;

  double C_cont = 8.0 * RGAS/AVOGADRO * temp / (3.0*mu) * AVOGADRO;
  double lambda = 3.0*mu/rho*sqrt(PI*Wmix/(8.0*RGAS*temp)) / pow(6.0*fMolarMassSoot/(PI*fSootDensity*AVOGADRO),1.0/3.0);

#ifdef AGG
  //Free Molecular
  if (i == 0 && j == 0)
    fm = - GetPsiSL(0.0, 0.0, 0.0, 0.0, moments, temp);
  else if (i < 2.01 && i > 1.99 && j == 0)
    fm = 2.0 * GetPsiSL(1.0, 0.0, 1.0, 0.0, moments, temp);
  else if (i < 1.01 && i > 0.99 && j < 1.01 && j > 0.99)
    fm = GetPsiSL(0.0, 1.0, 1.0, 0.0, moments, temp) + GetPsiSL(1.0, 0.0, 0.0, 1.0, moments, temp);
  else if (i == 0 && j < 2.01 && j > 1.99)
    fm = 2.0 * GetPsiSL(0.0, 1.0, 0.0, 1.0, moments, temp);
  else if (i < 3.01 && i > 2.99 && j == 0)
    fm = 3.0 * GetPsiSL(2.0, 0.0, 1.0, 0.0, moments, temp) + 3.0 * GetPsiSL(1.0, 0.0, 2.0, 0.0, moments, temp);
  else if (i < 2.01 && i > 1.99 && j < 1.01 && j > 0.99)
    fm = GetPsiSL(2.0, 0.0, 0.0, 1.0, moments, temp) + 2.0 * GetPsiSL(1.0, 1.0, 1.0, 0.0, moments, temp) + 2.0 * GetPsiSL(1.0, 0.0, 1.0, 1.0, moments, temp) + GetPsiSL(0.0, 1.0, 2.0, 0.0, moments, temp);
  else if (i < 1.01 && i > 0.99 && j < 2.01 && j > 1.99)
    fm = GetPsiSL(0.0, 2.0, 1.0, 0.0, moments, temp) + 2.0 * GetPsiSL(1.0, 1.0, 0.0, 1.0, moments, temp) + 2.0 * GetPsiSL(0.0, 1.0, 1.0, 1.0, moments, temp) + GetPsiSL(1.0, 0.0, 0.0, 2.0, moments, temp);
  else if (i == 0 && j < 3.01 && j > 2.99)
    fm = 3.0 * GetPsiSL(0.0, 2.0, 0.0, 1.0, moments, temp) + 3.0 * GetPsiSL(0.0, 1.0, 0.0, 2.0, moments, temp);
  else if (i < 4.01 && i > 3.99 && j == 0)
    fm = 4.0 * GetPsiSL(3.0, 0.0, 1.0, 0.0, moments, temp) + 6.0 * GetPsiSL(2.0, 0.0, 2.0, 0.0, moments, temp) + 4.0 * GetPsiSL(1.0, 0.0, 3.0, 0.0, moments, temp);
  else if (i < 3.01 && i > 2.99 && j < 1.01 && j > 0.99)
    fm = 3.0 * GetPsiSL(2.0, 1.0, 1.0, 0.0, moments, temp) + 3.0 * GetPsiSL(1.0, 1.0, 2.0, 0.0, moments, temp) + GetPsiSL(0.0, 1.0, 3.0, 0.0, moments, temp) + GetPsiSL(3.0, 0.0, 0.0, 1.0, moments, temp) + 3.0 * GetPsiSL(2.0, 0.0, 1.0, 1.0, moments, temp) + 3.0 * GetPsiSL(1.0, 0.0, 2.0, 1.0, moments, temp);
  else if (i < 2.01 && i > 1.99 && j < 2.01 && j > 1.99)
    fm = 2.0 * GetPsiSL(2.0, 1.0, 0.0, 1.0, moments, temp) + GetPsiSL(2.0, 0.0, 0.0, 2.0, moments, temp) + 2.0 * GetPsiSL(1.0, 2.0, 1.0, 0.0, moments, temp) + 4.0 * GetPsiSL(1.0, 1.0, 1.0, 1.0, moments, temp) + 2.0 * GetPsiSL(1.0, 0.0, 1.0, 2.0, moments, temp) + GetPsiSL(0.0, 2.0, 2.0, 0.0, moments, temp) + 2.0 * GetPsiSL(0.0, 1.0, 2.0, 1.0, moments, temp);
  else if (i < 1.01 && i > 0.99 && j < 3.01 && j > 2.99)
    fm = 3.0 * GetPsiSL(1.0, 2.0, 0.0, 1.0, moments, temp) + 3.0 * GetPsiSL(1.0, 1.0, 0.0, 2.0, moments, temp) + GetPsiSL(1.0, 0.0, 0.0, 3.0, moments, temp) + GetPsiSL(0.0, 3.0, 1.0, 0.0, moments, temp) + 3.0 * GetPsiSL(0.0, 2.0, 1.0, 1.0, moments, temp) + 3.0 * GetPsiSL(0.0, 1.0, 1.0, 2.0, moments, temp);
  else if (i == 0 && j < 4.01 && j > 3.99)
    fm = 4.0 * GetPsiSL(0.0, 3.0, 0.0, 1.0, moments, temp) + 6.0 * GetPsiSL(0.0, 2.0, 0.0, 2.0, moments, temp) + 4.0 * GetPsiSL(0.0, 1.0, 0.0, 3.0, moments, temp);
  else
    fm = GetPhiSL_Agg(i, j, moments, temp) - GetPsiSL(i, j, 0.0, 0.0, moments, temp) - GetPsiSL(0.0, 0.0, i, j, moments, temp);

  //Continuum
  if (i == 0 && j == 0)
    cont = -C_cont * (2.000*                        moments[fNSootMoments-1]*FracMomLarge( 0.0   , 0.0   ,moments) + 
		                   pow(V0, 1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-    av,-    as,moments) +
		                   pow(V0,-1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av,     as,moments) +
		      1.257*lambda*pow(V0,-1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge( 0.0   , 0.0   ,moments) +
		      1.257*lambda*                 moments[fNSootMoments-1]*FracMomLarge(-    av,-    as,moments) +
		      1.257*lambda*pow(V0,-2.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av,     as,moments) +
		      1.257*lambda*pow(V0, 1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-2.0*av,-2.0*as,moments));
  else if (i < 2.01 && i > 1.99 && j == 0)
    cont = 2.0 * C_cont * (2.000*           V0         *moments[fNSootMoments-1]*FracMomLarge(        1.0, 0.0   ,moments) + 
	                                pow(V0,4.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-    av+1.0,-    as,moments) +
	                                pow(V0,2.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av+1.0,     as,moments) +
	                   1.257*lambda*pow(V0,2.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(        1.0, 0.0   ,moments) +
	                   1.257*lambda*    V0         *moments[fNSootMoments-1]*FracMomLarge(-    av+1.0,-    as,moments) +
	                   1.257*lambda*pow(V0,1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av+1.0,     as,moments) +
		           1.257*lambda*pow(V0,4.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-2.0*av+1.0,-2.0*as,moments));
  else if (i < 1.01 && i > 0.99 && j < 1.01 && j > 0.99)
    cont = C_cont * (2.000*           V0         *moments[fNSootMoments-1]*FracMomLarge(        0.0,        1.0,moments) + 
	                          pow(V0,4.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-    av    ,-    as+1.0,moments) +
	                          pow(V0,2.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av    ,     as+1.0,moments) +
	             1.257*lambda*pow(V0,2.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(        0.0,        1.0,moments) +
	             1.257*lambda*    V0         *moments[fNSootMoments-1]*FracMomLarge(-    av    ,-    as+1.0,moments) +
	             1.257*lambda*pow(V0,1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av    ,     as+1.0,moments) +
		     1.257*lambda*pow(V0,4.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-2.0*av    ,-2.0*as+1.0,moments) +
		     2.000*       pow(V0,2.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(        1.0,        0.0,moments) + 
	                              V0         *moments[fNSootMoments-1]*FracMomLarge(-    av+1.0,-    as    ,moments) +
	                          pow(V0,1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av+1.0,     as    ,moments) +
	             1.257*lambda*pow(V0,1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(        1.0,        0.0,moments) +
		     1.257*lambda*pow(V0,2.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-    av+1.0,-    as    ,moments) +
	             1.257*lambda*                moments[fNSootMoments-1]*FracMomLarge(     av+1.0,     as    ,moments) +
		     1.257*lambda*    V0         *moments[fNSootMoments-1]*FracMomLarge(-2.0*av+1.0,-2.0*as    ,moments));
  else if (i == 0 && j < 2.01 && j > 1.99)
    cont = 2.0 * C_cont * (2.000*       pow(V0,2.0/3.0)*moments[fNSootMoments-1]*FracMomLarge( 0.0   ,        1.0,moments) + 
		                            V0         *moments[fNSootMoments-1]*FracMomLarge(-    av,-    as+1.0,moments) +
	                                pow(V0,1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av,     as+1.0,moments) +
	                   1.257*lambda*pow(V0,1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge( 0.0   ,        1.0,moments) +
	                   1.257*lambda*pow(V0,2.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-    av,-    as+1.0,moments) +
	                   1.257*lambda*                moments[fNSootMoments-1]*FracMomLarge(     av,     as+1.0,moments) +
		           1.257*lambda*    V0         *moments[fNSootMoments-1]*FracMomLarge(-2.0*av,-2.0*as+1.0,moments));
  else if (i < 3.01 && i > 2.99 && j == 0)
    cont = C_cont * (
		     3.0 * (2.000*pow(V0, 6.0/3.0)*moments[fNSootMoments-1]*FracMomLarge( 1.0   , 0.0   ,moments) + 
		                  pow(V0, 7.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-    av+1.0,-    as,moments) +
		                  pow(V0, 5.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av+1.0,     as,moments) +
		      1.257*lambda*pow(V0, 5.0/3.0)*moments[fNSootMoments-1]*FracMomLarge( 1.0   , 0.0   ,moments) +
		      1.257*lambda*pow(V0, 6.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-    av+1.0,-    as,moments) +
		      1.257*lambda*pow(V0, 4.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av+1.0,     as,moments) +
		      1.257*lambda*pow(V0, 7.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-2.0*av+1.0,-2.0*as,moments)
		     ) + 3.0 * ( 
			     2.000*pow(V0, 3.0/3.0)*moments[fNSootMoments-1]*FracMomLarge( 2.0   , 0.0   ,moments) + 
		                   pow(V0, 4.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-    av+2.0,-    as,moments) +
		                   pow(V0, 2.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av+2.0,     as,moments) +
		      1.257*lambda*pow(V0, 2.0/3.0)*moments[fNSootMoments-1]*FracMomLarge( 2.0   , 0.0   ,moments) +
		      1.257*lambda*pow(V0, 3.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-    av+2.0,-    as,moments) +
		      1.257*lambda*pow(V0, 1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av+2.0,     as,moments) +
		      1.257*lambda*pow(V0, 4.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-2.0*av+2.0,-2.0*as,moments)
		     ));
  else if (i < 2.01 && i > 1.99 && j < 1.01 && j > 0.99)
    cont = C_cont * (2.000*                        moments[fNSootMoments-1]*FracMomLarge( 0.0   , 0.0   ,moments) + 
		                   pow(V0, 1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-    av,-    as,moments) +
		                   pow(V0,-1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av,     as,moments) +
		      1.257*lambda*pow(V0,-1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge( 0.0   , 0.0   ,moments) +
		      1.257*lambda*                 moments[fNSootMoments-1]*FracMomLarge(-    av,-    as,moments) +
		      1.257*lambda*pow(V0,-2.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av,     as,moments) +
		      1.257*lambda*pow(V0, 1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-2.0*av,-2.0*as,moments));
  else if (i < 1.01 && i > 0.99 && j < 2.01 && j > 1.99)
    cont = C_cont * (2.000*                        moments[fNSootMoments-1]*FracMomLarge( 0.0   , 0.0   ,moments) + 
		                   pow(V0, 1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-    av,-    as,moments) +
		                   pow(V0,-1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av,     as,moments) +
		      1.257*lambda*pow(V0,-1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge( 0.0   , 0.0   ,moments) +
		      1.257*lambda*                 moments[fNSootMoments-1]*FracMomLarge(-    av,-    as,moments) +
		      1.257*lambda*pow(V0,-2.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av,     as,moments) +
		      1.257*lambda*pow(V0, 1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-2.0*av,-2.0*as,moments));
  else if (i == 0 && j < 3.01 && j > 2.99)
    cont = C_cont * (
	       3.0 * (2.000*       pow(V0, 4.0/3.0)*moments[fNSootMoments-1]*FracMomLarge( 0.0   , 1.0   ,moments) + 
		                   pow(V0, 5.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-    av,-    as+1.0,moments) +
		                   pow(V0, 3.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av,     as+1.0,moments) +
		      1.257*lambda*pow(V0, 3.0/3.0)*moments[fNSootMoments-1]*FracMomLarge( 0.0   , 1.0   ,moments) +
		      1.257*lambda*pow(V0, 4.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-    av,-    as+1.0,moments) +
		      1.257*lambda*pow(V0, 2.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av,     as+1.0,moments) +
		      1.257*lambda*pow(V0, 5.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-2.0*av,-2.0*as+1.0,moments)
		      ) + 3.0 * (
		      2.000*       pow(V0, 2.0/3.0)*moments[fNSootMoments-1]*FracMomLarge( 0.0   , 2.0   ,moments) + 
		                   pow(V0, 3.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-    av,-    as+2.0,moments) +
		                   pow(V0, 1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av,     as+2.0,moments) +
		      1.257*lambda*pow(V0, 1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge( 0.0   , 2.0   ,moments) +
		      1.257*lambda*pow(V0, 2.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-    av,-    as+2.0,moments) +
		      1.257*lambda*pow(V0, 0.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av,     as+2.0,moments) +
		      1.257*lambda*pow(V0, 3.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-2.0*av,-2.0*as+2.0,moments)
		     ));
  else
    cont = C_cont * ( (1.0 + 1.257 * lambda * pow(V0, -1.0/3.0)) * pow(V0, -1.0/3.0) + (1.0 + 1.257 * lambda * FracMomLarge(0,0,moments) / FracMomLarge(av,as,moments)) * FracMomLarge(0,0,moments) / FracMomLarge(av,as,moments)) * (pow(V0, 1.0/3.0) + FracMomLarge(0,0,moments) / FracMomLarge(av,as,moments)) * 
      (pow(V0 + FracMomLarge(1,0,moments) / FracMomLarge(0,0,moments), i) * pow(pow(V0,2.0/3.0) + FracMomLarge(0,1,moments) / FracMomLarge(0,0,moments), j) -
       pow(V0, i+2.0/3.0*j) - pow(FracMomLarge(1,0,moments)/FracMomLarge(0,0,moments),i) * pow(FracMomLarge(0,1,moments)/FracMomLarge(0,0,moments), j)) * 
      moments[fNSootMoments-1] * FracMomLarge(0,0,moments);
#endif

#ifdef BLEND
  //Free Molecular
  if (i == 0 && j == 0)
    fm = - GetPsiSL(0.0, 0.0, 0.0, 0.0, moments, temp);
  else if (i == 0 && j < 1.01 && j > 0.99)
    fm =   FitC * GetPsiSL(-2.0*FitExp-1.0, 3.0*FitExp+1.0, 1.0, 0.0, moments, temp) 
         -        GetPsiSL(            0.0,            0.0, 0.0, 1.0, moments, temp);
  else if (i < 2.01 && i > 1.99 && j == 0)
    fm = 2.0 * GetPsiSL(1.0, 0.0, 1.0, 0.0, moments, temp);
  else if (i < 1.01 && i > 0.99 && j < 1.01 && j > 0.99)
    fm =   FitC * GetPsiSL(-2.0*FitExp    , 3.0*FitExp+1.0, 1.0, 0.0, moments, temp)
         +        GetPsiSL(            0.0,            1.0, 1.0, 0.0, moments, temp)
         + FitC * GetPsiSL(-2.0*FitExp-1.0, 3.0*FitExp+1.0, 2.0, 0.0, moments, temp)
         -        GetPsiSL(            0.0,            0.0, 1.0, 1.0, moments, temp);
  else if (i == 0 && j < 2.01 && j > 1.99)
    fm =   2.0 * FitC      * GetPsiSL(-2.0*FitExp-1.0, 3.0*FitExp+2.0, 1.0, 0.0, moments, temp)
         +       FitC*FitC * GetPsiSL(-4.0*FitExp-2.0, 6.0*FitExp+2.0, 2.0, 0.0, moments, temp)
         -                   GetPsiSL(            0.0,            0.0, 0.0, 2.0, moments, temp);

  //Continuum
  if (i == 0 && j == 0)
    cont = -C_cont * (2.000*                        moments[fNSootMoments-1]*FracMomLarge( 0.0   , 0.0   ,moments) + 
		                   pow(V0, 1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-    av,-    as,moments) +
		                   pow(V0,-1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av,     as,moments) +
		      1.257*lambda*pow(V0,-1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge( 0.0   , 0.0   ,moments) +
		      1.257*lambda*                 moments[fNSootMoments-1]*FracMomLarge(-    av,-    as,moments) +
		      1.257*lambda*pow(V0,-2.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av,     as,moments) +
		      1.257*lambda*pow(V0, 1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-2.0*av,-2.0*as,moments));
  else if (i == 0 && j < 1.01 && j > 0.99)
    cont = C_cont * (
		      FitC * (2.000*           V0          *moments[fNSootMoments-1]*FracMomLarge(       -2.0*FitExp-1.0,        3.0*FitExp+1.0,moments) + 
		                           pow(V0, 4.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-    av-2.0*FitExp-1.0,-    as+3.0*FitExp+1.0,moments) +
		                           pow(V0, 2.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av-2.0*FitExp-1.0,     as+3.0*FitExp+1.0,moments) +
		              1.257*lambda*pow(V0, 2.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(       -2.0*FitExp-1.0,        3.0*FitExp+1.0,moments) +
		              1.257*lambda*    V0          *moments[fNSootMoments-1]*FracMomLarge(-    av-2.0*FitExp-1.0,-    as+3.0*FitExp+1.0,moments) +
		              1.257*lambda*pow(V0, 1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av-2.0*FitExp-1.0,     as+3.0*FitExp+1.0,moments) +
		              1.257*lambda*pow(V0, 4.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-2.0*av-2.0*FitExp-1.0,-2.0*as+3.0*FitExp+1.0,moments)
		      ) - (
			      2.000*       pow(V0, 2.0/3.0)*moments[fNSootMoments-1]*FracMomLarge( 0.0   , 0.0   ,moments) + 
		                           pow(V0, 3.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-    av,-    as,moments) +
		                           pow(V0, 1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av,     as,moments) +
		              1.257*lambda*pow(V0, 1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge( 0.0   , 0.0   ,moments) +
		              1.257*lambda*pow(V0, 2.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-    av,-    as,moments) +
		              1.257*lambda*pow(V0, 0.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av,     as,moments) +
		              1.257*lambda*pow(V0, 3.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-2.0*av,-2.0*as,moments)
		     ));
  else if (i < 2.01 && i > 1.99 && j == 0)
    cont = 2.0 * C_cont * (2.000*           V0         *moments[fNSootMoments-1]*FracMomLarge(        1.0, 0.0   ,moments) + 
	                                pow(V0,4.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-    av+1.0,-    as,moments) +
	                                pow(V0,2.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av+1.0,     as,moments) +
	                   1.257*lambda*pow(V0,2.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(        1.0, 0.0   ,moments) +
	                   1.257*lambda*    V0         *moments[fNSootMoments-1]*FracMomLarge(-    av+1.0,-    as,moments) +
	                   1.257*lambda*pow(V0,1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av+1.0,     as,moments) +
		           1.257*lambda*pow(V0,4.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-2.0*av+1.0,-2.0*as,moments));
  if (i < 1.01 && i > 0.99 && j < 1.01 && j > 0.99)
    cont = C_cont * (
		      FitC * (2.000*           V0          *moments[fNSootMoments-1]*FracMomLarge(       -2.0*FitExp    ,        3.0*FitExp+1.0,moments) + 
		                           pow(V0, 4.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-    av-2.0*FitExp    ,-    as+3.0*FitExp+1.0,moments) +
		                           pow(V0, 2.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av-2.0*FitExp    ,     as+3.0*FitExp+1.0,moments) +
		              1.257*lambda*pow(V0, 2.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(       -2.0*FitExp    ,        3.0*FitExp+1.0,moments) +
			      1.257*lambda*    V0          *moments[fNSootMoments-1]*FracMomLarge(-    av-2.0*FitExp    ,-    as+3.0*FitExp+1.0,moments) +
			      1.257*lambda*pow(V0, 1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av-2.0*FitExp    ,     as+3.0*FitExp+1.0,moments) +
			      1.257*lambda*pow(V0, 4.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-2.0*av-2.0*FitExp    ,-2.0*as+3.0*FitExp+1.0,moments)
		      ) + (
			      2.000*           V0          *moments[fNSootMoments-1]*FracMomLarge( 0.0   ,        1.0,moments) + 
		                           pow(V0, 4.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-    av,-    as+1.0,moments) +
		                           pow(V0, 2.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av,     as+1.0,moments) +
			      1.257*lambda*pow(V0, 2.0/3.0)*moments[fNSootMoments-1]*FracMomLarge( 0.0   ,        1.0,moments) +
			      1.257*lambda*    V0          *moments[fNSootMoments-1]*FracMomLarge(-    av,-    as+1.0,moments) +
			      1.257*lambda*pow(V0, 1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av,     as+1.0,moments) +
			      1.257*lambda*pow(V0, 4.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-2.0*av,-2.0*as+1.0,moments)
		      ) + FitC * (
			      2.000*       pow(V0, 6.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(       -2.0*FitExp-1.0,        3.0*FitExp+1.0,moments) + 
		                           pow(V0, 7.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-    av-2.0*FitExp-1.0,-    as+3.0*FitExp+1.0,moments) +
		                           pow(V0, 5.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av-2.0*FitExp-1.0,     as+3.0*FitExp+1.0,moments) +
			      1.257*lambda*pow(V0, 5.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(       -2.0*FitExp-1.0,        3.0*FitExp+1.0,moments) +
			      1.257*lambda*pow(V0, 6.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-    av-2.0*FitExp-1.0,-    as+3.0*FitExp+1.0,moments) +
			      1.257*lambda*pow(V0, 4.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av-2.0*FitExp-1.0,     as+3.0*FitExp+1.0,moments) +
			      1.257*lambda*pow(V0, 7.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-2.0*av-2.0*FitExp-1.0,-2.0*as+3.0*FitExp+1.0,moments)
		      ) - (
			      2.000*       pow(V0, 5.0/3.0)*moments[fNSootMoments-1]*FracMomLarge( 0.0   , 0.0   ,moments) + 
		                           pow(V0, 6.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-    av,-    as,moments) +
		                           pow(V0, 4.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av,     as,moments) +
			      1.257*lambda*pow(V0, 4.0/3.0)*moments[fNSootMoments-1]*FracMomLarge( 0.0   , 0.0   ,moments) +
			      1.257*lambda*pow(V0, 5.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-    av,-    as,moments) +
			      1.257*lambda*pow(V0, 3.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av,     as,moments) +
			      1.257*lambda*pow(V0, 6.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-2.0*av,-2.0*as,moments)
		     ));
  else if (i == 0 && j < 2.01 && j > 1.99)
    cont = 2.0 * C_cont * (
			   2.0 * FitC * (2.000*           V0          *moments[fNSootMoments-1]*FracMomLarge(       -2.0*FitExp-1.0,        3.0*FitExp+2.0,moments) + 
		                                      pow(V0, 4.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-    av-2.0*FitExp-1.0,-    as+3.0*FitExp+2.0,moments) +
					              pow(V0, 2.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av-2.0*FitExp-1.0,     as+3.0*FitExp+2.0,moments) +
		                         1.257*lambda*pow(V0, 2.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(       -2.0*FitExp-1.0,        3.0*FitExp+2.0,moments) +
					 1.257*lambda*    V0          *moments[fNSootMoments-1]*FracMomLarge(-    av-2.0*FitExp-1.0,-    as+3.0*FitExp+2.0,moments) +
					 1.257*lambda*pow(V0, 1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av-2.0*FitExp-1.0,     as+3.0*FitExp+2.0,moments) +
					 1.257*lambda*pow(V0, 4.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-2.0*av-2.0*FitExp-1.0,-2.0*as+3.0*FitExp+2.0,moments)
			   ) + FitC * FitC * (
					 2.000*       pow(V0, 6.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(       -4.0*FitExp-2.0,        6.0*FitExp+2.0,moments) + 
		                                      pow(V0, 7.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-    av-4.0*FitExp-2.0,-    as+6.0*FitExp+2.0,moments) +
		                                      pow(V0, 5.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av-4.0*FitExp-2.0,     as+6.0*FitExp+2.0,moments) +
		                         1.257*lambda*pow(V0, 5.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(       -4.0*FitExp-2.0,        6.0*FitExp+2.0,moments) +
					 1.257*lambda*pow(V0, 6.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-    av-4.0*FitExp-2.0,-    as+6.0*FitExp+2.0,moments) +
					 1.257*lambda*pow(V0, 4.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av-4.0*FitExp-2.0,     as+6.0*FitExp+2.0,moments) +
					 1.257*lambda*pow(V0, 7.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-2.0*av-4.0*FitExp-2.0,-2.0*as+6.0*FitExp+2.0,moments)
			   ) - (
				         2.000*       pow(V0, 4.0/3.0)*moments[fNSootMoments-1]*FracMomLarge( 0.0   , 0.0   ,moments) + 
					              pow(V0, 5.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-    av,-    as,moments) +
					              pow(V0, 3.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av,     as,moments) +
		                         1.257*lambda*pow(V0, 3.0/3.0)*moments[fNSootMoments-1]*FracMomLarge( 0.0   , 0.0   ,moments) +
					 1.257*lambda*pow(V0, 4.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-    av,-    as,moments) +
					 1.257*lambda*pow(V0, 2.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av,     as,moments) +
					 1.257*lambda*pow(V0, 5.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-2.0*av,-2.0*as,moments)
			  ));
#endif

  if (fm*cont == 0.0)
    return 0.0;
  else
    return (cont * fm) / (cont + fm);
}

double TSoot::GetPhiSL_Agg(double x, double y, double * moments, double temp)
{
  double C = GetC(temp);

  double phi1, phi2, phi3, phi4, phi5, phi6, phi7, phi8, phi9, phi10;
  phi1 = phi2 = phi3 = phi4 = phi5 = phi6 = phi7 = phi8 = phi9 = phi10 = 0;

  double V0 = 2.0 * dimer_nbrC2;

  double Df = 1.8;
  double av = 1.0 - (2.0 / Df);
  double as = (3.0 / Df) - 1.0;

  //phi1  = Phi(-1.0/2.0, 0)
  //phi2  = Phi( 1.0/2.0, 0)
  //phi3  = Phi(-1.0/2.0, 1)
  //phi4  = Phi( 3.0/2.0, 0)
  //phi5  = Phi( 1.0/2.0, 1)
  //phi6  = Phi(-1.0/2.0, 2)
  //phi7  = Phi( 5.0/2.0, 0)
  //phi8  = Phi( 3.0/2.0, 1)
  //phi9  = Phi( 1.0/2.0, 2)
  //phi10 = Phi(-1.0/2.0, 3)

  /*int index;

  if (x+y <= 0.5)
    index = 0;
  else if (x+y <= 1.5)
    if (x >= 0.5)
      index = 1;
    else if (y < 1.0)
      index = 2;
    else
      index = 3;
  else if ((x+y) <= 2.5)
    if (x >= 1.5)
      index = 4;
    else if (y < 1.0)
      index = 5;
    else if (x >= 0.5)
      index = 6;
    else if (y < 2.0)
      index = 7;
    else
      index = 8;
  else
    cerr << "Coagulation has not yet been implemented for more than 6+1 moments.\n";*/

  //if (index == 0)
    phi1 = (
	          pow(V0,   1.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(      -1.0/2.0,        0.0, moments) +
	    2.0 * pow(V0, - 1.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(    av-1.0/2.0,     as    , moments) +
	          pow(V0, - 1.0/2.0) * moments[fNSootMoments-1] * FracMomLarge(2.0*av-1.0/2.0, 2.0*as    , moments)
	   );

  //if (index == 0 || index == 1 || index == 2)
    phi2 = (
	          pow(V0,   7.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(      -1.0/2.0,        0.0, moments) +
	    2.0 * pow(V0,   5.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(    av-1.0/2.0,     as    , moments) +
	          pow(V0,   1.0/2.0) * moments[fNSootMoments-1] * FracMomLarge(2.0*av-1.0/2.0, 2.0*as    , moments)
	   ) + (
		  pow(V0,   1.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(       1.0/2.0,        0.0, moments) +
	    2.0 * pow(V0, - 1.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(    av+1.0/2.0,     as    , moments) +
		  pow(V0, - 1.0/2.0) * moments[fNSootMoments-1] * FracMomLarge(2.0*av+1.0/2.0, 2.0*as    , moments)
	   );

  //if (index == 0 || index == 2 || index == 3)
    phi3 = (
	          pow(V0,   5.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(      -1.0/2.0,        0.0, moments) +
	    2.0 * pow(V0,   1.0/2.0) * moments[fNSootMoments-1] * FracMomLarge(    av-1.0/2.0,     as    , moments) +
	          pow(V0,   1.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(2.0*av-1.0/2.0, 2.0*as    , moments)
	   ) + (
	          pow(V0,   1.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(      -1.0/2.0,        1.0, moments) +
	    2.0 * pow(V0, - 1.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(    av-1.0/2.0,     as+1.0, moments) +
		  pow(V0, - 1.0/2.0) * moments[fNSootMoments-1] * FracMomLarge(2.0*av-1.0/2.0, 2.0*as+1.0, moments)
	    );

  //if (index == 1 || index == 4 || index == 5)
    phi4 = (
	          pow(V0,  13.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(      -1.0/2.0,        0.0, moments) +
	    2.0 * pow(V0,  11.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(    av-1.0/2.0,     as    , moments) +
	          pow(V0,   3.0/2.0) * moments[fNSootMoments-1] * FracMomLarge(2.0*av-1.0/2.0, 2.0*as    , moments)
	   ) + 2.0 * (
		  pow(V0,   7.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(       1.0/2.0,        0.0, moments) +
	    2.0 * pow(V0,   5.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(    av+1.0/2.0,     as    , moments) +
		  pow(V0,   1.0/2.0) * moments[fNSootMoments-1] * FracMomLarge(2.0*av+1.0/2.0, 2.0*as    , moments)
	   ) + (
		  pow(V0,   1.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(       3.0/2.0,        0.0, moments) +
	    2.0 * pow(V0, - 1.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(    av+3.0/2.0,     as    , moments) +
		  pow(V0, - 1.0/2.0) * moments[fNSootMoments-1] * FracMomLarge(2.0*av+3.0/2.0, 2.0*as    , moments)
           );

  //if (index == 1 || index == 2 || index == 3 || index == 5 || index == 6 || index == 7)
    phi5 = (
	          pow(V0,  11.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(      -1.0/2.0,        0.0, moments) +
	    2.0 * pow(V0,   3.0/2.0) * moments[fNSootMoments-1] * FracMomLarge(    av-1.0/2.0,     as    , moments) +
	          pow(V0,   7.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(2.0*av-1.0/2.0, 2.0*as    , moments)
	   ) + (
	          pow(V0,   7.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(      -1.0/2.0,        1.0, moments) +
	    2.0 * pow(V0,   5.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(    av-1.0/2.0,     as+1.0, moments) +
		  pow(V0,   1.0/2.0) * moments[fNSootMoments-1] * FracMomLarge(2.0*av-1.0/2.0, 2.0*as+1.0, moments)
	   ) + (
		  pow(V0,   5.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(       1.0/2.0,        0.0, moments) +
	    2.0 * pow(V0,   1.0/2.0) * moments[fNSootMoments-1] * FracMomLarge(    av+1.0/2.0,     as    , moments) +
		  pow(V0,   1.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(2.0*av+1.0/2.0, 2.0*as    , moments)
	   ) + (
		  pow(V0,   1.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(       1.0/2.0,        1.0, moments) +
	    2.0 * pow(V0, - 1.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(    av+1.0/2.0,     as+1.0, moments) +
		  pow(V0, - 1.0/2.0) * moments[fNSootMoments-1] * FracMomLarge(2.0*av+1.0/2.0, 2.0*as+1.0, moments)
	   );

  //if (index == 3 || index == 7 || index == 8)
    phi6 = (
	          pow(V0,   3.0/2.0) * moments[fNSootMoments-1] * FracMomLarge(      -1.0/2.0,        0.0, moments) +
	    2.0 * pow(V0,   7.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(    av-1.0/2.0,     as    , moments) +
	          pow(V0,   5.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(2.0*av-1.0/2.0, 2.0*as    , moments)
	   ) + 2.0 * (
		  pow(V0,   5.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(      -1.0/2.0,        1.0, moments) +
	    2.0 * pow(V0,   1.0/2.0) * moments[fNSootMoments-1] * FracMomLarge(    av-1.0/2.0,     as+1.0, moments) +
		  pow(V0,   1.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(2.0*av-1.0/2.0, 2.0*as+1.0, moments)
	   ) + (
		  pow(V0,   1.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(      -1.0/2.0,        2.0, moments) +
	    2.0 * pow(V0, - 1.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(    av-1.0/2.0,     as+2.0, moments) +
		  pow(V0, - 1.0/2.0) * moments[fNSootMoments-1] * FracMomLarge(2.0*av-1.0/2.0, 2.0*as+2.0, moments)
	   );

  //if (index == 4)
    phi7 = (
	          pow(V0,  19.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(      -1.0/2.0,        0.0, moments) +
            2.0 * pow(V0,  17.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(    av-1.0/2.0,     as    , moments) +
                  pow(V0,   5.0/2.0) * moments[fNSootMoments-1] * FracMomLarge(2.0*av-1.0/2.0, 2.0*as    , moments)
	   ) + 3.0 * (
		  pow(V0,  13.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(       1.0/2.0,        0.0, moments) +
	    2.0 * pow(V0,  11.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(    av+1.0/2.0,     as    , moments) +
		  pow(V0,   3.0/2.0) * moments[fNSootMoments-1] * FracMomLarge(2.0*av+1.0/2.0, 2.0*as    , moments)
	   ) + 3.0 * (
		  pow(V0,   7.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(       3.0/2.0,        0.0, moments) +
	    2.0 * pow(V0,   5.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(    av+3.0/2.0,     as    , moments) +
		  pow(V0,   1.0/2.0) * moments[fNSootMoments-1] * FracMomLarge(2.0*av+3.0/2.0, 2.0*as    , moments)
	   ) + (
		  pow(V0,   1.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(       5.0/2.0,        0.0, moments) +
	    2.0 * pow(V0,  -1.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(    av+5.0/2.0,     as    , moments) +
		  pow(V0,  -1.0/2.0) * moments[fNSootMoments-1] * FracMomLarge(2.0*av+5.0/2.0, 2.0*as    , moments)
	   );

  //if (index == 4 || index == 5 || index == 6)
    phi8  = (
	          pow(V0,  17.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(      -1.0/2.0,        0.0, moments) +
	    2.0 * pow(V0,   5.0/2.0) * moments[fNSootMoments-1] * FracMomLarge(    av-1.0/2.0,     as    , moments) +
		  pow(V0,  13.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(2.0*av-1.0/2.0, 2.0*as    , moments)
	    ) + 2.0 * (
		  pow(V0,  11.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(       1.0/2.0,        0.0, moments) +
	    2.0 * pow(V0,   3.0/2.0) * moments[fNSootMoments-1] * FracMomLarge(    av+1.0/2.0,     as    , moments) +
		  pow(V0,   7.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(2.0*av+1.0/2.0, 2.0*as    , moments)
	    ) + (
		  pow(V0,   5.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(       3.0/2.0,        0.0, moments) +
	    2.0 * pow(V0,   1.0/2.0) * moments[fNSootMoments-1] * FracMomLarge(    av+3.0/2.0,     as    , moments) +
		  pow(V0,   1.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(2.0*av+3.0/2.0, 2.0*as    , moments)
	    ) + (
		  pow(V0,  13.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(      -1.0/2.0,        1.0, moments) +
	    2.0 * pow(V0,  11.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(    av-1.0/2.0,     as+1.0, moments) +
		  pow(V0,   3.0/2.0) * moments[fNSootMoments-1] * FracMomLarge(2.0*av-1.0/2.0, 2.0*as+1.0, moments)
	    ) + 2.0 * (
		  pow(V0,   7.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(       1.0/2.0,        1.0, moments) +
	    2.0 * pow(V0,   5.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(    av+1.0/2.0,     as+1.0, moments) +
		  pow(V0,   1.0/2.0) * moments[fNSootMoments-1] * FracMomLarge(2.0*av+1.0/2.0, 2.0*as+1.0, moments)
	    ) + (
		  pow(V0,   1.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(       3.0/2.0,        1.0, moments) +
	    2.0 * pow(V0, - 1.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(    av+3.0/2.0,     as+1.0, moments) +
		  pow(V0, - 1.0/2.0) * moments[fNSootMoments-1] * FracMomLarge(2.0*av+3.0/2.0, 2.0*as+1.0, moments)
	    );

  //if (index == 6 || index == 7 || index == 8)
    phi9  = (
	           pow(V0,   5.0/2.0) * moments[fNSootMoments-1] * FracMomLarge(      -1.0/2.0,        0.0, moments) +
	     2.0 * pow(V0,  13.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(    av-1.0/2.0,     as    , moments) +
	 	   pow(V0,  11.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(2.0*av-1.0/2.0, 2.0*as    , moments)
	    ) + (
		   pow(V0,   3.0/2.0) * moments[fNSootMoments-1] * FracMomLarge(       1.0/2.0,        0.0, moments) +
	     2.0 * pow(V0,   7.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(    av+1.0/2.0,     as    , moments) +
		   pow(V0,   5.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(2.0*av+1.0/2.0, 2.0*as    , moments)
 	    ) + 2.0 * (
		   pow(V0,  11.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(      -1.0/2.0,        1.0, moments) +
	     2.0 * pow(V0,   3.0/2.0) * moments[fNSootMoments-1] * FracMomLarge(    av-1.0/2.0,     as+1.0, moments) +
		   pow(V0,   7.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(2.0*av-1.0/2.0, 2.0*as+1.0, moments)
	    ) + 2.0 * (
		   pow(V0,   5.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(       1.0/2.0,        1.0, moments) +
	     2.0 * pow(V0,   1.0/2.0) * moments[fNSootMoments-1] * FracMomLarge(    av+1.0/2.0,     as+1.0, moments) +
		   pow(V0,   1.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(2.0*av+1.0/2.0, 2.0*as+1.0, moments)
	    ) + (
	           pow(V0,   7.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(      -1.0/2.0,        2.0, moments) +
	     2.0 * pow(V0,   5.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(    av-1.0/2.0,     as+2.0, moments) +
		   pow(V0,   1.0/2.0) * moments[fNSootMoments-1] * FracMomLarge(2.0*av-1.0/2.0, 2.0*as+2.0, moments)
	    ) + (
		   pow(V0,   1.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(       1.0/2.0,        2.0, moments) +
	     2.0 * pow(V0, - 1.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(    av+1.0/2.0,     as+2.0, moments) +
		   pow(V0, - 1.0/2.0) * moments[fNSootMoments-1] * FracMomLarge(2.0*av+1.0/2.0, 2.0*as+2.0, moments)
	    );

  //if (index == 8)
    phi10 = (
	           pow(V0,  13.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(      -1.0/2.0,        0.0, moments) +
	     2.0 * pow(V0,  11.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(    av-1.0/2.0,     as    , moments) +
	           pow(V0,   3.0/2.0) * moments[fNSootMoments-1] * FracMomLarge(2.0*av-1.0/2.0, 2.0*as    , moments)
	    ) + 3.0 * (
		   pow(V0,   3.0/2.0) * moments[fNSootMoments-1] * FracMomLarge(      -1.0/2.0,        1.0, moments) +
	     2.0 * pow(V0,   7.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(    av-1.0/2.0,     as+1.0, moments) +
		   pow(V0,   5.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(2.0*av-1.0/2.0, 2.0*as+1.0, moments)
	    ) + 3.0 * (
		   pow(V0,   5.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(      -1.0/2.0,        2.0, moments) +
	     2.0 * pow(V0,   1.0/2.0) * moments[fNSootMoments-1] * FracMomLarge(    av-1.0/2.0,     as+2.0, moments) +
		   pow(V0,   1.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(2.0*av-1.0/2.0, 2.0*as+2.0, moments)
	    ) + (
		   pow(V0,   1.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(      -1.0/2.0,        3.0, moments) +
	     2.0 * pow(V0, - 1.0/6.0) * moments[fNSootMoments-1] * FracMomLarge(    av-1.0/2.0,     as+3.0, moments) +
		   pow(V0, - 1.0/2.0) * moments[fNSootMoments-1] * FracMomLarge(2.0*av-1.0/2.0, 2.0*as+3.0, moments)
	    );

  /*switch (index)
  {
    case 0: return C * pow(phi1, (1.0/2.0)-x-y) * pow(phi2, x+(1.0/2.0)) * pow(phi3 , y);
    case 1: return C * pow(phi2, (3.0/2.0)-x-y) * pow(phi4, x-(1.0/2.0)) * pow(phi5 , y);
    case 2: return C * pow(phi5, x+y-(1.0/2.0)) * pow(phi3, (1.0/2.0)-x) * pow(phi2 , 1.0-y);
    case 3: return C * pow(phi3, (3.0/2.0)-x-y) * pow(phi5, x+(1.0/2.0)) * pow(phi6 , y-1.0);
    case 4: return C * pow(phi4, (5.0/2.0)-x-y) * pow(phi7, x-(3.0/2.0)) * pow(phi8 , y);
    case 5: return C * pow(phi8, x+y-(3.0/2.0)) * pow(phi5, (3.0/2.0)-x) * pow(phi4 , 1.0-y);
    case 6: return C * pow(phi5, (5.0/2.0)-x-y) * pow(phi8, x-(1.0/2.0)) * pow(phi9 , y-1.0);
    case 7: return C * pow(phi9, x+y-(3.0/2.0)) * pow(phi6, (1.0/2.0)-x) * pow(phi5 , 2.0-y);
    case 8: return C * pow(phi6, (5.0/2.0)-x-y) * pow(phi9, x+(1.0/2.0)) * pow(phi10, y-2.0);
  }*/

  double a = pow(phi1 ,   5.0/16.0)*pow(phi2 ,  15.0/16.0)*pow(phi4 , - 5.0/16.0)*pow(phi7 ,   1.0/16.0);
  double b = pow(phi1 , -23.0/24.0)*pow(phi2 ,   7.0/ 8.0)*pow(phi4 ,   1.0/ 8.0)*pow(phi7 , - 1.0/24.0);
  double c = pow(phi1 , -23.0/24.0)/    phi2              *pow(phi3 ,  15.0/ 8.0)*pow(phi4 ,   1.0/ 8.0)*pow(phi5 ,   5.0/ 4.0)*pow(phi6 , - 5.0/ 4.0)*pow(phi8 , - 1.0/ 8.0)*pow(phi9 , - 1.0/ 4.0)*pow(phi10,   1.0/ 3.0);
  double d = pow(phi1 ,   3.0/ 4.0)*pow(phi2 , - 7.0/ 4.0)*pow(phi4 ,   5.0/ 4.0)*pow(phi7 , - 1.0/ 4.0);
  double e = pow(phi1 ,   3.0/ 2.0)*pow(phi2 , - 3.0/ 2.0)*pow(phi3 , - 2.0     )*pow(phi5 ,   2.0     )*pow(phi6 ,   1.0/ 2.0)*pow(phi9 , - 1.0/ 2.0);
  double f = pow(phi1 ,   3.0/ 4.0)*pow(phi2 ,   1.0/ 4.0)*pow(phi3 , - 2.0     )*pow(phi5 , - 1.0/ 2.0)*pow(phi6 ,   7.0/ 4.0)*pow(phi9 ,   1.0/ 4.0)*pow(phi10, - 1.0/ 2.0);
  double g = pow(phi1 , - 1.0/ 6.0)*pow(phi2 ,   1.0/ 2.0)*pow(phi4 , - 1.0/ 2.0)*pow(phi7 ,   1.0/ 6.0);
  double h = pow(phi1 , - 1.0/ 2.0)*    phi2              *pow(phi3 ,   1.0/ 2.0)*pow(phi4 , - 1.0/ 2.0)*pow(phi5 , - 1.0     )*pow(phi8 ,   1.0/ 2.0);
  double i = pow(phi1 , - 1.0/ 2.0)*pow(phi2 ,   1.0/ 2.0)*    phi3              *pow(phi5 , - 1.0     )*pow(phi6 , - 1.0/ 2.0)*pow(phi9 ,   1.0/ 2.0);
  double j = pow(phi1 , - 1.0/ 6.0)*pow(phi3 ,   1.0/ 2.0)*pow(phi6 , - 1.0/ 2.0)*pow(phi10,   1.0/ 6.0);

  return C*a*pow(b,x)*pow(c,y)*pow(d,x*x)*pow(e,x*y)*pow(f,y*y)*pow(g,x*x*x)*pow(h,x*x*y)*pow(i,x*y*y)*pow(j,y*y*y);
}

double TSoot::GetPsiSL(double x, double y, double a, double b, double * moments, double temp)
{
  double C = GetC(temp);

  double psi1, psi2, psi3, psi4, psi5, psi6;

  double V0 = 2.0 * dimer_nbrC2;

  double Df = 1.8;
  double av = 1.0 - (2.0 / Df);
  double as = (3.0 / Df) - 1.0;

  psi1 = (
	        pow(V0,  1.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(      -0.5+x,        y, moments) +
	  2.0 * pow(V0, -1.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(    av-0.5+x,     as+y, moments) +
	        pow(V0, -3.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(2.0*av-0.5+x, 2.0*as+y, moments)
	 );

  psi2 = (
	        pow(V0,  7.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(      -0.5+x,        y, moments) +
	  2.0 * pow(V0,  5.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(    av-0.5+x,     as+y, moments) +
	        pow(V0,  3.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(2.0*av-0.5+x, 2.0*as+y, moments)
	 ) + (
		pow(V0,  1.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(       0.5+x,        y, moments) +
	  2.0 * pow(V0, -1.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(    av+0.5+x,     as+y, moments) +
		pow(V0, -3.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(2.0*av+0.5+x, 2.0*as+y, moments)
	 );

  psi3 = (
	        pow(V0, 13.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(      -0.5+x,        y, moments) +
	  2.0 * pow(V0, 11.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(    av-0.5+x,     as+y, moments) +
	        pow(V0,  9.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(2.0*av-0.5+x, 2.0*as+y, moments)
	 ) + 2.0 * (
		pow(V0,  7.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(       0.5+x,        y, moments) +
	  2.0 * pow(V0,  5.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(    av+0.5+x,     as+y, moments) +
		pow(V0,  3.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(2.0*av+0.5+x, 2.0*as+y, moments)
	 ) + (
		pow(V0,  1.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(       1.5+x,        y, moments) +
	  2.0 * pow(V0, -1.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(    av+1.5+x,     as+y, moments) +
		pow(V0, -3.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(2.0*av+1.5+x, 2.0*as+y, moments)
	 );

  if (fNSootMoments > 4)
  psi4 = (
	        pow(V0, 19.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(      -0.5+x,        y, moments) +
	  2.0 * pow(V0, 17.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(    av-0.5+x,     as+y, moments) +
	        pow(V0, 15.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(2.0*av-0.5+x, 2.0*as+y, moments)
	 ) + 3.0 * (
		pow(V0, 13.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(       0.5+x,        y, moments) +
	  2.0 * pow(V0, 11.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(    av+0.5+x,     as+y, moments) +
		pow(V0,  9.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(2.0*av+0.5+x, 2.0*as+y, moments)
	 ) + 3.0 * (
		pow(V0,  7.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(       1.5+x,        y, moments) +
	  2.0 * pow(V0,  5.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(    av+1.5+x,     as+y, moments) +
		pow(V0,  3.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(2.0*av+1.5+x, 2.0*as+y, moments)
	 ) + (
		pow(V0,  1.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(       2.5+x,        y, moments) +
	  2.0 * pow(V0, -1.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(    av+2.5+x,     as+y, moments) +
		pow(V0, -3.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(2.0*av+2.5+x, 2.0*as+y, moments)
	 );

  if (fNSootMoments > 7)
  psi5 = (
	        pow(V0, 25.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(      -0.5+x,        y, moments) +
	  2.0 * pow(V0, 23.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(    av-0.5+x,     as+y, moments) +
	        pow(V0, 21.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(2.0*av-0.5+x, 2.0*as+y, moments)
	 ) + 4.0 * (
		pow(V0, 19.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(       0.5+x,        y, moments) +
	  2.0 * pow(V0, 17.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(    av+0.5+x,     as+y, moments) +
		pow(V0, 15.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(2.0*av+0.5+x, 2.0*as+y, moments)
	 ) + 6.0 * (
		pow(V0, 13.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(       1.5+x,        y, moments) +
	  2.0 * pow(V0, 11.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(    av+1.5+x,     as+y, moments) +
		pow(V0,  9.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(2.0*av+1.5+x, 2.0*as+y, moments)
	 ) + 4.0 * (
		pow(V0,  7.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(       2.5+x,        y, moments) +
	  2.0 * pow(V0,  5.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(    av+2.5+x,     as+y, moments) +
		pow(V0,  3.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(2.0*av+2.5+x, 2.0*as+y, moments)
	 ) + (
		pow(V0,  1.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(       3.5+x,        y, moments) +
	  2.0 * pow(V0, -1.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(    av+3.5+x,     as+y, moments) +
		pow(V0, -3.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(2.0*av+3.5+x, 2.0*as+y, moments)
	 );

  if (fNSootMoments > 11)
  psi6 = (
	        pow(V0, 31.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(      -0.5+x,        y, moments) +
	  2.0 * pow(V0, 29.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(    av-0.5+x,     as+y, moments) +
	        pow(V0, 27.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(2.0*av-0.5+x, 2.0*as+y, moments)
	 ) + 5.0 * (
		pow(V0, 25.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(       0.5+x,        y, moments) +
	  2.0 * pow(V0, 23.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(    av+0.5+x,     as+y, moments) +
		pow(V0, 21.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(2.0*av+0.5+x, 2.0*as+y, moments)
	 ) + 10.0 * (
		pow(V0, 19.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(       1.5+x,        y, moments) +
	  2.0 * pow(V0, 17.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(    av+1.5+x,     as+y, moments) +
		pow(V0, 15.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(2.0*av+1.5+x, 2.0*as+y, moments)
	 ) + 10.0 * (
		pow(V0, 13.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(       2.5+x,        y, moments) +
	  2.0 * pow(V0, 11.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(    av+2.5+x,     as+y, moments) +
		pow(V0,  9.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(2.0*av+2.5+x, 2.0*as+y, moments)
	 ) + 5.0 * (
		pow(V0,  7.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(       3.5+x,        y, moments) +
	  2.0 * pow(V0,  5.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(    av+3.5+x,     as+y, moments) +
		pow(V0,  3.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(2.0*av+3.5+x, 2.0*as+y, moments)
	 ) + (
		pow(V0,  1.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(       4.5+x,        y, moments) +
	  2.0 * pow(V0, -1.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(    av+4.5+x,     as+y, moments) +
		pow(V0, -3.0/6.0+a+(2.0/3.0)*b) * moments[fNSootMoments-1] * FracMomLarge(2.0*av+4.5+x, 2.0*as+y, moments)
	 );

//   cerr << "Psi Small-Large\t" << a << "\t" << b << "\t" << x << "\t" << y << endl;
//   cerr << "Two Terms: " << C * sqrt(psi1*psi2) << endl;
//   cerr << "Three Terms: " << C * pow(psi1,3.0/8.0)*pow(psi2,3.0/4.0)*pow(psi3,-1.0/8.0) << endl;
//   cerr << "Four Terms: " << C * pow(psi1,5.0/16.0)*pow(psi2,15.0/16.0)*pow(psi3,-5.0/16.0)*pow(psi4,1.0/16.0) << endl;

  return C * sqrt(psi1 * psi2);
  if (fNSootMoments == 4)
    return C * pow(psi1,3.0/8.0)*pow(psi2,3.0/4.0)*pow(psi3,-1.0/8.0);
  else if (fNSootMoments == 7)
    return C * pow(psi1,5.0/16.0)*pow(psi2,15.0/16.0)*pow(psi3,-5.0/16.0)*pow(psi4,1.0/16.0);
  else if (fNSootMoments == 11)
    return C * pow(psi1,35.0/128.0)*pow(psi2,35.0/32.0)*pow(psi3,-35.0/64.0)*pow(psi4,7.0/32.0)*pow(psi5,-5.0/128.0);
  else if (fNSootMoments == 16)
    return C * pow(psi1,63.0/256.0)*pow(psi2,315.0/256.0)*pow(psi3,-105.0/128.0)*pow(psi4,63.0/128.0)*pow(psi5,-45.0/256.0)*pow(psi6,7.0/256.0);

  return 0.0;
}

double TSoot::GetPhiLL(double x, double y, double * moments, double temp)
{
  double C = GetC(temp);

  double phi1, phi2, phi3, phi4, phi5, phi6, phi7, phi8, phi9, phi10;
  phi1 = phi2 = phi3 = phi4 = phi5 = phi6 = phi7 = phi8 = phi9 = phi10 = 0;

  double Df = 1.8;
  double av = 1.0 - (2.0 / Df);
  double as = (3.0 / Df) - 1.0;

  //phi1  = Phi(-1.0/2.0, 0)
  //phi2  = Phi( 1.0/2.0, 0)
  //phi3  = Phi(-1.0/2.0, 1)
  //phi4  = Phi( 3.0/2.0, 0)
  //phi5  = Phi( 1.0/2.0, 1)
  //phi6  = Phi(-1.0/2.0, 2)
  //phi7  = Phi( 5.0/2.0, 0)
  //phi8  = Phi( 3.0/2.0, 1)
  //phi9  = Phi( 1.0/2.0, 2)
  //phi10 = Phi(-1.0/2.0, 3)

  /*int index;

  if (x+y <= 0.5)
    index = 0;
  else if (x+y <= 1.5)
    if (x >= 0.5)
      index = 1;
    else if (y < 1.0)
      index = 2;
    else
      index = 3;
  else if ((x+y) <= 2.5)
    if (x >= 1.5)
      index = 4;
    else if (y < 1.0)
      index = 5;
    else if (x >= 0.5)
      index = 6;
    else if (y < 2.0)
      index = 7;
    else
      index = 8;
  else
    cerr << "Coagulation has not yet been implemented for more than 6+1 moments.\n";*/

  //if (index == 0)
    phi1  = (
	           FracMomLarge(2.0*av-1.0/2.0, 2.0*as    , moments) * FracMomLarge(      -1.0/2.0,        0.0, moments) +
	     2.0 * FracMomLarge(    av-1.0/2.0,     as    , moments) * FracMomLarge(    av-1.0/2.0,     as    , moments) +
	           FracMomLarge(      -1.0/2.0,        0.0, moments) * FracMomLarge(2.0*av-1.0/2.0, 2.0*as    , moments)
	    );

  //if (index == 0 || index == 1 || index == 2)
    phi2  = (
	           FracMomLarge(2.0*av+1.0/2.0, 2.0*as    , moments) * FracMomLarge(      -1.0/2.0,        0.0, moments) +
	     2.0 * FracMomLarge(    av+1.0/2.0,     as    , moments) * FracMomLarge(    av-1.0/2.0,     as    , moments) +
	           FracMomLarge(       1.0/2.0,        0.0, moments) * FracMomLarge(2.0*av-1.0/2.0, 2.0*as    , moments)
	    ) + (
	 	   FracMomLarge(2.0*av-1.0/2.0, 2.0*as    , moments) * FracMomLarge(       1.0/2.0,        0.0, moments) +
	     2.0 * FracMomLarge(    av-1.0/2.0,     as    , moments) * FracMomLarge(    av+1.0/2.0,     as    , moments) +
		   FracMomLarge(      -1.0/2.0,        0.0, moments) * FracMomLarge(2.0*av+1.0/2.0, 2.0*as    , moments)
	    );

  //if (index == 0 || index == 2 || index == 3)
    phi3  = (
	           FracMomLarge(2.0*av-1.0/2.0, 2.0*as+1.0, moments) * FracMomLarge(      -1.0/2.0,        0.0, moments) +
	     2.0 * FracMomLarge(    av-1.0/2.0,     as+1.0, moments) * FracMomLarge(    av-1.0/2.0,     as    , moments) +
	           FracMomLarge(      -1.0/2.0,        1.0, moments) * FracMomLarge(2.0*av-1.0/2.0, 2.0*as    , moments)
	    ) + (
		   FracMomLarge(2.0*av-1.0/2.0, 2.0*as    , moments) * FracMomLarge(      -1.0/2.0,        1.0, moments) +
	     2.0 * FracMomLarge(    av-1.0/2.0,     as    , moments) * FracMomLarge(    av-1.0/2.0,     as+1.0, moments) +
                   FracMomLarge(      -1.0/2.0,        0.0, moments) * FracMomLarge(2.0*av-1.0/2.0, 2.0*as+1.0, moments)
	    );

  //if (index == 1 || index == 4 || index == 5)
    phi4  = (
	           FracMomLarge(2.0*av+3.0/2.0, 2.0*as    , moments) * FracMomLarge(      -1.0/2.0,        0.0, moments) +
	     2.0 * FracMomLarge(    av+3.0/2.0,     as    , moments) * FracMomLarge(    av-1.0/2.0,     as    , moments) +
	           FracMomLarge(       3.0/2.0,        0.0, moments) * FracMomLarge(2.0*av-1.0/2.0, 2.0*as    , moments)
	    ) + 2.0 * (
		   FracMomLarge(2.0*av+1.0/2.0, 2.0*as    , moments) * FracMomLarge(       1.0/2.0,        0.0, moments) +
	     2.0 * FracMomLarge(    av+1.0/2.0,     as    , moments) * FracMomLarge(    av+1.0/2.0,     as    , moments) +
		   FracMomLarge(       1.0/2.0,        0.0, moments) * FracMomLarge(2.0*av+1.0/2.0, 2.0*as    , moments)
	    ) + (
		   FracMomLarge(2.0*av-1.0/2.0, 2.0*as    , moments) * FracMomLarge(       3.0/2.0,        0.0, moments) +
	     2.0 * FracMomLarge(    av-1.0/2.0,     as    , moments) * FracMomLarge(    av+3.0/2.0,     as    , moments) +
		   FracMomLarge(      -1.0/2.0,        0.0, moments) * FracMomLarge(2.0*av+3.0/2.0, 2.0*as    , moments)
	    );

  //if (index == 1 || index == 2 || index == 3 || index == 5 || index == 6 || index == 7)
    phi5  = (
	           FracMomLarge(2.0*av+1.0/2.0, 2.0*as+1.0, moments) * FracMomLarge(      -1.0/2.0,        0.0, moments) +
	     2.0 * FracMomLarge(    av+1.0/2.0,     as+1.0, moments) * FracMomLarge(    av-1.0/2.0,     as    , moments) +
	           FracMomLarge(        1.0/2.0,       1.0, moments) * FracMomLarge(2.0*av-1.0/2.0, 2.0*as    , moments)
	    ) + (
	           FracMomLarge(2.0*av+1.0/2.0, 2.0*as    , moments) * FracMomLarge(      -1.0/2.0,        1.0, moments) +
	     2.0 * FracMomLarge(    av+1.0/2.0,     as    , moments) * FracMomLarge(    av-1.0/2.0,     as+1.0, moments) +
		   FracMomLarge(       1.0/2.0,        0.0, moments) * FracMomLarge(2.0*av-1.0/2.0, 2.0*as+1.0, moments)
	    ) + (
		   FracMomLarge(2.0*av-1.0/2.0, 2.0*as+1.0, moments) * FracMomLarge(       1.0/2.0,        0.0, moments) +
	     2.0 * FracMomLarge(    av-1.0/2.0,     as+1.0, moments) * FracMomLarge(    av+1.0/2.0,     as    , moments) +
		   FracMomLarge(      -1.0/2.0,        1.0, moments) * FracMomLarge(2.0*av+1.0/2.0, 2.0*as    , moments)
	    ) + (
		   FracMomLarge(2.0*av-1.0/2.0, 2.0*as    , moments) * FracMomLarge(       1.0/2.0,        1.0, moments) +
	     2.0 * FracMomLarge(    av-1.0/2.0,     as    , moments) * FracMomLarge(    av+1.0/2.0,     as+1.0, moments) +
		   FracMomLarge(      -1.0/2.0,        0.0, moments) * FracMomLarge(2.0*av+1.0/2.0, 2.0*as+1.0, moments)
	    );

  //if (index == 3 || index == 7 || index == 8)
    phi6  = (
		   FracMomLarge(2.0*av-1.0/2.0, 2.0*as+2.0, moments) * FracMomLarge(      -1.0/2.0,        0.0, moments) +
	     2.0 * FracMomLarge(    av-1.0/2.0,     as+2.0, moments) * FracMomLarge(    av-1.0/2.0,     as    , moments) +
		   FracMomLarge(      -1.0/2.0,        2.0, moments) * FracMomLarge(2.0*av-1.0/2.0, 2.0*as    , moments)
	    ) + 2.0 * (
		   FracMomLarge(2.0*av-1.0/2.0, 2.0*as+1.0, moments) * FracMomLarge(      -1.0/2.0,        1.0, moments) +
	     2.0 * FracMomLarge(    av-1.0/2.0,     as+1.0, moments) * FracMomLarge(    av-1.0/2.0,     as+1.0, moments) +
		   FracMomLarge(      -1.0/2.0,        1.0, moments) * FracMomLarge(2.0*av-1.0/2.0, 2.0*as+1.0, moments)
	    ) + (
		   FracMomLarge(2.0*av-1.0/2.0, 2.0*as    , moments) * FracMomLarge(      -1.0/2.0,        2.0, moments) +
	     2.0 * FracMomLarge(    av-1.0/2.0,     as    , moments) * FracMomLarge(    av-1.0/2.0,     as+2.0, moments) +
		   FracMomLarge(      -1.0/2.0,        0.0, moments) * FracMomLarge(2.0*av-1.0/2.0, 2.0*as+2.0, moments)
	    );

  //if (index == 4)
    phi7  = (
		   FracMomLarge(2.0*av+5.0/2.0, 2.0*as    , moments) * FracMomLarge(      -1.0/2.0,        0.0, moments) +
	     2.0 * FracMomLarge(    av+5.0/2.0,     as    , moments) * FracMomLarge(    av-1.0/2.0,     as    , moments) +
		   FracMomLarge(       5.0/2.0,        0.0, moments) * FracMomLarge(2.0*av-1.0/2.0, 2.0*as    , moments)
	    ) + 3.0 * (
		   FracMomLarge(2.0*av+3.0/2.0, 2.0*as    , moments) * FracMomLarge(       1.0/2.0,        0.0, moments) +
	     2.0 * FracMomLarge(    av+3.0/2.0,     as    , moments) * FracMomLarge(    av+1.0/2.0,     as    , moments) +
		   FracMomLarge(       3.0/2.0,        0.0, moments) * FracMomLarge(2.0*av+1.0/2.0, 2.0*as    , moments)
	    ) + 3.0 * (
		   FracMomLarge(2.0*av+1.0/2.0, 2.0*as    , moments) * FracMomLarge(       3.0/2.0,        0.0, moments) +
	     2.0 * FracMomLarge(    av+1.0/2.0,     as    , moments) * FracMomLarge(    av+3.0/2.0,     as    , moments) +
		   FracMomLarge(       1.0/2.0,        0.0, moments) * FracMomLarge(2.0*av+3.0/2.0, 2.0*as    , moments)
	    ) + (
		   FracMomLarge(2.0*av-1.0/2.0, 2.0*as    , moments) * FracMomLarge(       5.0/2.0,        0.0, moments) +
	     2.0 * FracMomLarge(    av-1.0/2.0,     as    , moments) * FracMomLarge(    av+5.0/2.0,     as    , moments) +
		   FracMomLarge(      -1.0/2.0,        0.0, moments) * FracMomLarge(2.0*av+5.0/2.0, 2.0*as    , moments)
	    );

  //if (index == 4 || index == 5 || index == 6)
    phi8  = (
		   FracMomLarge(2.0*av+3.0/2.0, 2.0*as+1.0, moments) * FracMomLarge(      -1.0/2.0,        0.0, moments) +
	     2.0 * FracMomLarge(    av+3.0/2.0,     as+1.0, moments) * FracMomLarge(    av-1.0/2.0,     as    , moments) +
		   FracMomLarge(       3.0/2.0,        1.0, moments) * FracMomLarge(2.0*av-1.0/2.0, 2.0*as    , moments)
	    ) + 2.0 * (
		   FracMomLarge(2.0*av+1.0/2.0, 2.0*as+1.0, moments) * FracMomLarge(       1.0/2.0,        0.0, moments) +
	     2.0 * FracMomLarge(    av+1.0/2.0,     as+1.0, moments) * FracMomLarge(    av+1.0/2.0,     as    , moments) +
		   FracMomLarge(       1.0/2.0,        1.0, moments) * FracMomLarge(2.0*av+1.0/2.0, 2.0*as    , moments)
	    ) + (
		   FracMomLarge(2.0*av-1.0/2.0, 2.0*as+1.0, moments) * FracMomLarge(       3.0/2.0,        0.0, moments) +
	     2.0 * FracMomLarge(    av-1.0/2.0,     as+1.0, moments) * FracMomLarge(    av+3.0/2.0,     as    , moments) +
		   FracMomLarge(      -1.0/2.0,        1.0, moments) * FracMomLarge(2.0*av+3.0/2.0, 2.0*as    , moments)
	    ) + (
		   FracMomLarge(2.0*av+3.0/2.0, 2.0*as    , moments) * FracMomLarge(      -1.0/2.0,        1.0, moments) +
	     2.0 * FracMomLarge(    av+3.0/2.0,     as    , moments) * FracMomLarge(    av-1.0/2.0,     as+1.0, moments) +
		   FracMomLarge(       3.0/2.0,        0.0, moments) * FracMomLarge(2.0*av-1.0/2.0, 2.0*as+1.0, moments)
	    ) + 2.0 * (
		   FracMomLarge(2.0*av+1.0/2.0, 2.0*as    , moments) * FracMomLarge(       1.0/2.0,        1.0, moments) +
	     2.0 * FracMomLarge(    av+1.0/2.0,     as    , moments) * FracMomLarge(    av+1.0/2.0,     as+1.0, moments) +
		   FracMomLarge(       1.0/2.0,        0.0, moments) * FracMomLarge(2.0*av+1.0/2.0, 2.0*as+1.0, moments)
	    ) + (
		   FracMomLarge(2.0*av-1.0/2.0, 2.0*as    , moments) * FracMomLarge(       3.0/2.0,        1.0, moments) +
	     2.0 * FracMomLarge(    av-1.0/2.0,     as    , moments) * FracMomLarge(    av+3.0/2.0,     as+1.0, moments) +
		   FracMomLarge(      -1.0/2.0,        0.0, moments) * FracMomLarge(2.0*av+3.0/2.0, 2.0*as+1.0, moments)
	    );

  //if (index == 6 || index == 7 || index == 8)
    phi9  = (
		   FracMomLarge(2.0*av+1.0/2.0, 2.0*as+2.0, moments) * FracMomLarge(      -1.0/2.0,        0.0, moments) +
	     2.0 * FracMomLarge(    av+1.0/2.0,     as+2.0, moments) * FracMomLarge(    av-1.0/2.0,     as    , moments) +
		   FracMomLarge(       1.0/2.0,        2.0, moments) * FracMomLarge(2.0*av-1.0/2.0, 2.0*as    , moments)
	    ) + (
		   FracMomLarge(2.0*av-1.0/2.0, 2.0*as+2.0, moments) * FracMomLarge(       1.0/2.0,        0.0, moments) +
	     2.0 * FracMomLarge(    av-1.0/2.0,     as+2.0, moments) * FracMomLarge(    av+1.0/2.0,     as    , moments) +
		   FracMomLarge(      -1.0/2.0,        2.0, moments) * FracMomLarge(2.0*av+1.0/2.0, 2.0*as    , moments)
	    ) + 2.0 * (
		   FracMomLarge(2.0*av+1.0/2.0, 2.0*as+1.0, moments) * FracMomLarge(      -1.0/2.0,        1.0, moments) +
	     2.0 * FracMomLarge(    av+1.0/2.0,     as+1.0, moments) * FracMomLarge(    av-1.0/2.0,     as+1.0, moments) +
		   FracMomLarge(       1.0/2.0,        1.0, moments) * FracMomLarge(2.0*av-1.0/2.0, 2.0*as+1.0, moments)
	    ) + 2.0 * (
		   FracMomLarge(2.0*av-1.0/2.0, 2.0*as+1.0, moments) * FracMomLarge(       1.0/2.0,        1.0, moments) +
	     2.0 * FracMomLarge(    av-1.0/2.0,     as+1.0, moments) * FracMomLarge(    av+1.0/2.0,     as+1.0, moments) +
		   FracMomLarge(      -1.0/2.0,        1.0, moments) * FracMomLarge(2.0*av+1.0/2.0, 2.0*as+1.0, moments)
	    ) + (
		   FracMomLarge(2.0*av+1.0/2.0, 2.0*as    , moments) * FracMomLarge(      -1.0/2.0,        2.0, moments) +
	     2.0 * FracMomLarge(    av+1.0/2.0,     as    , moments) * FracMomLarge(    av-1.0/2.0,     as+2.0, moments) +
		   FracMomLarge(       1.0/2.0,        0.0, moments) * FracMomLarge(2.0*av-1.0/2.0, 2.0*as+2.0, moments)
	    ) + (
		   FracMomLarge(2.0*av-1.0/2.0, 2.0*as    , moments) * FracMomLarge(       1.0/2.0,        2.0, moments) +
	     2.0 * FracMomLarge(    av-1.0/2.0,     as    , moments) * FracMomLarge(    av+1.0/2.0,     as+2.0, moments) +
		   FracMomLarge(      -1.0/2.0,        0.0, moments) * FracMomLarge(2.0*av+1.0/2.0, 2.0*as+2.0, moments)
	    );

  //if (index == 8)
    phi10 = (
		   FracMomLarge(2.0*av-1.0/2.0, 2.0*as+3.0, moments) * FracMomLarge(      -1.0/2.0,        0.0, moments) +
	     2.0 * FracMomLarge(    av-1.0/2.0,     as+3.0, moments) * FracMomLarge(    av-1.0/2.0,     as    , moments) +
		   FracMomLarge(      -1.0/2.0,        3.0, moments) * FracMomLarge(2.0*av-1.0/2.0, 2.0*as    , moments)
	    ) + 3.0 * (
		   FracMomLarge(2.0*av-1.0/2.0, 2.0*as+2.0, moments) * FracMomLarge(      -1.0/2.0,        1.0, moments) +
	     2.0 * FracMomLarge(    av-1.0/2.0,     as+2.0, moments) * FracMomLarge(    av-1.0/2.0,     as+1.0, moments) +
		   FracMomLarge(      -1.0/2.0,        2.0, moments) * FracMomLarge(2.0*av-1.0/2.0, 2.0*as+1.0, moments)
	    ) + 3.0 * (
		   FracMomLarge(2.0*av-1.0/2.0, 2.0*as+1.0, moments) * FracMomLarge(      -1.0/2.0,        2.0, moments) +
	     2.0 * FracMomLarge(    av-1.0/2.0,     as+1.0, moments) * FracMomLarge(    av-1.0/2.0,     as+2.0, moments) +
		   FracMomLarge(      -1.0/2.0,        1.0, moments) * FracMomLarge(2.0*av-1.0/2.0, 2.0*as+2.0, moments)
	    ) + (
		   FracMomLarge(2.0*av-1.0/2.0, 2.0*as    , moments) * FracMomLarge(      -1.0/2.0,        3.0, moments) +
	     2.0 * FracMomLarge(    av-1.0/2.0,     as    , moments) * FracMomLarge(    av-1.0/2.0,     as+3.0, moments) +
		   FracMomLarge(      -1.0/2.0,        0.0, moments) * FracMomLarge(2.0*av-1.0/2.0, 2.0*as+3.0, moments)
	    );

  /*switch (index)
  {
    case 0: return C * pow(phi1, (1.0/2.0)-x-y) * pow(phi2, x+(1.0/2.0)) * pow(phi3 , y);
    case 1: return C * pow(phi2, (3.0/2.0)-x-y) * pow(phi4, x-(1.0/2.0)) * pow(phi5 , y);
    case 2: return C * pow(phi5, x+y-(1.0/2.0)) * pow(phi3, (1.0/2.0)-x) * pow(phi2 , 1.0-y);
    case 3: return C * pow(phi3, (3.0/2.0)-x-y) * pow(phi5, x+(1.0/2.0)) * pow(phi6 , y-1.0);
    case 4: return C * pow(phi4, (5.0/2.0)-x-y) * pow(phi7, x-(3.0/2.0)) * pow(phi8 , y);
    case 5: return C * pow(phi8, x+y-(3.0/2.0)) * pow(phi5, (3.0/2.0)-x) * pow(phi4 , 1.0-y);
    case 6: return C * pow(phi5, (5.0/2.0)-x-y) * pow(phi8, x-(1.0/2.0)) * pow(phi9 , y-1.0);
    case 7: return C * pow(phi9, x+y-(3.0/2.0)) * pow(phi6, (1.0/2.0)-x) * pow(phi5 , 2.0-y);
    case 8: return C * pow(phi6, (5.0/2.0)-x-y) * pow(phi9, x+(1.0/2.0)) * pow(phi10, y-2.0);
  }*/

  double a = pow(phi1 ,   5.0/16.0)*pow(phi2 ,  15.0/16.0)*pow(phi4 , - 5.0/16.0)*pow(phi7 ,   1.0/16.0);
  double b = pow(phi1 , -23.0/24.0)*pow(phi2 ,   7.0/ 8.0)*pow(phi4 ,   1.0/ 8.0)*pow(phi7 , - 1.0/24.0);
  double c = pow(phi1 , -23.0/24.0)/    phi2              *pow(phi3 ,  15.0/ 8.0)*pow(phi4 ,   1.0/ 8.0)*pow(phi5 ,   5.0/ 4.0)*pow(phi6 , - 5.0/ 4.0)*pow(phi8 , - 1.0/ 8.0)*pow(phi9 , - 1.0/ 4.0)*pow(phi10,   1.0/ 3.0);
  double d = pow(phi1 ,   3.0/ 4.0)*pow(phi2 , - 7.0/ 4.0)*pow(phi4 ,   5.0/ 4.0)*pow(phi7 , - 1.0/ 4.0);
  double e = pow(phi1 ,   3.0/ 2.0)*pow(phi2 , - 3.0/ 2.0)*pow(phi3 , - 2.0     )*pow(phi5 ,   2.0     )*pow(phi6 ,   1.0/ 2.0)*pow(phi9 , - 1.0/ 2.0);
  double f = pow(phi1 ,   3.0/ 4.0)*pow(phi2 ,   1.0/ 4.0)*pow(phi3 , - 2.0     )*pow(phi5 , - 1.0/ 2.0)*pow(phi6 ,   7.0/ 4.0)*pow(phi9 ,   1.0/ 4.0)*pow(phi10, - 1.0/ 2.0);
  double g = pow(phi1 , - 1.0/ 6.0)*pow(phi2 ,   1.0/ 2.0)*pow(phi4 , - 1.0/ 2.0)*pow(phi7 ,   1.0/ 6.0);
  double h = pow(phi1 , - 1.0/ 2.0)*    phi2              *pow(phi3 ,   1.0/ 2.0)*pow(phi4 , - 1.0/ 2.0)*pow(phi5 , - 1.0     )*pow(phi8 ,   1.0/ 2.0);
  double i = pow(phi1 , - 1.0/ 2.0)*pow(phi2 ,   1.0/ 2.0)*    phi3              *pow(phi5 , - 1.0     )*pow(phi6 , - 1.0/ 2.0)*pow(phi9 ,   1.0/ 2.0);
  double j = pow(phi1 , - 1.0/ 6.0)*pow(phi3 ,   1.0/ 2.0)*pow(phi6 , - 1.0/ 2.0)*pow(phi10,   1.0/ 6.0);

  return C*a*pow(b,x)*pow(c,y)*pow(d,x*x)*pow(e,x*y)*pow(f,y*y)*pow(g,x*x*x)*pow(h,x*x*y)*pow(i,x*y*y)*pow(j,y*y*y);
}

double TSoot::GetPsiLL(double x, double y, double a, double b, double * moments, double temp)
{
  double C = GetC(temp);

  double psi1, psi2, psi3, psi4, psi5, psi6;

  double Df = 1.8;
  double av = 1.0 - (2.0 / Df);
  double as = (3.0 / Df) - 1.0;

  psi1 = (
	        FracMomLarge(2.0*av-1.0/2.0+x, 2.0*as+y, moments) * FracMomLarge(      -1.0/2.0+a,        b, moments) +
	  2.0 * FracMomLarge(    av-1.0/2.0+x,     as+y, moments) * FracMomLarge(    av-1.0/2.0+a,     as+b, moments) +
	        FracMomLarge(      -1.0/2.0+x,        y, moments) * FracMomLarge(2.0*av-1.0/2.0+a, 2.0*as+b, moments)
	 );

  psi2 = (
	        FracMomLarge(2.0*av+1.0/2.0+x, 2.0*as+y, moments) * FracMomLarge(      -1.0/2.0+a,        b, moments) +
	  2.0 * FracMomLarge(    av+1.0/2.0+x,     as+y, moments) * FracMomLarge(    av-1.0/2.0+a,     as+b, moments) +
	        FracMomLarge(       1.0/2.0+x,        y, moments) * FracMomLarge(2.0*av-1.0/2.0+a, 2.0*as+b, moments)
	 ) + (
	        FracMomLarge(2.0*av-1.0/2.0+x, 2.0*as+y, moments) * FracMomLarge(       1.0/2.0+a,        b, moments) +
	  2.0 * FracMomLarge(    av-1.0/2.0+x,     as+y, moments) * FracMomLarge(    av+1.0/2.0+a,     as+b, moments) +
	        FracMomLarge(      -1.0/2.0+x,        y, moments) * FracMomLarge(2.0*av+1.0/2.0+a, 2.0*as+b, moments)
	 );

  psi3 = (
	        FracMomLarge(2.0*av+3.0/2.0+x, 2.0*as+y, moments) * FracMomLarge(      -1.0/2.0+a,        b, moments) +
	  2.0 * FracMomLarge(    av+3.0/2.0+x,     as+y, moments) * FracMomLarge(    av-1.0/2.0+a,     as+b, moments) +
	        FracMomLarge(       3.0/2.0+x,        y, moments) * FracMomLarge(2.0*av-1.0/2.0+a, 2.0*as+b, moments)
	 ) + 2.0 * (
	        FracMomLarge(2.0*av+1.0/2.0+x, 2.0*as+y, moments) * FracMomLarge(       1.0/2.0+a,        b, moments) +
	  2.0 * FracMomLarge(    av+1.0/2.0+x,     as+y, moments) * FracMomLarge(    av+1.0/2.0+a,     as+b, moments) +
	        FracMomLarge(       1.0/2.0+x,        y, moments) * FracMomLarge(2.0*av+1.0/2.0+a, 2.0*as+b, moments)
	 ) + (
	        FracMomLarge(2.0*av-1.0/2.0+x, 2.0*as+y, moments) * FracMomLarge(       3.0/2.0+a,        b, moments) +
	  2.0 * FracMomLarge(    av-1.0/2.0+x,     as+y, moments) * FracMomLarge(    av+3.0/2.0+a,     as+b, moments) +
	        FracMomLarge(      -1.0/2.0+x,        y, moments) * FracMomLarge(2.0*av+3.0/2.0+a, 2.0*as+b, moments)
	 );

  if (fNSootMoments > 4)
  psi4 = (
	        FracMomLarge(2.0*av+5.0/2.0+x, 2.0*as+y, moments) * FracMomLarge(      -1.0/2.0+a,        b, moments) +
	  2.0 * FracMomLarge(    av+5.0/2.0+x,     as+y, moments) * FracMomLarge(    av-1.0/2.0+a,     as+b, moments) +
	        FracMomLarge(       5.0/2.0+x,        y, moments) * FracMomLarge(2.0*av-1.0/2.0+a, 2.0*as+b, moments)
	 ) + 3.0 * (
	        FracMomLarge(2.0*av+3.0/2.0+x, 2.0*as+y, moments) * FracMomLarge(       1.0/2.0+a,        b, moments) +
	  2.0 * FracMomLarge(    av+3.0/2.0+x,     as+y, moments) * FracMomLarge(    av+1.0/2.0+a,     as+b, moments) +
	        FracMomLarge(       3.0/2.0+x,        y, moments) * FracMomLarge(2.0*av+1.0/2.0+a, 2.0*as+b, moments)
	 ) + 3.0 * (
	        FracMomLarge(2.0*av+1.0/2.0+x, 2.0*as+y, moments) * FracMomLarge(       3.0/2.0+a,        b, moments) +
	  2.0 * FracMomLarge(    av+1.0/2.0+x,     as+y, moments) * FracMomLarge(    av+3.0/2.0+a,     as+b, moments) +
	        FracMomLarge(       1.0/2.0+x,        y, moments) * FracMomLarge(2.0*av+3.0/2.0+a, 2.0*as+b, moments)
	 ) + (
	        FracMomLarge(2.0*av-1.0/2.0+x, 2.0*as+y, moments) * FracMomLarge(       5.0/2.0+a,        b, moments) +
	  2.0 * FracMomLarge(    av-1.0/2.0+x,     as+y, moments) * FracMomLarge(    av+5.0/2.0+a,     as+b, moments) +
	        FracMomLarge(      -1.0/2.0+x,        y, moments) * FracMomLarge(2.0*av+5.0/2.0+a, 2.0*as+b, moments)
	 );

  if (fNSootMoments > 7)
  psi5 = (
	        FracMomLarge(2.0*av+7.0/2.0+x, 2.0*as+y, moments) * FracMomLarge(      -1.0/2.0+a,        b, moments) +
	  2.0 * FracMomLarge(    av+7.0/2.0+x,     as+y, moments) * FracMomLarge(    av-1.0/2.0+a,     as+b, moments) +
	        FracMomLarge(       7.0/2.0+x,        y, moments) * FracMomLarge(2.0*av-1.0/2.0+a, 2.0*as+b, moments)
	 ) + 4.0 * (
	        FracMomLarge(2.0*av+5.0/2.0+x, 2.0*as+y, moments) * FracMomLarge(       1.0/2.0+a,        b, moments) +
	  2.0 * FracMomLarge(    av+5.0/2.0+x,     as+y, moments) * FracMomLarge(    av+1.0/2.0+a,     as+b, moments) +
	        FracMomLarge(       5.0/2.0+x,        y, moments) * FracMomLarge(2.0*av+1.0/2.0+a, 2.0*as+b, moments)
	 ) + 6.0 * (
	        FracMomLarge(2.0*av+3.0/2.0+x, 2.0*as+y, moments) * FracMomLarge(       3.0/2.0+a,        b, moments) +
	  2.0 * FracMomLarge(    av+3.0/2.0+x,     as+y, moments) * FracMomLarge(    av+3.0/2.0+a,     as+b, moments) +
	        FracMomLarge(       3.0/2.0+x,        y, moments) * FracMomLarge(2.0*av+3.0/2.0+a, 2.0*as+b, moments)
	 ) + 4.0 * (
	        FracMomLarge(2.0*av+1.0/2.0+x, 2.0*as+y, moments) * FracMomLarge(       5.0/2.0+a,        b, moments) +
	  2.0 * FracMomLarge(    av+1.0/2.0+x,     as+y, moments) * FracMomLarge(    av+5.0/2.0+a,     as+b, moments) +
	        FracMomLarge(       1.0/2.0+x,        y, moments) * FracMomLarge(2.0*av+5.0/2.0+a, 2.0*as+b, moments)
	 ) + (
	        FracMomLarge(2.0*av-1.0/2.0+x, 2.0*as+y, moments) * FracMomLarge(       7.0/2.0+a,        b, moments) +
	  2.0 * FracMomLarge(    av-1.0/2.0+x,     as+y, moments) * FracMomLarge(    av+7.0/2.0+a,     as+b, moments) +
	        FracMomLarge(      -1.0/2.0+x,        y, moments) * FracMomLarge(2.0*av+7.0/2.0+a, 2.0*as+b, moments)
	 );

  if (fNSootMoments > 11)
  psi6 = (
	        FracMomLarge(2.0*av+9.0/2.0+x, 2.0*as+y, moments) * FracMomLarge(      -1.0/2.0+a,        b, moments) +
	  2.0 * FracMomLarge(    av+9.0/2.0+x,     as+y, moments) * FracMomLarge(    av-1.0/2.0+a,     as+b, moments) +
	        FracMomLarge(       9.0/2.0+x,        y, moments) * FracMomLarge(2.0*av-1.0/2.0+a, 2.0*as+b, moments)
	 ) + 5.0 * (
	        FracMomLarge(2.0*av+7.0/2.0+x, 2.0*as+y, moments) * FracMomLarge(       1.0/2.0+a,        b, moments) +
	  2.0 * FracMomLarge(    av+7.0/2.0+x,     as+y, moments) * FracMomLarge(    av+1.0/2.0+a,     as+b, moments) +
	        FracMomLarge(       7.0/2.0+x,        y, moments) * FracMomLarge(2.0*av+1.0/2.0+a, 2.0*as+b, moments)
	 ) + 10.0 * (
	        FracMomLarge(2.0*av+5.0/2.0+x, 2.0*as+y, moments) * FracMomLarge(       3.0/2.0+a,        b, moments) +
	  2.0 * FracMomLarge(    av+5.0/2.0+x,     as+y, moments) * FracMomLarge(    av+3.0/2.0+a,     as+b, moments) +
	        FracMomLarge(       5.0/2.0+x,        y, moments) * FracMomLarge(2.0*av+3.0/2.0+a, 2.0*as+b, moments)
	 ) + 10.0 * (
	        FracMomLarge(2.0*av+3.0/2.0+x, 2.0*as+y, moments) * FracMomLarge(       5.0/2.0+a,        b, moments) +
	  2.0 * FracMomLarge(    av+3.0/2.0+x,     as+y, moments) * FracMomLarge(    av+5.0/2.0+a,     as+b, moments) +
	        FracMomLarge(       3.0/2.0+x,        y, moments) * FracMomLarge(2.0*av+5.0/2.0+a, 2.0*as+b, moments)
	 ) + 5.0 * (
	        FracMomLarge(2.0*av+1.0/2.0+x, 2.0*as+y, moments) * FracMomLarge(       7.0/2.0+a,        b, moments) +
	  2.0 * FracMomLarge(    av+1.0/2.0+x,     as+y, moments) * FracMomLarge(    av+7.0/2.0+a,     as+b, moments) +
	        FracMomLarge(       1.0/2.0+x,        y, moments) * FracMomLarge(2.0*av+7.0/2.0+a, 2.0*as+b, moments)
	 ) + (
	        FracMomLarge(2.0*av-1.0/2.0+x, 2.0*as+y, moments) * FracMomLarge(       9.0/2.0+a,        b, moments) +
	  2.0 * FracMomLarge(    av-1.0/2.0+x,     as+y, moments) * FracMomLarge(    av+9.0/2.0+a,     as+b, moments) +
	        FracMomLarge(      -1.0/2.0+x,        y, moments) * FracMomLarge(2.0*av+9.0/2.0+a, 2.0*as+b, moments)
	 );

  /*cerr << "Psi Large-Large\t" << a << "\t" << b << "\t" << x << "\t" << y << endl;
  cerr << "Two Terms: " << C * sqrt(psi1*psi2) << endl;
  cerr << "Three Terms: " << C * pow(psi1,3.0/8.0)*pow(psi2,3.0/4.0)*pow(psi3,-1.0/8.0) << endl;
  cerr << "Four Terms: " << C * pow(psi1,5.0/16.0)*pow(psi2,15.0/16.0)*pow(psi3,-5.0/16.0)*pow(psi4,1.0/16.0) << endl;*/

  return C * sqrt(psi1 * psi2);
  if (fNSootMoments == 4)
    return C * pow(psi1,3.0/8.0)*pow(psi2,3.0/4.0)*pow(psi3,-1.0/8.0);
  else if (fNSootMoments == 7)
    return C * pow(psi1,5.0/16.0)*pow(psi2,15.0/16.0)*pow(psi3,-5.0/16.0)*pow(psi4,1.0/16.0);
  else if (fNSootMoments == 11)
    return C * pow(psi1,35.0/128.0)*pow(psi2,35.0/32.0)*pow(psi3,-35.0/64.0)*pow(psi4,7.0/32.0)*pow(psi5,-5.0/128.0);
  else if (fNSootMoments == 16)
    return C * pow(psi1,63.0/256.0)*pow(psi2,315.0/256.0)*pow(psi3,-105.0/128.0)*pow(psi4,63.0/128.0)*pow(psi5,-45.0/256.0)*pow(psi6,7.0/256.0);

  return 0.0;
}

double TSoot::CoagContSL(double i, double j, double k, double * moments, double temp)
{
#ifdef AGG
  if ((i < 1.01 && i > 0.99 && j == 0 && k == 0) || (i == 0 && j < 1.01 && j > 0.99 && k == 0) || (i == 0 && j == 0 && k < 1.01 && k > 0.99))
    return 0.0;
#endif

  double fm;
  double cont;

  fm = -GetPsiSL(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, moments, temp);

  double V0 = 2.0*dimer_nbrC2;

  double Df = 1.8;
  double av = 1.0 - (2.0 / Df);
  double as = (3.0 / Df) - 1.0;

  double lambda = 3.0*mu/rho*sqrt(PI*Wmix/(8.0*RGAS*temp)) / pow(6.0*fMolarMassSoot/(PI*fSootDensity*AVOGADRO),1.0/3.0);
  double C_cont = 8.0 * RGAS/AVOGADRO * temp / (3.0*mu) * AVOGADRO;

  cont = -C_cont * (2.000*                        moments[fNSootMoments-1]*FracMomLarge( 0.0   , 0.0   ,0.0,moments) + 
		                 pow(V0, 1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-    av,-    as,0.0,moments) +
		                 pow(V0,-1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av,     as,0.0,moments) +
		    1.257*lambda*pow(V0,-1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge( 0.0   , 0.0   ,0.0,moments) +
		    1.257*lambda*                 moments[fNSootMoments-1]*FracMomLarge(-    av,-    as,0.0,moments) +
		    1.257*lambda*pow(V0,-2.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av,     as,0.0,moments) +
		    1.257*lambda*pow(V0, 1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-2.0*av,-2.0*as,0.0,moments));

  return (fm * cont) / (fm + cont);
}

double TSoot::GetPsiSL(double u, double v, double w, double x, double y, double z, double * moments, double temp)
{
  double C = GetC(temp);
  double psi1, psi2;

  double V0 = 2.0*dimer_nbrC2;
  double H0 = 2.0*dimer_nbrH;

  double Df = 1.8;
  double av = 1.0 - (2.0 / Df);
  double as = (3.0 / Df) - 1.0;

  psi1 = (
                pow(V0, 1.0/6.0+u+(2.0/3.0)*v) * pow(H0, w) * moments[fNSootMoments-1] * FracMomLarge(      -1.0/2.0+x,        y, z, moments) +
          2.0 * pow(V0,-1.0/6.0+u+(2.0/3.0)*v) * pow(H0, w) * moments[fNSootMoments-1] * FracMomLarge(    av-1.0/2.0+x,     as+y, z, moments) +
                pow(V0,-3.0/6.0+u+(2.0/3.0)*v) * pow(H0, w) * moments[fNSootMoments-1] * FracMomLarge(2.0*av-1.0/2.0+x, 2.0*as+y, z, moments)
         );

  psi2 = (
                pow(V0, 7.0/6.0+u+(2.0/3.0)*v) * pow(H0, w) * moments[fNSootMoments-1] * FracMomLarge(      -1.0/2.0+x,        y, z, moments) +
          2.0 * pow(V0, 5.0/6.0+u+(2.0/3.0)*v) * pow(H0, w) * moments[fNSootMoments-1] * FracMomLarge(    av-1.0/2.0+x,     as+y, z, moments) +
                pow(V0, 3.0/6.0+u+(2.0/3.0)*v) * pow(H0, w) * moments[fNSootMoments-1] * FracMomLarge(2.0*av-1.0/2.0+x, 2.0*as+y, z, moments)
         ) + (
	        pow(V0, 1.0/6.0+u+(2.0/3.0)*v) * pow(H0, w) * moments[fNSootMoments-1] * FracMomLarge(       1.0/2.0+x,        y, z, moments) +
	  2.0 * pow(V0,-1.0/6.0+u+(2.0/3.0)*v) * pow(H0, w) * moments[fNSootMoments-1] * FracMomLarge(    av+1.0/2.0+x,     as+y, z, moments) +
	        pow(V0,-3.0/6.0+u+(2.0/3.0)*v) * pow(H0, w) * moments[fNSootMoments-1] * FracMomLarge(2.0*av+1.0/2.0+x, 2.0*as+y, z, moments)
         );

  return C * sqrt(psi1 * psi2);
}

double TSoot::GetPsiLL(double u, double v, double w, double x, double y, double z, double * moments, double temp)
{
  double C = GetC(temp);
  double psi1, psi2;

  double Df = 1.8;
  double av = 1.0 - (2.0 / Df);
  double as = (3.0 / Df) - 1.0;

  psi1 = (
                FracMomLarge(2.0*av-1.0/2.0+u, 2.0*as+v, w, moments) * FracMomLarge(      -1.0/2.0+x,        y, z, moments) +
          2.0 * FracMomLarge(    av-1.0/2.0+u,     as+v, w, moments) * FracMomLarge(    av-1.0/2.0+x,     as+y, z, moments) +
                FracMomLarge(      -1.0/2.0+u,        v, w, moments) * FracMomLarge(2.0*av-1.0/2.0+x, 2.0*as+y, z, moments)
         );

  psi2 = (
                FracMomLarge(2.0*av+1.0/2.0+u, 2.0*as+v, w, moments) * FracMomLarge(      -1.0/2.0+x,        y, z, moments) +
          2.0 * FracMomLarge(    av+1.0/2.0+u,     as+v, w, moments) * FracMomLarge(    av-1.0/2.0+x,     as+y, z, moments) +
                FracMomLarge(       1.0/2.0+u,        v, w, moments) * FracMomLarge(2.0*av-1.0/2.0+x, 2.0*as+y, z, moments)
         ) + (
                FracMomLarge(2.0*av-1.0/2.0+u, 2.0*as+v, w, moments) * FracMomLarge(       1.0/2.0+x,        y, z, moments) +
          2.0 * FracMomLarge(    av-1.0/2.0+u,     as+v, w, moments) * FracMomLarge(    av+1.0/2.0+x,     as+y, z, moments) +
                FracMomLarge(      -1.0/2.0+u,        v, w, moments) * FracMomLarge(2.0*av+1.0/2.0+x, 2.0*as+y, z, moments)
         );

  return C * sqrt(psi1 * psi2);
}

double TSoot::CoagContLL(double i, double j, double k, double * moments, double temp)
{
#ifdef AGG
  if ((i < 1.01 && i > 0.99 && j == 0 && k == 0) || (i == 0 && j < 1.01 && j > 0.99 && k == 0) || (i == 0 && j == 0 && k < 1.01 && k > 0.99))
    return 0.0;
#endif

  double fm;
  double cont;

  fm = -0.5 * GetPsiLL(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, moments, temp);

  double Df = 1.8;
  double av = 1.0 - (2.0 / Df);
  double as = (3.0 / Df) - 1.0;

  double lambda = 3.0*mu/rho*sqrt(PI*Wmix/(8.0*RGAS*temp)) / pow(6.0*fMolarMassSoot/(PI*fSootDensity*AVOGADRO),1.0/3.0);
  double C_cont = 8.0 * RGAS/AVOGADRO * temp / (3.0*mu) * AVOGADRO;

  cont = -0.5 * C_cont * (2.000*FracMomLarge(0,0,0,moments)*FracMomLarge(0,0,0,moments) + 
			  FracMomLarge(av,as,0,moments)*FracMomLarge(-av,-as,0,moments) +
			  FracMomLarge(-av,-as,0,moments)*FracMomLarge(av,as,0,moments) +
			  1.257*lambda*FracMomLarge(-av,-as,0,moments)*FracMomLarge(0,0,0,moments) +
			  1.257*lambda*FracMomLarge(0,0,0,moments)*FracMomLarge(-av,-as,0,moments) +
			  1.257*lambda*FracMomLarge(-2.0*av,-2.0*as,0,moments)*FracMomLarge(av,as,0,moments) +
			  1.257*lambda*FracMomLarge(av,as,0,moments)*FracMomLarge(-2.0*av,-2.0*as,0,moments));

  return (fm * cont) / (fm + cont);
}



Double TSoot::CoagulationSource(Double i, Double j, Double k, Double temp, Double * moments)
{
#ifdef HMOM
  double ss;
  double ssfm;
  double sscont;

  double sl;
  double ll;

  double V0 = 2.0 * dimer_nbrC2;
  double H0 = 2.0 * dimer_nbrH;

  double C_fm = GetC(temp);

  double C_cont = 8.0 * RGAS/AVOGADRO * temp / (3.0*mu) * AVOGADRO;
  double lambda = 3.0*mu/rho*sqrt(PI*Wmix/(8.0*RGAS*temp)) / pow(6.0*fMolarMassSoot/(PI*fSootDensity*AVOGADRO),1.0/3.0);

  double FitC = 2.0/3.0;
  double FitExp = -0.2043;

 #ifdef AGG
  if (i < 1.01 && i > 0.99 && j == 0 && k == 0)
    return 0.0;
  if (i == 0 && j < 1.01 && j > 0.99 && k == 0)
    return 0.0;
  if (i == 0 && j == 0 && k < 1.01 && k > 0.99)
    return 0.0;
 #endif

  //Contribution of source term due to collision of small particles with small particles
  ssfm = C_fm * pow(2.0, 2.5) * (pow(2.0, i+j+k-1.0) - 1.0) * 
    pow(V0, i+(2.0/3.0)*j+(1.0/6.0)) * pow(H0, k) * moments[fNSootMoments-1] * moments[fNSootMoments-1];
  sscont = 4.0 * C_cont * (1 + 1.257*lambda*pow(V0,-1.0/3.0)) * (pow(2.0, i+j+k-1.0) - 1.0) * 
    pow(V0, i+(2.0/3.0)*j) * pow(H0, k) * moments[fNSootMoments-1] * moments[fNSootMoments-1];

  //ss = (ssfm*sscont)/(ssfm+sscont);
  ss = ssfm;

  //Contribution of source term due to collision of small particles with large particles
  sl = CoagContSL(i,j,k,moments,temp);

  //Contribution of source term due to collision of large particles with large particles
  ll = CoagContLL(i,j,k,moments,temp);

  return ss + sl + ll;
#endif
}

//----------Condensation----------//
double TSoot::CondensationSource(double i, double temp, double * moments)
{
  double dV = dimer_nbrC2;
  double C = GetC(temp) / 2.2;

  double S0 =	    FracMom(i-1.0/3.0, moments) * pow(dV, - 1.0/2.0) +
              2.0 * FracMom(i-2.0/3.0, moments) * pow(dV, - 1.0/6.0) + 
		    FracMom(i-    1.0, moments) * pow(dV,   1.0/6.0) + 
	      0.5 * FracMom(i-4.0/3.0, moments) * pow(dV,   1.0/2.0) +
		    FracMom(i-5.0/3.0, moments) * pow(dV,   5.0/6.0) +
	      0.5 * FracMom(i-    2.0, moments) * pow(dV,   7.0/6.0);

  return i * dimer_conc * dV * C * S0;
}

double TSoot::CondensationSource(double i, double j, double temp, double * moments)
{
  if (i == 0 && j == 0)
    return 0.0;

  double C_fm = GetC(temp) / 2.2;
  double C_cont = 8.0 * RGAS/AVOGADRO * temp / (3.0*mu) * AVOGADRO;
  double lambda = 3.0*mu/rho*sqrt(PI*Wmix/(8.0*RGAS*temp)) / pow(6.0*fMolarMassSoot/(PI*fSootDensity*AVOGADRO),1.0/3.0);

  double dV = dimer_nbrC2;

  double Df = 1.8;
  double av = 1.0 - (2.0 / Df);
  double as = (3.0 / Df) - 1.0;

  double FitC = 2.0/3.0;
  double FitExp = -0.2043;

  double fm;
  double cont;

  double S0v =       FracMom(i+           2.0*av-1.0, j+           2.0*as, moments) * pow(dV,  3.0/6.0) +
	       2.0 * FracMom(i+               av-1.0, j+               as, moments) * pow(dV,  5.0/6.0) +
		     FracMom(i                  -1.0, j                  , moments) * pow(dV,  7.0/6.0) +
	       0.5 * FracMom(i+           2.0*av-2.0, j+           2.0*as, moments) * pow(dV,  9.0/6.0) +
		     FracMom(i+               av-2.0, j+               as, moments) * pow(dV, 11.0/6.0) +
               0.5 * FracMom(i                  -2.0, j                  , moments) * pow(dV, 13.0/6.0);

  double S0s =	     FracMom(i-2.0*FitExp+2.0*av-1.0, j+3.0*FitExp+2.0*as, moments) * pow(dV,  3.0/6.0) +
	       2.0 * FracMom(i-2.0*FitExp+    av-1.0, j+3.0*FitExp+    as, moments) * pow(dV,  5.0/6.0) +
		     FracMom(i-2.0*FitExp       -1.0, j+3.0*FitExp       , moments) * pow(dV,  7.0/6.0) +
	       0.5 * FracMom(i-2.0*FitExp+2.0*av-2.0, j+3.0*FitExp+2.0*as, moments) * pow(dV,  9.0/6.0) +
		     FracMom(i-2.0*FitExp+    av-2.0, j+3.0*FitExp+    as, moments) * pow(dV, 11.0/6.0) +
               0.5 * FracMom(i-2.0*FitExp       -2.0, j+3.0*FitExp       , moments) * pow(dV, 13.0/6.0);

  return C_fm * (i * S0v + FitC * j * S0s) * dimer_conc;

//   double S0v_cont = 2.000 *          FracMom(i                  -1.0, j                 ,  moments) * pow(dV, 3.0/3.0) +
// 		                     FracMom(i           +    av-1.0, j           +    as, moments) * pow(dV, 2.0/3.0) +
// 		                     FracMom(i           -    av-1.0, j           -    as, moments) * pow(dV, 4.0/3.0) +
// 		    1.257 * lambda * FracMom(i           -    av-1.0, j           -    as, moments) * pow(dV, 3.0/3.0) +
// 		    1.257 * lambda * FracMom(i                  -1.0, j                  , moments) * pow(dV, 2.0/3.0) +
// 		    1.257 * lambda * FracMom(i           -2.0*av-1.0, j           -2.0*as, moments) * pow(dV, 4.0/3.0) +
// 		    1.257 * lambda * FracMom(i           +    av-1.0, j           +    as, moments) * pow(dV, 1.0/3.0);

//   double S0s_cont = 2.000 *          FracMom(i-2.0*FitExp       -1.0, j+3.0*FitExp       , moments) * pow(dV, 3.0/3.0) +
// 		                     FracMom(i-2.0*FitExp+    av-1.0, j+3.0*FitExp+    as, moments) * pow(dV, 2.0/3.0) +
// 		                     FracMom(i-2.0*FitExp-    av-1.0, j+3.0*FitExp-    as, moments) * pow(dV, 4.0/3.0) +
// 		    1.257 * lambda * FracMom(i-2.0*FitExp-    av-1.0, j+3.0*FitExp-    as, moments) * pow(dV, 3.0/3.0) +
// 		    1.257 * lambda * FracMom(i-2.0*FitExp       -1.0, j+3.0*FitExp       , moments) * pow(dV, 2.0/3.0) +
// 		    1.257 * lambda * FracMom(i-2.0*FitExp-2.0*av-1.0, j+3.0*FitExp-2.0*as, moments) * pow(dV, 4.0/3.0) +
// 		    1.257 * lambda * FracMom(i-2.0*FitExp+    av-1.0, j+3.0*FitExp+    as, moments) * pow(dV, 1.0/3.0);

//   cont = C_cont * (i * S0v_cont + FitC * j * S0s_cont) * dimer_conc;

//   return (fm * cont) / (fm + cont);
}

double TSoot::CondensationSource(double i, double j, double k, double temp, double * moments)
{
  if (i == 0 && j == 0 && k == 0)
    return 0.0;

  double C_fm = GetC(temp) / 2.2;

  double dV = dimer_nbrC2;
  double dH = dimer_nbrH;

  double Df = 1.8;
  double av = 1.0 - (2.0 / Df);
  double as = (3.0 / Df) - 1.0;

  double FitC = 2.0/3.0;
  double FitExp = -0.2043;

  double S0v =       FracMom(i+           2.0*av-1.0, j+           2.0*as, k    , moments) * pow(dV,  3.0/6.0) +
               2.0 * FracMom(i+               av-1.0, j+               as, k    , moments) * pow(dV,  5.0/6.0) +
                     FracMom(i                  -1.0, j                  , k    , moments) * pow(dV,  7.0/6.0) +
               0.5 * FracMom(i+           2.0*av-2.0, j+           2.0*as, k    , moments) * pow(dV,  9.0/6.0) +
                     FracMom(i+               av-2.0, j+               as, k    , moments) * pow(dV, 11.0/6.0) +
               0.5 * FracMom(i                  -2.0, j                  , k    , moments) * pow(dV, 13.0/6.0);

  double S0s =	     FracMom(i-2.0*FitExp+2.0*av-1.0, j+3.0*FitExp+2.0*as, k    , moments) * pow(dV,  3.0/6.0) +
               2.0 * FracMom(i-2.0*FitExp+    av-1.0, j+3.0*FitExp+    as, k    , moments) * pow(dV,  5.0/6.0) +
                     FracMom(i-2.0*FitExp       -1.0, j+3.0*FitExp       , k    , moments) * pow(dV,  7.0/6.0) +
               0.5 * FracMom(i-2.0*FitExp+2.0*av-2.0, j+3.0*FitExp+2.0*as, k    , moments) * pow(dV,  9.0/6.0) +
                     FracMom(i-2.0*FitExp+    av-2.0, j+3.0*FitExp+    as, k    , moments) * pow(dV, 11.0/6.0) +
               0.5 * FracMom(i-2.0*FitExp       -2.0, j+3.0*FitExp       , k    , moments) * pow(dV, 13.0/6.0);

  double S0h =       FracMom(i+           2.0*av    , j+           2.0*as, k-1.0, moments) * pow(dV, -3.0/6.0) * dH +
               2.0 * FracMom(i+               av    , j+               as, k-1.0, moments) * pow(dV, -1.0/6.0) * dH +
		     FracMom(i                      , j                  , k-1.0, moments) * pow(dV,  1.0/6.0) * dH +
	       0.5 * FracMom(i+           2.0*av-1.0, j+           2.0*as, k-1.0, moments) * pow(dV,  3.0/6.0) * dH +
		     FracMom(i+               av-1.0, j+               as, k-1.0, moments) * pow(dV,  5.0/6.0) * dH +
               0.5 * FracMom(i                  -1.0, j                  , k-1.0, moments) * pow(dV,  7.0/6.0) * dH;

  return C_fm * (i * S0v + FitC * j * S0s + k * S0h) * dimer_conc;
}

//----------Surface Growth----------//
double TSoot::SurfaceGrowthSource(double i, double * moments, double * Y, double temp, double density, double * molarMass)
{
  double ksg;

  //if (i == 0.0)
  //  ksg = -SurfaceGrowthCoeffBack(Y, density, molarMass);
  //else
  ksg = SurfaceGrowthCoeffFor(Y, temp, density, molarMass) - SurfaceGrowthCoeffBack(Y, temp, density, molarMass);

  double Chi = fChi * pow(36.0*PI, 1.0/3.0) * pow(fMolarMassSoot / (AVOGADRO * fSootDensity), 2.0/3.0);

  return ksg * Chi * i * FracMom(i-(1.0/3.0), moments);
}

double TSoot::SurfaceGrowthSource(double i, double j, double * moments, double * Y, double temp, double density, double * molarMass)
{
  double FitC = 2.0/3.0;
  double FitExp = -0.2043;

  double ksg;

  //if (i == 0.0)
  //  ksg = -SurfaceGrowthCoeffBack(Y, density, molarMass);
  //else
  ksg = SurfaceGrowthCoeffFor(Y, temp, density, molarMass) - SurfaceGrowthCoeffBack(Y, temp, density, molarMass);

  double Chi = fChi * pow(36.0*PI, 1.0/3.0) * pow(fMolarMassSoot / (AVOGADRO * fSootDensity), 2.0/3.0);

  return ksg * Chi * (i * FracMom(i-1.0, j+1.0, moments) + j * FitC * FracMom(i-1.0-2.0*FitExp, j+1.0+3.0*FitExp, moments));
}

double TSoot::SurfaceGrowthSource(double i, double j, double k, double * moments, double * Y, double temp, double density, double * molarMass)
{
  double FitC = 2.0/3.0;
  double FitExp = -0.2043;

  double ksg;

  //if (i == 0.0)
  //  ksg = -SurfaceGrowthCoeffBack(Y, density, molarMass);
  //else
  ksg = SurfaceGrowthCoeffFor(Y, temp, density, molarMass) - SurfaceGrowthCoeffBack(Y, temp, density, molarMass);

  return ksg * (i * FracMom(i-1.0, j, k+1.0, moments) + j * FitC * FracMom(i-1.0-2.0*FitExp, j+3.0*FitExp, k+1.0, moments));
}

//----------Surface Oxidation----------//
double TSoot::OxidationSource(double i, double * moments, double * Y, double temp, double density, double * molarMass)
{
  //M_dot_x = k_ox * Chi * Sum([V_i-d^x - V_i^x] * V_i^2/3 * N_i)

  double kox = OxidationCoeff(Y, temp, density, molarMass);
  double Chi = fChi * pow(36.0*PI, 1.0/3.0) * pow(fMolarMassSoot / (AVOGADRO * fSootDensity), 2.0/3.0);

  return - i * kox * Chi * FracMom(i-(1.0/3.0), moments);
}

double TSoot::OxidationSource(double i, double j, double * moments, double * Y, double temp, double density, double * molarMass)
{
  //M_dot_xy = k_ox * Chi * Sum([V_i-d^x * S_k-d^y - V_i^x * S_k^y] * S_k * N_ik)

  double kox = OxidationCoeff(Y, temp, density, molarMass);
  double Chi = fChi * pow(36.0*PI, 1.0/3.0) * pow(fMolarMassSoot / (AVOGADRO * fSootDensity), 2.0/3.0);

#ifdef MOMIC
  return - (i+(2.0/3.0)*j) * kox * Chi * FracMom(i-1.0, j+1.0, moments);
#endif

#ifdef HMOM
  double V0 = 2.0 * dimer_nbrC2;

  double small = - kox * Chi * pow(V0, i-1.0+2.0/3.0*(j+1.0)) * moments[fNSootMoments-1];
  double large = 0.0;

  //OH Oxidation -- np is constant -- THIS IS GOOD
  double k_OH = OxidationCoeff_OH(Y, temp, density, molarMass);
  large += - (i+(2.0/3.0)*j) * kox * Chi * FracMomLarge(i-1.0, j+1.0, moments);

  //O2 Oxidation -- particles become less spherical -- now all -- THIS IS CRAP
  double FitC = 2.0/3.0;
  double FitExp = -0.2043;
  double k_O2 = OxidationCoeff_O2(Y, temp, density, molarMass);
  //large += -kox * Chi * (i * FracMomLarge(i-1.0, j+1.0, moments) + j * FitC * FracMomLarge(i-1.0-2.0*FitExp, j+1.0+3.0*FitExp, moments));

  return small + large;
#endif
}

double TSoot::OxidationRate_O2(double * moments, double * Y, double temp, double density, double * molarMass)
{
  //Retursn M_dot_10 / M_10

  double k_O2 = OxidationCoeff_O2(Y, temp, density, molarMass);
  double Chi = fChi * pow(36.0*PI, 1.0/3.0) * pow(fMolarMassSoot / (AVOGADRO * fSootDensity), 2.0/3.0);

  double V0 = 2.0 * dimer_nbrC2;

  double small = - k_O2 * Chi * pow(V0, 2.0/3.0) * moments[fNSootMoments-1];
  double large = - k_O2 * Chi * FracMomLarge(0.0, 1.0, moments);

  return (small+large) / moments[1];
}

double TSoot::OxidationRate_OH(double * moments, double * Y, double temp, double density, double * molarMass)
{
  //Retursn M_dot_10 / M_10

  double k_OH = OxidationCoeff_OH(Y, temp, density, molarMass);
  double Chi = fChi * pow(36.0*PI, 1.0/3.0) * pow(fMolarMassSoot / (AVOGADRO * fSootDensity), 2.0/3.0);

  double V0 = 2.0 * dimer_nbrC2;

  double small = - k_OH * Chi * pow(V0, 2.0/3.0) * moments[fNSootMoments-1];
  double large = - k_OH * Chi * FracMomLarge(0.0, 1.0, moments);

  return (small+large) / moments[1];
}

double TSoot::OxidationSource(double i, double j, double k, double * moments, double * Y, double temp, double density, double * molarMass)
{
  //M_dot_xy = k_ox * Sum([V_i-d^x * S_k-d^y * H_m-d^z - V_i^x * S_k^y * H_m^z] * H_m * N_ikm)
  
  double kox = OxidationCoeff(Y, temp, density, molarMass);
	
  return - (i+(2.0/3.0)*(j+k)) * kox * FracMom(i-1.0, j, k+1.0, moments);
}

//----------Aggregate Fragmentation----------//
double TSoot::FragmentationSource(double i, double j, double * moments)
{
  double * w = fSootReactionRate->vec;
  double dV = PI/12.0 * (6.0*fMolarMassSoot)/(PI*AVOGADRO*fSootDensity);
  double FragRate = w[ks5] * fMolarMassSoot / fSootDensity / FracMom(0.0, 1.0, moments) / (dV * AVOGADRO);

  return (pow(2.0, 1.0-i-j) - 1.0) * FragRate * FracMomLarge(i-1.0, j+1.0, moments);
}

//---------Molecular Diffusion--------------//
double T1DSoot::DiffusionSource(double i, double j, double * moments, T1DFlamePtr flame)
{
  TNewtonPtr bt = flame->GetSolver()->bt;
  TFlameNodePtr	flameNode = flame->GetFlameNode();
  NodeInfoPtr nodeInfo = bt->GetNodeInfo();

  double Df = 1.8;
  
  // omega_ij = 1.0
  static double pi = 4.0 * atan( 1.0 );
  static double fact = 1.5 * sqrt( 0.5 / pi );
  double *temp = flameNode->temp;
  double *rho  = flameNode->mixDensity;
  double *mixMolarMass = flameNode->mixMolarMass;

  double diffCoeff = fact*sqrt(RGAS) / (AVOGADRO*pow((6.0*fMolarMassSoot)/(PI*AVOGADRO*fSootDensity), 2.0/3.0))*pow(2.0*dimer_nbrC2,2.0/3.0);
  double diffSootP = diffCoeff * sqrt(temp[kPrev]*mixMolarMass[kPrev]);
  double diffSootC = diffCoeff * sqrt(temp[kCurr]*mixMolarMass[kCurr]);
  double diffSootN = diffCoeff * sqrt(temp[kNext]*mixMolarMass[kNext]);

  // Based on collision diameter
  double Index1 = i - 2.0*(1.0-2.0/Df); 
  double Index2 = j - 2.0*(3.0/Df-1.0); 
  
  double fracMomC = FracMom(Index1, Index2, moments) / rho[kCurr];
  double fracMomP = FracMom(Index1, Index2, &nodeInfo->yPrev[fOffsetMoments]);
  double fracMomN = FracMom(Index1, Index2, &nodeInfo->yNext[fOffsetMoments]);

  if (dt_diff==0.0) return 0.0;
  
  double alpha = 0.1 + 0.9*time_diff/abs(dt_diff);
  double hN,hP;
  if (dt_diff>0)
  {
    hN = nodeInfo->h + (1.0-alpha)*nodeInfo->hm;
    hP = alpha*nodeInfo->hm;
  }
  else
  {
    hN = alpha*nodeInfo->h;
    hP = nodeInfo->hm + (1.0-alpha)*nodeInfo->h;
  }
  double derN = 0.5*(diffSootN+diffSootC) * (fracMomN-fracMomC)/hN;
  double derP = 0.5*(diffSootC+diffSootP) * (fracMomC-fracMomP)/hP;
  
  double src = 2.0*(derN-derP)/(hN+hP);
  
  return src;
}

double T1DSoot::DiffusionSourceSmall(double * moments, T1DFlamePtr flame)
{
  TNewtonPtr bt = flame->GetSolver()->bt;
  TFlameNodePtr	flameNode = flame->GetFlameNode();
  NodeInfoPtr nodeInfo = bt->GetNodeInfo();

  double Df = 1.8;
  
  // omega_ij = 1.0
  static double pi = 4.0 * atan( 1.0 );
  static double fact = 1.5 * sqrt( 0.5 / pi );
  double *temp = flameNode->temp;
  double *rho  = flameNode->mixDensity;
  double *mixMolarMass = flameNode->mixMolarMass;
  
  double diffCoeff = fact*sqrt(RGAS) / (AVOGADRO*pow((6.0*fMolarMassSoot)/(PI*AVOGADRO*fSootDensity), 2.0/3.0));
  double diffSootP = diffCoeff * sqrt(temp[kPrev]*mixMolarMass[kPrev]);
  double diffSootC = diffCoeff * sqrt(temp[kCurr]*mixMolarMass[kCurr]);
  double diffSootN = diffCoeff * sqrt(temp[kNext]*mixMolarMass[kNext]);
  
  // Based on collision diameter
  double Index1 = 0.0 - 2.0*(1.0-2.0/Df); 
  double Index2 = 0.0 - 2.0*(3.0/Df-1.0); 
  
  double V0 = 2.0*dimer_nbrC2;

  double fracMomC = moments[fNSootMoments-1] * pow(V0, Index1+(2.0/3.0)*Index2) / rho[kCurr];
  double fracMomP = nodeInfo->yPrev[fOffsetMoments+fNSootMoments-1] * pow(V0, Index1+(2.0/3.0)*Index2);
  double fracMomN = nodeInfo->yNext[fOffsetMoments+fNSootMoments-1] * pow(V0, Index1+(2.0/3.0)*Index2);

  if (dt_diff==0.0) return 0.0;
  
  double alpha = 0.1 + 0.9*time_diff/abs(dt_diff);
  double hN,hP;
  if (dt_diff>0)
  {
    hN = nodeInfo->h + (1.0-alpha)*nodeInfo->hm;
    hP = alpha*nodeInfo->hm;
  }
  else
  {
    hN = alpha*nodeInfo->h;
    hP = nodeInfo->hm + (1.0-alpha)*nodeInfo->h;
  }
  double derN = 0.5*(diffSootN+diffSootC) * (fracMomN-fracMomC)/hN;
  double derP = 0.5*(diffSootC+diffSootP) * (fracMomC-fracMomP)/hP;
  
  double src = 2.0*(derN-derP)/(hN+hP);
  
  return src;
}


//Collision Kernels (Nucleation, Coagulation, and Condensation)
double TSoot::GetPhi(double x, double * moments, double temp)
{
  double C = GetC(temp);
  double phi1, phi2, phi3, phi4, phi5, phi6, phi7;

  //Phi1 = Phi(-1.0/2.0)
  //Phi2 = Phi(1.0/2.0)
  //Phi3 = Phi(3.0/2.0)
  //Phi4 = Phi(5.0/2.0)
  //Phi5 = Phi(7.0/2.0)
  //Phi6 = Phi(9.0/2.0)

  int index;

  if (x < 0.5)
    index = 0;
  else if (x < 1.5)
    index = 1;
  else if (x < 2.5)
    index = 2;
  else if (x < 3.5)
    index = 3;
  else if (x < 4.5)
    index = 4;
  else if (x < 5.5)
    index = 5;
  else
    cerr << "Coagulation not yet implement for more than 6 moments!" << endl;

  phi1 = (
                FracMom( 1.0/6.0, moments) * FracMom(-1.0/2.0, moments) + 
          2.0 * FracMom(-1.0/6.0, moments) * FracMom(-1.0/6.0, moments) +
		FracMom(-1.0/2.0, moments) * FracMom( 1.0/6.0, moments)
         );
  
  phi2 = (
                FracMom( 7.0/6.0, moments) * FracMom(-1.0/2.0, moments) +
	  2.0 * FracMom( 5.0/6.0, moments) * FracMom(-1.0/6.0, moments) +
		FracMom( 1.0/2.0, moments) * FracMom( 1.0/6.0, moments)
         ) + (
	        FracMom( 1.0/6.0, moments) * FracMom( 1.0/2.0, moments) +
	  2.0 * FracMom(-1.0/6.0, moments) * FracMom( 5.0/6.0, moments) +
		FracMom(-1.0/2.0, moments) * FracMom( 7.0/6.0, moments)
	 );
  
  if (index > 0)
    phi3 = (
                  FracMom(13.0/6.0, moments) * FracMom(-1.0/2.0, moments) +
	    2.0 * FracMom(11.0/6.0, moments) * FracMom(-1.0/6.0, moments) +
		  FracMom( 3.0/2.0, moments) * FracMom( 1.0/6.0, moments)
           ) + 2.0 * (
	          FracMom( 7.0/6.0, moments) * FracMom( 1.0/2.0, moments) +
	    2.0 * FracMom( 5.0/6.0, moments) * FracMom( 5.0/6.0, moments) +
		  FracMom( 1.0/2.0, moments) * FracMom( 7.0/6.0, moments)
	   ) + (
	          FracMom( 1.0/6.0, moments) * FracMom( 3.0/2.0, moments) +
	    2.0 * FracMom(-1.0/6.0, moments) * FracMom(11.0/6.0, moments) +
		  FracMom(-1.0/2.0, moments) * FracMom(13.0/6.0, moments)
	   );
  
  if (index > 1)
	phi4 = (
	              FracMom(19.0/6.0, moments) * FracMom(-1.0/2.0, moments) +
	        2.0 * FracMom(17.0/6.0, moments) * FracMom(-1.0/6.0, moments) +
		      FracMom( 5.0/2.0, moments) * FracMom( 1.0/6.0, moments)
	       ) + 3.0 * (
		      FracMom(13.0/6.0, moments) * FracMom( 1.0/2.0, moments) +
		2.0 * FracMom(11.0/6.0, moments) * FracMom( 5.0/6.0, moments) +
		      FracMom( 3.0/2.0, moments) * FracMom( 7.0/6.0, moments)
	       ) + 3.0 * (
		      FracMom( 7.0/6.0, moments) * FracMom( 3.0/2.0, moments) +
		2.0 * FracMom( 5.0/6.0, moments) * FracMom(11.0/6.0, moments) +
		      FracMom( 1.0/2.0, moments) * FracMom(13.0/6.0, moments)
	       ) + (
		      FracMom( 1.0/6.0, moments) * FracMom( 5.0/2.0, moments) +
		2.0 * FracMom(-1.0/6.0, moments) * FracMom(17.0/6.0, moments) +
		      FracMom(-1.0/2.0, moments) * FracMom(19.0/6.0, moments)
	       );

  if (index > 2)
    phi5 = (
                  FracMom(25.0/6.0, moments) * FracMom(-1.0/2.0, moments) +
	    2.0 * FracMom(23.0/6.0, moments) * FracMom(-1.0/6.0, moments) +
		  FracMom( 7.0/2.0, moments) * FracMom( 1.0/6.0, moments)
           ) + 4.0 * (
	          FracMom(19.0/6.0, moments) * FracMom( 1.0/2.0, moments) +
	    2.0 * FracMom(17.0/6.0, moments) * FracMom( 5.0/6.0, moments) +
		  FracMom( 5.0/2.0, moments) * FracMom( 7.0/6.0, moments)
	   ) + 6.0 * (
	          FracMom(13.0/6.0, moments) * FracMom( 3.0/2.0, moments) +
	    2.0 * FracMom(11.0/6.0, moments) * FracMom(11.0/6.0, moments) +
		  FracMom( 3.0/2.0, moments) * FracMom(13.0/6.0, moments)
	   ) + 4.0 * (
	          FracMom( 7.0/6.0, moments) * FracMom( 5.0/2.0, moments) +
	    2.0 * FracMom( 5.0/6.0, moments) * FracMom(17.0/6.0, moments) +
		  FracMom( 1.0/2.0, moments) * FracMom(19.0/6.0, moments)
	   ) + (
	          FracMom( 1.0/6.0, moments) * FracMom( 7.0/2.0, moments) +
	    2.0 * FracMom(-1.0/6.0, moments) * FracMom(23.0/6.0, moments) +
		  FracMom(-1.0/2.0, moments) * FracMom(25.0/6.0, moments)
	   );

  if (index > 3)
    phi6 = (
                  FracMom(31.0/6.0, moments) * FracMom(-1.0/2.0, moments) +
	    2.0 * FracMom(29.0/6.0, moments) * FracMom(-1.0/6.0, moments) +
		  FracMom( 9.0/2.0, moments) * FracMom( 1.0/6.0, moments)
           ) + 5.0 * (
	          FracMom(25.0/6.0, moments) * FracMom( 1.0/2.0, moments) +
	    2.0 * FracMom(23.0/6.0, moments) * FracMom( 5.0/6.0, moments) +
		  FracMom( 7.0/2.0, moments) * FracMom( 7.0/6.0, moments)
	   ) + 10.0 * (
	          FracMom(19.0/6.0, moments) * FracMom( 3.0/2.0, moments) +
	    2.0 * FracMom(17.0/6.0, moments) * FracMom(11.0/6.0, moments) +
		  FracMom( 5.0/2.0, moments) * FracMom(13.0/6.0, moments)
	   ) + 10.0 * (
	          FracMom(13.0/6.0, moments) * FracMom( 5.0/2.0, moments) +
	    2.0 * FracMom(11.0/6.0, moments) * FracMom(17.0/6.0, moments) +
		  FracMom( 3.0/2.0, moments) * FracMom(19.0/6.0, moments)
	   ) + 5.0 * (
	          FracMom( 7.0/6.0, moments) * FracMom( 7.0/2.0, moments) +
	    2.0 * FracMom( 5.0/6.0, moments) * FracMom(23.0/6.0, moments) +
		  FracMom( 1.0/2.0, moments) * FracMom(25.0/6.0, moments)
	   ) + (
		  FracMom( 1.0/6.0, moments) * FracMom( 9.0/2.0, moments) +
	    2.0 * FracMom(-1.0/6.0, moments) * FracMom(29.0/6.0, moments) +
		  FracMom(-1.0/2.0, moments) * FracMom(31.0/6.0, moments)
	   );
  
  if (index > 4)
    phi7 = (
                  FracMom(37.0/6.0, moments) * FracMom(-1.0/2.0, moments) +
            2.0 * FracMom(35.0/6.0, moments) * FracMom(-1.0/6.0, moments) +
		  FracMom(11.0/2.0, moments) * FracMom( 1.0/6.0, moments)
           ) + 6.0 * (
	          FracMom(31.0/6.0, moments) * FracMom( 1.0/2.0, moments) +
	    2.0 * FracMom(29.0/6.0, moments) * FracMom( 5.0/6.0, moments) +
		  FracMom( 9.0/2.0, moments) * FracMom( 7.0/6.0, moments)
	   ) + 15.0 * (
	          FracMom(25.0/6.0, moments) * FracMom( 3.0/2.0, moments) +
	    2.0 * FracMom(23.0/6.0, moments) * FracMom(11.0/6.0, moments) +
		  FracMom( 7.0/2.0, moments) * FracMom(13.0/6.0, moments)
	   ) + 20.0 * (
	          FracMom(19.0/6.0, moments) * FracMom( 5.0/2.0, moments) +
	    2.0 * FracMom(17.0/6.0, moments) * FracMom(17.0/6.0, moments) +
		  FracMom( 5.0/2.0, moments) * FracMom(19.0/6.0, moments)
	   ) + 15.0 * (
	          FracMom(13.0/6.0, moments) * FracMom( 7.0/2.0, moments) +
	    2.0 * FracMom(11.0/6.0, moments) * FracMom(23.0/6.0, moments) +
		  FracMom( 3.0/2.0, moments) * FracMom(25.0/6.0, moments)
	   ) + 6.0 * (
	          FracMom( 7.0/6.0, moments) * FracMom( 9.0/2.0, moments) +
	    2.0 * FracMom( 5.0/6.0, moments) * FracMom(29.0/6.0, moments) +
		  FracMom( 1.0/2.0, moments) * FracMom(31.0/6.0, moments)
	   ) + (
	          FracMom( 1.0/6.0, moments) * FracMom(11.0/2.0, moments) +
	    2.0 * FracMom(-1.0/6.0, moments) * FracMom(35.0/6.0, moments) +
		  FracMom(-1.0/2.0, moments) * FracMom(37.0/6.0, moments)
	   );

  switch (index)
  {
    case 0: return C * pow(phi1, ( 1.0/2.0)-x) * pow(phi2, x+( 1.0/2.0));
    case 1: return C * pow(phi2, ( 3.0/2.0)-x) * pow(phi3, x-( 1.0/2.0));
    case 2: return C * pow(phi3, ( 5.0/2.0)-x) * pow(phi4, x-( 3.0/2.0));
    case 3: return C * pow(phi4, ( 7.0/2.0)-x) * pow(phi5, x-( 5.0/2.0));
    case 4: return C * pow(phi5, ( 9.0/2.0)-x) * pow(phi6, x-( 7.0/2.0));
    case 5: return C * pow(phi6, (11.0/2.0)-x) * pow(phi7, x-( 9.0/2.0));
  }

  return 0.0;
}

double TSoot::GetPhi(double x, double y, double * moments, double temp)
{
	double C = GetC(temp);
	double phi1, phi2, phi3, phi4, phi5, phi6, phi7, phi8, phi9, phi10;

	double Df1 = 1.8;
	double Df2 = 1.8;

	double av1 = 1.0 - (2.0 / Df1);
	double as1 = (3.0 / Df1) - 1.0;

	double av2 = 1.0 - (2.0 / Df2);
	double as2 = (3.0 / Df2) - 1.0;

	//Phi1  = Phi(-1.0/2.0, 0)
	//Phi2  = Phi( 1.0/2.0, 0)
	//Phi3  = Phi(-1.0/2.0, 1)
	//Phi4  = Phi( 3.0/2.0, 0)
	//Phi5  = Phi( 1.0/2.0, 1)
	//Phi6  = Phi(-1.0/2.0, 2)
	//Phi7  = Phi( 5.0/2.0, 0)
	//Phi8  = Phi( 3.0/2.0, 1)
	//Phi9  = Phi( 1.0/2.0, 2)
	//Phi10 = Phi(-1.0/2.0, 3)

	int index;

	if ((x+y) < 0.5)
		index = 0;
	else if ((x+y) < 1.5)
	{
		if (x > 0.5)
			index = 1;
		else if (y < 1.0)
			index = 2;
		else
			index = 3;
	}
	else if ((x+y) <= 2.5)
	{
		if (x >= 1.5)
			index = 4;
		else if (y < 1.0)
			index = 5;
		else if (x >= 0.5)
			index = 6;
		else if (y < 2.0)
			index = 7;
		else
			index = 8;
	}
/*	else
	{
		if (y < 1.0)
			index = 9;
		if (x < 0.5)
			index = 10;
		else
			index = 11;
	}*/
	else
		cerr << "Coagulation is not yet implemented for more than six moments!" << endl;

	if (index == 0)
	phi1 =	(
				FracMom(2.0*av1-1.0/2.0, 2.0*as1, moments) * FracMom(-1.0/2.0, 0.0, moments) +
			2.0 *	FracMom(av1-1.0/2.0, as1, moments) * FracMom(av2-1.0/2.0, as2, moments) +
				FracMom(-1.0/2.0, 0.0, moments) * FracMom(2.0*av2-1.0/2.0, 2.0*as2, moments)
		);

	if (index == 0 || index == 1 || index == 2)
	phi2 =	(
				FracMom(2.0*av1+1.0/2.0, 2.0*as1, moments) * FracMom(-1.0/2.0, 0.0, moments) +
			2.0 *	FracMom(av1+1.0/2.0, as1, moments) * FracMom(av2-1.0/2.0, as2, moments) +
				FracMom(1.0/2.0, 0.0, moments) * FracMom(2.0*av2-1.0/2.0, 2.0*as2, moments)
		) + (
				FracMom(2.0*av1-1.0/2.0, 2.0*as1, moments) * FracMom(1.0/2.0, 0.0, moments) +
			2.0 *	FracMom(av1-1.0/2.0, as1, moments) * FracMom(av2+1.0/2.0, as2, moments) +
				FracMom(-1.0/2.0, 0.0, moments) * FracMom(2.0*av2+1.0/2.0, 2.0*as2, moments)
		);

	if (index == 0 || index == 2 || index == 3)
	phi3 =	(
				FracMom(2.0*av1-1.0/2.0, 2.0*as1+1.0, moments) * FracMom(-1.0/2.0, 0.0, moments) +
			2.0 *	FracMom(av1-1.0/2.0, as1+1.0, moments) * FracMom(av2-1.0/2.0, as2, moments) +
				FracMom(-1.0/2.0, 1.0, moments) * FracMom(2.0*av2-1.0/2.0, 2.0*as2, moments)
		) + (
				FracMom(2.0*av1-1.0/2.0, 2.0*as1, moments) * FracMom(-1.0/2.0, 1.0, moments) +
			2.0 *	FracMom(av1-1.0/2.0, as1, moments) * FracMom(av2-1.0/2.0, as2+1.0, moments) +
				FracMom(-1.0/2.0, 0.0, moments) * FracMom(2.0*av2-1.0/2.0, 2.0*as2+1.0, moments)
		);

	if (index == 1 || index == 4 || index == 5)
	phi4 = (
				FracMom(2.0*av1+3.0/2.0, 2.0*as1, moments) * FracMom(-1.0/2.0, 0.0, moments) +
			2.0 *	FracMom(av1+3.0/2.0, as1, moments) * FracMom(av2-1.0/2.0, as2, moments) +
				FracMom(3.0/2.0, 0.0, moments) * FracMom(2.0*av2-1.0/2.0, 2.0*as2, moments)
		) + 2.0 * (
				FracMom(2.0*av1+1.0/2.0, 2.0*as1, moments) * FracMom(1.0/2.0, 0.0, moments) +
			2.0 *	FracMom(av1+1.0/2.0, as1, moments) * FracMom(av2+1.0/2.0, as2, moments) +
				FracMom(1.0/2.0, 0.0, moments) * FracMom(2.0*av2+1.0/2.0, 2.0*as2, moments)
		) + (
				FracMom(2.0*av1-1.0/2.0, 2.0*as1, moments) * FracMom(3.0/2.0, 0.0, moments) +
			2.0 *	FracMom(av1-1.0/2.0, as1, moments) * FracMom(av2+3.0/2.0, as2, moments) +
				FracMom(-1.0/2.0, 0.0, moments) * FracMom(2.0*av2+3.0/2.0, 2.0*as2, moments)
		);

	if (index == 1 || index == 2 || index == 3 || index == 5 || index == 6 || index == 7)
	phi5 =	(
				FracMom(2.0*av1+1.0/2.0, 2.0*as1+1.0, moments) * FracMom(-1.0/2.0, 0.0, moments) +
			2.0 *	FracMom(av1+1.0/2.0, as1+1.0, moments) * FracMom(av2-1.0/2.0, as2, moments) +
				FracMom(1.0/2.0, 1.0, moments) * FracMom(2.0*av2-1.0/2.0, 2.0*as2, moments)
		) + (
				FracMom(2.0*av1+1.0/2.0, 2.0*as1, moments) * FracMom(-1.0/2.0, 1.0, moments) +
			2.0 *	FracMom(av1+1.0/2.0, as1, moments) * FracMom(av2-1.0/2.0, as2+1.0, moments) +
				FracMom(1.0/2.0, 0.0, moments) * FracMom(2.0*av2-1.0/2.0, 2.0*as2+1.0, moments)
		) + (
				FracMom(2.0*av1-1.0/2.0, 2.0*as1+1.0, moments) * FracMom(1.0/2.0, 0.0, moments) +
			2.0 *	FracMom(av1-1.0/2.0, as1+1.0, moments) * FracMom(av2+1.0/2.0, as2, moments) +
				FracMom(-1.0/2.0, 1.0, moments) * FracMom(2.0*av2+1.0/2.0, 2.0*as2, moments)
		) + (
				FracMom(2.0*av1-1.0/2.0, 2.0*as1, moments) * FracMom(1.0/2.0, 1.0, moments) +
			2.0 *	FracMom(av1-1.0/2.0, as1, moments) * FracMom(av2+1.0/2.0, as2+1.0, moments) +
				FracMom(-1.0/2.0, 0.0, moments) * FracMom(2.0*av2+1.0/2.0, 2.0*as2+1.0, moments)
		);

	if (index == 3 || index == 7 || index == 8)
	phi6 =	(
				FracMom(2.0*av1-1.0/2.0, 2.0*as1+2.0, moments) * FracMom(-1.0/2.0, 0.0, moments) +
			2.0 *	FracMom(av1-1.0/2.0, as1+2.0, moments) * FracMom(av2-1.0/2.0, as2, moments) +
				FracMom(-1.0/2.0, 2.0, moments) * FracMom(2.0*av2-1.0/2.0, 2.0*as2, moments)
		) + 2.0 * (
				FracMom(2.0*av1-1.0/2.0, 2.0*as1+1.0, moments) * FracMom(-1.0/2.0, 1.0, moments) +
			2.0 *	FracMom(av1-1.0/2.0, as1+1.0, moments) * FracMom(av2-1.0/2.0, as2+1.0, moments) +
				FracMom(-1.0/2.0, 1.0, moments) * FracMom(2.0*av2-1.0/2.0, 2.0*as2+1.0, moments)
		) + (
				FracMom(2.0*av1-1.0/2.0, 2.0*as1, moments) * FracMom(-1.0/2.0, 2.0, moments) +
			2.0 *	FracMom(av1-1.0/2.0, as1, moments) * FracMom(av2-1.0/2.0, as2+2.0, moments) +
				FracMom(-1.0/2.0, 0.0, moments) * FracMom(2.0*av2-1.0/2.0, 2.0*as2+2.0, moments)
		);

	if (index == 4)
	phi7 =	(
				FracMom(2.0*av1+5.0/2.0, 2.0*as1, moments) * FracMom(-1.0/2.0, 0.0, moments) +
			2.0 *	FracMom(av1+5.0/2.0, as1, moments) * FracMom(av2-1.0/2.0, as2, moments) +
				FracMom(5.0/2.0, 0.0, moments) * FracMom(2.0*av2-1.0/2.0, 2.0*as2, moments)
		) + 3.0 * (
				FracMom(2.0*av1+3.0/2.0, 2.0*as1, moments) * FracMom(1.0/2.0, 0.0, moments) +
			2.0 *	FracMom(av1+3.0/2.0, as1, moments) * FracMom(av2+1.0/2.0, as2, moments) +
				FracMom(3.0/2.0, 0.0, moments) * FracMom(2.0*av2+1.0/2.0, 2.0*as2, moments)
		) + 3.0 * (
				FracMom(2.0*av1+1.0/2.0, 2.0*as1, moments) * FracMom(3.0/2.0, 0.0, moments) +
			2.0 *	FracMom(av1+1.0/2.0, as1, moments) * FracMom(av2+3.0/2.0, as2, moments) +
				FracMom(1.0/2.0, 0.0, moments) * FracMom(2.0*av2+3.0/2.0, 2.0*as2, moments)
		) + (
				FracMom(2.0*av1-1.0/2.0, 2.0*as1, moments) * FracMom(5.0/2.0, 0.0, moments) +
			2.0 *	FracMom(av1-1.0/2.0, as1, moments) * FracMom(av2+5.0/2.0, as2, moments) +
				FracMom(-1.0/2.0, 0.0, moments) * FracMom(2.0*av2+5.0/2.0, 2.0*as2, moments)
		);

	if (index == 4 || index == 5 || index == 6)
	phi8 =	(
				FracMom(2.0*av1+3.0/2.0, 2.0*as1+1.0, moments) * FracMom(-1.0/2.0, 0.0, moments) +
			2.0 *	FracMom(av1+3.0/2.0, as1+1.0, moments) * FracMom(av2-1.0/2.0, as2, moments) +
				FracMom(3.0/2.0, 1.0, moments) * FracMom(2.0*av2-1.0/2.0, 2.0*as2, moments)
		) + 2.0 * (
				FracMom(2.0*av1+1.0/2.0, 2.0*as1+1.0, moments) * FracMom(1.0/2.0, 0.0, moments) +
			2.0 *	FracMom(av1+1.0/2.0, as1+1.0, moments) * FracMom(av2+1.0/2.0, as2, moments) +
				FracMom(1.0/2.0, 1.0, moments) * FracMom(2.0*av2+1.0/2.0, 2.0*as2, moments)
		) + (
				FracMom(2.0*av1-1.0/2.0, 2.0*as1+1.0, moments) * FracMom(3.0/2.0, 0.0, moments) +
			2.0 *	FracMom(av1-1.0/2.0, as1+1.0, moments) * FracMom(av2+3.0/2.0, as2, moments) +
				FracMom(-1.0/2.0, 1.0, moments) * FracMom(2.0*av2+3.0/2.0, 2.0*as2, moments)
		) + (
				FracMom(2.0*av1+3.0/2.0, 2.0*as1, moments) * FracMom(-1.0/2.0, 1.0, moments) +
			2.0 *	FracMom(av1+3.0/2.0, as1, moments) * FracMom(av2-1.0/2.0, as2+1.0, moments) +
				FracMom(3.0/2.0, 0.0, moments) * FracMom(2.0*av2-1.0/2.0, 2.0*as2+1.0, moments)
		) + 2.0 * (
				FracMom(2.0*av1+1.0/2.0, 2.0*as1, moments) * FracMom(1.0/2.0, 1.0, moments) +
			2.0 *	FracMom(av1+1.0/2.0, as1, moments) * FracMom(av2+1.0/2.0, as2+1.0, moments) +
				FracMom(1.0/2.0, 0.0, moments) * FracMom(2.0*av2+1.0/2.0, 2.0*as2+1.0, moments)
		) + (
				FracMom(2.0*av1-1.0/2.0, 2.0*as1, moments) * FracMom(3.0/2.0, 1.0, moments) +
			2.0 *	FracMom(av1-1.0/2.0, as1, moments) * FracMom(av2+3.0/2.0, as2+1.0, moments) +
				FracMom(-1.0/2.0, 0.0, moments) * FracMom(2.0*av2+3.0/2.0, 2.0*as2+1.0, moments)
		);

	if (index == 6 || index == 7 || index == 8)
	phi9 =	(
				FracMom(2.0*av1+1.0/2.0, 2.0*as1+2.0, moments) * FracMom(-1.0/2.0, 0.0, moments) +
			2.0 *	FracMom(av1+1.0/2.0, as1+2.0, moments) * FracMom(av2-1.0/2.0, as2, moments) +
				FracMom(1.0/2.0, 2.0, moments) * FracMom(2.0*av2-1.0/2.0, 2.0*as2, moments)
		) + (
				FracMom(2.0*av1-1.0/2.0, 2.0*as1+2.0, moments) * FracMom(1.0/2.0, 0.0, moments) +
			2.0 *	FracMom(av1-1.0/2.0, as1+2.0, moments) * FracMom(av2+1.0/2.0, as2, moments) +
				FracMom(-1.0/2.0, 2.0, moments) * FracMom(2.0*av2+1.0/2.0, 2.0*as2, moments)
		) + 2.0 * (
				FracMom(2.0*av1+1.0/2.0, 2.0*as1+1.0, moments) * FracMom(-1.0/2.0, 1.0, moments) +
			2.0 *	FracMom(av1+1.0/2.0, as1+1.0, moments) * FracMom(av2-1.0/2.0, as2+1.0, moments) +
				FracMom(1.0/2.0, 1.0, moments) * FracMom(2.0*av2-1.0/2.0, 2.0*as2+1.0, moments)
		) + 2.0 * (
				FracMom(2.0*av1-1.0/2.0, 2.0*as1+1.0, moments) * FracMom(1.0/2.0, 1.0, moments) +
			2.0 *	FracMom(av1-1.0/2.0, as1+1.0, moments) * FracMom(av2+1.0/2.0, as2+1.0, moments) +
				FracMom(-1.0/2.0, 1.0, moments) * FracMom(2.0*av2+1.0/2.0, 2.0*as2+1.0, moments)
		) + (
				FracMom(2.0*av1+1.0/2.0, 2.0*as1, moments) * FracMom(-1.0/2.0, 2.0, moments) +
			2.0 *	FracMom(av1+1.0/2.0, as1, moments) * FracMom(av2-1.0/2.0, as2+2.0, moments) +
				FracMom(1.0/2.0, 0.0, moments) * FracMom(2.0*av2-1.0/2.0, 2.0*as2+2.0, moments)
		) + (
				FracMom(2.0*av1-1.0/2.0, 2.0*as1, moments) * FracMom(1.0/2.0, 2.0, moments) +
			2.0 *	FracMom(av1-1.0/2.0, as1, moments) * FracMom(av2+1.0/2.0, as2+2.0, moments) +
				FracMom(-1.0/2.0, 0.0, moments) * FracMom(2.0*av2+1.0/2.0, 2.0*as2+2.0, moments)
		);

	if (index == 8)
	phi10 =	(
				FracMom(2.0*av1-1.0/2.0, 2.0*as1+3.0, moments) * FracMom(-1.0/2.0, 0.0, moments) +
			2.0 *	FracMom(av1-1.0/2.0, as1+3.0, moments) * FracMom(av2-1.0/2.0, as2, moments) +
				FracMom(-1.0/2.0, 3.0, moments) * FracMom(2.0*av2-1.0/2.0, 2.0*as2, moments)
		) + 3.0 * (
				FracMom(2.0*av1-1.0/2.0, 2.0*as1+2.0, moments) * FracMom(-1.0/2.0, 1.0, moments) +
			2.0 *	FracMom(av1-1.0/2.0, as1+2.0, moments) * FracMom(av2-1.0/2.0, as2+1.0, moments) +
				FracMom(-1.0/2.0, 2.0, moments) * FracMom(2.0*av2-1.0/2.0, 2.0*as2+1.0, moments)
		) + 3.0 * (
				FracMom(2.0*av1-1.0/2.0, 2.0*as1+1.0, moments) * FracMom(-1.0/2.0, 2.0, moments) +
			2.0 *	FracMom(av1-1.0/2.0, as1+1.0, moments) * FracMom(av2-1.0/2.0, as2+2.0, moments) +
				FracMom(-1.0/2.0, 1.0, moments) * FracMom(2.0*av2-1.0/2.0, 2.0*as2+2.0, moments)
		) + (
				FracMom(2.0*av1-1.0/2.0, 2.0*as1, moments) * FracMom(-1.0/2.0, 3.0, moments) +
			2.0 *	FracMom(av1-1.0/2.0, as1, moments) * FracMom(av2-1.0/2.0, as2+3.0, moments) +
				FracMom(-1.0/2.0, 0.0, moments) * FracMom(2.0*av2-1.0/2.0, 2.0*as2+3.0, moments)
		);

	switch (index)
	{
		case 0:	return C * pow(phi1, (1.0/2.0)-x-y) * pow(phi2, x+(1.0/2.0)) * pow(phi3 , y);
		case 1: return C * pow(phi2, (3.0/2.0)-x-y) * pow(phi4, x-(1.0/2.0)) * pow(phi5 , y);
		case 2: return C * pow(phi5, x+y-(1.0/2.0)) * pow(phi3, (1.0/2.0)-x) * pow(phi2 , 1.0-y);
		case 3: return C * pow(phi3, (3.0/2.0)-x-y) * pow(phi5, x+(1.0/2.0)) * pow(phi6 , y-1.0);
		case 4: return C * pow(phi4, (5.0/2.0)-x-y) * pow(phi7, x-(3.0/2.0)) * pow(phi8 , y);
		case 5: return C * pow(phi8, x+y-(3.0/2.0)) * pow(phi5, (3.0/2.0)-x) * pow(phi4 , 1.0-y);
		case 6: return C * pow(phi5, (5.0/2.0)-x-y) * pow(phi8, x-(1.0/2.0)) * pow(phi9 , y-1.0);
		case 7: return C * pow(phi9, x+y-(3.0/2.0)) * pow(phi6, (1.0/2.0)-x) * pow(phi5 , 2.0-y);
		case 8: return C * pow(phi6, (5.0/2.0)-x-y) * pow(phi9, x+(1.0/2.0)) * pow(phi10, y-2.0);
	}

	return 0.0;
}

double TSoot::GetPhi(double x, double y, double z, double * moments, double temp)
{
	Double C = GetC(temp);
	Double phi1, phi2, phi3, phi4;

	Double Df1 = 3.0;
	Double Df2 = 3.0;

	Double av1 = 1.0 - (2.0 / Df1);
	Double as1 = (3.0 / Df1) - 1.0;

	Double av2 = 1.0 - (2.0 / Df2);
	Double as2 = (3.0 / Df2) - 1.0;

	//Phi1 = Phi(-1.0/2.0, 0, 0) & Phi2 = Phi(1.0/2.0, 0, 0) & Phi3 = Phi(-1.0/2.0, 1, 0) & Phi4 = Phi(-1.0/2.0, 0, 1)

	phi1 =	(
				FracMom(2.0*av1-1.0/2.0, 2.0*as1, 0.0, moments) * FracMom(-1.0/2.0, 0.0, 0.0, moments) +
			2.0 *	FracMom(av1-1.0/2.0, as1, 0.0, moments) * FracMom(av2-1.0/2.0, as2, 0.0, moments) +
				FracMom(-1.0/2.0, 0.0, 0.0, moments) * FracMom(2.0*av2-1.0/2.0, 2.0*as2, 0.0, moments)
		);

	phi2 =	(		FracMom(2.0*av1+1.0/2.0, 2.0*as1, 0.0, moments) * FracMom(-1.0/2.0, 0.0, 0.0, moments) +
			2.0 *	FracMom(av1+1.0/2.0, as1, 0.0, moments) * FracMom(av2-1.0/2.0, as2, 0.0, moments) +
				FracMom(1.0/2.0, 0.0, 0.0, moments) * FracMom(2.0*av2-1.0/2.0, 2.0*as2, 0.0, moments)
		) + (
				FracMom(2.0*av1-1.0/2.0, 2.0*as1, 0.0, moments) * FracMom(1.0/2.0, 0.0, 0.0, moments) +
			2.0 *	FracMom(av1-1.0/2.0, as1, 0.0, moments) * FracMom(av2+1.0/2.0, as2, 0.0, moments) +
				FracMom(-1.0/2.0, 0.0, 0.0, moments) * FracMom(2.0*av2+1.0/2.0, 2.0*as2, 0.0, moments)
		);

	phi3 =	(
				FracMom(2.0*av1-1.0/2.0, 2.0*as1+1.0, 0.0, moments) * FracMom(-1.0/2.0, 0.0,0.0, moments) +
			2.0 *	FracMom(av1-1.0/2.0, as1+1.0, 0.0, moments) * FracMom(av2-1.0/2.0, as2, 0.0, moments) +
				FracMom(-1.0/2.0, 1.0, 0.0, moments) * FracMom(2.0*av2-1.0/2.0, 2.0*as2, 0.0, moments)
		) + (
				FracMom(2.0*av1-1.0/2.0, 2.0*as1, 0.0, moments) * FracMom(-1.0/2.0, 1.0, 0.0, moments) +
			2.0 *	FracMom(av1-1.0/2.0, as1, 0.0, moments) * FracMom(av2-1.0/2.0, as2+1.0, 0.0, moments) +
				FracMom(-1.0/2.0, 0.0, 0.0, moments) * FracMom(2.0*av2-1.0/2.0, 2.0*as2+1.0, 0.0, moments)
		);

	phi4 =	(
				FracMom(2.0*av1-1.0/2.0, 2.0*as1, 1.0, moments) * FracMom(-1.0/2.0, 0.0, 0.0, moments) +
			2.0 *	FracMom(av1-1.0/2.0, as1, 1.0, moments) * FracMom(av2-1.0/2.0, as2, 0.0, moments) +
				FracMom(-1.0/2.0, 0.0, 1.0, moments) * FracMom(2.0*av2-1.0/2.0, 2.0*as2, 0.0, moments)
		) + (
				FracMom(2.0*av1-1.0/2.0, 2.0*as1, 0.0, moments) * FracMom(-1.0/2.0, 0.0, 1.0, moments) +
			2.0 *	FracMom(av1-1.0/2.0, as1, 0.0, moments) * FracMom(av2-1.0/2.0, as2, 1.0, moments) +
				FracMom(-1.0/2.0, 0.0, 0.0, moments) * FracMom(2.0*av2-1.0/2.0, 2.0*as2, 1.0, moments)
		);

	return C * pow(phi1, (1.0/2.0)-x-y-z) * pow(phi2, (1.0/2.0)+x) * pow(phi3, y) * pow(phi4, z);
}

double TSoot::GetPsi(double u, double x, double * moments, double temp)
{
  double C = GetC(temp);
  double psi1, psi2, psi3, psi4, psi5, psi6;

  psi1 = (
	        FracMom( 1.0/6.0+u, moments) * FracMom(-3.0/6.0+x, moments) +
	  2.0 * FracMom(-1.0/6.0+u, moments) * FracMom(-1.0/6.0+x, moments) +
	        FracMom(-3.0/6.0+u, moments) * FracMom( 1.0/6.0+x, moments)
	 );

  psi2 = (
	        FracMom( 7.0/6.0+u, moments) * FracMom(-3.0/6.0+x, moments) +
	  2.0 * FracMom( 5.0/6.0+u, moments) * FracMom(-1.0/6.0+x, moments) +
	        FracMom( 3.0/6.0+u, moments) * FracMom( 1.0/6.0+x, moments)
	 ) + (
	        FracMom( 1.0/6.0+u, moments) * FracMom( 3.0/6.0+x, moments) +
	  2.0 * FracMom(-1.0/6.0+u, moments) * FracMom( 5.0/6.0+x, moments) +
	        FracMom(-3.0/6.0+u, moments) * FracMom( 7.0/6.0+x, moments)
	 );

  psi3 = (
	        FracMom(13.0/6.0+u, moments) * FracMom(-3.0/6.0+x, moments) +
	  2.0 * FracMom(11.0/6.0+u, moments) * FracMom(-1.0/6.0+x, moments) +
	        FracMom( 9.0/6.0+u, moments) * FracMom( 1.0/6.0+x, moments)
	 ) + 2.0 * (
	        FracMom( 7.0/6.0+u, moments) * FracMom( 3.0/6.0+x, moments) +
	  2.0 * FracMom( 5.0/6.0+u, moments) * FracMom( 5.0/6.0+x, moments) +
	        FracMom( 3.0/6.0+u, moments) * FracMom( 7.0/6.0+x, moments)
	 ) + (
	        FracMom( 1.0/6.0+u, moments) * FracMom( 9.0/6.0+x, moments) +
	  2.0 * FracMom(-1.0/6.0+u, moments) * FracMom(11.0/6.0+x, moments) +
	        FracMom(-3.0/6.0+u, moments) * FracMom(13.0/6.0+x, moments)
	 );

  if (fNSootMoments > 2)
  psi4 = (
	        FracMom(19.0/6.0+u, moments) * FracMom(-3.0/6.0+x, moments) +
	  2.0 * FracMom(17.0/6.0+u, moments) * FracMom(-1.0/6.0+x, moments) +
	        FracMom(15.0/6.0+u, moments) * FracMom( 1.0/6.0+x, moments)
	 ) + 3.0 * (
	        FracMom(13.0/6.0+u, moments) * FracMom( 3.0/6.0+x, moments) +
	  2.0 * FracMom(11.0/6.0+u, moments) * FracMom( 5.0/6.0+x, moments) +
	        FracMom( 9.0/6.0+u, moments) * FracMom( 7.0/6.0+x, moments)
	 ) + 3.0 * (
	        FracMom( 7.0/6.0+u, moments) * FracMom( 9.0/6.0+x, moments) +
	  2.0 * FracMom( 5.0/6.0+u, moments) * FracMom(11.0/6.0+x, moments) +
	        FracMom( 3.0/6.0+u, moments) * FracMom(13.0/6.0+x, moments)
	 ) + (
	        FracMom( 1.0/6.0+u, moments) * FracMom(15.0/6.0+x, moments) +
	  2.0 * FracMom(-1.0/6.0+u, moments) * FracMom(17.0/6.0+x, moments) +
	        FracMom(-3.0/6.0+u, moments) * FracMom(19.0/6.0+x, moments)
	 );

//   if (fNSootMoments > 3)
//   psi5 = (
// 	        FracMom(25.0/6.0+u, moments) * FracMom(-3.0/6.0+x, moments) +
// 	  2.0 * FracMom(23.0/6.0+u, moments) * FracMom(-1.0/6.0+x, moments) +
// 	        FracMom(21.0/6.0+u, moments) * FracMom( 1.0/6.0+x, moments)
// 	 ) + 4.0 * (
// 	        FracMom(19.0/6.0+u, moments) * FracMom( 3.0/6.0+x, moments) +
// 	  2.0 * FracMom(17.0/6.0+u, moments) * FracMom( 5.0/6.0+x, moments) +
// 	        FracMom(15.0/6.0+u, moments) * FracMom( 7.0/6.0+x, moments)
// 	 ) + 6.0 * (
// 	        FracMom(13.0/6.0+u, moments) * FracMom( 9.0/6.0+x, moments) +
// 	  2.0 * FracMom(11.0/6.0+u, moments) * FracMom(11.0/6.0+x, moments) +
// 	        FracMom( 9.0/6.0+u, moments) * FracMom(13.0/6.0+x, moments)
// 	 ) + 4.0 * (
// 	        FracMom( 7.0/6.0+u, moments) * FracMom(15.0/6.0+x, moments) +
// 	  2.0 * FracMom( 5.0/6.0+u, moments) * FracMom(17.0/6.0+x, moments) +
// 	        FracMom( 3.0/6.0+u, moments) * FracMom(19.0/6.0+x, moments)
// 	 ) + (
// 	        FracMom( 1.0/6.0+u, moments) * FracMom(21.0/6.0+x, moments) +
// 	  2.0 * FracMom(-1.0/6.0+u, moments) * FracMom(23.0/6.0+x, moments) +
// 	        FracMom(-3.0/6.0+u, moments) * FracMom(25.0/6.0+x, moments)
// 	 );

//   if (fNSootMoments > 4)
//   psi6 = (
// 	        FracMom(31.0/6.0+u, moments) * FracMom(-3.0/6.0+x, moments) +
// 	  2.0 * FracMom(29.0/6.0+u, moments) * FracMom(-1.0/6.0+x, moments) +
// 	        FracMom(27.0/6.0+u, moments) * FracMom( 1.0/6.0+x, moments)
// 	 ) + 5.0 * (
// 	        FracMom(25.0/6.0+u, moments) * FracMom( 3.0/6.0+x, moments) +
// 	  2.0 * FracMom(23.0/6.0+u, moments) * FracMom( 5.0/6.0+x, moments) +
// 	        FracMom(21.0/6.0+u, moments) * FracMom( 7.0/6.0+x, moments)
// 	 ) + 10.0 * (
// 	        FracMom(19.0/6.0+u, moments) * FracMom( 9.0/6.0+x, moments) +
// 	  2.0 * FracMom(17.0/6.0+u, moments) * FracMom(11.0/6.0+x, moments) +
// 	        FracMom(15.0/6.0+u, moments) * FracMom(13.0/6.0+x, moments)
// 	 ) + 10.0 * (
// 	        FracMom(13.0/6.0+u, moments) * FracMom(15.0/6.0+x, moments) +
// 	  2.0 * FracMom(11.0/6.0+u, moments) * FracMom(17.0/6.0+x, moments) +
// 	        FracMom( 9.0/6.0+u, moments) * FracMom(19.0/6.0+x, moments)
// 	 ) + 5.0 * (
// 	        FracMom( 7.0/6.0+u, moments) * FracMom(21.0/6.0+x, moments) +
// 	  2.0 * FracMom( 5.0/6.0+u, moments) * FracMom(23.0/6.0+x, moments) +
// 	        FracMom( 3.0/6.0+u, moments) * FracMom(25.0/6.0+x, moments)
// 	 ) + (
// 	        FracMom( 1.0/6.0+u, moments) * FracMom(27.0/6.0+x, moments) +
// 	  2.0 * FracMom(-1.0/6.0+u, moments) * FracMom(29.0/6.0+x, moments) +
// 	        FracMom(-3.0/6.0+u, moments) * FracMom(31.0/6.0+x, moments)
// 	 );

  //if (fNSootMoments > 1)
  //  return C * pow(psi1,3.0/8.0)*pow(psi2,3.0/4.0)*pow(psi3,-1.0/8.0);
  //else if (fNSootMoments == 3)
  //  return C * pow(psi1,5.0/16.0)*pow(psi2,15.0/16.0)*pow(psi3,-5.0/16.0)*pow(psi4,1.0/16.0);
  //else if (fNSootMoments == 4)
  //  return C * pow(psi1,35.0/128.0)*pow(psi2,35.0/32.0)*pow(psi3,-35.0/64.0)*pow(psi4,7.0/32.0)*pow(psi5,-5.0/128.0);
  //else if (fNSootMoments == 5)
  //  return C * pow(psi1,63.0/256.0)*pow(psi2,315.0/256.0)*pow(psi3,-105.0/128.0)*pow(psi4,63.0/128.0)*pow(psi5,-45.0/256.0)*pow(psi6,7.0/256.0);

  if (u < 2.99 && x < 2.99)
    return C * pow(psi1,5.0/16.0)*pow(psi2,15.0/16.0)*pow(psi3,-5.0/16.0)*pow(psi4,1.0/16.0);
  else if (u < 3.99 && x < 3.99)
    return C * pow(psi1,3.0/8.0)*pow(psi2,3.0/4.0)*pow(psi3,-1.0/8.0);
  else
    return C * pow(psi1,1.0/2.0)*pow(psi2,1.0/2.0);
}

double TSoot::GetPsiSMALL(double * moments, double temp)
{
  double C = GetC(temp);

  double psi1, psi2, psi3, psi4;

  double V0 = 2.0 * dimer_nbrC2;


  psi1 = (
	        pow(V0, 1.0/6.0) * FracMom(-3.0/6.0, moments) +
	  2.0 * pow(V0,-1.0/6.0) * FracMom(-1.0/6.0, moments) +
	        pow(V0,-3.0/6.0) * FracMom( 1.0/6.0, moments)
	 );

  psi2 = (
	        pow(V0, 7.0/6.0) * FracMom(-3.0/6.0, moments) +
	  2.0 * pow(V0, 5.0/6.0) * FracMom(-1.0/6.0, moments) +
	        pow(V0, 3.0/6.0) * FracMom( 1.0/6.0, moments)
	 ) + (
	        pow(V0, 1.0/6.0) * FracMom( 3.0/6.0, moments) +
	  2.0 * pow(V0,-1.0/6.0) * FracMom( 5.0/6.0, moments) +
	        pow(V0,-3.0/6.0) * FracMom( 7.0/6.0, moments)
	 );

  psi3 = (
	        pow(V0,13.0/6.0) * FracMom(-3.0/6.0, moments) +
	  2.0 * pow(V0,11.0/6.0) * FracMom(-1.0/6.0, moments) +
	        pow(V0, 9.0/6.0) * FracMom( 1.0/6.0, moments)
	 ) + 2.0 * (
	        pow(V0, 7.0/6.0) * FracMom( 3.0/6.0, moments) +
	  2.0 * pow(V0, 5.0/6.0) * FracMom( 5.0/6.0, moments) +
	        pow(V0, 3.0/6.0) * FracMom( 7.0/6.0, moments)
	 ) + (
	        pow(V0, 1.0/6.0) * FracMom( 9.0/6.0, moments) +
	  2.0 * pow(V0,-1.0/6.0) * FracMom(11.0/6.0, moments) +
	        pow(V0,-3.0/6.0) * FracMom(13.0/6.0, moments)
	 );

  psi4 = (
	        pow(V0,19.0/6.0) * FracMom(-3.0/6.0, moments) +
	  2.0 * pow(V0,17.0/6.0) * FracMom(-1.0/6.0, moments) +
	        pow(V0,15.0/6.0) * FracMom( 1.0/6.0, moments)
	 ) + 3.0 * (
	        pow(V0,13.0/6.0) * FracMom( 3.0/6.0, moments) +
	  2.0 * pow(V0,11.0/6.0) * FracMom( 5.0/6.0, moments) +
	        pow(V0, 9.0/6.0) * FracMom( 7.0/6.0, moments)
	 ) + 3.0 * (
	        pow(V0, 7.0/6.0) * FracMom( 9.0/6.0, moments) +
	  2.0 * pow(V0, 5.0/6.0) * FracMom(11.0/6.0, moments) +
	        pow(V0, 3.0/6.0) * FracMom(13.0/6.0, moments)
	 ) + (
	        pow(V0, 1.0/6.0) * FracMom(15.0/6.0, moments) +
	  2.0 * pow(V0,-1.0/6.0) * FracMom(17.0/6.0, moments) +
	        pow(V0,-3.0/6.0) * FracMom(19.0/6.0, moments)
	 );

  return C * moments[fNSootMoments-1] * pow(psi1,5.0/16.0)*pow(psi2,15.0/16.0)*pow(psi3,-5.0/16.0)*pow(psi4,1.0/16.0);
}

Double TSoot::GetPsiSL(double u, double x, double * moments, double temp)
{
  double C = GetC(temp);
  double psi1, psi2, psi3, psi4, psi5, psi6;

  double V0 = 2.0 * dimer_nbrC2;

  psi1 = (
	        pow(V0,  1.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge(-3.0/6.0+x, moments) +
	  2.0 * pow(V0, -1.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge(-1.0/6.0+x, moments) +
	        pow(V0, -3.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge( 1.0/6.0+x, moments)
	 );

  psi2 = (
	        pow(V0,  7.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge(-3.0/6.0+x, moments) +
	  2.0 * pow(V0,  5.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge(-1.0/6.0+x, moments) +
	        pow(V0,  3.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge( 1.0/6.0+x, moments)
	 ) + (
	        pow(V0,  1.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge( 3.0/6.0+x, moments) +
	  2.0 * pow(V0, -1.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge( 5.0/6.0+x, moments) +
	        pow(V0, -3.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge( 7.0/6.0+x, moments)
	 );

  psi3 = (
	        pow(V0, 13.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge(-3.0/6.0+x, moments) +
	  2.0 * pow(V0, 11.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge(-1.0/6.0+x, moments) +
	        pow(V0,  9.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge( 1.0/6.0+x, moments)
	 ) + 2.0 * (
	        pow(V0,  7.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge( 3.0/6.0+x, moments) +
	  2.0 * pow(V0,  5.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge( 5.0/6.0+x, moments) +
	        pow(V0,  3.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge( 7.0/6.0+x, moments)
	 ) + (
	        pow(V0,  1.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge( 9.0/6.0+x, moments) +
	  2.0 * pow(V0, -1.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge(11.0/6.0+x, moments) +
	        pow(V0, -3.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge(13.0/6.0+x, moments)
	 );

  if (fNSootMoments > 3)
  psi4 = (
	        pow(V0, 19.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge(-3.0/6.0+x, moments) +
	  2.0 * pow(V0, 17.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge(-1.0/6.0+x, moments) +
	        pow(V0, 15.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge( 1.0/6.0+x, moments)
	 ) + 3.0 * (
	        pow(V0, 13.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge( 3.0/6.0+x, moments) +
	  2.0 * pow(V0, 11.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge( 5.0/6.0+x, moments) +
	        pow(V0,  9.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge( 7.0/6.0+x, moments)
	 ) + 3.0 * (
	        pow(V0,  7.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge( 9.0/6.0+x, moments) +
	  2.0 * pow(V0,  5.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge(11.0/6.0+x, moments) +
	        pow(V0,  3.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge(13.0/6.0+x, moments)
	 ) + (
	        pow(V0,  1.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge(15.0/6.0+x, moments) +
	  2.0 * pow(V0, -1.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge(17.0/6.0+x, moments) +
	        pow(V0, -3.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge(19.0/6.0+x, moments)
	 );

//   if (fNSootMoments > 4)
//   psi5 = (
// 	        pow(V0, 25.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge(-3.0/6.0+x, moments) +
// 	  2.0 * pow(V0, 23.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge(-1.0/6.0+x, moments) +
// 	        pow(V0, 21.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge( 1.0/6.0+x, moments)
// 	 ) + 4.0 * (
// 	        pow(V0, 19.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge( 3.0/6.0+x, moments) +
// 	  2.0 * pow(V0, 17.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge( 5.0/6.0+x, moments) +
// 	        pow(V0, 15.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge( 7.0/6.0+x, moments)
// 	 ) + 6.0 * (
// 	        pow(V0, 13.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge( 9.0/6.0+x, moments) +
// 	  2.0 * pow(V0, 11.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge(11.0/6.0+x, moments) +
// 	        pow(V0,  9.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge(13.0/6.0+x, moments)
// 	 ) + 4.0 * (
// 	        pow(V0,  7.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge(15.0/6.0+x, moments) +
// 	  2.0 * pow(V0,  5.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge(17.0/6.0+x, moments) +
// 	        pow(V0,  3.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge(19.0/6.0+x, moments)
// 	 ) + (
// 	        pow(V0,  1.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge(21.0/6.0+x, moments) +
// 	  2.0 * pow(V0, -1.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge(23.0/6.0+x, moments) +
// 	        pow(V0, -3.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge(25.0/6.0+x, moments)
// 	 );

//   if (fNSootMoments > 5)
//   psi6 = (
// 	        pow(V0, 31.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge(-3.0/6.0+x, moments) +
// 	  2.0 * pow(V0, 29.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge(-1.0/6.0+x, moments) +
// 	        pow(V0, 27.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge( 1.0/6.0+x, moments)
// 	 ) + 5.0 * (
// 	        pow(V0, 25.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge( 3.0/6.0+x, moments) +
// 	  2.0 * pow(V0, 23.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge( 5.0/6.0+x, moments) +
// 	        pow(V0, 21.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge( 7.0/6.0+x, moments)
// 	 ) + 10.0 * (
// 	        pow(V0, 19.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge( 9.0/6.0+x, moments) +
// 	  2.0 * pow(V0, 17.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge(11.0/6.0+x, moments) +
// 	        pow(V0, 15.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge(13.0/6.0+x, moments)
// 	 ) + 10.0 * (
// 	        pow(V0, 13.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge(15.0/6.0+x, moments) +
// 	  2.0 * pow(V0, 11.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge(17.0/6.0+x, moments) +
// 	        pow(V0,  9.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge(19.0/6.0+x, moments)
// 	 ) + 5.0 * (
// 	        pow(V0,  7.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge(21.0/6.0+x, moments) +
// 	  2.0 * pow(V0,  5.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge(23.0/6.0+x, moments) +
// 	        pow(V0,  3.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge(25.0/6.0+x, moments)
// 	 ) + (
// 	        pow(V0,  1.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge(27.0/6.0+x, moments) +
// 	  2.0 * pow(V0, -1.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge(29.0/6.0+x, moments) +
// 	        pow(V0, -3.0/6.0+u) * moments[fNSootMoments-1] * FracMomLarge(31.0/6.0+x, moments)
// 	 );

//   if (fNSootMoments > 2)
//     return C * pow(psi1,3.0/8.0)*pow(psi2,3.0/4.0)*pow(psi3,-1.0/8.0);
//   else if (fNSootMoments == 4)
//     return C * pow(psi1,5.0/16.0)*pow(psi2,15.0/16.0)*pow(psi3,-5.0/16.0)*pow(psi4,1.0/16.0);
//   else if (fNSootMoments == 5)
//     return C * pow(psi1,35.0/128.0)*pow(psi2,35.0/32.0)*pow(psi3,-35.0/64.0)*pow(psi4,7.0/32.0)*pow(psi5,-5.0/128.0);
//   else if (fNSootMoments == 6)
//     return C * pow(psi1,63.0/256.0)*pow(psi2,315.0/256.0)*pow(psi3,-105.0/128.0)*pow(psi4,63.0/128.0)*pow(psi5,-45.0/256.0)*pow(psi6,7.0/256.0);

  if (u < 2.99 && x < 2.99)
    return C * pow(psi1,5.0/16.0)*pow(psi2,15.0/16.0)*pow(psi3,-5.0/16.0)*pow(psi4,1.0/16.0);
  else if (u < 3.99 && x < 3.99)
    return C * pow(psi1,3.0/8.0)*pow(psi2,3.0/4.0)*pow(psi3,-1.0/8.0);
  else
    return C * pow(psi1,1.0/2.0)*pow(psi2,1.0/2.0);
}

Double TSoot::GetPsiLL(double u, double x, double * moments, double temp)
{
  double C = GetC(temp);
  double psi1, psi2, psi3, psi4, psi5, psi6;

  psi1 = (
	        FracMomLarge( 1.0/6.0+u, moments) * FracMomLarge(-3.0/6.0+x, moments) +
	  2.0 * FracMomLarge(-1.0/6.0+u, moments) * FracMomLarge(-1.0/6.0+x, moments) +
	        FracMomLarge(-3.0/6.0+u, moments) * FracMomLarge( 1.0/6.0+x, moments)
	 );

  psi2 = (
	        FracMomLarge( 7.0/6.0+u, moments) * FracMomLarge(-3.0/6.0+x, moments) +
	  2.0 * FracMomLarge( 5.0/6.0+u, moments) * FracMomLarge(-1.0/6.0+x, moments) +
	        FracMomLarge( 3.0/6.0+u, moments) * FracMomLarge( 1.0/6.0+x, moments)
	 ) + (
	        FracMomLarge( 1.0/6.0+u, moments) * FracMomLarge( 3.0/6.0+x, moments) +
	  2.0 * FracMomLarge(-1.0/6.0+u, moments) * FracMomLarge( 5.0/6.0+x, moments) +
	        FracMomLarge(-3.0/6.0+u, moments) * FracMomLarge( 7.0/6.0+x, moments)
	 );

  psi3 = (
	        FracMomLarge(13.0/6.0+u, moments) * FracMomLarge(-3.0/6.0+x, moments) +
	  2.0 * FracMomLarge(11.0/6.0+u, moments) * FracMomLarge(-1.0/6.0+x, moments) +
	        FracMomLarge( 9.0/6.0+u, moments) * FracMomLarge( 1.0/6.0+x, moments)
	 ) + 2.0 * (
	        FracMomLarge( 7.0/6.0+u, moments) * FracMomLarge( 3.0/6.0+x, moments) +
	  2.0 * FracMomLarge( 5.0/6.0+u, moments) * FracMomLarge( 5.0/6.0+x, moments) +
	        FracMomLarge( 3.0/6.0+u, moments) * FracMomLarge( 7.0/6.0+x, moments)
	 ) + (
	        FracMomLarge( 1.0/6.0+u, moments) * FracMomLarge( 9.0/6.0+x, moments) +
	  2.0 * FracMomLarge(-1.0/6.0+u, moments) * FracMomLarge(11.0/6.0+x, moments) +
	        FracMomLarge(-3.0/6.0+u, moments) * FracMomLarge(13.0/6.0+x, moments)
	 );

  if (fNSootMoments > 3)
  psi4 = (
	        FracMomLarge(19.0/6.0+u, moments) * FracMomLarge(-3.0/6.0+x, moments) +
	  2.0 * FracMomLarge(17.0/6.0+u, moments) * FracMomLarge(-1.0/6.0+x, moments) +
	        FracMomLarge(15.0/6.0+u, moments) * FracMomLarge( 1.0/6.0+x, moments)
	 ) + 3.0 * (
	        FracMomLarge(13.0/6.0+u, moments) * FracMomLarge( 3.0/6.0+x, moments) +
	  2.0 * FracMomLarge(11.0/6.0+u, moments) * FracMomLarge( 5.0/6.0+x, moments) +
	        FracMomLarge( 9.0/6.0+u, moments) * FracMomLarge( 7.0/6.0+x, moments)
	 ) + 3.0 * (
	        FracMomLarge( 7.0/6.0+u, moments) * FracMomLarge( 9.0/6.0+x, moments) +
	  2.0 * FracMomLarge( 5.0/6.0+u, moments) * FracMomLarge(11.0/6.0+x, moments) +
	        FracMomLarge( 3.0/6.0+u, moments) * FracMomLarge(13.0/6.0+x, moments)
	 ) + (
	        FracMomLarge( 1.0/6.0+u, moments) * FracMomLarge(15.0/6.0+x, moments) +
	  2.0 * FracMomLarge(-1.0/6.0+u, moments) * FracMomLarge(17.0/6.0+x, moments) +
	        FracMomLarge(-3.0/6.0+u, moments) * FracMomLarge(19.0/6.0+x, moments)
	 );

//   if (fNSootMoments > 4)
//   psi5 = (
// 	        FracMomLarge(25.0/6.0+u, moments) * FracMomLarge(-3.0/6.0+x, moments) +
// 	  2.0 * FracMomLarge(23.0/6.0+u, moments) * FracMomLarge(-1.0/6.0+x, moments) +
// 	        FracMomLarge(21.0/6.0+u, moments) * FracMomLarge( 1.0/6.0+x, moments)
// 	 ) + 4.0 * (
// 	        FracMomLarge(19.0/6.0+u, moments) * FracMomLarge( 3.0/6.0+x, moments) +
// 	  2.0 * FracMomLarge(17.0/6.0+u, moments) * FracMomLarge( 5.0/6.0+x, moments) +
// 	        FracMomLarge(15.0/6.0+u, moments) * FracMomLarge( 7.0/6.0+x, moments)
// 	 ) + 6.0 * (
// 	        FracMomLarge(13.0/6.0+u, moments) * FracMomLarge( 9.0/6.0+x, moments) +
// 	  2.0 * FracMomLarge(11.0/6.0+u, moments) * FracMomLarge(11.0/6.0+x, moments) +
// 	        FracMomLarge( 9.0/6.0+u, moments) * FracMomLarge(13.0/6.0+x, moments)
// 	 ) + 4.0 * (
// 	        FracMomLarge( 7.0/6.0+u, moments) * FracMomLarge(15.0/6.0+x, moments) +
// 	  2.0 * FracMomLarge( 5.0/6.0+u, moments) * FracMomLarge(17.0/6.0+x, moments) +
// 	        FracMomLarge( 3.0/6.0+u, moments) * FracMomLarge(19.0/6.0+x, moments)
// 	 ) + (
// 	        FracMomLarge( 1.0/6.0+u, moments) * FracMomLarge(21.0/6.0+x, moments) +
// 	  2.0 * FracMomLarge(-1.0/6.0+u, moments) * FracMomLarge(23.0/6.0+x, moments) +
// 	        FracMomLarge(-3.0/6.0+u, moments) * FracMomLarge(25.0/6.0+x, moments)
// 	 );

//   if (fNSootMoments > 5)
//   psi6 = (
// 	        FracMomLarge(31.0/6.0+u, moments) * FracMomLarge(-3.0/6.0+x, moments) +
// 	  2.0 * FracMomLarge(29.0/6.0+u, moments) * FracMomLarge(-1.0/6.0+x, moments) +
// 	        FracMomLarge(27.0/6.0+u, moments) * FracMomLarge( 1.0/6.0+x, moments)
// 	 ) + 5.0 * (
// 	        FracMomLarge(25.0/6.0+u, moments) * FracMomLarge( 3.0/6.0+x, moments) +
// 	  2.0 * FracMomLarge(23.0/6.0+u, moments) * FracMomLarge( 5.0/6.0+x, moments) +
// 	        FracMomLarge(21.0/6.0+u, moments) * FracMomLarge( 7.0/6.0+x, moments)
// 	 ) + 10.0 * (
// 	        FracMomLarge(19.0/6.0+u, moments) * FracMomLarge( 9.0/6.0+x, moments) +
// 	  2.0 * FracMomLarge(17.0/6.0+u, moments) * FracMomLarge(11.0/6.0+x, moments) +
// 	        FracMomLarge(15.0/6.0+u, moments) * FracMomLarge(13.0/6.0+x, moments)
// 	 ) + 10.0 * (
// 	        FracMomLarge(13.0/6.0+u, moments) * FracMomLarge(15.0/6.0+x, moments) +
// 	  2.0 * FracMomLarge(11.0/6.0+u, moments) * FracMomLarge(17.0/6.0+x, moments) +
// 	        FracMomLarge( 9.0/6.0+u, moments) * FracMomLarge(19.0/6.0+x, moments)
// 	 ) + 5.0 * (
// 	        FracMomLarge( 7.0/6.0+u, moments) * FracMomLarge(21.0/6.0+x, moments) +
// 	  2.0 * FracMomLarge( 5.0/6.0+u, moments) * FracMomLarge(23.0/6.0+x, moments) +
// 	        FracMomLarge( 3.0/6.0+u, moments) * FracMomLarge(25.0/6.0+x, moments)
// 	 ) + (
// 	        FracMomLarge( 1.0/6.0+u, moments) * FracMomLarge(27.0/6.0+x, moments) +
// 	  2.0 * FracMomLarge(-1.0/6.0+u, moments) * FracMomLarge(29.0/6.0+x, moments) +
// 	        FracMomLarge(-3.0/6.0+u, moments) * FracMomLarge(31.0/6.0+x, moments)
// 	 );

//   if (fNSootMoments > 2)
//     return C * pow(psi1,3.0/8.0)*pow(psi2,3.0/4.0)*pow(psi3,-1.0/8.0);
//   else if (fNSootMoments == 4)
//     return C * pow(psi1,5.0/16.0)*pow(psi2,15.0/16.0)*pow(psi3,-5.0/16.0)*pow(psi4,1.0/16.0);
//   else if (fNSootMoments == 5)
//     return C * pow(psi1,35.0/128.0)*pow(psi2,35.0/32.0)*pow(psi3,-35.0/64.0)*pow(psi4,7.0/32.0)*pow(psi5,-5.0/128.0);
//   else if (fNSootMoments == 6)
//     return C * pow(psi1,63.0/256.0)*pow(psi2,315.0/256.0)*pow(psi3,-105.0/128.0)*pow(psi4,63.0/128.0)*pow(psi5,-45.0/256.0)*pow(psi6,7.0/256.0);

  if (u < 2.99 && x < 2.99)
    return C * pow(psi1,5.0/16.0)*pow(psi2,15.0/16.0)*pow(psi3,-5.0/16.0)*pow(psi4,1.0/16.0);
  else if (u < 3.99 && x < 3.99)
    return C * pow(psi1,3.0/8.0)*pow(psi2,3.0/4.0)*pow(psi3,-1.0/8.0);
  else
    return C * pow(psi1,1.0/2.0)*pow(psi2,1.0/2.0);
}

Double TSoot::GetPsiCond(Double x, Double * moments, Double temp)
{
	Double C = GetC(temp);
	Double psi1, psi2;

	psi1 =	(
				FracMom( 1.0/6.0+x, moments) * pow(dimer_nbrC2, -1.0/2.0) +
			2.0 *	FracMom(-1.0/6.0+x, moments) * pow(dimer_nbrC2, -1.0/6.0) +
				FracMom(-1.0/2.0+x, moments) * pow(dimer_nbrC2,  1.0/6.0)
		);

	psi2 =	(
				FracMom( 7.0/6.0+x, moments) * pow(dimer_nbrC2, -1.0/2.0) +
			2.0 *	FracMom( 5.0/6.0+x, moments) * pow(dimer_nbrC2, -1.0/6.0) +
				FracMom( 1.0/2.0+x, moments) * pow(dimer_nbrC2,  1.0/6.0)
		) + (
				FracMom( 1.0/6.0+x, moments) * pow(dimer_nbrC2,  1.0/2.0) +
			2.0 *	FracMom(-1.0/6.0+x, moments) * pow(dimer_nbrC2,  5.0/6.0) +
				FracMom(-1.0/2.0+x, moments) * pow(dimer_nbrC2,  7.0/6.0)
		);

	return C * sqrt(psi1 * psi2);
}

double TSoot::GetPsi(double u, double v, double x, double y, double * moments, double temp)
{
  double C = GetC(temp);
  double psi1, psi2, psi3, psi4, psi5, psi6;

  double Df1 = 1.8;
  double Df2 = 1.8;

  double av1 = 1.0 - (2.0 / Df1);
  double as1 = (3.0 / Df1) - 1.0;

  double av2 = 1.0 - (2.0 / Df2);
  double as2 = (3.0 / Df2) - 1.0;

  psi1 = (
	        FracMom(2.0*av1-1.0/2.0+u, 2.0*as1+v, moments) * FracMom(       -1.0/2.0+x,         y, moments) +
	  2.0 * FracMom(    av1-1.0/2.0+u,     as1+v, moments) * FracMom(    av2-1.0/2.0+x,     as2+y, moments) +
	        FracMom(       -1.0/2.0+u,         v, moments) * FracMom(2.0*av2-1.0/2.0+x, 2.0*as2+y, moments)
	 );

  psi2 = (
	        FracMom(2.0*av1+1.0/2.0+u, 2.0*as1+v, moments) * FracMom(       -1.0/2.0+x,         y, moments) +
	  2.0 * FracMom(    av1+1.0/2.0+u,     as1+v, moments) * FracMom(    av2-1.0/2.0+x,     as2+y, moments) +
	        FracMom(        1.0/2.0+u,         v, moments) * FracMom(2.0*av2-1.0/2.0+x, 2.0*as2+y, moments)
	 ) + (
	        FracMom(2.0*av1-1.0/2.0+u, 2.0*as1+v, moments) * FracMom(        1.0/2.0+x,         y, moments) +
	  2.0 * FracMom(    av1-1.0/2.0+u,     as1+v, moments) * FracMom(    av2+1.0/2.0+x,     as2+y, moments) +
	        FracMom(       -1.0/2.0+u,         v, moments) * FracMom(2.0*av2+1.0/2.0+x, 2.0*as2+y, moments)
	 );

  psi3 = (
	        FracMom(2.0*av1+3.0/2.0+u, 2.0*as1+v, moments) * FracMom(       -1.0/2.0+x,         y, moments) +
	  2.0 * FracMom(    av1+3.0/2.0+u,     as1+v, moments) * FracMom(    av2-1.0/2.0+x,     as2+y, moments) +
	        FracMom(        3.0/2.0+u,         v, moments) * FracMom(2.0*av2-1.0/2.0+x, 2.0*as2+y, moments)
	 ) + 2.0 * (
		FracMom(2.0*av1+1.0/2.0+u, 2.0*as1+v, moments) * FracMom(        1.0/2.0+x,         y, moments) +
	  2.0 * FracMom(    av1+1.0/2.0+u,     as1+v, moments) * FracMom(    av2+1.0/2.0+x,     as2+y, moments) +
	        FracMom(        1.0/2.0+u,         v, moments) * FracMom(2.0*av2+1.0/2.0+x, 2.0*as2+y, moments)
	 ) + (
	        FracMom(2.0*av1-1.0/2.0+u, 2.0*as1+v, moments) * FracMom(        3.0/2.0+x,         y, moments) +
	  2.0 * FracMom(    av1-1.0/2.0+u,     as1+v, moments) * FracMom(    av2+3.0/2.0+x,     as2+y, moments) +
	        FracMom(       -1.0/2.0+u,         v, moments) * FracMom(2.0*av2+3.0/2.0+x, 2.0*as2+y, moments)
	 );

  if (fNSootMoments > 3)
  psi4 = (
	        FracMom(2.0*av1+5.0/2.0+u, 2.0*as1+v, moments) * FracMom(       -1.0/2.0+x,         y, moments) +
	  2.0 * FracMom(    av1+5.0/2.0+u,     as1+v, moments) * FracMom(    av2-1.0/2.0+x,     as2+y, moments) +
	        FracMom(        5.0/2.0+u,         v, moments) * FracMom(2.0*av2-1.0/2.0+x, 2.0*as2+y, moments)
	 ) + 3.0 * (
	        FracMom(2.0*av1+3.0/2.0+u, 2.0*as1+v, moments) * FracMom(        1.0/2.0+x,         y, moments) +
	  2.0 * FracMom(    av1+3.0/2.0+u,     as1+v, moments) * FracMom(    av2+1.0/2.0+x,     as2+y, moments) +
	        FracMom(        3.0/2.0+u,         v, moments) * FracMom(2.0*av2+1.0/2.0+x, 2.0*as2+y, moments)
	 ) + 3.0 * (
	        FracMom(2.0*av1+1.0/2.0+u, 2.0*as1+v, moments) * FracMom(        3.0/2.0+x,         y, moments) +
	  2.0 * FracMom(    av1+1.0/2.0+u,     as1+v, moments) * FracMom(    av2+3.0/2.0+x,     as2+y, moments) +
	        FracMom(        1.0/2.0+u,         v, moments) * FracMom(2.0*av2+3.0/2.0+x, 2.0*as2+y, moments)
	 ) + (
	        FracMom(2.0*av1-1.0/2.0+u, 2.0*as1+v, moments) * FracMom(        5.0/2.0+x,         y, moments) +
	  2.0 * FracMom(    av1-1.0/2.0+u,     as1+v, moments) * FracMom(    av2+5.0/2.0+x,     as2+y, moments) +
	        FracMom(       -1.0/2.0+u,         v, moments) * FracMom(2.0*av2+5.0/2.0+x, 2.0*as2+y, moments)
	 );

  if (fNSootMoments > 6)
  psi5 = (
	        FracMom(2.0*av1+7.0/2.0+u, 2.0*as1+v, moments) * FracMom(       -1.0/2.0+x,         y, moments) +
	  2.0 * FracMom(    av1+7.0/2.0+u,     as1+v, moments) * FracMom(    av2-1.0/2.0+x,     as2+y, moments) +
	        FracMom(        7.0/2.0+u,         v, moments) * FracMom(2.0*av2-1.0/2.0+x, 2.0*as2+y, moments)
	 ) + 4.0 * (
	        FracMom(2.0*av1+5.0/2.0+u, 2.0*as1+v, moments) * FracMom(        1.0/2.0+x,         y, moments) +
	  2.0 * FracMom(    av1+5.0/2.0+u,     as1+v, moments) * FracMom(    av2+1.0/2.0+x,     as2+y, moments) +
	        FracMom(        5.0/2.0+u,         v, moments) * FracMom(2.0*av2+1.0/2.0+x, 2.0*as2+y, moments)
	 ) + 6.0 * (
	        FracMom(2.0*av1+3.0/2.0+u, 2.0*as1+v, moments) * FracMom(        3.0/2.0+x,         y, moments) +
	  2.0 * FracMom(    av1+3.0/2.0+u,     as1+v, moments) * FracMom(    av2+3.0/2.0+x,     as2+y, moments) +
	        FracMom(        3.0/2.0+u,         v, moments) * FracMom(2.0*av2+3.0/2.0+x, 2.0*as2+y, moments)
	 ) + 4.0 * (
	        FracMom(2.0*av1+1.0/2.0+u, 2.0*as1+v, moments) * FracMom(        5.0/2.0+x,         y, moments) +
	  2.0 * FracMom(    av1+1.0/2.0+u,     as1+v, moments) * FracMom(    av2+5.0/2.0+x,     as2+y, moments) +
	        FracMom(        1.0/2.0+u,         v, moments) * FracMom(2.0*av2+5.0/2.0+x, 2.0*as2+y, moments)
	 ) + (
	        FracMom(2.0*av1-1.0/2.0+u, 2.0*as1+v, moments) * FracMom(        7.0/2.0+x,         y, moments) +
	  2.0 * FracMom(    av1-1.0/2.0+u,     as1+v, moments) * FracMom(    av2+7.0/2.0+x,     as2+y, moments) +
	        FracMom(       -1.0/2.0+u,         v, moments) * FracMom(2.0*av2+7.0/2.0+x, 2.0*as2+y, moments)
	 );

  if (fNSootMoments > 10)
  psi6 = (
	        FracMom(2.0*av1+9.0/2.0+u, 2.0*as1+v, moments) * FracMom(       -1.0/2.0+x,         y, moments) +
	  2.0 * FracMom(    av1+9.0/2.0+u,     as1+v, moments) * FracMom(    av2-1.0/2.0+x,     as2+y, moments) +
	        FracMom(        9.0/2.0+u,         v, moments) * FracMom(2.0*av2-1.0/2.0+x, 2.0*as2+y, moments)
	 ) + 5.0 * (
	        FracMom(2.0*av1+7.0/2.0+u, 2.0*as1+v, moments) * FracMom(        1.0/2.0+x,         y, moments) +
	  2.0 * FracMom(    av1+7.0/2.0+u,     as1+v, moments) * FracMom(    av2+1.0/2.0+x,     as2+y, moments) +
	        FracMom(        7.0/2.0+u,         v, moments) * FracMom(2.0*av2+1.0/2.0+x, 2.0*as2+y, moments)
	 ) + 10.0 * (
	        FracMom(2.0*av1+5.0/2.0+u, 2.0*as1+v, moments) * FracMom(        3.0/2.0+x,         y, moments) +
	  2.0 * FracMom(    av1+5.0/2.0+u,     as1+v, moments) * FracMom(    av2+3.0/2.0+x,     as2+y, moments) +
	        FracMom(        5.0/2.0+u,         v, moments) * FracMom(2.0*av2+3.0/2.0+x, 2.0*as2+y, moments)
	 ) + 10.0 * (
	        FracMom(2.0*av1+3.0/2.0+u, 2.0*as1+v, moments) * FracMom(        5.0/2.0+x,         y, moments) +
	  2.0 * FracMom(    av1+3.0/2.0+u,     as1+v, moments) * FracMom(    av2+5.0/2.0+x,     as2+y, moments) +
	        FracMom(        3.0/2.0+u,         v, moments) * FracMom(2.0*av2+5.0/2.0+x, 2.0*as2+y, moments)
	 ) + 5.0 * (
	        FracMom(2.0*av1+1.0/2.0+u, 2.0*as1+v, moments) * FracMom(        7.0/2.0+x,         y, moments) +
	  2.0 * FracMom(    av1+1.0/2.0+u,     as1+v, moments) * FracMom(    av2+7.0/2.0+x,     as2+y, moments) +
	        FracMom(        1.0/2.0+u,         v, moments) * FracMom(2.0*av2+7.0/2.0+x, 2.0*as2+y, moments)
	 ) + (
	        FracMom(2.0*av1-1.0/2.0+u, 2.0*as1+v, moments) * FracMom(        9.0/2.0+x,         y, moments) +
	  2.0 * FracMom(    av1-1.0/2.0+u,     as1+v, moments) * FracMom(    av2+9.0/2.0+x,     as2+y, moments) +
	        FracMom(       -1.0/2.0+u,         v, moments) * FracMom(2.0*av2+9.0/2.0+x, 2.0*as2+y, moments)
	 );

  //return C * sqrt(psi1 * psi2);
  if (fNSootMoments == 3)
    return C * pow(psi1,3.0/8.0)*pow(psi2,3.0/4.0)*pow(psi3,-1.0/8.0);
  else if (fNSootMoments == 6)
    return C * pow(psi1,5.0/16.0)*pow(psi2,15.0/16.0)*pow(psi3,-5.0/16.0)*pow(psi4,1.0/16.0);
  else if (fNSootMoments == 10)
    return C * pow(psi1,35.0/128.0)*pow(psi2,35.0/32.0)*pow(psi3,-35.0/64.0)*pow(psi4,7.0/32.0)*pow(psi5,-5.0/128.0);
  else if (fNSootMoments == 15)
    return C * pow(psi1,63.0/256.0)*pow(psi2,315.0/256.0)*pow(psi3,-105.0/128.0)*pow(psi4,63.0/128.0)*pow(psi5,-45.0/256.0)*pow(psi6,7.0/256.0);

  return 0.0;
}

double TSoot::GetPsi(double u, double v, double w, double x, double y, double z, double * moments, double temp)
{
	Double C = GetC(temp);
	Double psi1, psi2;

	Double Df1 = 3.0;
	Double Df2 = 3.0;

	Double av1 = 1.0 - (2.0 / Df1);
	Double as1 = (3.0 / Df1) - 1.0;

	Double av2 = 1.0 - (2.0 / Df2);
	Double as2 = (3.0 / Df2) - 1.0;

	psi1 =	(
				FracMom(2.0*av1-1.0/2.0+u, 2.0*as1+v, w, moments) * FracMom(-1.0/2.0+x, y, z, moments) +
			2.0 *	FracMom(av1-1.0/2.0+u, as1+v, w, moments) * FracMom(av2-1.0/2.0+x, as2+y, z, moments) +
				FracMom(-1.0/2.0+u, v, w, moments) * FracMom(2.0*av2-1.0/2.0+x, 2.0*as2+y, z, moments)
		);

	psi2 =	(
				FracMom(2.0*av1+1.0/2.0+u, 2.0*as1+v, w, moments) * FracMom(-1.0/2.0+x, y, z, moments) +
			2.0 *	FracMom(av1+1.0/2.0+u, as1+v, w, moments) * FracMom(av2-1.0/2.0+x, as2+y, z, moments) +
				FracMom(1.0/2.0+u, v, w, moments) * FracMom(2.0*av2-1.0/2.0+x, 2.0*as2+y, z, moments)
		) + (
				FracMom(2.0*av1-1.0/2.0+u, 2.0*as1+v, w, moments) * FracMom(1.0/2.0+x, y, z, moments) +
			2.0 *	FracMom(av1-1.0/2.0+u, as1+v, w, moments) * FracMom(av2+1.0/2.0+x, as2+y, z, moments) +
				FracMom(-1.0/2.0+u, v, w, moments) * FracMom(2.0*av2+1.0/2.0+x, 2.0*as2+y, z, moments)
		);

	return C * sqrt(psi1 * psi2);
}


//Reaction Rates (Surface Growth and Oxidation)
double TSoot::SurfaceGrowthCoeffFor(double * Y, double temp, double density, double * molarMass)
{
#ifndef GLOBAL_SOOT
  double C_C2H2 = density * Y[f_C2H2] / molarMass[f_C2H2];
  double C_H = density * Y[f_H] / molarMass[f_H];
  //double locSG = (fSootRateCoeffs->vec[ks10f] * fThirdBodyConc * C_C2H2 * fCSootStar + fSootRateCoeffs->vec[ks13b] * C_H) * fFk10;
  double locSG = fSootRateCoeffs->vec[ks4] * C_C2H2 * fCSootStar; // + fSootRateCoeffs->vec[ks13b] * C_H) * fFk10;

  return locSG;
#else
  //Concentrations in SI units (mol/m^3)
  double C_C2H2 = density * Y[f_C2H2] / molarMass[f_C2H2] * 1.0e3;
  return 3.3e13*pow(C_C2H2,0.725)*exp(-3.5e4/temp);
#endif
}

double TSoot::CalcSurfaceGrowthCoeff(TReactionPtr reaction, double * Y, double temp, double density, double * molarMass, double mixMolarMass)
{
  double C_C2H2 = density * Y[f_C2H2] / molarMass[f_C2H2];

  double * k = fSootRateCoeffs->vec;

  ComputeSootRateCoeffs(k, temp, reaction);
  ComputeCSootStar(k, Y, density, molarMass, mixMolarMass);

  double locSG = k[ks4] * C_C2H2 * fCSootStar;

  return locSG;
}

double TSoot::SurfaceGrowthCoeffBack(double * Y, double temp, double density, double * molarMass)
{
  double C_H = density * Y[f_H] / molarMass[f_H];
  //double locSG = fSootRateCoeffs->vec[ks13b] * C_H;

  //return locSG;
  return 0.0;
}

double TSoot::OxidationCoeff(double * Y, double temp, double density, double * molarMass)
{
#ifndef GLOBAL_SOOT
  double C_OH = density * Y[f_OH] / molarMass[f_OH];
  double C_O2 = density * Y[f_O2] / molarMass[f_O2];
  double * k = fSootRateCoeffs->vec;

  double Chi = fChi;

  return k[ks5] * C_O2 * fCSootStar + (1.0/ (2.0 * Chi)) * k[ks6] * C_OH;
#else
  //Concentrations in SI units (mol/m^3)
  double C_O2 = density * Y[f_O2] / molarMass[f_O2] * 1.0e3;
  return 7.2e5*pow(C_O2,1.429)*exp(-7.3e3/temp) + 2.9e7*pow(C_O2,1.177)*exp(-1.0e4/temp);
#endif
}

double TSoot::OxidationCoeff_O2(double * Y, double temp, double density, double * molarMass)
{
#ifndef GLOBAL_SOOT
  double C_O2 = density * Y[f_O2] / molarMass[f_O2];
  double * k = fSootRateCoeffs->vec;

  return k[ks5] * C_O2 * fCSootStar;
#else
  //Concentrations in SI units (mol/m^3)
  double C_O2 = density * Y[f_O2] / molarMass[f_O2] * 1.0e3;
  return 7.2e5*pow(C_O2,1.429)*exp(-7.3e3/temp);
#endif
}

double TSoot::OxidationCoeff_OH(double * Y, double temp, double density, double * molarMass)
{
#ifndef GLOBAL_SOOT
  double C_OH = density * Y[f_OH] / molarMass[f_OH];
  double * k = fSootRateCoeffs->vec;

  double Chi = fChi;

  return (1.0/ (2.0 * Chi)) * k[ks6] * C_OH;
#else
  //Concentrations in SI units (mol/m^3)
  double C_O2 = density * Y[f_O2] / molarMass[f_O2] * 1.0e3;
  return 2.9e7*pow(C_O2,1.177)*exp(-1.0e4/temp);
#endif
}

double TSoot::CalcOxidationCoeff(TReactionPtr reaction, double * Y, double temp, double density, double * molarMass, double mixMolarMass)
{
  double C_OH = density * Y[f_OH] / molarMass[f_OH];
  double C_O2 = density * Y[f_O2] / molarMass[f_O2];

  double * k = fSootRateCoeffs->vec;

  ComputeSootRateCoeffs(k, temp, reaction);
  ComputeCSootStar(k, Y, density, molarMass, mixMolarMass);

  double Chi = fChi;

  return k[ks5] * C_O2 * fCSootStar + (1.0/ (2.0 * Chi)) * k[ks6] * C_OH;
}

double TSoot::CalcOxidationCoeff_O2(TReactionPtr reaction, double * Y, double temp, double density, double * molarMass, double mixMolarMass)
{
  double C_O2 = density * Y[f_O2] / molarMass[f_O2];

  double * k = fSootRateCoeffs->vec;

  ComputeSootRateCoeffs(k, temp, reaction);
  ComputeCSootStar(k, Y, density, molarMass, mixMolarMass);

  double Chi = fChi;

  return k[ks5] * C_O2 * fCSootStar;
}



//**************************************//
//****************Closure***************//
//**************************************//
/*
double fn1(double u, double s, double * moments)
{
  return exp(u+1.0/2.0*s*s)*(1.0+erf((s*s+u-log(20.0))/(sqrt(2.0)*s)))/(1.0+erf((u-log(20.0))/(sqrt(2.0)*s))) - (moments[1]/moments[0]);
}

double fn2(double u, double s, double * moments)
{
  return exp(2.0*u+2.0*s*s)*(1.0+erf((2.0*s*s+u-log(20.0))/(sqrt(2.0)*s)))/(1.0+erf((u-log(20.0))/(sqrt(2.0)*s))) - (moments[2]/moments[0]);
}
*/
Double TSoot::FracMom(Double Index, Double * moments)
{
#ifdef MOMIC
double temp;
  double fracmom = moments[0];

  for (int i = 1; i < fNSootMoments; i++)
  {
    temp = 1.0;

    for (int j = 0; j < fNSootMoments; j++)
    {
      if (i != j)
	temp *= (Index - (double(Geti(j))/6.0)) / ((double(Geti(i))/6.0) - (double(Geti(j))/6.0));
    }
 
    fracmom *= pow(moments[i]/moments[0], temp);
  }

  return fracmom;
#endif
#ifdef FRENKLACH
  double temp;
  double fracmom = moments[0];

  if (Index >= 0.0)
  {
    for (int i = 1; i < fNSootMoments-1; i++)
    {
      temp = 1.0;

      for (int j = 0; j < fNSootMoments-1; j++)
      {
	if (i != j)
	  temp *= (Index - (double(Geti(j))/6.0)) / ((double(Geti(i))/6.0) - (double(Geti(j))/6.0));
      }

      fracmom *= pow(moments[i]/moments[0], temp);
    }
  }
  else
  {
    fracmom = moments[fNSootMoments-1] * 
      pow(pow(moments[fNSootMoments-1],-3.0/2.0)*pow(moments[0], 2.0)*pow(moments[1],-1.0/2.0),pow(2.0,Index)) * 
      pow(pow(moments[fNSootMoments-1], 1.0/2.0)*pow(moments[0],-1.0)*pow(moments[1], 1.0/2.0),pow(4.0,Index));
}

  return fracmom;
#endif

#ifdef HMOM
  double temp;

  double * momentslarge = new double[fNSootMoments-1];

  for (int i = 0; i < fNSootMoments-1; i++)
  {
    momentslarge[i] = moments[i] - moments[fNSootMoments-1] * pow(2.0*dimer_nbrC2, double(Geti(i))/6.0);
    momentslarge[i] = MAX(momentslarge[i], 1.0e-60);
  }

  double fracmom = momentslarge[0];

  for (int i = 1; i < fNSootMoments-1; i++)
  {
    temp = 1.0;

    for (int j = 0; j < fNSootMoments-1; j++)
    {
      if (i != j)
	temp *= (Index - (double(Geti(j))/6.0)) / ((double(Geti(i))/6.0) - (double(Geti(j))/6.0));
    }
 
    fracmom *= pow(momentslarge[i] / momentslarge[0], temp);
  }
  /*
  double mu = 2.0 * log(momentslarge[1]) - 3.0/2.0 * log(momentslarge[0]) - 1.0/2.0 * log(momentslarge[2]);
  double sigma2 = log(momentslarge[2]) + log(momentslarge[0]) - 2.0 * log(momentslarge[1]);

  if (sigma2 <= 0.5)
    ;
  else
  {
    //fracmom *= 1.0/2.0 * (1.0 + erf((Index*sigma2 + mu - log(2.0*dimer_nbrC2)) / (sqrt(2.0 * sigma2))));

    double uold = 0;
    double u = 2.0*log(momentslarge[1]) - 3.0/2.0*log(momentslarge[0]) - 1.0/2.0*log(momentslarge[2]);

    double sold = 0;
    double s = sqrt(log(momentslarge[2]) + log(momentslarge[0]) - 2.0*log(momentslarge[1]));

    double del = 1.0e-6;

    double J11, J12, J21, J22, Ji11, Ji12, Ji21, Ji22, f1, f2, det;

    while (abs(uold-u)/abs(u) > 0.0001 && abs(sold-s)/abs(s) > 0.0001)
    {
      uold = u;
      sold = s;

      J11 = (fn1(u+del,s,momentslarge) - fn1(u-del,s,momentslarge)) / (2.0*del);
      J12 = (fn1(u,s+del,momentslarge) - fn1(u,s-del,momentslarge)) / (2.0*del);
      J21 = (fn2(u+del,s,momentslarge) - fn2(u-del,s,momentslarge)) / (2.0*del);
      J22 = (fn2(u,s+del,momentslarge) - fn2(u,s-del,momentslarge)) / (2.0*del);

      det = J11*J22 - J12*J21;

      f1 = fn1(u,s,momentslarge);
      f2 = fn2(u,s,momentslarge);

      Ji11 =  J22 / det;
      Ji12 = -J21 / det;
      Ji21 = -J12 / det;
      Ji22 =  J11 / det;

      u = uold - Ji11*f1 - Ji12*f2;
      s = sold - Ji21*f1 - Ji22*f2;
    }

    fracmom = momentslarge[0]*exp(Index*u+1.0/2.0*Index*Index*s*s)*(1.0+erf((Index*s*s+u-log(2.0*dimer_nbrC2))/(sqrt(2.0)*s)))/(1.0+erf((u-log(2.0*dimer_nbrC2))/(sqrt(2.0)*s)));
  }
  */
  fracmom += moments[fNSootMoments-1] * pow(2.0*dimer_nbrC2, Index);

  delete [] momentslarge;

  return fracmom;
#endif
}

double TSoot::FracMomLarge(double Index, double * moments)
{
  double temp = FracMom(Index, moments) - moments[fNSootMoments-1] * pow(2.0*dimer_nbrC2, Index);

  if (temp <= 0)
    return pow(2.0*dimer_nbrC2, Index) * MAX(moments[0]-moments[fNSootMoments-1],1.0e-60);

  return temp;
}

Double TSoot::FracMom(Double Index1, Double Index2, Double * moments)
{
  if (fNSootMoments == 3)
    return pow(moments[0], 1.0-Index1-Index2) * pow(moments[1], Index1) * pow(moments[2], Index2);

  if (fNSootMoments == 4) 
  {
    double M00 = moments[0] - moments[3];
    double M10 = moments[1] - 2.0*dimer_nbrC2*moments[3];
    double M01 = moments[2] - pow(2.0*dimer_nbrC2, 2.0/3.0)*moments[3];

    double fracmom;
    double firstpeak;
    double bothpeaks;

    firstpeak = pow(moments[0], 1.0-Index1-Index2) * pow(moments[1], Index1) * pow(moments[2], Index2);
    bothpeaks = moments[3] * pow(2.0*dimer_nbrC2, Index1+(2.0/3.0)*Index2) + pow(M00, 1.0-Index1-Index2) * pow(M10, Index1) * pow(M01, Index2);

    //Check for a delta function
    if (M00 < 1.0e-20 || M10 < 1.0e-20 || M01 < 1.0e-20)
      fracmom = firstpeak;
    else
      fracmom = bothpeaks;
    
    if (fracmom != fracmom)
    {
      cerr << "The fractional moment is NaN...\n";
      cerr << "M00:" << "\t" << moments[0] << endl;
      cerr << "M10:" << "\t" << moments[1] << endl;
      cerr << "M01:" << "\t" << moments[2] << endl;
      cerr << "N0:" << "\t" << moments[3] << endl;
      cerr << "M'00:" << "\t" << M00 << endl;
      cerr << "M'10:" << "\t" << M10 << endl;
      cerr << "M'01:" << "\t" << M01 << endl;
      cerr << "First Peak: " << "\t" << firstpeak << endl;
      cerr << "Both Peaks: " << "\t"<< bothpeaks << endl;
      cerr << "Frac Mom(" << Index1 << "," << Index2 << "):" << "\t" << fracmom << endl;
    }

    return fracmom;
  }

  int deltaextrap = 0;

  if (fNSootMoments == 6)
  {
    if (!deltaextrap)
    {
      //Analytical Inversion: (0,0), (1,0), (0,1), (2,0), (1,1), (0,2)
      double fact1 = moments[0];
      double fact2 = pow(pow(moments[1], 2.0) * pow(moments[0], -1.5) * pow(moments[3], -0.5), Index1);
      double fact3 = pow(pow(moments[2], 2.0) * pow(moments[0], -1.5) * pow(moments[5], -0.5), Index2);
      double fact4 = pow(pow(moments[3], 0.5) * pow(moments[0], 0.5) * pow(moments[1], -1.0), Index1*Index1);
      double fact5 = pow((moments[4] * moments[0]) / (moments[1] * moments[2]), Index1*Index2);
      double fact6 = pow(pow(moments[5], 0.5) * pow(moments[0], 0.5) * pow(moments[2], -1.0), Index2*Index2);

      //Analytical Inversion: (0,0), (1,0), (0,1), (h,0), (h,h), (0,h)
      //double fact1 = moments[0];
      //double fact2 = pow(pow(moments[3], 4.0) * pow(moments[0], -3.0) / moments[1], Index1);
      //double fact3 = pow(pow(moments[5], 4.0) * pow(moments[0], -3.0) / moments[2], Index2);
      //double fact4 = pow(pow(moments[0], 2.0) * pow(moments[1], 2.0) * pow(moments[3], -4.0), Index1*Index1);
      //double fact5 = pow((moments[0] * moments[4]) / (moments[3] * moments[5]), 4.0*Index1*Index2);
      //double fact6 = pow(pow(moments[0], 2.0) * pow(moments[2], 2.0) * pow(moments[5], -4.0), Index2*Index2);

      return fact1 * fact2 * fact3 * fact4 * fact5 * fact6;
    }
    else
    {
      //Analytical Inversion: (0,0), (1,0), (0,1), (2,0), (1,1), (0,2)
      if ((Index1+Index2) <= 2.0)
      {
	double fact1 = moments[0];
	double fact2 = pow(pow(moments[1], 2.0) * pow(moments[0], -1.5) * pow(moments[3], -0.5), Index1);
	double fact3 = pow(pow(moments[2], 2.0) * pow(moments[0], -1.5) * pow(moments[5], -0.5), Index2);
	double fact4 = pow(pow(moments[3], 0.5) * pow(moments[0], 0.5) * pow(moments[1], -1.0), Index1*Index1);
	double fact5 = pow((moments[4] * moments[0]) / (moments[1] * moments[2]), Index1*Index2);
	double fact6 = pow(pow(moments[5], 0.5) * pow(moments[0], 0.5) * pow(moments[2], -1.0), Index2*Index2);

	return fact1 * fact2 * fact3 * fact4 * fact5 * fact6;
      }
      else
      {
	if (Index1 <= 1.0)
	  return pow(moments[2], 2.0-Index1-Index2) * pow(moments[4], Index1) * pow(moments[5], Index2-1.0);
	else if (Index2 <= 1.0)
	  return pow(moments[1], 2.0-Index1-Index2) * pow(moments[3], Index1-1.0) * pow(moments[4], Index2);
	else
	{
	  double momleft = pow(moments[2], 1.0-Index2) * moments[4] * pow(moments[5], Index2-1.0);
	  double momdown = pow(moments[1], 1.0-Index1) * pow(moments[3], Index1-1.0) * moments[4];

	  return momleft * momdown / moments[4];
	}
      }
    }
  }

  if (fNSootMoments == 7)
  {
    double M00 = moments[0] - moments[6];
    double M10 = moments[1] - 2.0*dimer_nbrC2*moments[6];
    double M01 = moments[2] - pow(2.0*dimer_nbrC2, 2.0/3.0)*moments[6];

    double M20 = moments[3] - pow(2.0*dimer_nbrC2, 2.0)*moments[6]; //M20
    double M11 = moments[4] - pow(2.0*dimer_nbrC2, 5.0/3.0)*moments[6]; //M11
    double M02 = moments[5] - pow(2.0*dimer_nbrC2, 4.0/3.0)*moments[6]; //M02

    //double M20 = moments[3] - pow(2.0*dimer_nbrC2, 1.0/2.0)*moments[6]; //Mh0
    //double M11 = moments[4] - pow(2.0*dimer_nbrC2, 5.0/6.0)*moments[6]; //Mhh
    //double M02 = moments[5] - pow(2.0*dimer_nbrC2, 1.0/3.0)*moments[6]; //M0h

    double fracmom;
    double firstpeak;
    double bothpeaks;

    if (!deltaextrap)
    {
      //Analytical Inversion: (0,0), (1,0), (0,1), (2,0), (1,1), (0,2)
      double fact1 = moments[0];
      double fact2 = pow(pow(moments[1], 2.0) * pow(moments[0], -1.5) * pow(moments[3], -0.5), Index1);
      double fact3 = pow(pow(moments[2], 2.0) * pow(moments[0], -1.5) * pow(moments[5], -0.5), Index2);
      double fact4 = pow(pow(moments[3], 0.5) * pow(moments[0], 0.5) * pow(moments[1], -1.0), Index1*Index1);
      double fact5 = pow((moments[4] * moments[0]) / (moments[1] * moments[2]), Index1*Index2);
      double fact6 = pow(pow(moments[5], 0.5) * pow(moments[0], 0.5) * pow(moments[2], -1.0), Index2*Index2);

      //Analytical Inversion: (0,0), (1,0), (0,1), (h,0), (h,h), (0,h)
      //double fact1 = moments[0];
      //double fact2 = pow(pow(moments[3], 4.0) * pow(moments[0], -3.0) / moments[1], Index1);
      //double fact3 = pow(pow(moments[5], 4.0) * pow(moments[0], -3.0) / moments[2], Index2);
      //double fact4 = pow(pow(moments[0], 2.0) * pow(moments[1], 2.0) * pow(moments[3], -4.0), Index1*Index1);
      //double fact5 = pow((moments[0] * moments[4]) / (moments[3] * moments[5]), 4.0*Index1*Index2);
      //double fact6 = pow(pow(moments[0], 2.0) * pow(moments[2], 2.0) * pow(moments[5], -4.0), Index2*Index2);

      firstpeak = fact1 * fact2 * fact3 * fact4 * fact5 * fact6;
    }
    else
    {
      //if ((Index1+Index2) <= 2.0)
      if ((Index1+Index2) <= 1.0)
      {
	//Analytical Inversion: (0,0), (1,0), (0,1), (2,0), (1,1), (0,2)
	double fact1 = moments[0];
	double fact2 = pow(pow(moments[1], 2.0) * pow(moments[0], -1.5) * pow(moments[3], -0.5), Index1);
	double fact3 = pow(pow(moments[2], 2.0) * pow(moments[0], -1.5) * pow(moments[5], -0.5), Index2);
	double fact4 = pow(pow(moments[3], 0.5) * pow(moments[0], 0.5) * pow(moments[1], -1.0), Index1*Index1);
	double fact5 = pow((moments[4] * moments[0]) / (moments[1] * moments[2]), Index1*Index2);
	double fact6 = pow(pow(moments[5], 0.5) * pow(moments[0], 0.5) * pow(moments[2], -1.0), Index2*Index2);

	//Analytical Inversion: (0,0), (1,0), (0,1), (h,0), (h,h), (0,h)
	//double fact1 = moments[0];
	//double fact2 = pow(pow(moments[3], 4.0) * pow(moments[0], -3.0) / moments[1], Index1);
	//double fact3 = pow(pow(moments[5], 4.0) * pow(moments[0], -3.0) / moments[2], Index2);
	//double fact4 = pow(pow(moments[0], 2.0) * pow(moments[1], 2.0) * pow(moments[3], -4.0), Index1*Index1);
	//double fact5 = pow((moments[0] * moments[4]) / (moments[3] * moments[5]), 4.0*Index1*Index2);
	//double fact6 = pow(pow(moments[0], 2.0) * pow(moments[2], 2.0) * pow(moments[5], -4.0), Index2*Index2);

	firstpeak = fact1 * fact2 * fact3 * fact4 * fact5 * fact6;
      }
      else
      {
	if (Index1 <= 1.0) //M2
	  firstpeak = pow(moments[2], 2.0-Index1-Index2) * pow(moments[4], Index1) * pow(moments[5], Index2-1.0); //M2
	//if (Index1 <= 0.5)
	//  firstpeak = pow(moments[5], 2.0-2.0*Index1-2.0*Index2) * pow(moments[4], 2.0*Index1) * pow(moments[2], 2.0*Index2-1.0);
	else if (Index2 <= 1.0) //M2
	  firstpeak = pow(moments[1], 2.0-Index1-Index2) * pow(moments[3], Index1-1.0) * pow(moments[4], Index2); //M2
	//else if (Index2 <= 0.5)
	//  firstpeak = pow(moments[3], 2.0-2.0*Index1-2.0*Index2) * pow(moments[1], 2.0*Index1-1.0) * pow(moments[4], 2.0*Index2);
	else
	{
	  double momleft = pow(moments[2], 1.0-Index2) * moments[4] * pow(moments[5], Index2-1.0); //M2
	  //double momleft = pow(moments[5], 1.0-2.0*Index2) * moments[4] * pow(moments[2], 2.0*Index2-1.0);
	  double momdown = pow(moments[1], 1.0-Index1) * pow(moments[3], Index1-1.0) * moments[4]; //M2
	  //double momdown = pow(moments[3], 1.0-2.0*Index1) * pow(moments[1], 2.0*Index1-1.0) * moments[4];

	  firstpeak = momleft * momdown / moments[4];
	}
      }
    }

    if (!deltaextrap)
    {
      //Analytical Inversion: (0,0), (1,0), (0,1), (2,0), (1,1), (0,2)
      double fact1 = M00;
      double fact2 = pow(pow(M10, 2.0) * pow(M00, -1.5) * pow(M20, -0.5), Index1);
      double fact3 = pow(pow(M01, 2.0) * pow(M00, -1.5) * pow(M02, -0.5), Index2);
      double fact4 = pow(pow(M20, 0.5) * pow(M00,  0.5) * pow(M10, -1.0), Index1*Index1);
      double fact5 = pow((M11 * M00) / (M10 * M01), Index1*Index2);
      double fact6 = pow(pow(M02, 0.5) * pow(M00,  0.5) * pow(M01, -1.0), Index2*Index2);

      //Analytical Inversion: (0,0), (1,0), (0,1), (h,0), (h,h), (0,h)
      //double fact1 = M00;
      //double fact2 = pow(pow(M20, 4.0) * pow(M00, -3.0) / M10, Index1);
      //double fact3 = pow(pow(M02, 4.0) * pow(M00, -3.0) / M01, Index2);
      //double fact4 = pow(pow(M00, 2.0) * pow(M10, 2.0) * pow(M20, -4.0), Index1*Index1);
      //double fact5 = pow((M00 * M11) / (M20 * M02), 4.0*Index1*Index2);
      //double fact6 = pow(pow(M00, 2.0) * pow(M01, 2.0) * pow(M02, -4.0), Index2*Index2);

      bothpeaks = moments[6] * pow(2.0*dimer_nbrC2, Index1+(2.0/3.0)*Index2) + fact1 * fact2 * fact3 * fact4 * fact5 * fact6;
    }
    else
    {
      if ((Index1+Index2) <= 2.0)
      {
	//Analytical Inversion: (0,0), (1,0), (0,1), (2,0), (1,1), (0,2)
	double fact1 = M00;
	double fact2 = pow(pow(M10, 2.0) * pow(M00, -1.5) * pow(M20, -0.5), Index1);
	double fact3 = pow(pow(M01, 2.0) * pow(M00, -1.5) * pow(M02, -0.5), Index2);
	double fact4 = pow(pow(M20, 0.5) * pow(M00,  0.5) * pow(M10, -1.0), Index1*Index1);
	double fact5 = pow((M11 * M00) / (M10 * M01), Index1*Index2);
	double fact6 = pow(pow(M02, 0.5) * pow(M00,  0.5) * pow(M01, -1.0), Index2*Index2);

	//Analytical Inversion: (0,0), (1,0), (0,1), (h,0), (h,h), (0,h)
	//double fact1 = M00;
	//double fact2 = pow(pow(M20, 4.0) * pow(M00, -3.0) / M10, Index1);
	//double fact3 = pow(pow(M02, 4.0) * pow(M00, -3.0) / M01, Index2);
	//double fact4 = pow(pow(M00, 2.0) * pow(M10, 2.0) * pow(M20, -4.0), Index1*Index1);
	//double fact5 = pow((M00 * M11) / (M20 * M02), 4.0*Index1*Index2);
	//double fact6 = pow(pow(M00, 2.0) * pow(M01, 2.0) * pow(M02, -4.0), Index2*Index2);

	bothpeaks = moments[6] * pow(2.0*dimer_nbrC2, Index1+(2.0/3.0)*Index2) + fact1 * fact2 * fact3 * fact4 * fact5 * fact6;
      }
      else
      {
	if (Index1 <= 1.0) //M2
	  bothpeaks = moments[6] * pow(2.0*dimer_nbrC2, Index1+(2.0/3.0)*Index2) + pow(M01, 2.0-Index1-Index2) * pow(M11, Index1) * pow(M02, Index2-1.0); //M2
	//if (Index1 <= 0.5)
	//  firstpeak = moments[6] * pow(2.0*dimer_nbrC2, Index1+(2.0/3.0)*Index2) + pow(M02, 2.0-2.0*Index1-2.0*Index2) * pow(M11, 2.0*Index1) * pow(M01, 2.0*Index2-1.0);
	else if (Index2 <= 1.0) //M2
	  bothpeaks = moments[6] * pow(2.0*dimer_nbrC2, Index1+(2.0/3.0)*Index2) + pow(M10, 2.0-Index1-Index2) * pow(M20, Index1-1.0) * pow(M11, Index2); //M2
	//else if (Index2 <= 0.5)
	//  firstpeak = moments[6] * pow(2.0*dimer_nbrC2, Index1+(2.0/3.0)*Index2) + pow(M20, 2.0-2.0*Index1-2.0*Index2) * pow(M10, 2.0*Index1-1.0) * pow(M11, 2.0*Index2);
	else
	{
	  double momleft = pow(M01, 1.0-Index2) * M11 * pow(M02, Index2-1.0); //M2
	  //double momleft = pow(M02, 1.0-2.0*Index2) * M11 * pow(M01, 2.0*Index2-1.0);
	  double momdown = pow(M10, 1.0-Index1) * pow(M20, Index1-1.0) * M11; //M2
	  //double momdown = pow(M20, 1.0-2.0*Index1) * pow(M10, 2.0*Index1-1.0) * M11;

	  bothpeaks = moments[6] * pow(2.0*dimer_nbrC2, Index1+(2.0/3.0)*Index2) + momleft * momdown / moments[4];
	}
      }
    }

    //Check for a delta function
    if (M00 < 1.0e-30 || M10 < 1.0e-30 || M01 < 1.0e-30 || M20 < 1.0e-30 || M11 < 1.0e-30 || M02 < 1.0e-30)
      fracmom = firstpeak;
    else
      fracmom = bothpeaks;

    if (fracmom != fracmom)
    {
      cerr << "The fractional moment is NaN...\n";
      cerr << "M00:" << "\t" << moments[0] << endl;
      cerr << "M10:" << "\t" << moments[1] << endl;
      cerr << "M01:" << "\t" << moments[2] << endl;
      cerr << "M20:" << "\t" << moments[3] << endl;
      cerr << "M11:" << "\t" << moments[4] << endl;
      cerr << "M02:" << "\t" << moments[5] << endl;
      cerr << "N0:" << "\t" << moments[6] << endl;
      cerr << "M'00:" << "\t" << M00 << endl;
      cerr << "M'10:" << "\t" << M10 << endl;
      cerr << "M'01:" << "\t" << M01 << endl;
      cerr << "M'20:" << "\t" << M20 << endl;
      cerr << "M'11:" << "\t" << M11 << endl;
      cerr << "M'02:" << "\t" << M02 << endl;
      cerr << "First Peak: " << "\t" << firstpeak << endl;
      cerr << "Both Peaks: " << "\t" << bothpeaks << endl;
      cerr << "Frac Mom(" << Index1 << "," << Index2 << "):" << "\t" << fracmom << endl;
    }

    return fracmom;
  }

  if (fNSootMoments == 10)
  {
    //Analytical Inversion: (0,0), (1,0), (0,1), (2,0), (1,1), (0,2), (3,0), (2,1), (1,2), (0,3)
    double fact01 = moments[0];
    double fact02 = pow(pow(moments[0], -11.0/6.0) * pow(moments[1], 3.0) * pow(moments[3], -3.0/2.0) * pow(moments[6], 1.0/3.0), Index1);
    double fact03 = pow(pow(moments[0], -11.0/6.0) * pow(moments[2], 3.0) * pow(moments[5], -3.0/2.0) * pow(moments[9], 1.0/3.0), Index2);
    double fact04 = pow(moments[0] * pow(moments[1], -5.0/2.0) * pow(moments[3], 2.0) * pow(moments[6], -1.0/2.0), Index1*Index1);
    double fact05 = pow(pow(moments[0], 2.0) * pow(moments[1], -5.0/2.0) * pow(moments[2], -5.0/2.0) * pow(moments[3], 1.0/2.0) * pow(moments[4], 3.0) * pow(moments[5], 1.0/2.0) * pow(moments[7], -1.0/2.0) * pow(moments[8], -1.0/2.0), Index1*Index2);
    double fact06 = pow(moments[0] * pow(moments[2], -5.0/2.0) * pow(moments[5], 2.0) * pow(moments[9], -1.0/2.0), Index2*Index2);
    double fact07 = pow(pow(moments[0], -1.0/6.0) * pow(moments[1], 1.0/2.0) * pow(moments[3], -1.0/2.0) * pow(moments[6], 1.0/6.0), Index1*Index1*Index1);
    double fact08 = pow(pow(moments[0], -1.0/2.0) * moments[1] * pow(moments[2], 1.0/2.0) * pow(moments[3], -1.0/2.0) / moments[4] * pow(moments[7], 1.0/2.0), Index1*Index1*Index2);
    double fact09 = pow(pow(moments[0], -1.0/2.0) * pow(moments[1], 1.0/2.0) * moments[2] / moments[4] * pow(moments[5], -1.0/2.0) * pow(moments[8], 1.0/2.0) , Index1*Index2*Index2);
    double fact10 = pow(pow(moments[0], -1.0/6.0) * pow(moments[2], 1.0/2.0) * pow(moments[5], -1.0/2.0) * pow(moments[9], 1.0/6.0), Index2*Index2*Index2);

    return fact01 * fact02 * fact03 * fact04 * fact05 * fact06 * fact07 * fact08 * fact09 * fact10;
  }

  if (fNSootMoments == 11)
  {
    double M00 = moments[0] -                               moments[10];
    double M10 = moments[1] -     2.0*dimer_nbrC2          *moments[10];
    double M01 = moments[2] - pow(2.0*dimer_nbrC2, 2.0/3.0)*moments[10];
    double M20 = moments[3] - pow(2.0*dimer_nbrC2, 6.0/3.0)*moments[10];
    double M11 = moments[4] - pow(2.0*dimer_nbrC2, 5.0/3.0)*moments[10];
    double M02 = moments[5] - pow(2.0*dimer_nbrC2, 4.0/3.0)*moments[10];
    double M30 = moments[6] - pow(2.0*dimer_nbrC2, 9.0/3.0)*moments[10];
    double M21 = moments[7] - pow(2.0*dimer_nbrC2, 8.0/3.0)*moments[10];
    double M12 = moments[8] - pow(2.0*dimer_nbrC2, 7.0/3.0)*moments[10];
    double M03 = moments[9] - pow(2.0*dimer_nbrC2, 6.0/3.0)*moments[10];

    double fracmom;
    double firstpeak;
    double bothpeaks;

    //Analytical Inversion: (0,0), (1,0), (0,1), (2,0), (1,1), (0,2), (3,0), (2,1), (1,2), (0,3)
    double ffact01 = moments[0];
    double ffact02 = pow(pow(moments[0], -11.0/6.0) * pow(moments[1], 3.0) * pow(moments[3], -3.0/2.0) * pow(moments[6], 1.0/3.0), Index1);
    double ffact03 = pow(pow(moments[0], -11.0/6.0) * pow(moments[2], 3.0) * pow(moments[5], -3.0/2.0) * pow(moments[9], 1.0/3.0), Index2);
    double ffact04 = pow(moments[0] * pow(moments[1], -5.0/2.0) * pow(moments[3], 2.0) * pow(moments[6], -1.0/2.0), Index1*Index1);
    double ffact05 = pow(pow(moments[0], 2.0) * pow(moments[1], -5.0/2.0) * pow(moments[2], -5.0/2.0) * pow(moments[3], 1.0/2.0) * pow(moments[4], 3.0) * pow(moments[5], 1.0/2.0) * pow(moments[7], -1.0/2.0) * pow(moments[8], -1.0/2.0), Index1*Index2);
    double ffact06 = pow(moments[0] * pow(moments[2], -5.0/2.0) * pow(moments[5], 2.0) * pow(moments[9], -1.0/2.0), Index2*Index2);
    double ffact07 = pow(pow(moments[0], -1.0/6.0) * pow(moments[1], 1.0/2.0) * pow(moments[3], -1.0/2.0) * pow(moments[6], 1.0/6.0), Index1*Index1*Index1);
    double ffact08 = pow(pow(moments[0], -1.0/2.0) * moments[1] * pow(moments[2], 1.0/2.0) * pow(moments[3], -1.0/2.0) / moments[4] * pow(moments[7], 1.0/2.0), Index1*Index1*Index2);
    double ffact09 = pow(pow(moments[0], -1.0/2.0) * pow(moments[1], 1.0/2.0) * moments[2] / moments[4] * pow(moments[5], -1.0/2.0) * pow(moments[8], 1.0/2.0) , Index1*Index2*Index2);
    double ffact10 = pow(pow(moments[0], -1.0/6.0) * pow(moments[2], 1.0/2.0) * pow(moments[5], -1.0/2.0) * pow(moments[9], 1.0/6.0), Index2*Index2*Index2);

    firstpeak = ffact01 * ffact02 * ffact03 * ffact04 * ffact05 * ffact06 * ffact07 * ffact08 * ffact09 * ffact10;

    //Analytical Inversion: (0,0), (1,0), (0,1), (2,0), (1,1), (0,2), (3,0), (2,1), (1,2), (0,3)
    double fact01 = M00;
    double fact02 = pow(pow(M00, -11.0/6.0) * pow(M10, 3.0) * pow(M20, -3.0/2.0) * pow(M30, 1.0/3.0), Index1);
    double fact03 = pow(pow(M00, -11.0/6.0) * pow(M01, 3.0) * pow(M02, -3.0/2.0) * pow(M03, 1.0/3.0), Index2);
    double fact04 = pow(M00 * pow(M10, -5.0/2.0) * pow(M20, 2.0) * pow(M30, -1.0/2.0), Index1*Index1);
    double fact05 = pow(pow(M00, 2.0) * pow(M10, -5.0/2.0) * pow(M01, -5.0/2.0) * pow(M20, 1.0/2.0) * pow(M11, 3.0) * pow(M02, 1.0/2.0) * pow(M21, -1.0/2.0) * pow(M12, -1.0/2.0), Index1*Index2);
    double fact06 = pow(M00 * pow(M01, -5.0/2.0) * pow(M02, 2.0) * pow(M03, -1.0/2.0), Index2*Index2);
    double fact07 = pow(pow(M00, -1.0/6.0) * pow(M10, 1.0/2.0) * pow(M20, -1.0/2.0) * pow(M30, 1.0/6.0), Index1*Index1*Index1);
    double fact08 = pow(pow(M00, -1.0/2.0) * M10 * pow(M01, 1.0/2.0) * pow(M20, -1.0/2.0) / M11 * pow(M21, 1.0/2.0), Index1*Index1*Index2);
    double fact09 = pow(pow(M00, -1.0/2.0) * pow(M10, 1.0/2.0) * M01 / M11 * pow(M02, -1.0/2.0) * pow(M12, 1.0/2.0) , Index1*Index2*Index2);
    double fact10 = pow(pow(M00, -1.0/6.0) * pow(M01, 1.0/2.0) * pow(M02, -1.0/2.0) * pow(M03, 1.0/6.0), Index2*Index2*Index2);

    bothpeaks = moments[10] * pow(2.0*dimer_nbrC2, Index1+(2.0/3.0)*Index2) + fact01 * fact02 * fact03 * fact04 * fact05 * fact06 * fact07 * fact08 * fact09 * fact10;

    //Check for a delta function
    if (M00 < 1.0e-30 || M10 < 1.0e-30 || M01 < 1.0e-30 || M20 < 1.0e-30 || M11 < 1.0e-30 || M02 < 1.0e-30 || M30 < 1.0e-30 || M21 < 1.0e-30 || M12 < 1.0e-30 || M03 < 1.0e-30)
      fracmom = firstpeak;
    else
      fracmom = bothpeaks;

    if (fracmom != fracmom)
    {
      cerr << "The fractional moment is NaN...\n";
      cerr << "M00:" << "\t" << moments[0] << endl;
      cerr << "M10:" << "\t" << moments[1] << endl;
      cerr << "M01:" << "\t" << moments[2] << endl;
      cerr << "M20:" << "\t" << moments[3] << endl;
      cerr << "M11:" << "\t" << moments[4] << endl;
      cerr << "M02:" << "\t" << moments[5] << endl;
      cerr << "M30:" << "\t" << moments[6] << endl;
      cerr << "M21:" << "\t" << moments[7] << endl;
      cerr << "M12:" << "\t" << moments[8] << endl;
      cerr << "M03:" << "\t" << moments[9] << endl;
      cerr << "N0:" << "\t" << moments[10] << endl;
      cerr << "M'00:" << "\t" << M00 << endl;
      cerr << "M'10:" << "\t" << M10 << endl;
      cerr << "M'01:" << "\t" << M01 << endl;
      cerr << "M'20:" << "\t" << M20 << endl;
      cerr << "M'11:" << "\t" << M11 << endl;
      cerr << "M'02:" << "\t" << M02 << endl;
      cerr << "M'30:" << "\t" << M30 << endl;
      cerr << "M'21:" << "\t" << M21 << endl;
      cerr << "M'12:" << "\t" << M12 << endl;
      cerr << "M'03:" << "\t" << M03 << endl;
      cerr << "First Peak: " << "\t" << firstpeak << endl;
      cerr << "Both Peaks: " << "\t" << bothpeaks << endl;
      cerr << "Frac Mom(" << Index1 << "," << Index2 << "):" << "\t" << fracmom << endl;
    }

    return fracmom;
  }

  if (fNSootMoments == 15)
  {
    //Analytical Inversion: (0,0), (1,0), (0,1), (2,0), (1,1), (0,2), (3,0), (2,1), (1,2), (0,3), (4,0), (3,1), (2,2), (1,3), (0,4)
    double fact01 = moments[0];
    double fact02 = pow(pow(moments[0],-25.0/12.0)*pow(moments[1],4.0)*pow(moments[3],-3.0)*pow(moments[6],4.0/3.0)*pow(moments[10],-1.0/4.0),Index1);
    double fact03 = pow(pow(moments[0],-25.0/12.0)*pow(moments[2],4.0)*pow(moments[5],-3.0)*pow(moments[9],4.0/3.0)*pow(moments[14],-1.0/4.0),Index2);
    double fact04 = pow(pow(moments[0],35.0/24.0)*pow(moments[1],-13.0/3.0)*pow(moments[3],19.0/4.0)*pow(moments[6],-7.0/3.0)*pow(moments[10],11.0/24.0),Index1*Index1);
    double fact05 = pow(pow(moments[0],35.0/12.0)*pow(moments[1],-13.0/3.0)*pow(moments[2],-13.0/3.0)*pow(moments[3],7.0/4.0)*pow(moments[4],6.0)*pow(moments[5],7.0/4.0)*pow(moments[6],-1.0/3.0)*pow(moments[7],-2.0)*pow(moments[8],-2.0)*pow(moments[9],-1.0/3.0)*pow(moments[11],1.0/3.0)*pow(moments[12],1.0/4.0)*pow(moments[13],1.0/3.0),Index1*Index2);
    double fact06 = pow(pow(moments[0],35.0/24.0)*pow(moments[2],-13.0/3.0)*pow(moments[5],19.0/4.0)*pow(moments[9],-7.0/3.0)*pow(moments[14],11.0/24.0),Index2*Index2);
    double fact07 = pow(pow(moments[0],-5.0/12.0)*pow(moments[1],3.0/2.0)*pow(moments[3],-2.0)*pow(moments[6],7.0/6.0)*pow(moments[10],-1.0/4.0),Index1*Index1*Index1);
    double fact08 = pow(pow(moments[0],-5.0/4.0)*pow(moments[1],3.0)*pow(moments[2],3.0/2.0)*pow(moments[3],-9.0/4.0)*pow(moments[4],-7.0/2.0)*pow(moments[5],-1.0/4.0)*pow(moments[6],1.0/2.0)*pow(moments[7],5.0/2.0)*pow(moments[8],1.0/2.0)*pow(moments[11],-1.0/2.0)*pow(moments[12],-1.0/4.0),Index1*Index1*Index2);
    double fact09 = pow(pow(moments[0],-5.0/4.0)*pow(moments[1],3.0/2.0)*pow(moments[2],3.0)*pow(moments[3],-1.0/4.0)*pow(moments[4],-7.0/2.0)*pow(moments[5],-9.0/4.0)*pow(moments[7],1.0/2.0)*pow(moments[8],5.0/2.0)*pow(moments[9],1.0/2.0)*pow(moments[12],-1.0/4.0)*pow(moments[13],-1.0/2.0),Index1*Index2*Index2);
    double fact10 = pow(pow(moments[0],-5.0/12.0)*pow(moments[2],3.0/2.0)*pow(moments[5],-2.0)*pow(moments[9],7.0/6.0)*pow(moments[14],-1.0/4.0),Index2*Index2*Index2);
    double fact11 = pow(pow(moments[0],1.0/24.0)*pow(moments[1],-1.0/6.0)*pow(moments[3],1.0/4.0)*pow(moments[6],-1.0/6.0)*pow(moments[10],1.0/24.0),Index1*Index1*Index1*Index1);
    double fact12 = pow(pow(moments[0],1.0/6.0)*pow(moments[1],-1.0/2.0)*pow(moments[2],-1.0/6.0)*pow(moments[3],1.0/2.0)*pow(moments[4],1.0/2.0)*pow(moments[6],-1.0/6.0)*pow(moments[7],-1.0/2.0)*pow(moments[11],1.0/6.0),Index1*Index1*Index1*Index2);
    double fact13 = pow(pow(moments[0],1.0/4.0)*pow(moments[1],-1.0/2.0)*pow(moments[2],-1.0/2.0)*pow(moments[3],1.0/4.0)*moments[4]*pow(moments[5],1.0/4.0)*pow(moments[7],-1.0/2.0)*pow(moments[8],-1.0/2.0)*pow(moments[12],1.0/4.0),Index1*Index1*Index2*Index2);
    double fact14 = pow(pow(moments[0],1.0/6.0)*pow(moments[1],-1.0/6.0)*pow(moments[2],-1.0/2.0)*pow(moments[4],1.0/2.0)*pow(moments[5],1.0/2.0)*pow(moments[8],-1.0/2.0)*pow(moments[9],-1.0/6.0)*pow(moments[13],1.0/6.0),Index1*Index2*Index2*Index2);
    double fact15 = pow(pow(moments[0],1.0/24.0)*pow(moments[2],-1.0/6.0)*pow(moments[5],1.0/4.0)*pow(moments[9],-1.0/6.0)*pow(moments[14],1.0/24.0),Index2*Index2*Index2*Index2);

    return fact01 * fact02 * fact03 * fact04 * fact05 * fact06 * fact07 * fact08 * fact09 * fact10 * fact11 * fact12 * fact13 * fact14 * fact15;
  }

  if (fNSootMoments == 16)
  {
    double M00 = moments[ 0] -                                moments[15];
    double M10 = moments[ 1] -     2.0*dimer_nbrC2           *moments[15];
    double M01 = moments[ 2] - pow(2.0*dimer_nbrC2,  2.0/3.0)*moments[15];
    double M20 = moments[ 3] - pow(2.0*dimer_nbrC2,  6.0/3.0)*moments[15];
    double M11 = moments[ 4] - pow(2.0*dimer_nbrC2,  5.0/3.0)*moments[15];
    double M02 = moments[ 5] - pow(2.0*dimer_nbrC2,  4.0/3.0)*moments[15];
    double M30 = moments[ 6] - pow(2.0*dimer_nbrC2,  9.0/3.0)*moments[15];
    double M21 = moments[ 7] - pow(2.0*dimer_nbrC2,  8.0/3.0)*moments[15];
    double M12 = moments[ 8] - pow(2.0*dimer_nbrC2,  7.0/3.0)*moments[15];
    double M03 = moments[ 9] - pow(2.0*dimer_nbrC2,  6.0/3.0)*moments[15];
    double M40 = moments[10] - pow(2.0*dimer_nbrC2, 12.0/3.0)*moments[15];
    double M31 = moments[11] - pow(2.0*dimer_nbrC2, 11.0/3.0)*moments[15];
    double M22 = moments[12] - pow(2.0*dimer_nbrC2, 10.0/3.0)*moments[15];
    double M13 = moments[13] - pow(2.0*dimer_nbrC2,  9.0/3.0)*moments[15];
    double M04 = moments[14] - pow(2.0*dimer_nbrC2,  8.0/3.0)*moments[15];

    double fracmom;
    double firstpeak;
    double bothpeaks;

    //Analytical Inversion: (0,0), (1,0), (0,1), (2,0), (1,1), (0,2), (3,0), (2,1), (1,2), (0,3), (4,0), (3,1), (2,2), (1,3), (0,4)
    double ffact01 = moments[0];
    double ffact02 = pow(pow(moments[0],-25.0/12.0)*pow(moments[1],4.0)*pow(moments[3],-3.0)*pow(moments[6],4.0/3.0)*pow(moments[10],-1.0/4.0),Index1);
    double ffact03 = pow(pow(moments[0],-25.0/12.0)*pow(moments[2],4.0)*pow(moments[5],-3.0)*pow(moments[9],4.0/3.0)*pow(moments[14],-1.0/4.0),Index2);
    double ffact04 = pow(pow(moments[0],35.0/24.0)*pow(moments[1],-13.0/3.0)*pow(moments[3],19.0/4.0)*pow(moments[6],-7.0/3.0)*pow(moments[10],11.0/24.0),Index1*Index1);
    double ffact05 = pow(pow(moments[0],35.0/12.0)*pow(moments[1],-13.0/3.0)*pow(moments[2],-13.0/3.0)*pow(moments[3],7.0/4.0)*pow(moments[4],6.0)*pow(moments[5],7.0/4.0)*pow(moments[6],-1.0/3.0)*pow(moments[7],-2.0)*pow(moments[8],-2.0)*pow(moments[9],-1.0/3.0)*pow(moments[11],1.0/3.0)*pow(moments[12],1.0/4.0)*pow(moments[13],1.0/3.0),Index1*Index2);
    double ffact06 = pow(pow(moments[0],35.0/24.0)*pow(moments[2],-13.0/3.0)*pow(moments[5],19.0/4.0)*pow(moments[9],-7.0/3.0)*pow(moments[14],11.0/24.0),Index2*Index2);
    double ffact07 = pow(pow(moments[0],-5.0/12.0)*pow(moments[1],3.0/2.0)*pow(moments[3],-2.0)*pow(moments[6],7.0/6.0)*pow(moments[10],-1.0/4.0),Index1*Index1*Index1);
    double ffact08 = pow(pow(moments[0],-5.0/4.0)*pow(moments[1],3.0)*pow(moments[2],3.0/2.0)*pow(moments[3],-9.0/4.0)*pow(moments[4],-7.0/2.0)*pow(moments[5],-1.0/4.0)*pow(moments[6],1.0/2.0)*pow(moments[7],5.0/2.0)*pow(moments[8],1.0/2.0)*pow(moments[11],-1.0/2.0)*pow(moments[12],-1.0/4.0),Index1*Index1*Index2);
    double ffact09 = pow(pow(moments[0],-5.0/4.0)*pow(moments[1],3.0/2.0)*pow(moments[2],3.0)*pow(moments[3],-1.0/4.0)*pow(moments[4],-7.0/2.0)*pow(moments[5],-9.0/4.0)*pow(moments[7],1.0/2.0)*pow(moments[8],5.0/2.0)*pow(moments[9],1.0/2.0)*pow(moments[12],-1.0/4.0)*pow(moments[13],-1.0/2.0),Index1*Index2*Index2);
    double ffact10 = pow(pow(moments[0],-5.0/12.0)*pow(moments[2],3.0/2.0)*pow(moments[5],-2.0)*pow(moments[9],7.0/6.0)*pow(moments[14],-1.0/4.0),Index2*Index2*Index2);
    double ffact11 = pow(pow(moments[0],1.0/24.0)*pow(moments[1],-1.0/6.0)*pow(moments[3],1.0/4.0)*pow(moments[6],-1.0/6.0)*pow(moments[10],1.0/24.0),Index1*Index1*Index1*Index1);
    double ffact12 = pow(pow(moments[0],1.0/6.0)*pow(moments[1],-1.0/2.0)*pow(moments[2],-1.0/6.0)*pow(moments[3],1.0/2.0)*pow(moments[4],1.0/2.0)*pow(moments[6],-1.0/6.0)*pow(moments[7],-1.0/2.0)*pow(moments[11],1.0/6.0),Index1*Index1*Index1*Index2);
    double ffact13 = pow(pow(moments[0],1.0/4.0)*pow(moments[1],-1.0/2.0)*pow(moments[2],-1.0/2.0)*pow(moments[3],1.0/4.0)*moments[4]*pow(moments[5],1.0/4.0)*pow(moments[7],-1.0/2.0)*pow(moments[8],-1.0/2.0)*pow(moments[12],1.0/4.0),Index1*Index1*Index2*Index2);
    double ffact14 = pow(pow(moments[0],1.0/6.0)*pow(moments[1],-1.0/6.0)*pow(moments[2],-1.0/2.0)*pow(moments[4],1.0/2.0)*pow(moments[5],1.0/2.0)*pow(moments[8],-1.0/2.0)*pow(moments[9],-1.0/6.0)*pow(moments[13],1.0/6.0),Index1*Index2*Index2*Index2);
    double ffact15 = pow(pow(moments[0],1.0/24.0)*pow(moments[2],-1.0/6.0)*pow(moments[5],1.0/4.0)*pow(moments[9],-1.0/6.0)*pow(moments[14],1.0/24.0),Index2*Index2*Index2*Index2);

    firstpeak = ffact01 * ffact02 * ffact03 * ffact04 * ffact05 * ffact06 * ffact07 * ffact08 * ffact09 * ffact10 * ffact11 * ffact12 * ffact13 * ffact14 *ffact15;

    //Analytical Inversion: (0,0), (1,0), (0,1), (2,0), (1,1), (0,2), (3,0), (2,1), (1,2), (0,3), (4,0), (3,1), (2,2), (1,3), (0,4)
    double fact01 = M00;
    double fact02 = pow(pow(M00,-25.0/12.0)*pow(M10,4.0)*pow(M20,-3.0)*pow(M30,4.0/3.0)*pow(M40,-1.0/4.0),Index1);
    double fact03 = pow(pow(M00,-25.0/12.0)*pow(M01,4.0)*pow(M02,-3.0)*pow(M03,4.0/3.0)*pow(M04,-1.0/4.0),Index2);
    double fact04 = pow(pow(M00,35.0/24.0)*pow(M10,-13.0/3.0)*pow(M20,19.0/4.0)*pow(M30,-7.0/3.0)*pow(M40,11.0/24.0),Index1*Index1);
    double fact05 = pow(pow(M00,35.0/12.0)*pow(moments[1],-13.0/3.0)*pow(moments[2],-13.0/3.0)*pow(moments[3],7.0/4.0)*pow(moments[4],6.0)*pow(moments[5],7.0/4.0)*pow(moments[6],-1.0/3.0)*pow(moments[7],-2.0)*pow(moments[8],-2.0)*pow(moments[9],-1.0/3.0)*pow(moments[11],1.0/3.0)*pow(moments[12],1.0/4.0)*pow(moments[13],1.0/3.0),Index1*Index2);
    double fact06 = pow(pow(M00,35.0/24.0)*pow(M01,-13.0/3.0)*pow(M02,19.0/4.0)*pow(M03,-7.0/3.0)*pow(M04,11.0/24.0),Index2*Index2);
    double fact07 = pow(pow(M00,-5.0/12.0)*pow(M10,3.0/2.0)*pow(M20,-2.0)*pow(M30,7.0/6.0)*pow(M40,-1.0/4.0),Index1*Index1*Index1);
    double fact08 = pow(pow(M00,-5.0/4.0)*pow(moments[1],3.0)*pow(moments[2],3.0/2.0)*pow(moments[3],-9.0/4.0)*pow(moments[4],-7.0/2.0)*pow(moments[5],-1.0/4.0)*pow(moments[6],1.0/2.0)*pow(moments[7],5.0/2.0)*pow(moments[8],1.0/2.0)*pow(moments[11],-1.0/2.0)*pow(moments[12],-1.0/4.0),Index1*Index1*Index2);
    double fact09 = pow(pow(M00,-5.0/4.0)*pow(moments[1],3.0/2.0)*pow(moments[2],3.0)*pow(moments[3],-1.0/4.0)*pow(moments[4],-7.0/2.0)*pow(moments[5],-9.0/4.0)*pow(moments[7],1.0/2.0)*pow(moments[8],5.0/2.0)*pow(moments[9],1.0/2.0)*pow(moments[12],-1.0/4.0)*pow(moments[13],-1.0/2.0),Index1*Index2*Index2);
    double fact10 = pow(pow(M00,-5.0/12.0)*pow(M01,3.0/2.0)*pow(M02,-2.0)*pow(M03,7.0/6.0)*pow(M04,-1.0/4.0),Index2*Index2*Index2);
    double fact11 = pow(pow(M00,1.0/24.0)*pow(M10,-1.0/6.0)*pow(M20,1.0/4.0)*pow(M30,-1.0/6.0)*pow(M40,1.0/24.0),Index1*Index1*Index1*Index1);
    double fact12 = pow(pow(M00,1.0/6.0)*pow(moments[1],-1.0/2.0)*pow(moments[2],-1.0/6.0)*pow(moments[3],1.0/2.0)*pow(moments[4],1.0/2.0)*pow(moments[6],-1.0/6.0)*pow(moments[7],-1.0/2.0)*pow(moments[11],1.0/6.0),Index1*Index1*Index1*Index2);
    double fact13 = pow(pow(M00,1.0/4.0)*pow(moments[1],-1.0/2.0)*pow(moments[2],-1.0/2.0)*pow(moments[3],1.0/4.0)*moments[4]*pow(moments[5],1.0/4.0)*pow(moments[7],-1.0/2.0)*pow(moments[8],-1.0/2.0)*pow(moments[12],1.0/4.0),Index1*Index1*Index2*Index2);
    double fact14 = pow(pow(M00,1.0/6.0)*pow(moments[1],-1.0/6.0)*pow(moments[2],-1.0/2.0)*pow(moments[4],1.0/2.0)*pow(moments[5],1.0/2.0)*pow(moments[8],-1.0/2.0)*pow(moments[9],-1.0/6.0)*pow(moments[13],1.0/6.0),Index1*Index2*Index2*Index2);
    double fact15 = pow(pow(M00,1.0/24.0)*pow(M01,-1.0/6.0)*pow(M02,1.0/4.0)*pow(M03,-1.0/6.0)*pow(M04,1.0/24.0),Index2*Index2*Index2*Index2);

    bothpeaks = moments[10] * pow(2.0*dimer_nbrC2, Index1+(2.0/3.0)*Index2) + fact01 * fact02 * fact03 * fact04 * fact05 * fact06 * fact07 * fact08 * fact09 * fact10 * fact11 * fact12 * fact13 * fact14 * fact15;

    //Check for a delta function
    if (M00 < 1.0e-30 || M10 < 1.0e-30 || M01 < 1.0e-30 || M20 < 1.0e-30 || M11 < 1.0e-30 || M02 < 1.0e-30 || M30 < 1.0e-30 || M21 < 1.0e-30 || M12 < 1.0e-30 || M03 < 1.0e-30 || M40 < 1.0e-30 || M31 < 1.0e-30 || M22 < 1.0e-30 || M13 < 1.0e-30 || M04 < 1.0e-30)
      fracmom = firstpeak;
    else
      fracmom = bothpeaks;

    if (fracmom != fracmom)
    {
      cerr << "The fractional moment is NaN...\n";
      cerr << "M00:" << "\t" << moments[0] << endl;
      cerr << "M10:" << "\t" << moments[1] << endl;
      cerr << "M01:" << "\t" << moments[2] << endl;
      cerr << "M20:" << "\t" << moments[3] << endl;
      cerr << "M11:" << "\t" << moments[4] << endl;
      cerr << "M02:" << "\t" << moments[5] << endl;
      cerr << "M30:" << "\t" << moments[6] << endl;
      cerr << "M21:" << "\t" << moments[7] << endl;
      cerr << "M12:" << "\t" << moments[8] << endl;
      cerr << "M03:" << "\t" << moments[9] << endl;
      cerr << "M40:" << "\t" << moments[10] << endl;
      cerr << "M31:" << "\t" << moments[11] << endl;
      cerr << "M22:" << "\t" << moments[12] << endl;
      cerr << "M13:" << "\t" << moments[13] << endl;
      cerr << "M04:" << "\t" << moments[14] << endl;
      cerr << "N0:" << "\t" << moments[15] << endl;
      cerr << "M'00:" << "\t" << M00 << endl;
      cerr << "M'10:" << "\t" << M10 << endl;
      cerr << "M'01:" << "\t" << M01 << endl;
      cerr << "M'20:" << "\t" << M20 << endl;
      cerr << "M'11:" << "\t" << M11 << endl;
      cerr << "M'02:" << "\t" << M02 << endl;
      cerr << "M'30:" << "\t" << M30 << endl;
      cerr << "M'21:" << "\t" << M21 << endl;
      cerr << "M'12:" << "\t" << M12 << endl;
      cerr << "M'03:" << "\t" << M03 << endl;
      cerr << "M'40:" << "\t" << M40 << endl;
      cerr << "M'31:" << "\t" << M31 << endl;
      cerr << "M'22:" << "\t" << M22 << endl;
      cerr << "M'13:" << "\t" << M13 << endl;
      cerr << "M'04:" << "\t" << M04 << endl;
      cerr << "First Peak: " << "\t" << firstpeak << endl;
      cerr << "Both Peaks: " << "\t" << bothpeaks << endl;
      cerr << "Frac Mom(" << Index1 << "," << Index2 << "):" << "\t" << fracmom << endl;
    }

    return fracmom;
  }

  return 0.0;
}

Double TSoot::FracMomLarge(double Index1, double Index2, double * moments)
{
  double temp = FracMom(Index1, Index2, moments) - moments[fNSootMoments-1] * pow(2.0*dimer_nbrC2, Index1+(2.0/3.0)*Index2); //

  if (temp <= 0)
    return pow(2.0*dimer_nbrC2, Index1+(2.0/3.0)*Index2) * MAX(moments[0]-moments[fNSootMoments-1],1.0e-60); //

  return temp;
}

double TSoot::FracMom(double Index1, double Index2, double Index3, double * moments)
{
  double M000 = moments[0] - moments[4];
  double M100 = moments[1] - 2.0*dimer_nbrC2*moments[4];
  double M010 = moments[2] - pow(2.0*dimer_nbrC2, 2.0/3.0)*moments[4];
  double M001 = moments[3] - 2.0*dimer_nbrH*moments[4]; 

  double fracmom;
  double firstpeak;
  double bothpeaks;

  firstpeak = pow(moments[0], 1.0-Index1-Index2-Index3) * pow(moments[1], Index1) * pow(moments[2], Index2) * pow(moments[3], Index3);
  bothpeaks = moments[4] * pow(2.0*dimer_nbrC2, Index1+(2.0/3.0)*Index2) * pow(2.0*dimer_nbrH, Index3) +
    pow(M000, 1.0-Index1-Index2-Index3) * pow(M100, Index1) * pow(M010, Index2) * pow(M001, Index3);

  //Check for a delta function
  if (M000 < 1.0e-20 || M100 < 1.0e-20 || M010 < 1.0e-20 || M001 < 1.0e-20)
    fracmom = firstpeak;
  else
    fracmom = bothpeaks;

  if (fracmom != fracmom)
  {
    cerr << "The fractional moment is NaN...\n";
    cerr << "M000:" << "\t" << moments[0] << endl;
    cerr << "M100:" << "\t" << moments[1] << endl;
    cerr << "M010:" << "\t" << moments[2] << endl;
    cerr << "M001:" << "\t" << moments[3] << endl;
    cerr << "N0:" << "\t" << moments[4] << endl;
    cerr << "M'000:" << "\t" << M000 << endl;
    cerr << "M'100:" << "\t" << M100 << endl;
    cerr << "M'010:" << "\t" << M010 << endl;
    cerr << "M'001:" << "\t" << M001 << endl;
    cerr << "First Peak: " << "\t" << firstpeak << endl;
    cerr << "Both Peaks: " << "\t" << bothpeaks << endl;
    cerr << "Frac Mom(" << Index1 << "," << Index2 << "," << Index3 << "):" << "\t" << fracmom << endl;
  }

  return fracmom;
}

double TSoot::FracMomLarge(double Index1, double Index2, double Index3, double * moments)
{
  double temp = FracMom(Index1, Index2, Index3, moments) - moments[fNSootMoments-1] * pow(2.0*dimer_nbrC2, Index1+2.0/3.0*Index2) * pow(2.0*dimer_nbrH, Index3);

  if (temp < 0)
    return 1.0e-60;

  return temp;
}

// LINDSTEDT MODEL //
double TSoot::LindstedtNucleation(double temp, double * Y, double density, double * molarMass)
{
  double C_C2H2 = density * Y[f_C2H2] / molarMass[f_C2H2];
  double k = 1.0e4 * exp(-21100.0/temp);

  return k * C_C2H2;
}

double TSoot::LindstedtGrowth(double * moments, double temp, double * Y, double density, double * molarMass)
{
  double C_C2H2 = density * Y[f_C2H2] / molarMass[f_C2H2];
  double k = 6.0e3 * exp(-12100/temp);

  return k * C_C2H2 * sqrt(PI*pow(6.0*fMolarMassSoot/(PI*fSootDensity),2.0/3.0)) * pow(moments[1], 1.0/3.0) * pow(moments[0]*AVOGADRO, 1.0/6.0);
}

double TSoot::LindstedtOxidation(double * moments, double temp, double * Y, double density, double * molarMass)
{
  double C_O2 = density * Y[f_O2] / molarMass[f_O2];
  double k = 1.0e4 * sqrt(temp) * exp(-19680/temp);

  double surf = PI*pow(6.0*moments[1]*fMolarMassSoot/(PI*fSootDensity*moments[0]*AVOGADRO),2.0/3.0)*moments[0]*AVOGADRO;

  return -0.5 * k * surf * C_O2;
}

double TSoot::LindstedtCoagulation(double * moments, double temp)
{
  double Ca = 9.0;

  double dp = pow(6.0*moments[1]*fMolarMassSoot/(PI*fSootDensity*moments[0]*AVOGADRO),1.0/3.0);

  return -2.0*Ca*sqrt(dp)*sqrt(6.0*RGAS*temp/(AVOGADRO*fSootDensity))*(moments[0]*AVOGADRO)*(moments[0]*AVOGADRO)/AVOGADRO;
}
  

#ifdef VSH
void TSoot::CheckSolution(Double * moments, double density)
{
  double x, y, z;

#ifdef MOMIC
  for (int i = 0; i < fNSootMoments; ++i)
#endif
#ifdef HMOM
  for (int i = 0; i < fNSootMoments-1; ++i)
#endif
  {
    x = double(Geti(i))/6.0;
    y = double(Getj(i))/6.0;
    z = double(Getk(i))/6.0;

    switch (i)
    {
      case 0:
	moments[i] = MAX(moments[0]*density, 1.0e-20) / density;
	break;
      case 1:
	moments[i] = MAX(moments[0] * 2.0*nucl_nbrC2*density, moments[i]*density) / density;
	break;
      case 2:
	moments[i] = MAX(moments[0] * pow(2.0*nucl_nbrC2, 2.0/3.0)*density, moments[i]*density) / density;
	break;
      case 3:
	moments[i] = MAX(moments[0] * 2.0*nucl_nbrH*density, moments[i]*density) / density;
      default:
	moments[i] = MAX(pow(moments[0]*density, 1.0-x-y) * pow(moments[1]*density, x) * pow(moments[2]*density, y) * pow(moments[3]*density, z), moments[i]*density) / density;
	break;
    }
  }

#ifdef HMOM
	moments[fNSootMoments-1] = MAX(moments[fNSootMoments-1]*density, 1.0e-20) / density;
#endif
}
#endif
#ifdef VS
void TSoot::CheckSolution(Double * mom_rho, double density)
{
  double x, y;

  if (density == 0.0) density = 1.0;

#ifdef MOMIC
  for (int i = 0; i < fNSootMoments; ++i)
#endif
#ifdef HMOM
  for (int i = 0; i < fNSootMoments-1; ++i)
#endif
  {
    x = double(Geti(i))/6.0;
    y = double(Getj(i))/6.0;

    if (i==0)
    {
#ifndef SOLVELOG
      mom_rho[i] = MAX(mom_rho[0]*density, 1.0e-20) / density;
#else
      mom_rho[i] = log(MAX(exp(mom_rho[0])*density, 1.0e-20) / density);
#endif
    }
    else if (i==1)
    {
#ifndef SOLVELOG
      mom_rho[i] = MAX(mom_rho[1]*density, mom_rho[0]*density*2.0*dimer_nbrC2) / density;
#else
      mom_rho[i] = log(MAX(exp(mom_rho[1])*density, exp(mom_rho[0])*density*2.0*dimer_nbrC2) / density);
#endif
    }
    else if (i==2)
    {
#ifndef SOLVELOG
      mom_rho[i] = MAX(mom_rho[2]*density, mom_rho[0]*density*pow(2.0*dimer_nbrC2, 2.0/3.0)) / density;
#else
      mom_rho[i] = log(MAX(exp(mom_rho[2])*density, exp(mom_rho[0])*density*pow(2.0*dimer_nbrC2, 2.0/3.0)) / density);
#endif
    }
    else
    {
#ifndef SOLVELOG
      mom_rho[i] = MAX(pow(mom_rho[0]*density, 1.0-x-y) * pow(mom_rho[1]*density, x) * pow(mom_rho[2]*density, y), mom_rho[i]*density) / density;
#else
      mom_rho[i] = log(MAX(pow(exp(mom_rho[0])*density, 1.0-x-y) * pow(exp(mom_rho[1])*density, x) * pow(exp(mom_rho[2])*density, y), exp(mom_rho[i])*density) / density);
#endif
    }
  }

#ifdef HMOM
#ifndef SOLVELOG
  mom_rho[fNSootMoments-1] = MAX(mom_rho[fNSootMoments-1]*density, 1.0e-20) / density;
#else
  mom_rho[fNSootMoments-1] = log(MAX(exp(mom_rho[fNSootMoments-1])*density, 1.0e-20) / density);
#endif
#endif
}
#endif
#ifdef V
void TSoot::CheckSolution(Double * moments, double density)
{
  double x;

  for (int i = 0; i < fNSootMoments; ++i)
  {
    x = double(Geti(i))/6.0;
    int index = 0;
    for (int j=1; j<fNSootMoments; j++)
      if ((Geti(i)>Geti(j)) && Geti(j)>Geti(index)) index = j;
    if (Geti(index)>Geti(i)) index = -1;
    
    switch (i)
    {
      case 0:
	moments[i] = MAX(moments[0], 1.0e-20) / density;
	break;
      default:
	moments[i] = MAX(moments[0] * pow(2.0*nucl_nbrC2, x), moments[i]) / density;
	if (index!=-1)
	  moments[i] = MAX(moments[i],moments[index]) / density;
	break;
    }
  }
}
#endif

void TSoot::ComputeCSootStar(Double * k, Double * Y, Double density, Double * molarMass, Double mixMolarMass)
{
  double C_H = density * MAX(Y[f_H],1.0e-60) / molarMass[f_H];
  double C_OH = density * MAX(Y[f_OH],1.0e-60) / molarMass[f_OH];
  double C_H2 = density * MAX(Y[f_H2],1.0e-60) / molarMass[f_H2];
  double C_H2O = density * MAX(Y[f_H2O],1.0e-60) / molarMass[f_H2O];
  double C_O2 = density * MAX(Y[f_O2],1.0e-60) / molarMass[f_O2];
  double C_C2H2 = density * MAX(Y[f_C2H2],1.0e-60) / molarMass[f_C2H2];

  fCSootStar = (k[ks1f] * C_OH + k[ks2f] * C_H + k[ks3f]) / (k[ks1b] * C_H2O + k[ks2b] * C_H2 + k[ks3b] * C_H + k[ks4] * C_C2H2);
  fCSootStar = fCSootStar / (1.0 + fCSootStar);	
}

void TSoot::ComputeSootReactionRates(Double * w, Double * Y, Double * moments, Double density, Double * molarMass)
{
  //Use Chi
  double SphereConst = pow(36.0 * PI, 1.0/3.0);
  double MassConst = pow(fMolarMassSoot / (AVOGADRO * fSootDensity), 2.0 / 3.0);
  double SGConst = fChi * SphereConst * MassConst;
  double partsurffact = PI * pow((6.0*fMolarMassSoot)/(PI*AVOGADRO*fSootDensity), 2.0/3.0);

#ifdef V
  double C_Soot = SGConst * FracMom(2.0/3.0, moments);
#endif
#ifdef VS
  double C_Soot = SGConst * FracMom(0.0, 1.0, moments);
#endif
#ifdef VSH
  double C_Soot = FracMom(0.0, 0.0, 1.0, moments);
#endif
  double C_SootStar = fCSootStar * C_Soot;
  C_Soot = (1.0-fCSootStar)*C_Soot;
  double C_C2H2 = density * Y[f_C2H2] / molarMass[f_C2H2];
  double C_H = density * Y[f_H] / molarMass[f_H];
  double C_O2 = density * Y[f_O2] / molarMass[f_O2];
  double * k = fSootRateCoeffs->vec;
  double C_H2 = density * Y[f_H2] / molarMass[f_H2];
  double C_OH = density * Y[f_OH] / molarMass[f_OH];

  double C_H2O = density * Y[f_H2O] / molarMass[f_H2O];
  double CSootC2H2 = 0.0;

  w[ks1f] = k[ks1f] * C_Soot * C_OH;
  w[ks1b] = k[ks1b] * C_SootStar * C_H2O;
  w[ks2f] = k[ks2f] * C_Soot * C_H;
  w[ks2b] = k[ks2b] * C_SootStar * C_H2;
  w[ks3f] = k[ks3f] * C_Soot;
  w[ks3b] = k[ks3b] * C_SootStar * C_H;
  w[ks4] = k[ks4] * C_SootStar * C_C2H2;

  w[ks5] = 0.0;
  w[ks6] = 0.0;  
  if (fSurfaceOxidation)
  {
    // Surface Oxidation by O2 - Kazakov, Wang & Frenklach 1995
    // Total number of sites
    // -> Given by number of Hydrogenated carbon surface sites
    w[ks5]  = k[ks5] * C_O2 * C_SootStar;

    //cerr << C_O2 << "\t" << fCSootStar << "\t" << C_SootStar << endl;
    
    // Surface Oxidation by OH - Neoh & Sarofim 1981
    // Collision efficiency : 0.13
    // -> Total surface
#ifdef V
    double Surf = partsurffact * FracMom(2.0/3.0, moments);
#endif
#ifdef VS
    double Surf = partsurffact * FracMom(0.0, 1.0, moments);
#endif
#ifdef VSH
    double Surf = partsurffact * FracMom(0.0, 1.0, 0.0, moments);
#endif
    static const double GammaOH = 0.13;
    w[ks6] = k[ks6] * C_OH * Surf;
  }
}

double TSoot::GetSootRadiationCoeff(double temp)
{
  // alphas taken from:
  //   Hubbard, G. L. , Tien C. L: 
  //   Infrared Mean Absorption Coefficient
  //   of Luminous Flames and Smoke
  //   Journal of Heat Transfer, vol 100, p. 235ff, 1978
  double alphas = -3.75e5 + 1735.0 * temp; //[m^-1] 
	
  //MM-Add fake absorption term to keep temperature above 298K
  //return MAX(MAGIC * 4.0 * alphas * STEFBOLTZ * pow(temp, 4.0), 0.0);
  return MAX(4.0 * alphas * STEFBOLTZ * (pow(temp, 4.0)-pow(298.0,4.0)), 0.0);
  //return MAX(4.0 * alphas * STEFBOLTZ * (pow(temp, 4.0)-pow(812.0,4.0)), 0.0);
  //return MAX(4.0 * alphas * STEFBOLTZ * (pow(temp, 4.0)-pow(897.0,4.0)), 0.0);
}

//Mueller (7/23/07)
#ifdef VS
Double TSoot::GetSootRadiation(Double temp, Double * moments)
{
  const double	vol1 = 2.46e-10 * 2.46e-10 * 3.51e-10; //[m^3]
  //double		fvo = vol1 * AVOGADRO * moments[1];
  double fvo = moments[1] * fMolarMassSoot / fSootDensity;
		
  // alphas taken from:
  //		Hubbard, G. L. , Tien C. L: 
  //		Infrared Mean Absorption Coefficient
  //		of Luminous Flames and Smoke
  //		Journal of Heat Transfer, vol 100, p. 235ff, 1978
  double		alphas = -3.75e5 + 1735.0 * temp; //[m^-1] 

  /*
  if (temp > 1400.0){
  cout<<"fv = "<<fvo<<" T = "<<temp<<" Qrad = "<<MAX(MAGIC * 4.0 * alphas * fvo * STEFBOLTZ * (pow(temp, 4.0)-pow(298.0,4.0)), 0.0)<<"\n";
  }
  */
	
  //MM-Add fake absorption term to keep temperature above 298K
  //return MAX(MAGIC * 4.0 * alphas * fvo * STEFBOLTZ * pow(temp, 4.0), 0.0);
  return MAX(MAGIC * 4.0 * alphas * fvo * STEFBOLTZ * (pow(temp, 4.0)-pow(298.0,4.0)), 0.0);
  //return MAX(MAGIC * 4.0 * alphas * fvo * STEFBOLTZ * (pow(temp, 4.0)-pow(812.0,4.0)), 0.0);
  //return MAX(MAGIC * 4.0 * alphas * fvo * STEFBOLTZ * (pow(temp, 4.0)-pow(897.0,4.0)), 0.0);
}
#else
Double TSoot::GetSootRadiation(Double temp, Double * moments)
{
  const double	vol1 = 2.46e-10 * 2.46e-10 * 3.51e-10; //[m^3]
  double		fvo = vol1 * AVOGADRO * moments[1];
		
  // alphas taken from:
  //		Hubbard, G. L. , Tien C. L: 
  //		Infrared Mean Absorption Coefficient
  //		of Luminous Flames and Smoke
  //		Journal of Heat Transfer, vol 100, p. 235ff, 1978
  double		alphas = -3.75e5 + 1735.0 * temp; //[m^-1] 
	
  return MAX(MAGIC * 4.0 * alphas * fvo * STEFBOLTZ * pow(temp, 4.0), 0.0);
}
#endif

//Mueller (7/23/07)
#ifdef VS
Double TSoot::GetSootRadRossCoeff(Double temp, Double * moments, Double cutoff)
{
  const double	vol1 = 2.46e-10 * 2.46e-10 * 3.51e-10; //[m^3]
  double		fvo = vol1 * AVOGADRO * moments[1];

  // alphas taken from:
  //		Hubbard, G. L. , Tien C. L: 
  //		Infrared Mean Absorption Coefficient
  //		of Luminous Flames and Smoke
  //		Journal of Heat Transfer, vol 100, p. 235ff, 1978

  // Rosseland model is del q_R = del(-16/3 sigma T^3/alpha del(T))
  // equal to del(coeff del(T)) = coeff del^2 T + del(coeff) del(T)
  // The coefficient returned here is     coeff = -16/3 sigma T^3/alpha

  double		alphas = -3.75e5 + 1735.0 * MAX(273.0, temp); //[m^-1] 
	
  cutoff = MAX(1.0e-30, cutoff);

  return -MAX(MAGIC * 16.0/3.0 * STEFBOLTZ * pow(temp, 3.0) / (alphas * MAX(cutoff, fvo)), 0.0);
}
#else
Double TSoot::GetSootRadRossCoeff(Double temp, Double * moments, Double cutoff)
{
  const double	vol1 = 2.46e-10 * 2.46e-10 * 3.51e-10; //[m^3]
  double		fvo = vol1 * AVOGADRO * moments[1];

  // alphas taken from:
  //		Hubbard, G. L. , Tien C. L: 
  //		Infrared Mean Absorption Coefficient
  //		of Luminous Flames and Smoke
  //		Journal of Heat Transfer, vol 100, p. 235ff, 1978
  
  // Rosseland model is del q_R = del(-16/3 sigma T^3/alpha del(T))
  // equal to del(coeff del(T)) = coeff del^2 T + del(coeff) del(T)
  // The coefficient returned here is     coeff = -16/3 sigma T^3/alpha

  double		alphas = -3.75e5 + 1735.0 * MAX(273.0, temp); //[m^-1] 
	
  cutoff = MAX(1.0e-30, cutoff);

  return -MAX(MAGIC * 16.0/3.0 * STEFBOLTZ * pow(temp, 3.0) / (alphas * MAX(cutoff, fvo)), 0.0);
}
#endif

void TSoot::UpdateSoot(TReactionPtr reaction, TSpeciesPtr species, Double * moments, Double temp, Double * Y, Double density, Double mixMolarMass)
{
  double *	kSoot = fSootRateCoeffs->vec;
  double *	wSoot = fSootReactionRate->vec;
  double *	molarMass = species->GetMolarMass()->vec;

  ComputeDimerParticules(species, temp, Y, density, 1.0, mixMolarMass, moments);
  ComputeSootRateCoeffs(kSoot, temp, reaction);
  ComputeCSootStar(kSoot, Y, density, molarMass, mixMolarMass);
  ComputeSootReactionRates(wSoot, Y, moments, density, molarMass);
}

void TSoot::MOverRhoToM(Double * momSource, Double * momDest, int nMom, Double * Y, Double temp, Double pressure, Double * molarmass, int nSpeciesIn, TPropertiesPtr props)
{
  double	mixMolarMass;
  double	rho;
	
  props->ComputeMixtureMolarMass(mixMolarMass, Y, molarmass, nSpeciesIn);
  rho = pressure * mixMolarMass / (RGAS * temp);

  CheckSolution(momSource, rho);
	
  for (int i = 0; i < nMom; ++i)
  {
#ifdef SOLVELOG
    momDest[i] = exp(momSource[i]) * rho;
#else
    momDest[i] = momSource[i] * rho;
#endif
  }
}

void TSoot::MOverRhoToM(Double * momSource, Double * momDest, int nMom, Double rho)
{
  for (int i = 0; i < nMom; ++i)
  {
    momDest[i] = momSource[i] * rho;
  }
}

//#ifndef ZEROD
void T1DSoot::UpdateDimensions(int len)
{
  fMoments->cols = len;
}

void T1DSoot::UpdateSolution(T1DFlamePtr flame, Double * y, int gridPoint)
{
  int		nSpeciesInSystem = flame->GetSpecies()->GetNSpeciesInSystem();
  int		firstSpeciesOff = flame->GetOffsetFirstSpecies();
  double *	moments = fMoments->mat[gridPoint];
  double		mixMM, density;

  flame->GetProperties()->ComputeMixtureMolarMass(mixMM, &y[firstSpeciesOff], flame->GetSpecies()->GetMolarMass()->vec, nSpeciesInSystem);
  density = flame->GetPressure() * mixMM / (RGAS * y[flame->GetOffsetTemperature()]);
  
  for (int i = 0; i < fNSootMoments; ++i)
  {
#ifdef SOLVELOG
    moments[i] = exp(y[fOffsetMoments+i]) * density;
#else
    moments[i] = y[fOffsetMoments+i] * density;
#endif
  }
}

void T1DSoot::UpdateSolution(T1DFlamePtr flame, MatrixPtr yMat, VectorPtr yLeftVec, VectorPtr yRightVec)
{
  T1DPropertiesPtr prop = flame->GetProperties();
  int nGridPoints = yMat->cols;
  int nSpeciesInSystem = flame->GetSpecies()->GetNSpeciesInSystem();
  int firstSpeciesOff = flame->GetOffsetFirstSpecies();
  int tempOff = flame->GetOffsetTemperature();
  double ** y = yMat->mat;
  double * yLeft = yLeftVec->vec;
  double * yRight = yRightVec->vec;
  double ** moments = fMoments->mat;
  double * molarMass = flame->GetSpecies()->GetMolarMass()->vec;
  double mixMM, density;
  double pressure = flame->GetPressure();

  //Set boundary values at left side
  prop->ComputeMixtureMolarMass(mixMM, &yLeft[firstSpeciesOff], molarMass, nSpeciesInSystem);
  density = pressure * mixMM / (RGAS * yLeft[tempOff]);
  
  for (int i = 0; i < fNSootMoments; ++i)
  {
#ifdef SOLVELOG
    moments[kPrev][i] = exp(yLeft[fOffsetMoments+i]) * density;
#else
    moments[kPrev][i] = yLeft[fOffsetMoments+i] * density;
#endif
  }

  //Set boundary value at right side
  prop->ComputeMixtureMolarMass(mixMM, &yRight[firstSpeciesOff], molarMass, nSpeciesInSystem);
  density = pressure * mixMM / (RGAS * yRight[tempOff]);
  
  for (int i = 0; i < fNSootMoments; ++i)
  {
#ifdef SOLVELOG
    moments[nGridPoints][i] = exp(yRight[fOffsetMoments+i]) * density;
#else
    moments[nGridPoints][i] = yRight[fOffsetMoments+i] * density;
#endif
  }

//   CheckSolution(nGridPoints, y);

  //Set inner values
  for (int k = 0; k < nGridPoints; ++k)
  {
    prop->ComputeMixtureMolarMass(mixMM, &y[k][firstSpeciesOff], molarMass, nSpeciesInSystem);
    density = pressure * mixMM / (RGAS * y[k][tempOff]);
    
    for (int i = 0; i < fNSootMoments; ++i)
    {
#ifdef SOLVELOG
      moments[k][i] = exp(y[k][fOffsetMoments+i]) * density;
#else
      moments[k][i] = y[k][fOffsetMoments+i] * density;
#endif
    }
  }
}

void T1DSoot::SolutionToSolver(Double ** y)
{
  int		nGridPoints = fMoments->cols;
  double **	moments = fMoments->mat;
  
  fprintf(stderr, "check if 'T1DSoot::SolutionToSolver' is correct\n");
  
  for (int k = 0; k < nGridPoints; ++k)
  {
    for (int i = 0; i < fNSootMoments; ++i)
    {
      y[k][fOffsetMoments+i] = moments[k][i];
    }
  }
}

void T1DSoot::SaveSolution(void)
{
  int		len = fMoments->cols;
  double **	moments = fMoments->mat;
  double **	savedMoments = fSavedMoments->mat;

//   fprintf(stderr, "check if 'T1DSoot::SaveSolution' is correct\n");

  fSavedMoments->cols = fMoments->cols;
  fSavedMoments->rows = fMoments->rows;
  
  for (int k = -1; k <= len; ++k)
  {
    for (int i = 0; i < fNSootMoments; ++i)
    {
      savedMoments[k][i] = moments[k][i];
    }
  }
}

void T1DSoot::RestoreSolution(void)
{
  int		len = fSavedMoments->cols;
  double **	moments = fMoments->mat;
  double **	savedMoments = fSavedMoments->mat;

  fprintf(stderr, "check if 'T1DSoot::RestoreSolution' is correct\n");
  
  for (int k = -1; k <= len; ++k)
  {
    for (int i = 0; i < fNSootMoments; ++i)
    {
      moments[k][i] = savedMoments[k][i];
    }
  }
}

#ifdef VSH
void T1DSoot::FillRHS(T1DFlamePtr flame, NodeInfoPtr nodeInfo, CoordType coordinate)
{
  if (!fSootDASSL)
  {
    int i, j, k, l, ioff;
    TFlameNodePtr flameNode = flame->GetFlameNode();
    double * molarMass = flame->GetSpecies()->GetMolarMass()->vec;
    double temp = flameNode->temp[kCurr];
    double density = flameNode->mixDensity[kCurr];
    double * rhs = nodeInfo->rhs;
    double * moments = flameNode->moments;
    double * Y = flameNode->Y[kCurr];
    double * kSoot = fSootRateCoeffs->vec;
    double * wSoot = fSootReactionRate->vec;
    double * theSource = New1DArray(fNSootMoments);
    double rhodot = flameNode->rhodot[kCurr];

    ComputeSootRateCoeffs(kSoot, temp, flame->GetReaction());
    ComputeCSootStar(kSoot, Y, density, molarMass, flameNode->mixMolarMass[kCurr]);
    ComputeSootReactionRates(wSoot, Y, moments, density, molarMass);

    if (coordinate == kPhysical)
    {
      int fVVelocity = flame->GetOffsetVVelocity();
      for (l = 0; l < fNSootMoments; ++l)
      {
	ioff = fOffsetMoments + l;
	i = Geti(l);
	j = Getj(l);
	k = Getk(l);

	//Convection
	rhs[ioff] += NonlinearConvectUpwind(nodeInfo->y[fVVelocity] , nodeInfo->yPrev[ioff], nodeInfo->y[ioff] , nodeInfo->yNext[ioff], nodeInfo->hm, nodeInfo->h);

	//Density Correction
	rhs[ioff] += rhodot * nodeInfo->y[ioff];

	//Diffusion
	if (fSizeDepDiff)
	{
//#ifndef NODIFF
	  //rhs[ioff] -= SootDiffusion(i, j, k, kPhysical, flame, nodeInfo);
//#endif
	}
	else
	{
	  //rhs[ioff] -= SootDiffusionNew(i, j, k, kPhysical, flame, nodeInfo);
	}

	//Thermophoresis
	if (fThermoPhoresis)
	{
	  rhs[ioff] -= SootThermoPhoresis(l, kPhysical, flame, nodeInfo);
	}

	//Source Terms
	rhs[ioff] -= flameNode->sootSource[l];
      }
    }
    else
    {
      fprintf(stderr, "#error: invalid coordinate type %d\n", coordinate);
      exit(2);
    }

    Free1DArray(theSource);
  }
  else
  {
    TFlameNodePtr flameNode = flame->GetFlameNode();
    double * rhs = nodeInfo->rhs;
    int i, ioff;
  
    //If using DASSL, automatically satisfy the Newton solver.
    for (i = 0; i < fNSootMoments; ++i)
    {
      ioff = i + fOffsetMoments;
      rhs[ioff] = nodeInfo->y[ioff] - ySaved->mat[nodeInfo->gridPoint][i];
    }
  }
}
#endif
#ifdef VS
void T1DSoot::FillRHS(T1DFlamePtr flame, NodeInfoPtr nodeInfo, CoordType coordinate)
{
  if (!fSootDASSL)
  {
    int k, ioff, koff;
    double i, j;
    TFlameNodePtr flameNode = flame->GetFlameNode();
    double * molarMass = flame->GetSpecies()->GetMolarMass()->vec;
    double temp = flameNode->temp[kCurr];
    double density = flameNode->mixDensity[kCurr];
    double * rhs = nodeInfo->rhs;
    double * moments = flameNode->moments;
    double * Y = flameNode->Y[kCurr];
    double * kSoot = fSootRateCoeffs->vec;
    double * theSource = New1DArray(fNSootMoments);
    double rhodot = flameNode->rhodot[kCurr];
 
    rho = flameNode->mixDensity[kCurr];
    mu = flameNode->mixViscosity[kCurr];
    Wmix = flameNode->mixMolarMass[kCurr];   

    ComputeDimerParticules(flame, moments);
    ComputeSootRateCoeffs(kSoot, temp, flame->GetReaction());
    ComputeCSootStar(kSoot, Y, density, molarMass, flameNode->mixMolarMass[kCurr]);

    if (coordinate == kPhysical)
    {
      int fVVelocity = flame->GetOffsetVVelocity();

#ifdef MOMIC
      for (k = 0; k < fNSootMoments; ++k)
#endif
#ifdef HMOM
      for (k = 0; k < fNSootMoments-1; ++k)
#endif
      {
	ioff = fOffsetMoments + k;
	i = double(Geti(k))/6.0;
	j = double(Getj(k))/6.0;

	//Convection
	rhs[ioff] += NonlinearConvectUpwind(nodeInfo->y[fVVelocity] , nodeInfo->yPrev[ioff], nodeInfo->y[ioff] , nodeInfo->yNext[ioff], nodeInfo->hm, nodeInfo->h);

	//Density Correction
	#ifdef SOLVELOG
	rhs[ioff] += rhodot;
	#else
	rhs[ioff] += rhodot * nodeInfo->y[ioff];
	#endif

	//Diffusion
	if (fSizeDepDiff)
	{
	  rhs[ioff] -= SootDiffusion(i, j, kPhysical, flame, nodeInfo);
	}
	else
	{
	  rhs[ioff] -= SootDiffusionNew(k, kPhysical, flame, nodeInfo);
	}

	//Thermophoresis
	if (fThermoPhoresis)
	{
	  rhs[ioff] -= SootThermoPhoresis(k, kPhysical, flame, nodeInfo);
	}

	//Source Terms
	rhs[ioff] -= flameNode->sootSource[k];
      }

#ifdef HMOM
      ioff = fOffsetMoments + fNSootMoments - 1;

      // Convection
      rhs[ioff] += NonlinearConvectUpwind(nodeInfo->y[fVVelocity] , nodeInfo->yPrev[ioff], nodeInfo->y[ioff] , nodeInfo->yNext[ioff], nodeInfo->hm, nodeInfo->h);

      //Density Correction
      #ifdef SOLVELOG
      rhs[ioff] += rhodot;
      #else
      rhs[ioff] += rhodot * nodeInfo->y[ioff];
      #endif

      //Diffusion
      if (fSizeDepDiff)
      {
	rhs[ioff] -= SootDiffusionSmall(kPhysical, flame, nodeInfo);
      }
      else
      {
	rhs[ioff] -= SootDiffusionNew(fNSootMoments-1, kPhysical, flame, nodeInfo);
      }

      //Thermophoresis
      if (fThermoPhoresis)
	rhs[ioff] -= SootThermoPhoresis(fNSootMoments-1, kPhysical, flame, nodeInfo);
      
      //Source Terms
      rhs[ioff] -= flameNode->sootSource[fNSootMoments-1];
#endif HMOM
    }
    else if (coordinate == kMixtureFraction)
    {
      for (int k = 0; k < fNSootMoments; ++k)
      {
	koff = k + fOffsetMoments;
	i = double(Geti(k))/6.0;
	j = double(Getj(k))/6.0;

	// Soot Source Terms
	if (k != (fNSootMoments-1))
	{
	  if (fNucleation)
	    rhs[koff] += NucleationSource(i, j, temp) / (moments[k] / density);
	  if (fCoagulation)
	    rhs[koff] += CoagulationSource(i, j, temp, moments) / (moments[k] / density);
	  if (fCondensation)
	    rhs[koff] += CondensationSource(i, j, temp, moments) / (moments[k] / density);
	  if (fSurfaceGrowth)
	    rhs[koff] += SurfaceGrowthSource(i, j, moments, Y, temp, density, molarMass) / (moments[k] / density);
	  if (fSurfaceOxidation)
	    rhs[koff] += OxidationSource(i, j, moments, Y, temp, density, molarMass) / (moments[k] / density);
	}
	else
	{
	  if (fNucleation)
	    rhs[koff] += NucleationSource(0, 0, temp) / (moments[k] / density);
	  if (fCoagulation)
	    rhs[koff] += CoagulationSourceSmall(temp, moments) / (moments[k] / density);
	  if (fCondensation)
	    rhs[koff] += CondensationSourceSmall(temp, moments) / (moments[k] / density);
	  if (fSurfaceGrowth)
	    rhs[koff] += SurfaceGrowthSourceSmall(moments, Y, temp, density, molarMass) / (moments[k] / density);
	  if (fSurfaceOxidation)
	    rhs[koff] += OxidationSourceSmall(moments, Y, temp, density, molarMass) / (moments[k] / density);
	}
      }
    }
    else
    {
      fprintf(stderr, "#error: invalid coordinate type %d\n", coordinate);
      exit(2);
    }

    Free1DArray(theSource);
  }
  else
  {
    TFlameNodePtr flameNode = flame->GetFlameNode();
    double * rhs = nodeInfo->rhs;
    int i, ioff;
  
    //If using DASSL, automatically satisfy the Newton solver.
    for (i = 0; i < fNSootMoments; ++i)
    {
      ioff = i + fOffsetMoments;
      rhs[ioff] = nodeInfo->y[ioff] - ySaved->mat[nodeInfo->gridPoint][i];
    }
  }
}
#endif
#ifdef V
void T1DSoot::FillRHS(T1DFlamePtr flame, NodeInfoPtr nodeInfo, CoordType coordinate)
{
  if (!fSootDASSL)
  {
    int		i, j, k, ioff;
    TFlameNodePtr	flameNode = flame->GetFlameNode();
    double		*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
    double		temp = flameNode->temp[kCurr];
    double		density = flameNode->mixDensity[kCurr];
    double		*rhs = nodeInfo->rhs;
//     double		*pahMoments = flameNode->pahMoments;
    double		*moments = flameNode->moments;
    double		*Y = flameNode->Y[kCurr];
    double		*kSoot = fSootRateCoeffs->vec;
    double		*theSource = New1DArray(fNSootMoments);
    double		rhodot = flameNode->rhodot[kCurr];

    ComputeSootRateCoeffs(kSoot, temp, flame->GetReaction());
    ComputeCSootStar(kSoot, Y, density, molarMass, flameNode->mixMolarMass[kCurr]);

    if (coordinate == kPhysical)
    {
      int fVVelocity = flame->GetOffsetVVelocity();
      for (k = 0; k < fNSootMoments; ++k)
      {
	ioff = fOffsetMoments + k;
	i = Geti(k);
	
	//Convection
	rhs[ioff] += NonlinearConvectUpwind(nodeInfo->y[fVVelocity] , nodeInfo->yPrev[ioff], nodeInfo->y[ioff] , nodeInfo->yNext[ioff], nodeInfo->hm, nodeInfo->h);
	
	//Density Correction
	rhs[ioff] += rhodot * nodeInfo->y[ioff];
	
	//Diffusion
	if (fSizeDepDiff)
	{
#ifndef NODIFF
	  rhs[ioff] -= SootDiffusion(i, kPhysical, flame, nodeInfo);
#endif
	}
	else
	{
	  rhs[ioff] -= SootDiffusionNew(i, kPhysical, flame, nodeInfo);
	}
	
	//Thermophoresis
	if (fThermoPhoresis)
	{
	  rhs[ioff] -= SootThermoPhoresis(k, kPhysical, flame, nodeInfo);
	}

	//Source Terms
	rhs[ioff] -= flameNode->sootSource[k];
      }
    }
/*	else if (coordinate == kSimilarity)
	{
		int	fVVelocity = flame->GetOffsetVVelocity();
		Double	oneOverRhoMuRef = 1.0 / (flameNode->rhoInf * flameNode->viscosityInf);
		Double	oneOverRhoA = 1.0 / (flameNode->mixDensity[kCurr] * flame->GetStrainRate());

		for (k = 0; k < fNSootMoments; ++k)
		{
				ioff = fOffsetMoments + k;
				i = Geti(k);
				
				// convection
				rhs[ioff] += NonlinearConvectUpwind(nodeInfo->y[fVVelocity] ,nodeInfo->yPrev[ioff], nodeInfo->y[ioff], nodeInfo->yNext[ioff], nodeInfo->hm, nodeInfo->h, FALSE);

				// diffusion
				if ( fSizeDepDiff )
				{
					rhs[ioff] += oneOverRhoMuRef * SootDiffusion(i, j, kSimilarity, flame, nodeInfo);
				}
				else
				{
					rhs[ioff] += oneOverRhoMuRef * SootDiffusionNew(i, j, kSimilarity, flame, nodeInfo);
				}

				// thermophoresis
				if (fThermoPhoresis)
				{
					rhs[ioff] += oneOverRhoMuRef * SootThermoPhoresis(i, j, kSimilarity, flame, nodeInfo);
				}

				// nucleation
				if (fNucleation)
				{
					rhs[ioff] += oneOverRhoA * NucleationNew(i, j, temp, pahMoments, pahMoments);
				}

				// coagulation
				if (fCoagulation)
				{
					fprintf( stderr, "###Attention: old Coagulation used here\n" );
//					rhs[ioff] += oneOverRhoA * SourceCoagulation(i);
				}
		
				// condensation
				if (fCondensation)
				{
					rhs[ioff] += oneOverRhoA * SourceCondensationNew(i, j, temp, pahMoments, moments);
				}
	
				// surface growth
				if (fSurfaceGrowth)
				{
					rhs[ioff] += oneOverRhoA * SourceSurfGrowthNew(i, j, moments, Y, density, molarMass);
				}
			
				// surface oxidation
				if (fSurfaceOxidation)
				{
					rhs[ioff] += oneOverRhoA * SourceSootOxidationNew(i, j, moments, Y, density, molarMass);
				}
			}
		}
	}*/ //Mueller
/*	else if ( coordinate == kMixtureFraction ) { //Mueller
		for ( i = 0; i < fNSootMoments; ++i ) {
			ioff = i + fOffsetMoments;

			//	diffusion
			if ( fSizeDepDiff ) {
//				rhs[ioff] += SecondDeriv( nodeInfo->yPrev[ioff], nodeInfo->y[ioff]
//								, nodeInfo->yNext[ioff], nodeInfo->hm, nodeInfo->h );
//				rhs[ioff] += SootDiffusion( i, kMixtureFraction, flame, nodeInfo )
//							/ ( 2.0 * GetLewis1() );
			}
			else {
//				rhs[ioff] += SootDiffusionNew( i, kMixtureFraction, flame, nodeInfo )
//						/ ( 2.0 * GetLewis1() );
			}

			// nucleation
			if ( fNucleation ) {
//				fprintf( stderr, "Nuc = %g\n", NucleationNew( i, temp, pahMoments ) );
				rhs[ioff] += NucleationNew( i, temp, pahMoments ) / density;
			}

			//	coagulation
			if ( fCoagulation ) {
//				rhs[ioff] += SourceCoagulation( i ) / density;
				rhs[ioff] += SourceCoagulationNew( i, temp, moments ) / density;
			}

			//	condensation
			if ( fCondensation ) {
				rhs[ioff] += SourceCondensationNew( i, temp, pahMoments, moments, Y, density, molarMass ) / density;
			}

			//	surface growth
			if ( fSurfaceGrowth ) {
//				rhs[ioff] += SourceSurfGrowthNew( i, Y, density, molarMass ) / density;
				rhs[ioff] += SourceSurfGrowthNew( i, moments, Y, density, molarMass ) / density;
			}

			//	surface oxidation
			if ( fSurfaceOxidation ) {
				rhs[ioff] += SourceSootOxidationNew( i, moments, Y, density, molarMass ) / density;
			}
		}
	}*/ //Mueller
    else
    {
      fprintf(stderr, "#error: invalid coordinate type %d\n", coordinate);
      exit(2);
    }
    
    Free1DArray(theSource);
  }
  else
  {
    TFlameNodePtr flameNode = flame->GetFlameNode ();
    double * rhs = nodeInfo->rhs;
    int i, ioff;
  
    //If using DASSL, automatically satisfy the Newton solver.
    for (i = 0; i < fNSootMoments; ++i)
    {
      ioff = i + fOffsetMoments;
      rhs[ioff] = nodeInfo->y[ioff] - ySaved->mat[nodeInfo->gridPoint][i];
    }
  }
}
#endif

#ifdef VSH
void T1DSoot::FillSource(Double * source, T1DFlamePtr flame)
{
  double i, j, k;
  TFlameNodePtr	flameNode = flame->GetFlameNode();
  double * molarMass = flame->GetSpecies()->GetMolarMass()->vec;
  double temp = flameNode->temp[kCurr];
  double density = flameNode->mixDensity[kCurr];
  double * moments = flameNode->moments;
  double * Y = flameNode->Y[kCurr];
  double * kSoot = fSootRateCoeffs->vec;

  rho = flameNode->mixDensity[kCurr];
  mu = flameNode->mixViscosity[kCurr];
  Wmix = flameNode->mixMolarMass[kCurr];

  ComputeDimerParticules(flame, moments);
  ComputeSootRateCoeffs(kSoot, temp, flame->GetReaction());
  ComputeCSootStar(kSoot, Y, density, molarMass, flameNode->mixMolarMass[kCurr]);

#ifdef MOMIC
  for (int l = 0; l < fNSootMoments; ++l)
#endif
#ifdef HMOM
  for (int l = 0; l < fNSootMoments-1; ++l)
#endif  
  {
    i = double(Geti(l))/6.0;
    j = double(Getj(l))/6.0;
    k = double(Getk(l))/6.0;

    source[l] = 0.0;

    if (fNucleation)
      source[l] += NucleationSource(i, j, k, temp);
    
    if (fCoagulation)
      source[l] += CoagulationSource(i, j, k, temp, moments);
    
    if (fCondensation)
      source[l] += CondensationSource(i, j, k, temp, moments);
    
    if (fSurfaceGrowth)
      source[l] += SurfaceGrowthSource(i, j, k, moments, Y, temp, density, molarMass);
    
    if (fSurfaceOxidation)
      source[l] += OxidationSource(i, j, k, moments, Y, temp, density, molarMass);

#ifdef HMOM
  //Source terms for number of small particles
  int r = (fNSootMoments - 1);

  source[r] = 0.0;

  if (fNucleation)
    source[r] += NucleationSource(0, 0, 0, temp);
  
  if (fCoagulation)
    source[r] += CoagulationSourceSmall(temp, moments);
  
  if (fCondensation)
    source[r] += CondensationSourceSmall(temp, moments);
  
  if (fSurfaceGrowth)
    source[r] += SurfaceGrowthSourceSmall(moments, Y, temp, density, molarMass);
  
  if (fSurfaceOxidation)
    source[r] += OxidationSourceSmall(momets, Y, temp, density, molarMass);
#endif
  }
}
#endif
#ifdef VS
void T1DSoot::FillSource(double * source, T1DFlamePtr flame)
{
  double i, j;
  TFlameNodePtr	flameNode = flame->GetFlameNode();
  double * molarMass = flame->GetSpecies()->GetMolarMass()->vec;
  double temp = flameNode->temp[kCurr];
  double density = flameNode->mixDensity[kCurr];
  double * moments = flameNode->moments;
  double * Y = flameNode->Y[kCurr];
  double * kSoot = fSootRateCoeffs->vec;

  rho = flameNode->mixDensity[kCurr];
  mu = flameNode->mixViscosity[kCurr];
  Wmix = flameNode->mixMolarMass[kCurr];

  ComputeDimerParticules(flame, moments);
  ComputeSootRateCoeffs(kSoot, temp, flame->GetReaction());
  ComputeCSootStar(kSoot, Y, density, molarMass, flameNode->mixMolarMass[kCurr]);

#ifdef MOMIC
  for (int k = 0; k < fNSootMoments; ++k)
#endif
#ifdef HMOM
  for (int k = 0; k < fNSootMoments-1; ++k)
#endif
  {
    i = double(Geti(k))/6.0;
    j = double(Getj(k))/6.0;

    source[k] = 0.0;
    
    if (fNucleation)
      source[k] += NucleationSource(i, j, temp);

    if (fCoagulation)
      source[k] += 0.0 *CoagulationSource(i, j, temp, moments);
    
    if (fCondensation)
      source[k] += 0.0 *CondensationSource(i, j, temp, moments);

    if (fSurfaceGrowth)
      source[k] += 0.0 *SurfaceGrowthSource(i, j, moments, Y, temp, density, molarMass);
    
    if (fSurfaceOxidation)
      source[k] += OxidationSource(i, j, moments, Y, temp, density, molarMass);

    if (fFragmentation)
      source[k] += FragmentationSource(i, j, moments);
  }

#ifdef HMOM
  //Source terms for number of small particles
  int r = (fNSootMoments - 1);

  source[r] = 0.0;
  
    if (fNucleation)
    source[r] += NucleationSource(0, 0, temp);
  
  if (fCoagulation)
    source[r] += 0.0 *CoagulationSourceSmall(temp, moments);
  
  if (fCondensation)
    source[r] += 0.0 *CondensationSourceSmall(temp, moments);
  
  if (fSurfaceGrowth)
    source[r] += 0.0 *SurfaceGrowthSourceSmall(moments, Y, temp, density, molarMass);
  
  if (fSurfaceOxidation)
    source[r] += OxidationSourceSmall(moments, Y, temp, density, molarMass);

  if (fFragmentation)
    source[r] += FragmentationSourceSmall(moments);
#endif

#ifdef SOLVELOG
  for (int k = 0; k < fNSootMoments; ++k)
  {
    source[k] /= (moments[k] / density);
  }
#endif
}
#endif

double TSoot::CoagulationSourceSmall(double temp, double * moments)
{
#ifdef V
#ifdef FRENKLACH
  return -GetPsiSMALL(moments, temp);
#endif
#ifdef HMOM
  double C_fm = GetC(temp);

  double V0 = 2.0 * dimer_nbrC2;

  //Collision of the small particles with themselves
  double ssfm;
  ssfm = -C_fm * pow(2.0, 2.5) * pow(1.0/(2.0*dimer_nbrC2), 0.5) * pow(2.0*dimer_nbrC2, 2.0/3.0) * moments[fNSootMoments-1] * moments[fNSootMoments-1];

  //Collision of the small particles with the large particles
  double slfm;
  slfm = -GetPsiSL(0.0, 0.0, moments, temp);
  
  return ssfm + slfm;
#endif
#endif
#ifdef VS
  double C_fm = GetC(temp);
  double C_cont = 8.0 * RGAS/AVOGADRO * temp / (3.0*mu) * AVOGADRO;
  double lambda = 3.0*mu/rho*sqrt(PI*Wmix/(8.0*RGAS*temp)) / pow(6.0*fMolarMassSoot/(PI*fSootDensity*AVOGADRO),1.0/3.0);

  double V0 = 2.0 * dimer_nbrC2;

  double Df = 1.8;
  double av = 1.0 - (2.0 / Df);
  double as = (3.0 / Df) - 1.0;

  //Collision of the small particles with themselves
  double ss;
  double ssfm;
  double sscont;

  ssfm = -C_fm * pow(2.0, 2.5) * pow(1.0/(2.0*dimer_nbrC2), 0.5) * pow(2.0*dimer_nbrC2, 2.0/3.0) * moments[fNSootMoments-1] * moments[fNSootMoments-1];
  sscont = -4.0 * C_cont * (1 + 1.257*lambda*pow(V0,-1.0/3.0)) * moments[fNSootMoments-1] * moments[fNSootMoments-1];

  if (ssfm*sscont == 0.0)
    ss = 0.0;
  else
    ss = (ssfm * sscont) / (ssfm + sscont);

  //Collision of the small particles with the large particles
  double sl;
  double slfm;
  double slcont;

  slfm = -GetPsiSL(0.0,0.0,0.0,0.0,moments,temp);

  slcont = -C_cont * (2.000*                        moments[fNSootMoments-1]*FracMomLarge( 0.0   , 0.0   ,moments) + 
		                   pow(V0, 1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-    av,-    as,moments) +
		                   pow(V0,-1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av,     as,moments) +
		      1.257*lambda*pow(V0,-1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge( 0.0   , 0.0   ,moments) +
		      1.257*lambda*                 moments[fNSootMoments-1]*FracMomLarge(-    av,-    as,moments) +
		      1.257*lambda*pow(V0,-2.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av,     as,moments) +
		      1.257*lambda*pow(V0, 1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-2.0*av,-2.0*as,moments));

  if (slfm*slcont == 0.0)
    sl = 0.0;
  else
    sl = (slfm * slcont) / (slfm + slcont);
  
  return ss + sl;
#endif
#ifdef VSH
  double C_fm = GetC(temp);
  double C_cont = 8.0 * RGAS/AVOGADRO * temp / (3.0*mu) * AVOGADRO;
  double lambda = 3.0*mu/rho*sqrt(PI*Wmix/(8.0*RGAS*temp)) / pow(6.0*fMolarMassSoot/(PI*fSootDensity*AVOGADRO),1.0/3.0);

  double V0 = 2.0 * dimer_nbrC2;

  double Df = 1.8;
  double av = 1.0 - (2.0 / Df);
  double as = (3.0 / Df) - 1.0;

  //Collision of the small particles with themselves
  double ss;
  double ssfm;
  double sscont;

  ssfm = -C_fm * pow(2.0, 2.5) * pow(1.0/(2.0*dimer_nbrC2), 0.5) * pow(2.0*dimer_nbrC2, 2.0/3.0) * moments[fNSootMoments-1] * moments[fNSootMoments-1];
  sscont = -4.0 * C_cont * (1 + 1.257*lambda*pow(V0,-1.0/3.0)) * moments[fNSootMoments-1] * moments[fNSootMoments-1];

  //ss = (ssfm * sscont) / (ssfm + sscont);
  ss = ssfm;

  //Collision of the small particles with the large particles
  double sl;
  double slfm;
  double slcont;

  slfm = -GetPsiSL(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, moments, temp);

  slcont = -C_cont * (2.000*                        moments[fNSootMoments-1]*FracMomLarge( 0.0   , 0.0   , 0.0, moments) + 
		                   pow(V0, 1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-    av,-    as, 0.0, moments) +
		                   pow(V0,-1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av,     as, 0.0, moments) +
		      1.257*lambda*pow(V0,-1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge( 0.0   , 0.0   , 0.0, moments) +
		      1.257*lambda*                 moments[fNSootMoments-1]*FracMomLarge(-    av,-    as, 0.0, moments) +
		      1.257*lambda*pow(V0,-2.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(     av,     as, 0.0, moments) +
		      1.257*lambda*pow(V0, 1.0/3.0)*moments[fNSootMoments-1]*FracMomLarge(-2.0*av,-2.0*as, 0.0, moments));

  //sl = (slfm * slcont) / (slfm + slcont);
  sl = slfm;

  return ss + sl;
#endif
}

double TSoot::CondensationSourceSmall(double temp, double * moments)
{
  double C_fm = GetC(temp) / 2.2;
  double dV = dimer_nbrC2;
  double V0 = 2.0 * dimer_nbrC2;
  double lambda = 3.0*mu/rho*sqrt(PI*Wmix/(8.0*RGAS*temp)) / pow(6.0*fMolarMassSoot/(PI*fSootDensity*AVOGADRO),1.0/3.0);
  double C_cont = 8.0 * RGAS/AVOGADRO * temp / (3.0*mu) * AVOGADRO;

  double fm;
  double cont;

  return -C_fm * pow(((1/V0) + (1/dV)), 0.5) * pow((pow(V0, 1.0/3.0) + pow(dV, 1.0/3.0)), 2.0) * dimer_conc * moments[fNSootMoments-1];

//   double dc_dV = pow(dV,1.0/3.0);
//   double dc_V0 = pow(V0,1.0/3.0);

//   cont = -C_cont * (2.0 + dc_V0/dc_dV + dc_dV/dc_V0 + 1.257*lambda/dc_V0 + 1.257*lambda/dc_dV + 1.257*lambda*dc_dV/(dc_V0*dc_V0) + 1.257*lambda*dc_V0/(dc_dV*dc_dV)) * dimer_conc * moments[fNSootMoments-1];

//   return (fm * cont) / (fm + cont);
}

double TSoot::SurfaceGrowthSourceSmall(double * moments, double * Y, double temp, double density, double * molarMass)
{
  double ksg = SurfaceGrowthCoeffFor(Y, temp, density, molarMass) - SurfaceGrowthCoeffBack(Y, temp, density, molarMass);
  double Chi = fChi * pow(36.0*PI, 1.0/3.0) * pow(fMolarMassSoot / (AVOGADRO * fSootDensity), 2.0/3.0);
#ifndef VSH
  double H0 = Chi * pow(2.0 * dimer_nbrC2, 2.0/3.0);
#else
  double H0 = 2.0 * dimer_nbrH;
#endif

  return - ksg * H0 * moments[fNSootMoments-1];
}

double TSoot::OxidationSourceSmall(double * moments, double * Y, double temp, double density, double * molarMass)
{
  double kox = OxidationCoeff(Y, temp, density, molarMass);
  double Chi = fChi * pow(36.0*PI, 1.0/3.0) * pow(fMolarMassSoot / (AVOGADRO * fSootDensity), 2.0/3.0);
  double V0 = 2.0 * dimer_nbrC2;
  double dV = 1.0;

#ifndef VSH
  double H0 = Chi * pow(2.0 * dimer_nbrC2, 2.0/3.0);
#else
  double H0 = 2.0 * dimer_nbrH;
#endif

  double dV_V0 = dV / V0;

  // No model for moving particles to the first bin
  //return -kox * dV_V0 * H0 * moments[fNSootMoments-1];

  // Model for moving particles to the first bin
  double C_inter;
  C_inter = V0 / (FracMomLarge(1.0,0.0,moments) / FracMomLarge(0.0,0.0,moments)); //Revised CI Paper: Volume Ratio

  double small = -kox * dV_V0 * H0 * moments[fNSootMoments-1];
  double large = C_inter * kox * Chi * dV * FracMomLarge(-1.0,1.0,moments);
  //return small;
  return small + large;
}

double TSoot::FragmentationSourceSmall(double * moments)
{
  // No model for moving particles to the first bin
  //return 0.0;

  double V0 = 2.0*dimer_nbrC2;

  double C_inter;
  C_inter = V0 / (FracMomLarge(1.0,0.0,moments) / FracMomLarge(0.0,0.0,moments)); //Revised CI Paper: Volume Ratio

  return C_inter * FragmentationSource(0.0,0.0,moments);
}

#ifdef V
void T1DSoot::FillSource(double * source, T1DFlamePtr flame)
{
  double i;
  TFlameNodePtr	flameNode = flame->GetFlameNode();
  double * molarMass = flame->GetSpecies()->GetMolarMass()->vec;
  double temp = flameNode->temp[kCurr];
  double density = flameNode->mixDensity[kCurr];
  double * moments = flameNode->moments;
  double * Y = flameNode->Y[kCurr];
  double * kSoot = fSootRateCoeffs->vec;

  rho = flameNode->mixDensity[kCurr];
  mu = flameNode->mixViscosity[kCurr];
  Wmix = flameNode->mixMolarMass[kCurr];

  ComputeDimerParticules(flame, moments);
  ComputeSootRateCoeffs(kSoot, temp, flame->GetReaction());
  ComputeCSootStar(kSoot, Y, density, molarMass, flameNode->mixMolarMass[kCurr]);

#ifdef LINDSTEDT
  source[0] = 0.0;

  if (fNucleation)
    source[0] += (1.0/50.0) * LindstedtNucleation(temp, Y, density, molarMass);
  if (fCoagulation)
    source[0] += LindstedtCoagulation(moments, temp);

  source[1] = 0.0;

  if (fNucleation)
    source[1] += LindstedtNucleation(temp, Y, density, molarMass);
  if (fSurfaceGrowth)
    source[1] += LindstedtGrowth(moments, temp, Y, density, molarMass);
  if (fSurfaceOxidation)
    source[1] += LindstedtOxidation(moments, temp, Y, density, molarMass);

#else
#ifdef MOMIC
  for (int j = 0; j < fNSootMoments; ++j)
#endif
#ifdef FRENKLACH
  for (int j = 0; j < fNSootMoments-1; ++j)
#endif
#ifdef HMOM
  for (int j = 0; j < fNSootMoments-1; ++j)
#endif	
  {
    i = double(Geti(j))/6.0;
    
    source[j] = 0.0;
    
    if (fNucleation)
      source[j] += NucleationSource(i, temp);

    if (fCoagulation)
      source[j] += CoagulationSource(i, temp, moments);

    if (fCondensation)
      source[j] += CondensationSource(i, temp, moments);

    if (fSurfaceGrowth)
      source[j] += SurfaceGrowthSource(i, moments, Y, temp, density, molarMass);

    if (fSurfaceOxidation)
      source[j] += OxidationSource(i, moments, Y, temp, density, molarMass);

#ifdef FRENKLACH
  //Source terms for number of small particles
  int r = (fNSootMoments - 1);

  source[r] = 0.0;

  if (fNucleation)
    source[r] += NucleationSource(0, temp);
  
  if (fCoagulation)
    source[r] += CoagulationSourceSmall(temp, moments);
  
  if (fCondensation)
    source[r] += CondensationSourceSmall(temp, moments);
  
  if (fSurfaceGrowth)
    source[r] += SurfaceGrowthSourceSmall(moments, Y, temp, density, molarMass);
  
  if (fSurfaceOxidation)
    source[r] += OxidationSourceSmall(moments, Y, temp, density, molarMass);
#endif

#ifdef HMOM
  //Source terms for number of small particles
  int r = (fNSootMoments - 1);

  source[r] = 0.0;

  if (fNucleation)
    source[r] += NucleationSource(0, temp);
  
  if (fCoagulation)
    source[r] += CoagulationSourceSmall(temp, moments);
  
  if (fCondensation)
    source[r] += CondensationSourceSmall(temp, moments);
  
  if (fSurfaceGrowth)
    source[r] += SurfaceGrowthSourceSmall(moments, Y, temp, density, molarMass);
  
  if (fSurfaceOxidation)
    source[r] += OxidationSourceSmall(moments, Y, temp, density, molarMass);
#endif
  }
#endif
}
#endif

Double T1DSoot::SootConvection( int i, T1DFlamePtr flame, NodeInfoPtr nodeInfo, Flag velocityPositive )
{
//	velocityPositive should have the value FALSE, if 'V' has the negative 
//	direction of the physical velocity

// returns     V * d(M_i/rho)/dy

/*	int		iOff = i + fOffsetMoments;
	Double	*y = nodeInfo->y;
	Double	*yNext = nodeInfo->yNext;
	Double	*yPrev = nodeInfo->yPrev;
	Double	*rho = flame->GetFlameNode()->mixDensity;
	Double	V = y[flame->GetOffsetVVelocity()];

	if ( ( V > 0.0 && velocityPositive ) || ( V < 0.0 && !velocityPositive ) ) {
		return ( V * ( y[iOff]/rho[kCurr] - yPrev[iOff]/rho[kPrev] ) / nodeInfo->hm );
	}
	else {
		return ( V * ( yNext[iOff]/rho[kNext] - y[iOff]/rho[kCurr] ) / nodeInfo->h );
	}*/
  return 0.0;
}

//Mueller
double T1DSoot::SootDiffusion(double i, double j, CoordType coordinate, T1DFlamePtr flame, NodeInfoPtr nodeInfo)
{
  double Df = 1.8;
  double av = 1.0 - 2.0 / Df;
  double as = 3.0 / Df - 1.0;

  TFlameNodePtr	flameNode = flame->GetFlameNode();
  double diff = flameNode->diffSoot[kCurr]*pow(2.0*dimer_nbrC2,2.0/3.0);
  double diffNext = flameNode->diffSoot[kNext]*pow(2.0*dimer_nbrC2,2.0/3.0);
  double diffPrev = flameNode->diffSoot[kPrev]*pow(2.0*dimer_nbrC2,2.0/3.0);
  double * rho = flameNode->mixDensity;
  double Index1 = i - 2.0*av;
  double Index2 = j - 2.0*as;

  double fracMom = FracMom(Index1, Index2, &nodeInfo->y[fOffsetMoments]);
  double fracMomPrev = FracMom(Index1, Index2, &nodeInfo->yPrev[fOffsetMoments]);
  double fracMomNext = FracMom(Index1, Index2, &nodeInfo->yNext[fOffsetMoments]);
  double diffPlusHm, diffMinusH;

  if (coordinate == kPhysical)
  {
    diffPlusHm = nodeInfo->hm * (diff * rho[kCurr] + diffNext * rho[kNext]);
    diffMinusH = nodeInfo->h * (diffPrev * rho[kPrev] + diff * rho[kCurr]);
  }
  else if (coordinate == kSimilarity)
  {
    diffPlusHm = nodeInfo->hm * (diff * rho[kCurr] * rho[kCurr] + diffNext * rho[kNext] * rho[kNext]);
    diffMinusH = nodeInfo->h * (diffPrev * rho[kPrev] * rho[kPrev] + diff * rho[kCurr] * rho[kCurr]);
  }
  else if (coordinate == kMixtureFraction)
  {
    diffPlusHm = 2.0 * nodeInfo->hm;
    diffMinusH = 2.0 * nodeInfo->h;
  }
	
  return (diffPlusHm * (fracMomNext - fracMom) + diffMinusH * (fracMomPrev - fracMom)) / nodeInfo->hnenn;
}

double T1DSoot::SootDiffusionSmall(CoordType coordinate, T1DFlamePtr flame, NodeInfoPtr nodeInfo)
{
  double Df = 1.8;
  double av = 1.0 - 2.0 / Df;
  double as = 3.0 / Df - 1.0;

  TFlameNodePtr	flameNode = flame->GetFlameNode();
  double diff = flameNode->diffSoot[kCurr];
  double diffNext = flameNode->diffSoot[kNext];
  double diffPrev = flameNode->diffSoot[kPrev];
  double * rho = flameNode->mixDensity;
  double Index1 = - 2.0*av;
  double Index2 = - 2.0*as;

  double fracMom = nodeInfo->y[fOffsetMoments+fNSootMoments-1];
  double fracMomPrev = nodeInfo->y[fOffsetMoments+fNSootMoments-1];
  double fracMomNext = nodeInfo->y[fOffsetMoments+fNSootMoments-1];
  double diffPlusHm, diffMinusH;

  if (coordinate == kPhysical)
  {
    diffPlusHm = nodeInfo->hm * (diff * rho[kCurr] + diffNext * rho[kNext]);
    diffMinusH = nodeInfo->h * (diffPrev * rho[kPrev] + diff * rho[kCurr]);
  }
  else if (coordinate == kSimilarity)
  {
    diffPlusHm = nodeInfo->hm * (diff * rho[kCurr] * rho[kCurr] + diffNext * rho[kNext] * rho[kNext]);
    diffMinusH = nodeInfo->h * (diffPrev * rho[kPrev] * rho[kPrev] + diff * rho[kCurr] * rho[kCurr]);
  }
  else if (coordinate == kMixtureFraction)
  {
    diffPlusHm = 2.0 * nodeInfo->hm;
    diffMinusH = 2.0 * nodeInfo->h;
  }
	
  return (diffPlusHm * (fracMomNext - fracMom) + diffMinusH * (fracMomPrev - fracMom)) / nodeInfo->hnenn;
}

Double T1DSoot::SootDiffusion( int r, CoordType coordinate, T1DFlamePtr flame, NodeInfoPtr nodeInfo )
{
	TFlameNodePtr	flameNode = flame->GetFlameNode();
/*	Double		diff = flameNode->diffusivity[f_A1];
	Double		diffNext = flameNode->diffusivityNext[f_A1];
	Double		diffPrev = flameNode->diffusivityPrev[f_A1];*/
	Double		diff = flameNode->diffSoot[kCurr];
	Double		diffNext = flameNode->diffSoot[kNext];
	Double		diffPrev = flameNode->diffSoot[kPrev];
	Double		*rho = flameNode->mixDensity;
	Double		fracIndex = r - 2.0 / 3.0;
//	Double		fracMom = FracMom2( fracIndex, &nodeInfo->y[fOffsetMoments] );
	Double		fracMom = FracMom( fracIndex, &nodeInfo->y[fOffsetMoments] );
//	Double		fracMomPrev = FracMom2( fracIndex, &nodeInfo->yPrev[fOffsetMoments] );
	Double		fracMomPrev = FracMom( fracIndex, &nodeInfo->yPrev[fOffsetMoments] );
//	Double		fracMomNext = FracMom2( fracIndex, &nodeInfo->yNext[fOffsetMoments] );
	Double		fracMomNext = FracMom( fracIndex, &nodeInfo->yNext[fOffsetMoments] );
	Double		diffPlusHm, diffMinusH;
	
	if ( coordinate == kPhysical ) {
		diffPlusHm = nodeInfo->hm * ( diff * rho[kCurr]
											+ diffNext * rho[kNext] );
		diffMinusH = nodeInfo->h * ( diffPrev * rho[kPrev]
						+ diff * rho[kCurr] );
	}
	else if ( coordinate == kSimilarity ) {
		diffPlusHm = nodeInfo->hm * ( diff * rho[kCurr] * rho[kCurr]
											+ diffNext * rho[kNext] * rho[kNext] );
		diffMinusH = nodeInfo->h * ( diffPrev * rho[kPrev] * rho[kPrev]
						+ diff * rho[kCurr] * rho[kCurr] );
	}
	else if ( coordinate == kMixtureFraction ) {
		diffPlusHm = 2.0 * nodeInfo->hm;
		diffMinusH = 2.0 * nodeInfo->h;
	}
	
	return ( diffPlusHm * ( fracMomNext - fracMom ) 
			+ diffMinusH * ( fracMomPrev - fracMom ) ) 
				/ nodeInfo->hnenn;
}

//Mueller
Double T1DSoot::SootDiffusionNew(int i, int j, CoordType coordinate, T1DFlamePtr flame, NodeInfoPtr nodeInfo)
{
	int		k = i ? (2*i) : (2*j-1);
	int		iOff = fOffsetMoments + k;
	TFlameNodePtr	flameNode = flame->GetFlameNode();
	Double		diff = flameNode->diffSoot[kCurr];
	Double		diffNext = flameNode->diffSoot[kNext];
	Double		diffPrev = flameNode->diffSoot[kPrev];
	Double		*rho = flameNode->mixDensity;
	Double		mom = nodeInfo->y[iOff];
	Double		momPrev = nodeInfo->yPrev[iOff];
	Double		momNext = nodeInfo->yNext[iOff];
	Double		diffPlusHm, diffMinusH;
	
	if (coordinate == kPhysical)
	{
		diffPlusHm = nodeInfo->hm * (diff * rho[kCurr] + diffNext * rho[kNext]);
		diffMinusH = nodeInfo->h * (diffPrev * rho[kPrev] + diff * rho[kCurr]);
	}
	else if (coordinate == kSimilarity)
	{
		diffPlusHm = nodeInfo->hm * (diff * rho[kCurr] * rho[kCurr] + diffNext * rho[kNext] * rho[kNext]);
		diffMinusH = nodeInfo->h * (diffPrev * rho[kPrev] * rho[kPrev] + diff * rho[kCurr] * rho[kCurr]);
	}
	else if (coordinate == kMixtureFraction)
	{
		diffPlusHm = 2.0 * nodeInfo->hm;
		diffMinusH = 2.0 * nodeInfo->h;
	}

	return (diffPlusHm * (momNext - mom) + diffMinusH * (momPrev - mom)) / nodeInfo->hnenn;
}

Double T1DSoot::SootDiffusionNew( int r, CoordType coordinate, T1DFlamePtr flame, NodeInfoPtr nodeInfo )
{
  int rOff = r + fOffsetMoments;
  TFlameNodePtr	flameNode = flame->GetFlameNode();
  double diff = flameNode->diffSoot[kCurr];
  double diffNext = flameNode->diffSoot[kNext];
  double diffPrev = flameNode->diffSoot[kPrev];
  double * rho = flameNode->mixDensity;
  double mom = nodeInfo->y[rOff];
  double momPrev = nodeInfo->yPrev[rOff];
  double momNext = nodeInfo->yNext[rOff];
  double diffPlusHm, diffMinusH;

  #ifdef SOLVELOG
  mom = exp(mom);
  momPrev = exp(momPrev);
  momNext = exp(momNext);
  #endif
	
  if (coordinate == kPhysical)
  {
    diffPlusHm = nodeInfo->hm * (diff * rho[kCurr] + diffNext * rho[kNext]);
    diffMinusH = nodeInfo->h * (diffPrev * rho[kPrev] + diff * rho[kCurr]);
  }
  else if (coordinate == kSimilarity)
  {
    diffPlusHm = nodeInfo->hm * (diff * rho[kCurr] * rho[kCurr] + diffNext * rho[kNext] * rho[kNext]);
    diffMinusH = nodeInfo->h * (diffPrev * rho[kPrev] * rho[kPrev] + diff * rho[kCurr] * rho[kCurr]);
  }
  else if (coordinate == kMixtureFraction)
  {
    diffPlusHm = 2.0 * nodeInfo->hm;
    diffMinusH = 2.0 * nodeInfo->h;
  }

  #ifdef SOLVELOG
  return (diffPlusHm * (momNext - mom) + diffMinusH * (momPrev - mom)) / nodeInfo->hnenn / mom;
  #else
  return (diffPlusHm * (momNext - mom) + diffMinusH * (momPrev - mom)) / nodeInfo->hnenn;
  #endif
}

//Mueller
Double T1DSoot::SootThermoPhoresis(int i, int j, CoordType coordinate, T1DFlamePtr flame, NodeInfoPtr nodeInfo)
{
  // returns d/dy ( 0.55 * nu / T * M_r * dT/dy )	

  const Double	A = 0.9; // accommodation coefficient following
                         // R. J. Santoro: The Transport and Growth
			 // of Soot Particles ...
			 // Comb. Sci. Tech., 1987, Vol. 53 p. 89
  static Double	pi = 4.0 * atan( 1.0 );
  const Double	fact = 3.0 / ( 4.0 * ( 1.0 + pi * A / 8.0 ) );
  int		k = i ? (2*i) : (2*j-1);
  int		iOff = fOffsetMoments + k;
  TFlameNodePtr	flameNode = flame->GetFlameNode();
  Double		mu = flameNode->mixViscosity[kCurr];
  Double		muNext = flameNode->mixViscosity[kNext];
  Double		muPrev = flameNode->mixViscosity[kPrev];
  Double		M = nodeInfo->y[iOff];
  Double		MNext = nodeInfo->yNext[iOff];
  Double		MPrev = nodeInfo->yPrev[iOff];
  Double		*temp = flameNode->temp;
  Double		*rho = flameNode->mixDensity;
  Double		diffPlusHm, diffMinusH;
	
  if (coordinate == kPhysical)
  {
    diffPlusHm = nodeInfo->hm * fact * (mu * M / temp[kCurr] + muNext * MNext / temp[kNext]);
    diffMinusH = nodeInfo->h * fact * (mu * M / temp[kCurr] + muPrev * MPrev / temp[kPrev]);
  }
  else
  {
    diffPlusHm = nodeInfo->hm * fact * (rho[kCurr] * mu * M / temp[kCurr] + rho[kNext] * muNext * MNext / temp[kNext]);
    diffMinusH = nodeInfo->h * fact * (rho[kCurr] * mu * M / temp[kCurr] + rho[kPrev] * muPrev * MPrev /temp[kPrev]);
  }
  return (diffPlusHm * (temp[kNext] - temp[kCurr]) + diffMinusH * (temp[kPrev] - temp[kCurr])) / nodeInfo->hnenn;
}

Double T1DSoot::SootThermoPhoresis( int r, CoordType coordinate, T1DFlamePtr flame, NodeInfoPtr nodeInfo )
{
//returns d/dy ( 0.55 * nu / T * M_r * dT/dy )	

  if (!fSootDASSL)
  {
    const Double  A = 0.9; // accommodation coefficient following
			   // R. J. Santoro: The Transport and Growth
			   // of Soot Particles ...
			   // Comb. Sci. Tech., 1987, Vol. 53 p. 89
    static Double pi = 4.0 * atan( 1.0 );
    const Double  fact = 3.0 / ( 4.0 * ( 1.0 + pi * A / 8.0 ) );
    int		  rOff = r + fOffsetMoments;
    TFlameNodePtr flameNode = flame->GetFlameNode();
    Double	  mu = flameNode->mixViscosity[kCurr];
    Double	  muNext = flameNode->mixViscosity[kNext];
    Double	  muPrev = flameNode->mixViscosity[kPrev];
    Double	  M = nodeInfo->y[rOff];
    Double	  MNext = nodeInfo->yNext[rOff];
    Double	  MPrev = nodeInfo->yPrev[rOff];
    Double *      temp = flameNode->temp;
    Double *      rho = flameNode->mixDensity;
    Double	  diffPlusHm, diffMinusH;

  #ifdef SOLVELOG
    M = exp(M);
    MNext = exp(MNext);
    MPrev = exp(MPrev);
  #endif

    if ( coordinate == kPhysical )
    {
      diffPlusHm = nodeInfo->hm * fact * ( mu * M / temp[kCurr] + muNext * MNext / temp[kNext] );
      diffMinusH = nodeInfo->h  * fact * ( mu * M / temp[kCurr] + muPrev * MPrev / temp[kPrev] );
    }
    else
    {
      diffPlusHm = nodeInfo->hm * fact * ( rho[kCurr] * mu * M / temp[kCurr] + rho[kNext] * muNext * MNext / temp[kNext] );
      diffMinusH = nodeInfo->h  * fact * ( rho[kCurr] * mu * M / temp[kCurr] + rho[kPrev] * muPrev * MPrev / temp[kPrev] );
    }

  #ifdef SOLVELOG
    return ( diffPlusHm * ( temp[kNext] - temp[kCurr] ) + diffMinusH * ( temp[kPrev] - temp[kCurr] ) ) / nodeInfo->hnenn / M;
  #else
    return ( diffPlusHm * ( temp[kNext] - temp[kCurr] ) + diffMinusH * ( temp[kPrev] - temp[kCurr] ) ) / nodeInfo->hnenn;
  #endif
  }
  else
  {
    TFlameNodePtr	flameNode = flame->GetFlameNode();

    // accommodation coefficient following
    // R. J. Santoro: The Transport and Growth
    // of Soot Particles ...
    // Comb. Sci. Tech., 1987, Vol. 53 p. 89
    static Double	pi = 4.0 * atan( 1.0 );
    static const double A = 0.9;
    static const double fact = 3.0 / ( 4.0 * ( 1.0+pi*A/8.0 ) );

    double *mu   = flameNode->mixViscosity;
    double *rho  = flameNode->mixDensity;
    double *temp = flameNode->temp;
  
    double d2 = nodeInfo->hm * fact * ( mu[kCurr]/temp[kCurr] + mu[kNext]/temp[kNext] );
    double d1 = nodeInfo->h  * fact * ( mu[kCurr]/temp[kCurr] + mu[kPrev]/temp[kPrev] );
  
    double src = ( d2*(temp[kNext]-temp[kCurr]) - d1*(temp[kCurr]-temp[kPrev]) ) / nodeInfo->hnenn;
    src *= nodeInfo->y[fOffsetMoments+r] / rho[kCurr];
  
    return src;
  }
}

//Mueller (8/1/07)
//#ifdef VS
void T1DSoot::PostIter(T1DFlamePtr flame)
{
  cerr << "If this function is called, this is bad...\n";
}

void T1DSoot::PostIter(T1DFlamePtr flame, int Instruction)
{
  if (!fSootDASSL)
    PostIterNewton(flame);
  if (fSootDASSL)
    if (Instruction == 1)
      PostIterDASSL(flame);
}

void T1DSoot::PostIterNewton(T1DFlamePtr flame)
{
  TNewtonPtr bt = flame->GetSolver()->bt;
  TGridPtr currGrid = bt->GetGrid()->GetCurrentGrid();
  int nGridPoints = currGrid->GetNGridPoints();
  double * yLeft = currGrid->GetYLeft()->vec;
  double * yRight = currGrid->GetYRight()->vec;
  double ** y = currGrid->GetY()->mat;
  double * density = flame->GetProperties()->GetDensity()->vec;

  //CheckSolution(nGridPoints, y, density);

  //SAVE this solution for later use
  for (int i=0; i < fNSootMoments; i++)
  {
    ySaved->mat[0][i] = y[0][i+fOffsetMoments];
    ySavedPrev->mat[0][i] = yLeft[i+fOffsetMoments];
    ySavedNext->mat[0][i] = y[1][i+fOffsetMoments];
  }
  for (int k=1; k<nGridPoints-1; k++)
  {
    for (int i=0; i<fNSootMoments; i++)
    {
      ySaved->mat[k][i] = y[k][i+fOffsetMoments];
      ySavedPrev->mat[k][i] = y[k-1][i+fOffsetMoments];
      ySavedNext->mat[k][i] = y[k+1][i+fOffsetMoments];
    }
  }
  for (int i=0; i < fNSootMoments; i++)
  {
    ySaved->mat[nGridPoints-1][i] = y[nGridPoints-1][i+fOffsetMoments];
    ySavedPrev->mat[nGridPoints-1][i] = y[nGridPoints-2][i+fOffsetMoments];
    ySavedNext->mat[nGridPoints-1][i] = yRight[i+fOffsetMoments];
  }

  UpdateSolution(flame, currGrid->GetY(), currGrid->GetYLeft(), currGrid->GetYRight()); //Mueller
}

void T1DSoot::PostIterDASSL(T1DFlamePtr flame)
{
  TNewtonPtr bt = flame->GetSolver()->bt;
  TGridPtr currGrid = bt->GetGrid()->GetCurrentGrid();
  int nGridPoints = currGrid->GetNGridPoints();
  double * yLeft = currGrid->GetYLeft()->vec;
  double * yRight = currGrid->GetYRight()->vec;
  double ** y = currGrid->GetY()->mat;
  int fVVelocity = flame->GetOffsetVVelocity();

  int nSteps = 1;
  if (fSootDiffusion) nSteps = 2000;

  //If it is a problem in physical space
  if (fVVelocity >= 0)
  {
    //Do the integration in 'time' for points in between
    //Only on the fine grid
    if (bt->GetGrid()->IsFine())
      for (int step=0; step<nSteps; step++) {
	TimeIntegration(flame);}
  }
  else
  {
    //If not a problem in physical space
    //Just put some default values
    SootSetDefault(&yLeft [fOffsetMoments]);
    SootSetDefault(&yRight[fOffsetMoments]);
    
    for (int k = 0; k < nGridPoints; k++)
      SootSetDefault(&y[k][fOffsetMoments]);
  }
  
  //SAVE this solution for later use
  for (int k=0; k<nGridPoints; k++)
    for (int i=0; i<fNSootMoments; i++)
      ySaved->mat[k][i] = y[k][i+fOffsetMoments];
  
  //Update the solution
  UpdateSolution(flame, currGrid->GetY(), currGrid->GetYLeft(), currGrid->GetYRight());
}
//#endif

/*
//Old VS PostIter -- PreDASSL
void T1DSoot::PostIter(T1DFlamePtr flame)
{
	int		i, j, k;
	Double		ii, jj;
	TNewtonPtr	bt = flame->GetSolver()->bt;
	TGridPtr 	currGrid = bt->GetGrid()->GetCurrentGrid();
	int		nGridPoints = currGrid->GetNGridPoints();
	Double		*yLeft = currGrid->GetYLeft()->vec;
	Double		*yRight = currGrid->GetYRight()->vec;
	Double		**y = currGrid->GetY()->mat;

	CheckSolution(nGridPoints, y);

	UpdateSolution(flame, currGrid->GetY(), currGrid->GetYLeft(), currGrid->GetYRight()); //Mueller

#ifdef MOMIC
	for (k = 0; k < fNSootMoments; ++k)
#endif
#ifdef HMOM
	for (k = 0; k < fNSootMoments-1; ++k)
#endif
	{
			ii = double(Geti(k)/6.0);
			jj = double(Getj(k)/6.0);

			if (k==0)
			{
				yLeft[fOffsetMoments+k] = SMALLSOOT;
				//yRight[fOffsetMoments+k] = SMALLSOOT;
			}
			else
			{
				yLeft[fOffsetMoments+k] = pow(2.0*nucl_nbrC2, ii+jj*(2.0/3.0)) * yLeft[fOffsetMoments];
				//yRight[fOffsetMoments+k] = pow(2.0*nucl_nbrC2, ii+jj*(2.0/3.0)) * yRight[fOffsetMoments];
			}
	}
#ifdef HMOM
	yLeft[fOffsetMoments+fNSootMoments-1] = SMALLSOOT;
	//yRight[fOffsetMoments+fNSootMoments-1] = SMALLSOOT;
#endif
}
*/

//#ifdef DASSL
//##############################//
//####DASSL Time Integration####//
//##############################//

void SootResFunc1D(double * time, double * yS, double * yP, double * delta, int * iRes, double * rPar, int * iPar)
{ 
  T1DFlamePtr flame = (T1DFlamePtr) iPar;
  int neq = flame->GetSoot()->GetNSootMoments();

  //Store time for Diffusion
  flame->GetSoot()->SetTime(*time);
  
  //Compute the source terms for the equation
  flame->GetSoot()->GetYPrime(flame,yS,delta);

  //Compute the residual
  for (int i=0; i<neq; i++)
    delta[i] -= yP[i];
}

// ** For 1D simulations **
// After each Newton Iterations
// Operate a time integration of the equations
void T1DSoot::TimeIntegration(T1DFlamePtr flame)
{    
  TNewtonPtr bt = flame->GetSolver()->bt;
  TGridPtr currGrid = bt->GetGrid()->GetCurrentGrid();
  int nGridPoints = currGrid->GetNGridPoints();

  double * yLeft = currGrid->GetYLeft()->vec;
  double * yRight = currGrid->GetYRight()->vec;
  double ** y = currGrid->GetY()->mat;
  int offVel = flame->GetOffsetVVelocity();
  int k;
  double vel, vel1;
  
  //------------- FROM THE LEFT -------------
  
  if (yLeft[offVel] > 0)
  {
    //Set initial solution
    SootSetDefault(&yLeft [fOffsetMoments]);

    for (k = 0; k<nGridPoints; k++)
    {
      //Set grid point
      flame->SetFlameNode(k);
      TFlameNodePtr flameNode = flame->GetFlameNode();
      bt->SetNodeInfo(flame, k);
      NodeInfoPtr nodeInfo = bt->GetNodeInfo();
      double density = flameNode->mixDensity[kCurr];
      
      //Get full velocity
      vel = GetFullVelocity(flame,k);
      if (k < nGridPoints-1)
	vel1 = GetFullVelocity(flame,k+1);
      else
	vel1 = vel;
      
      //Advance it this way only if velocity is positive
      if (vel > 0.0)
      //if (vel1 > 0.0)
      {
	//cerr << nodeInfo->x[kCurr] << "\t" << flameNode->moments[1]*fMolarMassSoot/fSootDensity << "\t" << vel << endl; //this is not what i want to see

	//Set time step
	double dt = nodeInfo->hm/vel;
	dt_diff = dt;
	time_diff = 0.0;
	
	//Force new problem
	double time = 0.0;
	double newTime = dt;
	int idid = 0;
	info[0] = 0;
	
	//Get a consistent derivative to start DASSL
	GetYPrime(flame,ySol,yPrime);

	//Run DASSL
	do
	{
	  DDASSL(SootResFunc1D, (int*)&fNSootMoments, &time, ySol, yPrime, &newTime, info, &rtol, &atol, &idid, rwork, &lrw, iwork, &liw, NULL,(int*)flame,NULL);
	  info[0] = 1;
	}
	while (idid == -1);

	if (idid<0) exit(2);
	
	//Copy back new solution
	for (int i = 0; i < fNSootMoments; i++)
	  y[k][i + fOffsetMoments] = ySol[i];

	//TSoot::CheckSolution(&y[k][fOffsetMoments], density);
      }
    }
    
    //Set last point
    if (vel > 0 && k == nGridPoints)
      for (int i=0; i< fNSootMoments; i++)
	yRight[i + fOffsetMoments] = ySol[i];   
  }
  
  //------------- FROM THE RIGHT -------------
  
  if (yRight[offVel] < 0)
  {
    //Set initial solution
    SootSetDefault(&yRight[fOffsetMoments]);
    
    for (k = nGridPoints-1; k >= 0; k--)
    {
      //Set grid point
      flame->SetFlameNode (k);
      TFlameNodePtr flameNode = flame->GetFlameNode ();
      bt->SetNodeInfo (flame, k);
      NodeInfoPtr nodeInfo = bt->GetNodeInfo ();
      double density = flameNode->mixDensity[kCurr];

      //Get full velocity
      vel = GetFullVelocity(flame,k);
      if (k > 0) 
	vel1 = GetFullVelocity(flame,k-1);
      else
	vel1 = vel;
      
      //Advance it this way only if velocity is positive
      if (vel < 0.0)
      //if (vel1 < 0.0)
      {
	//cerr << nodeInfo->x[kCurr] << "\t" << flameNode->moments[1]*fMolarMassSoot/fSootDensity << "\t" << vel << endl;

	//Set time step
	double dt = -nodeInfo->h/vel;
	dt_diff = -dt;
	time_diff = 0.0;
	if (k == nGridPoints-1)
	  dt_diff = 0.0;
	
	//Force new problem
	double time    = 0.0;
	double newTime = dt;
	int idid = 0;
	info[0] = 0;
	
	//Get a consistent derivative to start DASSL
	GetYPrime(flame,ySol,yPrime);
	
	//Run DASSL
	do
	{
	  DDASSL(SootResFunc1D, (int*)&fNSootMoments, &time, ySol, yPrime, &newTime, info, &rtol, &atol, &idid, rwork, &lrw, iwork, &liw, NULL,(int*)flame,NULL);
	  info[0] = 1;
	}
	while (idid == -1);

	if (idid<0) exit(2);
	
	//Copy back new solution
	for (int i = 0; i < fNSootMoments; i++)
	  y[k][i + fOffsetMoments] = ySol[i];

	//TSoot::CheckSolution(&y[k][fOffsetMoments], density);
      }
    }
    
    //Set last point
    if (vel < 0 && k == -1)
      for (int i = 0; i < fNSootMoments; i++)
	yLeft[i + fOffsetMoments] = ySol[i];
  }
}

double T1DSoot::GetFullVelocity(T1DFlamePtr flame, int k)
{
  TNewtonPtr bt = flame->GetSolver ()->bt;
  TGridPtr currGrid = bt->GetGrid ()->GetCurrentGrid ();
  TFlameNodePtr	flameNode = flame->GetFlameNode();

  double mu = flameNode->mixViscosity[kCurr];
  double rho = flameNode->mixDensity[kCurr];
  double ** y = currGrid->GetY ()->mat;
  double * x = currGrid->GetX ()->vec;
  double * temp = flameNode->temp;
  int offVel = flame->GetOffsetVVelocity ();

  //Convective Velocity
  double vel = y[k][offVel] / rho;
  double vel0 = vel;

  //Thermophoretic Velocity Correction
  if (fThermoPhoresis)
  {
    // accommodation coefficient following
    // R. J. Santoro: The Transport and Growth
    // of Soot Particles ...
    // Comb. Sci. Tech., 1987, Vol. 53 p. 89
    static const double A = 0.9;
    static const double fact = 3.0 / ( 4.0 * ( 1.0+PI*A/8.0 ) );
    
    if (vel > 0)
      vel -= fact * mu/rho * (log(temp[kCurr]) - log(temp[kPrev])) / (x[k] - x[k-1]);
    else
      vel -= fact * mu/rho * (log(temp[kNext]) - log(temp[kCurr])) / (x[k+1] - x[k]);
  }
  
  //if (vel*vel0<0) vel = vel0;

  return vel;
}

//** For 1D simulations **
//Compute the residual of the soot evolution equations
void T1DSoot::GetYPrime(T1DFlamePtr flame, double * yS, double * yP)
{  
  TFlameNodePtr flameNode = flame->GetFlameNode();
  double density = flameNode->mixDensity[kCurr];
  double * moments = flameNode->moments;

  // Copy the variables from DASSL
  for (int i = 0; i< fNSootMoments; i++)
  {
    moments[i] = yS[i];
    #ifdef SOLVELOG
    moments[i] = exp(moments[i]) * density;
    #else
    moments[i] *= density;
    #endif
  }
  
  // Compute source terms
  FillSource(flame, moments, src);

  // Recompute the derivative in canonical variables
  for (int i = 0; i < fNSootMoments; i++)
    #ifdef SOLVELOG
    yP[i] = src[i] / moments[i];
    #else
    yP[i] = src[i] / density;
    #endif
}

#ifdef VSH
void T1DSoot::FillSource(T1DFlamePtr flame, double * moments, double * source)
{
  double i, j, k;
  TFlameNodePtr flameNode = flame->GetFlameNode();
  double * molarMass = flame->GetSpecies()->GetMolarMass()->vec;
  double temp = flameNode->temp[kCurr];
  double density = flameNode->mixDensity[kCurr];
  double * Y = flameNode->Y[kCurr];
  double * kSoot = fSootRateCoeffs->vec;

  rho = flameNode->mixDensity[kCurr];
  mu = flameNode->mixViscosity[kCurr];
  Wmix = flameNode->mixMolarMass[kCurr];

  ComputeDimerParticules(flame, moments);
  ComputeSootRateCoeffs(kSoot, temp, flame->GetReaction());
  ComputeCSootStar(kSoot, Y, density, molarMass, flameNode->mixMolarMass[kCurr]);

#ifdef MOMIC
  for (int l = 0; l < fNSootMoments; ++l)
#endif
#ifdef HMOM
  for (int l = 0; l < fNSootMoments-1; ++l)
#endif
  {
    i = double(Geti(l))/6.0;
    j = double(Getj(l))/6.0;
    k = double(Getk(l))/6.0;

    source[l] = 0.0;

    if (fNucleation)
      source[l] += NucleationSource(i, j, k, temp);

    if (fCoagulation)
      source[l] += CoagulationSource(i, j, k, temp, moments);

    if (fCondensation)
      source[l] += CondensationSource(i, j, k, temp, moments);

    if (fSurfaceGrowth)
      source[l] += SurfaceGrowthSource(i, j, k, moments, Y, temp, density, molarMass);

    if (fSurfaceOxidation)
      source[l] += OxidationSource(i, j, k, moments, Y, temp, density, molarMass);

    if (fThermoPhoresis)
      source[l] += SourceThermoPhoresis(i, j, k, moments, flame);
  }

#ifdef HMOM
//Source terms for number of small particles
  int r = (fNSootMoments - 1);

  source[r] = 0.0;

  if (fNucleation)
    source[r] += NucleationSource(0.0, 0.0, 0.0, temp);

  if (fCoagulation)
    source[r] += CoagulationSourceSmall(temp, moments);

  if (fCondensation)
    source[r] += CondensationSourceSmall(temp, moments);

  if (fSurfaceGrowth)
    source[r] += SurfaceGrowthSourceSmall(moments, Y, temp, density, molarMass);

  if (fSurfaceOxidation)
    source[r] += OxidationSourceSmall(moments, Y, temp, density, molarMass);

  if (fThermoPhoresis)
    source[r] += SourceThermoPhoresisSmall(moments, flame);
#endif
}
#endif
#ifdef VS
void T1DSoot::FillSource(T1DFlamePtr flame, double * moments, double * source)
{
  double i, j;
  TFlameNodePtr flameNode = flame->GetFlameNode();
  double * molarMass = flame->GetSpecies()->GetMolarMass()->vec;
  double temp = flameNode->temp[kCurr];
  double density = flameNode->mixDensity[kCurr];
  double * Y = flameNode->Y[kCurr];
  double * kSoot = fSootRateCoeffs->vec;
  double * wSoot = fSootReactionRate->vec;
  double rhodot = flameNode->rhodot[kCurr];

  rho = flameNode->mixDensity[kCurr];
  mu = flameNode->mixViscosity[kCurr];
  Wmix = flameNode->mixMolarMass[kCurr];

  ComputeDimerParticules(flame, moments);
  ComputeSootRateCoeffs(kSoot, temp, flame->GetReaction());
  ComputeCSootStar(kSoot, Y, density, molarMass, flameNode->mixMolarMass[kCurr]);
  ComputeSootReactionRates(wSoot, Y, moments, density, molarMass);

  NodeInfoPtr nodeInfo = flame->GetSolver()->bt->GetNodeInfo();

#ifdef MOMIC
  for (int k = 0; k < fNSootMoments; ++k)
#endif
#ifdef HMOM
  for (int k = 0; k < fNSootMoments-1; ++k)
#endif
  {
    i = double(Geti(k))/6.0;
    j = double(Getj(k))/6.0;

    source[k] = 0.0;
    source[k] -= rhodot * moments[k] / density;

    //if ((moments[1] > 1.0e-9 || temp < 1800.0) && nodeInfo->x[kCurr] < 0.015) { //DiffFlame: 1x is 8,1800,0.015; 2x is 9,1800,0.015; 3x is 10,1800,0.015 (don't solve for log)
    //if (moments[1] > 1.0e-9) { //Good for PremFlame 0.95, don't need anything for 1.15
    //if (moments[1] / moments[0] > (2.0*dimer_nbrC2+1.0)) {
    if (fNucleation)
      source[k] += NucleationSource(i, j, temp);
    if (fCoagulation)
      source[k] += CoagulationSource(i, j, temp, moments);
    if (fCondensation)
      source[k] += CondensationSource(i, j, temp, moments);
    if (fSurfaceGrowth)
      source[k] += SurfaceGrowthSource(i, j, moments, Y, temp, density, molarMass);
    if (fSurfaceOxidation)
    {
      //if (moments[1]/rho > 1.0e-3)
      //source[k] += exp(-1.0e-16/moments[0]) * OxidationSource(i, j, moments, Y, density, molarMass);
      //if (nodeInfo->x[kCurr] > 1.0e-5)
      source[k] += OxidationSource(i, j, moments, Y, temp, density, molarMass);
    }
    if (fThermoPhoresis)
      source[k] += SourceThermoPhoresis(i, j, moments, flame);
    if (fFragmentation)
    {
      //if (moments[0] > 1.0e-20)
      //if (moments[1] > 20.0e-20)
	source[k] += FragmentationSource(i, j, moments);
    }
    if (fSootDiffusion)
      source[k] += DiffusionSource(i, j, moments, flame);//}
  }

  //cerr << "Volume Fraction: " << moments[1] * 24.0/1800.0 << "\t" << moments[1] / moments[0] << endl;
  //cerr << "Location: " << flame->GetSolver()->bt->GetNodeInfo()->x[kCurr] << endl;

#ifdef HMOM
  //Source terms for number of small particles
  int r = (fNSootMoments - 1);

  source[r] = 0.0;
  source[r] -= rhodot * moments[r] / density;

  //Limit as very large negative moments is used...
  //if ((moments[1] > 1.0e-9 || temp < 1800.0) && nodeInfo->x[kCurr] < 0.015) { //DiffFlame: see above
  //if (moments[1] > 1.0e-9) { //Good for PremFlame 0.95, don't need anything for 1.15
  //if (moments[1] / moments[0] > (2.0*dimer_nbrC2+1.0)) {
  if (fNucleation)
    source[r] += NucleationSource(0, 0, temp);
  if (fCoagulation)
    source[r] += CoagulationSourceSmall(temp, moments);
  if (fCondensation)
    source[r] += CondensationSourceSmall(temp, moments);
  if (fSurfaceGrowth)
    source[r] += SurfaceGrowthSourceSmall(moments, Y, temp, density, molarMass);
  if (fSurfaceOxidation)
  {
    //if (moments[1]/rho > 1.0e-3)
    //source[r] += exp(-1.0e-16/moments[0]) * OxidationSourceSmall(moments, Y, density, molarMass);
    //if (nodeInfo->x[kCurr] > 1.0e-5)
    source[r] += OxidationSourceSmall(moments, Y, temp, density, molarMass);
  }
  if (fThermoPhoresis)
    source[r] += SourceThermoPhoresisSmall(moments, flame);
  if (fFragmentation)
  {
    //if (moments[0]/rho > 1.0e-20)
    //if (moments[1]/rho > 20.0e-20)
      source[r] += FragmentationSourceSmall(moments);
  }
  if (fSootDiffusion)
    source[r] += DiffusionSourceSmall(moments, flame);//}

  //General form for 3+1
//   double Mdh, Md10, Md01, Md00;
//   double V0 = 2.0 * dimer_nbrC2;
//   double S0 = pow(V0, 2.0/3.0);
//   double VL = FracMomLarge(1.0, 0.0, moments) / FracMomLarge(0.0, 0.0, moments);
//   double SL = FracMomLarge(0.0, 1.0, moments) / FracMomLarge(0.0, 0.0, moments);

//   double M00 = FracMomLarge(0.0, 0.0, moments);
//   double M10 = FracMomLarge(1.0, 0.0, moments);
//   double M01 = FracMomLarge(0.0, 1.0, moments);

//   double aaa = -1.0;
//   double bbb = -1.0;

//   if (fNucleation)
//   {
//     if (M00 < 1.0e-30 || M10 < 1.0e-30 || M01 < 1.0e-30)
//       source[r] += NucleationSource(0, 0, temp);
//     else
//     {
//       Mdh = NucleationSource(-aaa, -bbb, temp) / (pow(V0, -aaa) * pow(S0, -bbb));
//       Md00 = NucleationSource(0.0, 0.0, temp);
//       Md10 = NucleationSource(1.0, 0.0, temp);
//       Md01 = NucleationSource(0.0, 1.0, temp);

//       source[r] += (Mdh + pow(VL/V0, -aaa) * pow(SL/S0, -bbb) * (aaa*Md10/VL + bbb*Md01/SL - (1.0+aaa+bbb)*Md00)) /
// 	           (1.0 + pow(VL/V0, -aaa) * pow(SL/S0, -bbb) * (aaa*V0  /VL + bbb*  S0/SL - (1.0+aaa+bbb)     ));
//     }
//   }
//   if (fCoagulation)
//   {
//     if (M00 < 1.0e-30 || M10 < 1.0e-30 || M01 < 1.0e-30)
//       source[r] += CoagulationSourceSmall(temp, moments);
//     else
//     {
//       Mdh = CoagulationSource(-aaa, -bbb, temp, moments) / (pow(V0, -aaa) * pow(S0, -bbb));
//       Md00 = CoagulationSource(0.0, 0.0, temp, moments);
//       Md10 = CoagulationSource(1.0, 0.0, temp, moments);
//       Md01 = CoagulationSource(0.0, 1.0, temp, moments);

//       source[r] += (Mdh + pow(VL/V0, -aaa) * pow(SL/S0, -bbb) * (aaa*Md10/VL + bbb*Md01/SL - (1.0+aaa+bbb)*Md00)) /
//                    (1.0 + pow(VL/V0, -aaa) * pow(SL/S0, -bbb) * (aaa*V0  /VL + bbb*  S0/SL - (1.0+aaa+bbb)     ));
//     }
//   }
//   if (fCondensation)
//   {
//     if (M00 < 1.0e-30 || M10 < 1.0e-30 || M01 < 1.0e-30)
//       source[r] += CondensationSourceSmall(temp, moments);
//     else
//     {
//       Mdh = CondensationSource(-aaa, -bbb, temp, moments) / (pow(V0, -aaa) * pow(S0, -bbb));
//       Md00 = CondensationSource(0.0, 0.0, temp, moments);
//       Md10 = CondensationSource(1.0, 0.0, temp, moments);
//       Md01 = CondensationSource(0.0, 1.0, temp, moments);

//       source[r] += (Mdh + pow(VL/V0, -aaa) * pow(SL/S0, -bbb) * (aaa*Md10/VL + bbb*Md01/SL - (1.0+aaa+bbb)*Md00)) /
//                    (1.0 + pow(VL/V0, -aaa) * pow(SL/S0, -bbb) * (aaa*V0  /VL + bbb*  S0/SL - (1.0+aaa+bbb)     ));
//     }
//   }
//   if (fSurfaceGrowth)
//   {
//     if (M00 < 1.0e-30 || M10 < 1.0e-30 || M01 < 1.0e-30)
//       source[r] += SurfaceGrowthSourceSmall(moments, Y, density, molarMass);
//     else
//     {
//       Mdh = SurfaceGrowthSource(-aaa, -bbb, moments, Y, density, molarMass) / (pow(V0, -aaa) * pow(S0, -bbb));
//       Md00 = SurfaceGrowthSource(0.0, 0.0, moments, Y, density, molarMass);
//       Md10 = SurfaceGrowthSource(1.0, 0.0, moments, Y, density, molarMass);
//       Md01 = SurfaceGrowthSource(0.0, 1.0, moments, Y, density, molarMass);

//       source[r] += (Mdh + pow(VL/V0, -aaa) * pow(SL/S0, -bbb) * (aaa*Md10/VL + bbb*Md01/SL - (1.0+aaa+bbb)*Md00)) /
//                    (1.0 + pow(VL/V0, -aaa) * pow(SL/S0, -bbb) * (aaa*V0  /VL + bbb*  S0/SL - (1.0+aaa+bbb)     ));
//     }
//   }
//   if (fSurfaceOxidation)
//   {
//     if (M00 < 1.0e-30 || M10 < 1.0e-30 || M01 < 1.0e-30)
//       source[r] += OxidationSourceSmall(moments, Y, density, molarMass);
//     else
//     {
//       Mdh = OxidationSource(-aaa, -bbb, moments, Y, density, molarMass) / (pow(V0, -aaa) * pow(S0, -bbb));
//       Md00 = OxidationSource(0.0, 0.0, moments, Y, density, molarMass);
//       Md10 = OxidationSource(1.0, 0.0, moments, Y, density, molarMass);
//       Md01 = OxidationSource(0.0, 1.0, moments, Y, density, molarMass);

//       source[r] += (Mdh + pow(VL/V0, -aaa) * pow(SL/S0, -bbb) * (aaa*Md10/VL + bbb*Md01/SL - (1.0+aaa+bbb)*Md00)) /
//                    (1.0 + pow(VL/V0, -aaa) * pow(SL/S0, -bbb) * (aaa*V0  /VL + bbb*  S0/SL - (1.0+aaa+bbb)     ));
//     }
//   }
//   if (fThermoPhoresis)
//     source[r] += SourceThermoPhoresisSmall(moments, flame);

  //General form for 6+1
//   double Mdh, Md20, Md11, Md02, Md10, Md01, Md00;
//   double V0 = 2.0 * dimer_nbrC2;
//   double S0 = pow(V0, 2.0/3.0);
//   double VL = pow(FracMomLarge(1.0, 0.0, moments), 2.0) * pow(FracMomLarge(0.0, 0.0, moments), -3.0/2.0) * pow(FracMomLarge(2.0, 0.0, moments), -1.0/2.0);
//   double SL = pow(FracMomLarge(0.0, 1.0, moments), 2.0) * pow(FracMomLarge(0.0, 0.0, moments), -3.0/2.0) * pow(FracMomLarge(0.0, 2.0, moments), -1.0/2.0);
//   double SVV = pow(FracMomLarge(0.0, 0.0, moments), 1.0/2.0) * pow(FracMomLarge(2.0, 0.0, moments), 1.0/2.0) / FracMomLarge(1.0, 0.0, moments);
//   double SVS = FracMomLarge(0.0, 0.0, moments) * FracMomLarge(1.0, 1.0, moments) / (FracMomLarge(1.0, 0.0, moments) * FracMomLarge(0.0, 1.0, moments));
//   double SSS = pow(FracMomLarge(0.0, 0.0, moments), 1.0/2.0) * pow(FracMomLarge(0.0, 2.0, moments), 1.0/2.0) / FracMomLarge(0.0, 1.0, moments);

//   double M00 = FracMomLarge(0.0, 0.0, moments) / FracMom(0.0, 0.0, moments);
//   double M10 = FracMomLarge(1.0, 0.0, moments);
//   double M01 = FracMomLarge(0.0, 1.0, moments);
//   double M20 = FracMomLarge(2.0, 0.0, moments);
//   double M11 = FracMomLarge(1.0, 1.0, moments);
//   double M02 = FracMomLarge(0.0, 2.0, moments);

//   double aaa = -1.0;
//   double bbb = -2.0;

//   if (fNucleation)
//   {
//     //if (M00 < 1.0e-20 || M10 < 1.0e-20 || M01 < 1.0e-20 || M20 < 1.0e-20 || M11 < 1.0e-20 || M02 < 1.0e-20
//     if (M00 < 1.0e-5)
//       source[r] += NucleationSource(0, 0, temp);
//     else
//     {
//       Mdh = NucleationSource(-aaa, -bbb, temp) / (pow(V0, -aaa) * pow(S0, -bbb));
//       Md00 = NucleationSource(0.0, 0.0, temp);
//       Md10 = NucleationSource(1.0, 0.0, temp);
//       Md01 = NucleationSource(0.0, 1.0, temp);
//       Md20 = NucleationSource(2.0, 0.0, temp);
//       Md11 = NucleationSource(1.0, 1.0, temp);
//       Md02 = NucleationSource(0.0, 2.0, temp);

//       source[r] += (Mdh + pow(VL/V0, -aaa) * pow(SL/S0, -bbb) * pow(SVV, aaa*aaa) * pow(SVS, aaa*bbb) * pow(SSS, bbb*bbb) * 
// 		    (-(1.0/2.0)*aaa*(aaa+1.0)*Md20 /(VL*VL*SVV*SVV*SVV*SVV) - aaa*bbb*Md11 /(VL*SL*SVV*SVS*SSS) - (1.0/2.0)*bbb*(bbb+1.0)*Md02 /(SL*SL*SSS*SSS*SSS*SSS) +
// 		     aaa*(aaa+bbb+2.0)*Md10/(VL*SVV) + bbb*(aaa+bbb+2.0)*Md01/(SL*SSS) - (1.0/2.0)*(aaa+bbb+1.0)*(aaa+bbb+2.0)*Md00)) /
//                    (1.0 + pow(VL/V0, -aaa) * pow(SL/S0, -bbb) * pow(SVV, aaa*aaa) * pow(SVS, aaa*bbb) * pow(SSS, bbb*bbb) * 
// 		    (-(1.0/2.0)*aaa*(aaa+1.0)*V0*V0/(VL*VL*SVV*SVV*SVV*SVV) - aaa*bbb*V0*S0/(VL*SL*SVV*SVS*SSS) - (1.0/2.0)*bbb*(bbb+1.0)*S0*S0/(SL*SL*SSS*SSS*SSS*SSS) +
// 		     aaa*(aaa+bbb+2.0)*V0  /(VL*SVV) + bbb*(aaa+bbb+2.0)*S0  /(SL*SSS) - 1.0/2.0*(aaa+bbb+1.0)*(aaa+bbb+2.0)     ));
//     }
//   }
//   if (fCoagulation)
//   {
//     //if (M00 < 1.0e-20 || M10 < 1.0e-20 || M01 < 1.0e-20 || M20 < 1.0e-20 || M11 < 1.0e-20 || M02 < 1.0e-20)
//     if (M00 < 1.0e-5)
//       source[r] += CoagulationSourceSmall(temp, moments);
//     else
//     {
//       Mdh = CoagulationSource(-aaa, -bbb, temp, moments) / (pow(V0, -aaa) * pow(S0, -bbb));
//       Md00 = CoagulationSource(0.0, 0.0, temp, moments);
//       Md10 = CoagulationSource(1.0, 0.0, temp, moments);
//       Md01 = CoagulationSource(0.0, 1.0, temp, moments);
//       Md20 = CoagulationSource(2.0, 0.0, temp, moments);
//       Md11 = CoagulationSource(1.0, 1.0, temp, moments);
//       Md02 = CoagulationSource(0.0, 2.0, temp, moments);

//       source[r] += (Mdh + pow(VL/V0, -aaa) * pow(SL/S0, -bbb) * pow(SVV, aaa*aaa) * pow(SVS, aaa*bbb) * pow(SSS, bbb*bbb) * 
// 		    (-1.0/2.0*aaa*(aaa+1.0)*Md20 /(VL*VL*SVV*SVV*SVV*SVV) - aaa*bbb*Md11 /(VL*SL*SVV*SVS*SSS) - 1.0/2.0*bbb*(bbb+1.0)*Md02 /(SL*SL*SSS*SSS*SSS*SSS) +
// 		     aaa*(aaa+bbb+2.0)*Md10/(VL*SVV) + bbb*(aaa+bbb+2.0)*Md01/(SL*SSS) - 1.0/2.0*(aaa+bbb+1.0)*(aaa+bbb+2.0)*Md00)) /
//                    (1.0 + pow(VL/V0, -aaa) * pow(SL/S0, -bbb) * pow(SVV, aaa*aaa) * pow(SVS, aaa*bbb) * pow(SSS, bbb*bbb) * 
// 		    (-1.0/2.0*aaa*(aaa+1.0)*V0*V0/(VL*VL*SVV*SVV*SVV*SVV) - aaa*bbb*V0*S0/(VL*SL*SVV*SVS*SSS) - 1.0/2.0*bbb*(bbb+1.0)*S0*S0/(SL*SL*SSS*SSS*SSS*SSS) +
// 		     aaa*(aaa+bbb+2.0)*V0  /(VL*SVV) + bbb*(aaa+bbb+2.0)*S0  /(SL*SSS) - 1.0/2.0*(aaa+bbb+1.0)*(aaa+bbb+2.0)     ));
//     }
//   }
//   if (fCondensation)
//   {
//     //if (M00 < 1.0e-20 || M10 < 1.0e-20 || M01 < 1.0e-20 || M20 < 1.0e-20 || M11 < 1.0e-20 || M02 < 1.0e-20)
//     if (M00 < 1.0e-5)
//       source[r] += CondensationSourceSmall(temp, moments);
//     else
//     {
//       Mdh = CondensationSource(-aaa, -bbb, temp, moments) / (pow(V0, -aaa) * pow(S0, -bbb));
//       Md00 = CondensationSource(0.0, 0.0, temp, moments);
//       Md10 = CondensationSource(1.0, 0.0, temp, moments);
//       Md01 = CondensationSource(0.0, 1.0, temp, moments);
//       Md20 = CondensationSource(2.0, 0.0, temp, moments);
//       Md11 = CondensationSource(1.0, 1.0, temp, moments);
//       Md02 = CondensationSource(0.0, 2.0, temp, moments);

//       source[r] += (Mdh + pow(VL/V0, -aaa) * pow(SL/S0, -bbb) * pow(SVV, aaa*aaa) * pow(SVS, aaa*bbb) * pow(SSS, bbb*bbb) * 
// 		    (-(1.0/2.0)*aaa*(aaa+1.0)*Md20 /(VL*VL*SVV*SVV*SVV*SVV) - aaa*bbb*Md11 /(VL*SL*SVV*SVS*SSS) - (1.0/2.0)*bbb*(bbb+1.0)*Md02 /(SL*SL*SSS*SSS*SSS*SSS) +
// 		     aaa*(aaa+bbb+2.0)*Md10/(VL*SVV) + bbb*(aaa+bbb+2.0)*Md01/(SL*SSS) - (1.0/2.0)*(aaa+bbb+1.0)*(aaa+bbb+2.0)*Md00)) /
//                    (1.0 + pow(VL/V0, -aaa) * pow(SL/S0, -bbb) * pow(SVV, aaa*aaa) * pow(SVS, aaa*bbb) * pow(SSS, bbb*bbb) * 
// 		    (-(1.0/2.0)*aaa*(aaa+1.0)*V0*V0/(VL*VL*SVV*SVV*SVV*SVV) - aaa*bbb*V0*S0/(VL*SL*SVV*SVS*SSS) - (1.0/2.0)*bbb*(bbb+1.0)*S0*S0/(SL*SL*SSS*SSS*SSS*SSS) +
// 		     aaa*(aaa+bbb+2.0)*V0  /(VL*SVV) + bbb*(aaa+bbb+2.0)*S0  /(SL*SSS) - 1.0/2.0*(aaa+bbb+1.0)*(aaa+bbb+2.0)     ));
//     }
//   }
//   if (fSurfaceGrowth)
//   {
//     //if (M00 < 1.0e-20 || M10 < 1.0e-20 || M01 < 1.0e-20 || M20 < 1.0e-20 || M11 < 1.0e-20 || M02 < 1.0e-20)
//     if (M00 < 1.0e-5)
//       source[r] += SurfaceGrowthSourceSmall(moments, Y, density, molarMass);
//     else
//     {
//       Mdh = SurfaceGrowthSource(-aaa, -bbb, moments, Y, density, molarMass) / (pow(V0, -aaa) * pow(S0, -bbb));
//       Md00 = SurfaceGrowthSource(0.0, 0.0, moments, Y, density, molarMass);
//       Md10 = SurfaceGrowthSource(1.0, 0.0, moments, Y, density, molarMass);
//       Md01 = SurfaceGrowthSource(0.0, 1.0, moments, Y, density, molarMass);
//       Md20 = SurfaceGrowthSource(2.0, 0.0, moments, Y, density, molarMass);
//       Md11 = SurfaceGrowthSource(1.0, 1.0, moments, Y, density, molarMass);
//       Md02 = SurfaceGrowthSource(0.0, 2.0, moments, Y, density, molarMass);

//       source[r] += (Mdh + pow(VL/V0, -aaa) * pow(SL/S0, -bbb) * pow(SVV, aaa*aaa) * pow(SVS, aaa*bbb) * pow(SSS, bbb*bbb) * 
// 		    (-1.0/2.0*aaa*(aaa+1.0)*Md20 /(VL*VL*SVV*SVV*SVV*SVV) - aaa*bbb*Md11 /(VL*SL*SVV*SVS*SSS) - 1.0/2.0*bbb*(bbb+1.0)*Md02 /(SL*SL*SSS*SSS*SSS*SSS) +
// 		     aaa*(aaa+bbb+2.0)*Md10/(VL*SVV) + bbb*(aaa+bbb+2.0)*Md01/(SL*SSS) - 1.0/2.0*(aaa+bbb+1.0)*(aaa+bbb+2.0)*Md00)) /
//                    (1.0 + pow(VL/V0, -aaa) * pow(SL/S0, -bbb) * pow(SVV, aaa*aaa) * pow(SVS, aaa*bbb) * pow(SSS, bbb*bbb) * 
// 		    (-1.0/2.0*aaa*(aaa+1.0)*V0*V0/(VL*VL*SVV*SVV*SVV*SVV) - aaa*bbb*V0*S0/(VL*SL*SVV*SVS*SSS) - 1.0/2.0*bbb*(bbb+1.0)*S0*S0/(SL*SL*SSS*SSS*SSS*SSS) +
// 		     aaa*(aaa+bbb+2.0)*V0  /(VL*SVV) + bbb*(aaa+bbb+2.0)*S0  /(SL*SSS) - 1.0/2.0*(aaa+bbb+1.0)*(aaa+bbb+2.0)     ));
//     }
//   }
//   if (fSurfaceOxidation)
//   {
//     //if (M00 < 1.0e-20 || M10 < 1.0e-20 || M01 < 1.0e-20 || M20 < 1.0e-20 || M11 < 1.0e-20 || M02 < 1.0e-20)
//     if (M00 < 1.0e-5)
//       source[r] += OxidationSourceSmall(moments, Y, density, molarMass);
//     else
//     {
//       Mdh = OxidationSource(-aaa, -bbb, moments, Y, density, molarMass) / (pow(V0, -aaa) * pow(S0, -bbb));
//       Md00 = OxidationSource(0.0, 0.0, moments, Y, density, molarMass);
//       Md10 = OxidationSource(1.0, 0.0, moments, Y, density, molarMass);
//       Md01 = OxidationSource(0.0, 1.0, moments, Y, density, molarMass);
//       Md20 = OxidationSource(2.0, 0.0, moments, Y, density, molarMass);
//       Md11 = OxidationSource(1.0, 1.0, moments, Y, density, molarMass);
//       Md02 = OxidationSource(0.0, 2.0, moments, Y, density, molarMass);

//       source[r] += (Mdh + pow(VL/V0, -aaa) * pow(SL/S0, -bbb) * pow(SVV, aaa*aaa) * pow(SVS, aaa*bbb) * pow(SSS, bbb*bbb) * 
// 		    (-1.0/2.0*aaa*(aaa+1.0)*Md20 /(VL*VL*SVV*SVV*SVV*SVV) - aaa*bbb*Md11 /(VL*SL*SVV*SVS*SSS) - 1.0/2.0*bbb*(bbb+1.0)*Md02 /(SL*SL*SSS*SSS*SSS*SSS) +
// 		     aaa*(aaa+bbb+2.0)*Md10/(VL*SVV) + bbb*(aaa+bbb+2.0)*Md01/(SL*SSS) - 1.0/2.0*(aaa+bbb+1.0)*(aaa+bbb+2.0)*Md00)) /
//                    (1.0 + pow(VL/V0, -aaa) * pow(SL/S0, -bbb) * pow(SVV, aaa*aaa) * pow(SVS, aaa*bbb) * pow(SSS, bbb*bbb) * 
// 		    (-1.0/2.0*aaa*(aaa+1.0)*V0*V0/(VL*VL*SVV*SVV*SVV*SVV) - aaa*bbb*V0*S0/(VL*SL*SVV*SVS*SSS) - 1.0/2.0*bbb*(bbb+1.0)*S0*S0/(SL*SL*SSS*SSS*SSS*SSS) +
// 		     aaa*(aaa+bbb+2.0)*V0  /(VL*SVV) + bbb*(aaa+bbb+2.0)*S0  /(SL*SSS) - 1.0/2.0*(aaa+bbb+1.0)*(aaa+bbb+2.0)     ));
//     }
//   }
//   if (fThermoPhoresis)
//     source[r] += SourceThermoPhoresisSmall(moments, flame);
#endif
}
#endif
#ifdef V
void T1DSoot::FillSource(T1DFlamePtr flame, double * moments, double * source)
{
  double i;
  TFlameNodePtr flameNode = flame->GetFlameNode();
  double * molarMass = flame->GetSpecies()->GetMolarMass()->vec;
  double temp = flameNode->temp[kCurr];
  double density = flameNode->mixDensity[kCurr];
  double * Y = flameNode->Y[kCurr];
  double * kSoot = fSootRateCoeffs->vec;

  rho = flameNode->mixDensity[kCurr];
  mu = flameNode->mixViscosity[kCurr];
  Wmix = flameNode->mixMolarMass[kCurr];

  ComputeDimerParticules(flame, moments);
  ComputeSootRateCoeffs(kSoot, temp, flame->GetReaction());
  ComputeCSootStar(kSoot, Y, density, molarMass, flameNode->mixMolarMass[kCurr]);

#ifdef LINDSTEDT
  source[0] = 0.0;

  if (fNucleation)
    source[0] += (1.0/50.0) * LindstedtNucleation(temp, Y, density, molarMass);
  if (fCoagulation)
    source[0] += LindstedtCoagulation(moments, temp);
  if (fThermoPhoresis)
    source[0] += SourceThermoPhoresis(0, moments, flame);

  source[1] = 0.0;

  if (fNucleation)
    source[1] += LindstedtNucleation(temp, Y, density, molarMass);
  if (fSurfaceGrowth)
    source[1] += LindstedtGrowth(moments, temp, Y, density, molarMass);
  if (fSurfaceOxidation)
    source[1] += LindstedtOxidation(moments, temp, Y, density, molarMass);
  if (fThermoPhoresis)
    source[1] += SourceThermoPhoresis(1, moments, flame);

#else
#ifdef MOMIC
  for (int k = 0; k < fNSootMoments; ++k)
#endif
#ifdef FRENKLACH
  for (int k = 0; k < fNSootMoments-1; ++k)
#endif
#ifdef HMOM
  for (int k = 0; k < fNSootMoments-1; ++k)
#endif
  {
    i = double(Geti(k))/6.0;

    source[k] = 0.0;

    if (fNucleation)
      source[k] += NucleationSource(i, temp);

    if (fCoagulation)
      source[k] += CoagulationSource(i, temp, moments);

    if (fCondensation)
      source[k] += CondensationSource(i, temp, moments);

    if (fSurfaceGrowth)
      source[k] += SurfaceGrowthSource(i, moments, Y, temp, density, molarMass);

    if (fSurfaceOxidation)
      source[k] += OxidationSource(i, moments, Y, temp, density, molarMass);

    if (fThermoPhoresis)
      source[k] += SourceThermoPhoresis(i, moments, flame);
  }

#ifdef FRENKLACH
  //Source terms for number of small particles
  int r = (fNSootMoments - 1);

  source[r] = 0.0;

  if (fNucleation)
    source[r] += NucleationSource(0, temp);
  
  if (fCoagulation)
    source[r] += CoagulationSourceSmall(temp, moments);
  
  if (fCondensation)
    source[r] += CondensationSourceSmall(temp, moments);
  
  if (fSurfaceGrowth)
    source[r] += SurfaceGrowthSourceSmall(moments, Y, temp, density, molarMass);
  
  if (fSurfaceOxidation)
    source[r] += OxidationSourceSmall(moments, Y, temp, density, molarMass);

  if (fThermoPhoresis)
    source[r] += SourceThermoPhoresisSmall(moments, flame);
#endif

#ifdef HMOM
//Source terms for number of small particles
  int r = (fNSootMoments - 1);

  //double V0 = 2.0 * dimer_nbrC2;
  //double VL = pow(FracMomLarge(0.0, moments), -3.0/2.0) * pow(FracMomLarge(1.0, moments), 2.0) * pow(FracMomLarge(2.0, moments), -1.0/2.0);
  //double SL = FracMomLarge(0.0, moments) * pow(FracMomLarge(1.0, moments), -2.0) * FracMomLarge(2.0, moments);
  //double denom = V0*V0 - 2.0*V0*VL*SL*SL*SL + VL*VL*SL*SL*SL*SL;

  source[r] = 0.0;

  //double ML0 = moments[0] - moments[3];

  if (fNucleation)
    source[r] += NucleationSource(0, temp);
  if (fCoagulation)
  {
    //if (ML0 < 1.0e-25)
    source[r] += CoagulationSourceSmall(temp, moments);
    //else
    //source[r] += (CoagulationSource(2.0, temp, moments) - 2.0*VL*SL*SL*SL*CoagulationSource(1.0, temp, moments) + VL*VL*SL*SL*SL*SL*CoagulationSource(0.0, temp, moments)) / denom;
  }
  if (fCondensation)
  {
    //if (ML0 < 1.0e-25)
    source[r] += CondensationSourceSmall(temp, moments);
    //else
    //source[r] += (CondensationSource(2.0, temp, moments) - 2.0*VL*SL*SL*SL*CondensationSource(1.0, temp, moments) + VL*VL*SL*SL*SL*SL*CondensationSource(0.0, temp, moments)) / denom;
  }
  if (fSurfaceGrowth)
  {
    //if (ML0 < 1.0e-25)
    source[r] += SurfaceGrowthSourceSmall(moments, Y, temp, density, molarMass);
    //else
    //source[r] += (SurfaceGrowthSource(2.0, moments, Y, density, molarMass) - 2.0*VL*SL*SL*SL*SurfaceGrowthSource(1.0, moments, Y, density, molarMass) + VL*VL*SL*SL*SL*SL*SurfaceGrowthSource(0.0, moments, Y, density, molarMass)) / denom;
  }
  if (fSurfaceOxidation)
    source[r] += OxidationSourceSmall(moments, Y, temp, density, molarMass);
    //else
    //source[r] += (OxidationSource(2.0, moments, Y, density, molarMass) - 2.0*VL*SL*SL*SL*OxidationSource(1.0, moments, Y, density, molarMass) + VL*VL*SL*SL*SL*SL*OxidationSource(0.0, moments, Y, density, molarMass)) / denom;
  if (fThermoPhoresis)
    source[r] += SourceThermoPhoresisSmall(moments, flame);
#endif
#endif
}
#endif

double T1DSoot::SourceThermoPhoresis(double i, double j, double k, double * moments, T1DFlamePtr flame)
{
  TNewtonPtr bt = flame->GetSolver()->bt;
  TFlameNodePtr	flameNode = flame->GetFlameNode();
  NodeInfoPtr nodeInfo = bt->GetNodeInfo();

  // accommodation coefficient following
  // R. J. Santoro: The Transport and Growth
  // of Soot Particles ...
  // Comb. Sci. Tech., 1987, Vol. 53 p. 89
  static const double A = 0.9;
  static const double fact = 3.0 / ( 4.0 * ( 1.0+PI*A/8.0 ) );

  double * mu = flameNode->mixViscosity;
  double * rho = flameNode->mixDensity;
  double * temp = flameNode->temp;
  
  double d2 = nodeInfo->hm * fact * ( mu[kCurr]/temp[kCurr] + mu[kNext]/temp[kNext] );
  double d1 = nodeInfo->h  * fact * ( mu[kCurr]/temp[kCurr] + mu[kPrev]/temp[kPrev] );
  
  double src = ( d2*(temp[kNext]-temp[kCurr]) - d1*(temp[kCurr]-temp[kPrev]) ) / nodeInfo->hnenn;
  src *= FracMom(i, j, k, moments) / rho[kCurr];
  
  return src;
}

double T1DSoot::SourceThermoPhoresis(double i, double j, double * moments, T1DFlamePtr flame)
{
  TNewtonPtr bt = flame->GetSolver()->bt;
  TFlameNodePtr	flameNode = flame->GetFlameNode();
  NodeInfoPtr nodeInfo = bt->GetNodeInfo();

  // accommodation coefficient following
  // R. J. Santoro: The Transport and Growth
  // of Soot Particles ...
  // Comb. Sci. Tech., 1987, Vol. 53 p. 89
  static const double A = 0.9;
  static const double fact = 3.0 / ( 4.0 * ( 1.0+PI*A/8.0 ) );

  double * mu = flameNode->mixViscosity;
  double * rho = flameNode->mixDensity;
  double * temp = flameNode->temp;
  
  double d2 = nodeInfo->hm * fact * ( mu[kCurr]/temp[kCurr] + mu[kNext]/temp[kNext] );
  double d1 = nodeInfo->h  * fact * ( mu[kCurr]/temp[kCurr] + mu[kPrev]/temp[kPrev] );
  
  double src = ( d2*(temp[kNext]-temp[kCurr]) - d1*(temp[kCurr]-temp[kPrev]) ) / nodeInfo->hnenn;
  src *= FracMom(i, j, moments) / rho[kCurr];
  
  return src;
}

double T1DSoot::SourceThermoPhoresis(double i, double * moments, T1DFlamePtr flame)
{
  TNewtonPtr bt = flame->GetSolver()->bt;
  TFlameNodePtr	flameNode = flame->GetFlameNode();
  NodeInfoPtr nodeInfo = bt->GetNodeInfo();

  // accommodation coefficient following
  // R. J. Santoro: The Transport and Growth
  // of Soot Particles ...
  // Comb. Sci. Tech., 1987, Vol. 53 p. 89
  static const double A = 0.9;
  static const double fact = 3.0 / ( 4.0 * ( 1.0+PI*A/8.0 ) );

  double * mu = flameNode->mixViscosity;
  double * rho = flameNode->mixDensity;
  double * temp = flameNode->temp;
  
  double d2 = nodeInfo->hm * fact * ( mu[kCurr]/temp[kCurr] + mu[kNext]/temp[kNext] );
  double d1 = nodeInfo->h  * fact * ( mu[kCurr]/temp[kCurr] + mu[kPrev]/temp[kPrev] );
  
  double src = ( d2*(temp[kNext]-temp[kCurr]) - d1*(temp[kCurr]-temp[kPrev]) ) / nodeInfo->hnenn;
  src *= FracMom(i, moments) / rho[kCurr];
  
  return src;
}

double T1DSoot::SourceThermoPhoresisSmall(double * moments, T1DFlamePtr flame)
{
  TNewtonPtr bt = flame->GetSolver()->bt;
  TFlameNodePtr	flameNode = flame->GetFlameNode();
  NodeInfoPtr nodeInfo = bt->GetNodeInfo();

  // accommodation coefficient following
  // R. J. Santoro: The Transport and Growth
  // of Soot Particles ...
  // Comb. Sci. Tech., 1987, Vol. 53 p. 89
  static const double A = 0.9;
  static const double fact = 3.0 / ( 4.0 * ( 1.0+PI*A/8.0 ) );

  double * mu = flameNode->mixViscosity;
  double * rho = flameNode->mixDensity;
  double * temp = flameNode->temp;
  
  double d2 = nodeInfo->hm * fact * ( mu[kCurr]/temp[kCurr] + mu[kNext]/temp[kNext] );
  double d1 = nodeInfo->h  * fact * ( mu[kCurr]/temp[kCurr] + mu[kPrev]/temp[kPrev] );
  
  double src = ( d2*(temp[kNext]-temp[kCurr]) - d1*(temp[kCurr]-temp[kPrev]) ) / nodeInfo->hnenn;
  src *= moments[fNSootMoments-1] / rho[kCurr];
  
  return src;
}

#ifdef VSH
void TSoot::SootSetDefault(double * yIn)
{
  double i, j, k;

#ifdef MOMIC
  for (int l = 0; l < fNSootMoments; ++l)
#endif
#ifdef HMOM		
  for (int l = 0; l < fNSootMoments-1; ++l)
#endif
  {
    i = double(Geti(l)/6.0);
    j = double(Getj(l)/6.0);
    k = double(Getk(l)/6.0);
		
    if (l == 0)
      yIn[fOffsetMoments] = 1.0e-20;
    else
      yIn[l+fOffsetMoments] = pow(2.0*nucl_nbrC2, i+j*(2.0/3.0)) * pow(2.0*nucl_nbrH, k) * yIn[fOffsetMoments];

    ySol[l] = yIn[l+fOffsetMoments];
  }
#ifdef HMOM
  yIn[fOffsetMoments+fNSootMoments-1] = 1.0e-20;
  ySol[fNSootMoments-1] = yIn[fOffsetMoments+fNSootMoments-1];
#endif
}
#endif
#ifdef VS
void TSoot::SootSetDefault(double * yIn)
{
  double i, j;

#ifdef MOMIC
  for (int k = 0; k < fNSootMoments; ++k)
#endif
#ifdef HMOM		
  for (int k = 0; k < fNSootMoments-1; ++k)
#endif
  {
    i = double(Geti(k)/6.0);
    j = double(Getj(k)/6.0);
		
    if (k == 0)
    {
      #ifdef SOLVELOG
      yIn[k] = log(1.0e-20);
      #else
      yIn[k] = 1.0e-20;
      #endif
    }
    else
    {
      #ifdef SOLVELOG
      yIn[k] = log(pow(2.0*nucl_nbrC2, i+j*(2.0/3.0))) + yIn[0];
      #else
      yIn[k] = pow(2.0*nucl_nbrC2, i+j*(2.0/3.0)) * yIn[0];
      #endif
    }

    ySol[k] = yIn[k];
  }
#ifdef HMOM
  #ifdef SOLVELOG
  yIn[fNSootMoments-1] = log(1.0e-20);
  #else
  yIn[fNSootMoments-1] = 1.0e-20;
  #endif
  ySol[fNSootMoments-1] = yIn[fNSootMoments-1];
#endif

#ifdef FRAGTEST
  #ifdef SOLVELOG
//     #ifdef PHI_0_85
//       yIn[0] = ySol[0] = log(2.0804e-13 / 1.043797);
//       yIn[1] = ySol[1] = log(1.4667e-05 / 1.043797);
//       yIn[2] = ySol[2] = log(3.1108e-07 / 1.043797);
//       yIn[3] = ySol[3] = log(1.6877e+04 / 1.043797);
//       yIn[4] = ySol[4] = log(2.6058e+02 / 1.043797);
//       yIn[5] = ySol[5] = log(4.1942e+00 / 1.043797);
//       yIn[6] = ySol[6] = log(8.4762e-16 / 1.043797);
//     #endif
//     #ifdef PHI_0_95
//       yIn[0] = ySol[0] = log(2.1810e-13 / 1.035113);
//       yIn[1] = ySol[1] = log(1.5376e-05 / 1.035113);
//       yIn[2] = ySol[2] = log(3.2613e-07 / 1.035113);
//       yIn[3] = ySol[3] = log(1.7694e+04 / 1.035113);
//       yIn[4] = ySol[4] = log(2.7319e+02 / 1.035113);
//       yIn[5] = ySol[5] = log(4.3971e+00 / 1.035113);
//       yIn[6] = ySol[6] = log(8.8862e-16 / 1.035113);
//     #endif
//     #ifdef PHI_1_05
//       yIn[0] = ySol[0] = log(2.2552e-13 / 1.028150);
//       yIn[1] = ySol[1] = log(1.5899e-05 / 1.028150);
//       yIn[2] = ySol[2] = log(3.3722e-07 / 1.028150);
//       yIn[3] = ySol[3] = log(1.8295e+04 / 1.028150);
//       yIn[4] = ySol[4] = log(2.8247e+02 / 1.028150);
//       yIn[5] = ySol[5] = log(4.5466e+00 / 1.028150);
//       yIn[6] = ySol[6] = log(9.1884e-16 / 1.028150);
//     #endif
//     #ifdef PHI_1_10
//       yIn[0] = ySol[0] = log(2.3787e-13 / 1.020561);
//       yIn[1] = ySol[1] = log(1.6770e-05 / 1.020561);
//       yIn[2] = ySol[2] = log(3.5569e-07 / 1.020561);
//       yIn[3] = ySol[3] = log(1.9298e+04 / 1.020561);
//       yIn[4] = ySol[4] = log(2.9795e+02 / 1.020561);
//       yIn[5] = ySol[5] = log(4.7956e+00 / 1.020561);
//       yIn[6] = ySol[6] = log(9.6916e-16 / 1.020561);
//     #endif
//     #ifdef PHI_1_15
//       yIn[0] = ySol[0] = log(2.7136e-13 / 1.002392);
//       yIn[1] = ySol[1] = log(1.9131e-05 / 1.002392);
//       yIn[2] = ySol[2] = log(4.0577e-07 / 1.002392);
//       yIn[3] = ySol[3] = log(2.2014e+04 / 1.002392);
//       yIn[4] = ySol[4] = log(3.3989e+02 / 1.002392);
//       yIn[5] = ySol[5] = log(5.4708e+00 / 1.002392);
//       yIn[6] = ySol[6] = log(1.1056e-15 / 1.002392);
//     #endif
    #ifdef PHI_0_85
      yIn[0] = ySol[0] = log(1.6354e-12 / 1.050789);
      yIn[1] = ySol[1] = log(5.4311e-05 / 1.050789);
      yIn[2] = ySol[2] = log(5.8471e-07 / 1.050789);
      yIn[3] = ySol[3] = log(3.3149e+03 / 1.050789);
      yIn[4] = ySol[4] = log(3.4201e+01 / 1.050789);
      yIn[5] = ySol[5] = log(3.5481e-01 / 1.050789);
      yIn[6] = ySol[6] = log(2.2160e-13 / 1.050789);
    #endif
    #ifdef PHI_0_95
      yIn[0] = ySol[0] = log(1.7151e-12 / 1.042384);
      yIn[1] = ySol[1] = log(5.6957e-05 / 1.042384);
      yIn[2] = ySol[2] = log(6.1319e-07 / 1.042384);
      yIn[3] = ySol[3] = log(3.4763e+03 / 1.042384);
      yIn[4] = ySol[4] = log(3.5867e+01 / 1.042384);
      yIn[5] = ySol[5] = log(3.7209e-01 / 1.042384);
      yIn[6] = ySol[6] = log(2.3240e-13 / 1.042384);
    #endif
    #ifdef PHI_1_05
      yIn[0] = ySol[0] = log(1.7738e-12 / 1.035620);
      yIn[1] = ySol[1] = log(5.8907e-05 / 1.035620);
      yIn[2] = ySol[2] = log(6.3419e-07 / 1.035620);
      yIn[3] = ySol[3] = log(3.5954e+03 / 1.035620);
      yIn[4] = ySol[4] = log(3.7095e+01 / 1.035620);
      yIn[5] = ySol[5] = log(3.8484e-01 / 1.035620);
      yIn[6] = ySol[6] = log(2.4035e-13 / 1.035620);
    #endif
    #ifdef PHI_1_10
      yIn[0] = ySol[0] = log(1.8717e-12 / 1.028385);
      yIn[1] = ySol[1] = log(6.2158e-05 / 1.028385);
      yIn[2] = ySol[2] = log(6.6919e-07 / 1.028385);
      yIn[3] = ySol[3] = log(3.7938e+03 / 1.028385);
      yIn[4] = ySol[4] = log(3.9142e+01 / 1.028385);
      yIn[5] = ySol[5] = log(4.0608e-01 / 1.028385);
      yIn[6] = ySol[6] = log(2.5362e-13 / 1.028385);
    #endif
    #ifdef PHI_1_15
      yIn[0] = ySol[0] = log(2.1376e-12 / 1.011168);
      yIn[1] = ySol[1] = log(7.0986e-05 / 1.011168);
      yIn[2] = ySol[2] = log(7.6423e-07 / 1.011168);
      yIn[3] = ySol[3] = log(4.3326e+03 / 1.011168);
      yIn[4] = ySol[4] = log(4.4702e+01 / 1.011168);
      yIn[5] = ySol[5] = log(4.6375e-01 / 1.011168);
      yIn[6] = ySol[6] = log(2.8964e-13 / 1.011168);
    #endif
  #else
  cerr << "FragTest only works when solving for log...\n";
  #endif
#endif
}
#endif
#ifdef V
void TSoot::SootSetDefault(double * yIn)
{
  double i;

#ifdef MOMIC
  for (int k = 0; k < fNSootMoments; ++k)
#endif
#ifdef FRENKLACH
  for (int k = 0; k < fNSootMoments-1; ++k)
#endif
#ifdef HMOM		
  for (int k = 0; k < fNSootMoments-1; ++k)
#endif
  {
    i = double(Geti(k)/6.0);
		
    if (k==0)
      yIn[fOffsetMoments] = 1.0e-20;
    else
      yIn[k+fOffsetMoments] = pow(2.0*nucl_nbrC2, i) * yIn[fOffsetMoments];

    ySol[k] = yIn[k+fOffsetMoments];
  }
#ifdef FRENKLACH
  yIn[fOffsetMoments+fNSootMoments-1] = 1.0e-20;
  ySol[fNSootMoments-1] = yIn[fOffsetMoments+fNSootMoments-1];
#endif
#ifdef HMOM
  yIn[fOffsetMoments+fNSootMoments-1] = 1.0e-20;
  ySol[fNSootMoments-1] = yIn[fOffsetMoments+fNSootMoments-1];
#endif
}
#endif

/*
//Old Post Iter for V -- PreDASS
#ifdef V
void T1DSoot::PostIter( T1DFlamePtr flame )
{
	int		i, k;
	Double		ii;
	TNewtonPtr	bt = flame->GetSolver()->bt;
	TGridPtr 	currGrid = bt->GetGrid()->GetCurrentGrid();
	int			nGridPoints = currGrid->GetNGridPoints();
	Double		*yLeft = currGrid->GetYLeft()->vec;
	Double		*yRight = currGrid->GetYRight()->vec;
	Double		**y = currGrid->GetY()->mat;
#ifdef FLUXBC
	Double		*x = currGrid->GetX()->vec;
	Double		hFirst = x[0] - bt->GetLeft();
	Double		hLast = bt->GetRight() - x[nGridPoints-1];
	Double		*rho0 = &flame->GetProperties()->GetDensity()->vec[0];
	Double		*rhoN = &flame->GetProperties()->GetDensity()->vec[nGridPoints-1];
	Double		cLeft = rho0[kPrev] * GetSootDiff()->vec[kPrev] 
						/ ( yLeft[flame->GetOffsetVVelocity()] * hFirst );
	Double		cRight = rhoN[kNext] * GetSootDiff()->vec[nGridPoints] 
						/ ( yRight[flame->GetOffsetVVelocity()] * hLast );
	TFlameNodePtr	flameNode = flame->GetFlameNode();
	Double			pressure = flame->GetPressure();

	flame->SetFlameNode( 0 );
	flame->ComputeProperties( flameNode, flameNode->temp[kCurr], flameNode->Y[kCurr], pressure );
	flame->SetFlameNode( nGridPoints-1 );
	flame->ComputeProperties( flameNode, flameNode->temp[kCurr], flameNode->Y[kCurr], pressure );
#endif
	CheckSolution(nGridPoints, y);
#ifdef FLUXBC
	if ( fSizeDepDiff ) {
		int counter = 0;
		Double	normLeft, normRight, val;
		Double	*MLeft = new Double[fNSootMoments];
		Double	*MRight = new Double[fNSootMoments];
		for ( i = 0; i < fNSootMoments; ++i ) {
			yLeft[i+fOffsetMoments] = MAX( yLeft[i+fOffsetMoments], SMALLSOOT );
			yRight[i+fOffsetMoments] = MAX( yRight[i+fOffsetMoments], SMALLSOOT );
			MLeft[i] = yLeft[i+fOffsetMoments];
			MRight[i] = yRight[i+fOffsetMoments];
		}
		
		do {
			++counter;
			normLeft = normRight = 0.0;
			for ( i = 0; i < fNSootMoments; ++i ) {
// 				yLeft[i+fOffsetMoments] = cLeft * ( rho0[kPrev]/rho0[kCurr] 
// 						* FracMom2( i-2.0/3.0, &y[0][fOffsetMoments] ) 
// 						- FracMom2( i-2.0/3.0, &yLeft[fOffsetMoments] ) )
// 	//										/ ( 1.0 + cLeft )
// 											;
// 				yRight[i+fOffsetMoments] = cRight * (
// 							FracMom2( i-2.0/3.0, &yRight[fOffsetMoments] )
// 							- FracMom2( i-2.0/3.0, &y[nGridPoints-1][fOffsetMoments] ) 
// 								* rhoN[kNext]/rhoN[kCurr] ) 
// 	//										/ ( 1.0 - cRight )
// 											;
				Double	c1, c2, f, fPrime;
	
// 				c1 = cLeft / ( 1.0 + cLeft );
// 				c2 = rho0[kPrev]/rho0[kCurr] 
// 						* FracMom2( i-2.0/3.0, &y[0][fOffsetMoments] );
// 				f = yLeft[i+fOffsetMoments] 
// 						- c1 * ( c2 - FracMom2( i-2.0/3.0, &yLeft[fOffsetMoments] ) );
// 				fPrime = 1.0 + c1 * GetDerivFracMom( i, i-2.0/3.0, &yLeft[fOffsetMoments] );
// 				yLeft[i+fOffsetMoments] -= f / fPrime;
	
				c1 = cRight / ( 1.0 - cRight );
				c2 = FracMom2( i-2.0/3.0, &y[nGridPoints-1][fOffsetMoments] );
				f = yRight[i+fOffsetMoments] 
						- c1 * ( FracMom2( i-2.0/3.0, &yRight[fOffsetMoments] ) - c2 );
				fPrime = 1.0 - c1 * GetDerivFracMom( i, i-2.0/3.0, &yRight[fOffsetMoments] );
				yRight[i+fOffsetMoments] -= f / fPrime;
				
				val = ( yLeft[i+fOffsetMoments] - MLeft[i] ) / MLeft[i];
				normLeft += val * val;
				val = ( yRight[i+fOffsetMoments] - MRight[i] ) / MRight[i];
				normRight += val * val;
				if ( i == 1 ) {
// 					cerr << "i = " << i << NEWL << "left = " << yLeft[i+fOffsetMoments]
// 						<< TAB << "old = " << MLeft[i]
// 						<< TAB << "diff = " << yLeft[i+fOffsetMoments] - MLeft[i]
// 						<< TAB << "normLeft = " << sqrt( normLeft ) << NEWL;
					cerr << "right = " << yRight[i+fOffsetMoments] 
						<< TAB << "old = " << MRight[i]
						<< TAB << "diff = " << yRight[i+fOffsetMoments] - MRight[i]
						<< TAB << "normRight = " << sqrt( normRight ) << NEWL;
				}
				MLeft[i] = yLeft[i+fOffsetMoments];
				MRight[i] = yRight[i+fOffsetMoments];
			}
		} while ( ( sqrt( normLeft ) > 0.01 || sqrt( normRight ) > 0.01 ) && counter < 100 );
		cerr << "#counter = " << counter << NEWL;
		if ( counter >=100 ) {
			cerr << "#error: iteration of boundary values not converged" << NEWL;
		}
		for ( i = 0; i < fNSootMoments; ++i ) {
			switch( i ) {
				case 0:
					yLeft[i+fOffsetMoments] = SMALLSOOT;
	//				yRight[i+fOffsetMoments] = SMALLSOOT;
					break;
				default:
					yLeft[i+fOffsetMoments] = 9.0 * yLeft[i+fOffsetMoments-1];
	//				yRight[i+fOffsetMoments] = 9.0 * yRight[i+fOffsetMoments-1];
					break;
			}
		}
		delete MLeft;
		delete MRight;
	else {
		for ( i = 0; i < fNSootMoments; ++i ) {
			yLeft[i+fOffsetMoments] = cLeft * rho0[kPrev]/rho0[kCurr] * y[0][i+fOffsetMoments] 
										/ ( 1.0 + cLeft );
			yRight[i+fOffsetMoments] = - cRight * rhoN[kNext]/rhoN[kCurr] * y[nGridPoints-1][i+fOffsetMoments] 
										/ ( 1.0 - cRight );
		}
	}
#else
	for (k = 0; k < fNSootMoments; ++k)
	{
			ii = ((double) Geti(k))/6.0;

			if (k==0)
			{
				yLeft[fOffsetMoments+k] = SMALLSOOT;
				yRight[fOffsetMoments+k] = SMALLSOOT;
			}
			else
			{
				yLeft[fOffsetMoments+k] = pow(2.0*nucl_nbrC2, ii) * yLeft[fOffsetMoments];
				yRight[fOffsetMoments+k] = pow(2.0*nucl_nbrC2, ii) * yRight[fOffsetMoments];
			}
	}
#endif
}
#endif
*/

void T1DSoot::CheckSolution(int nGridPoints, Double ** y, double * density)
{
	for (int k = 0; k < nGridPoints; ++k)
	{
	  TSoot::CheckSolution(&y[k][fOffsetMoments], density[k]);
	}
}

#ifdef VSH
void T1DSoot::PrintFlameletFile(int gridPoints, T1DFlamePtr flame, FILE *fp)
{
  int i, j, k, l, m;
  double ii, jj, ll;
  double ** moments = fMoments->mat;
  double fvFact = fMolarMassSoot / fSootDensity;
  double numdensfact = AVOGADRO * 1.0e-6; //kmole/m^3 -> parts/cm^3
  double partdiamfact = pow((6.0*fMolarMassSoot)/(PI*AVOGADRO*fSootDensity), 1.0/3.0) * 1.0e9; //nm
  double chifact = pow((6.0*fMolarMassSoot)/(PI*AVOGADRO*fSootDensity), -2.0/3.0); //1/m^2
  double * molarMass = flame->GetSpecies()->GetMolarMass()->vec;
  double * density = flame->GetProperties()->GetDensity()->vec;
  double * mixMolarMass = flame->GetProperties()->GetMolarMass()->vec;
  double ** Y = flame->GetMassFracs()->mat;
  double * temp = flame->GetTemperature()->vec;
  double * kSoot = fSootRateCoeffs->vec;
  double * rhodot = fRhoDot->vec;
  T1DSpeciesPtr species = flame->GetSpecies();

  //Write Moments
  for (m = 0; m < fNSootMoments; ++m)
  {
    i = Geti(m);
    j = Getj(m);
    l = Getk(m);

    fprintf(fp, "conc-SOOT%d-%d-%d\n", i, j, l);
			
    for (k = 0; k < gridPoints+2; ++k)
    {
      fprintf(fp, "\t%-.6e", moments[k-1][m]);
      if ((k+1) % 5 == 0)
	fprintf(fp, "\n");
    }
    if (k % 5)
      fprintf(fp, "\n");
  }

  //Write Soot Mass Fraction
  fprintf(fp, "MassFrac_Soot\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    fprintf(fp, "\t%-.6e", moments[k-1][1] * fMolarMassSoot / density[k-1]);
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");

#ifdef HMOM
  //Alpha
  fprintf(fp, "alpha\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    fprintf(fp, "\t%-.6e", (moments[k-1][0]-moments[k-1][fNSootMoments-1]) / moments[k-1][0]);
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");
#endif

  //Write Soot Volume Fraction
  fprintf(fp, "fv\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    fprintf(fp, "\t%-.6e", moments[k-1][1] * fvFact);
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");

  //Write Number Density
  fprintf(fp, "numdens [cm^-3]\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    fprintf(fp, "\t%-.6e", moments[k-1][0] * numdensfact);
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");

  //Write Particle "Diameter" Based on Volume
  fprintf(fp, "partdiamV [nm]\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    fprintf(fp, "\t%-.6e", partdiamfact * FracMom(1.0/3.0, 0.0, 0.0, moments[k-1]) / MAX(moments[k-1][0], 1.0e-30));
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");

  //Write Particle "Diameter" Based on Surface
  fprintf(fp, "partdiamS [nm]\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    fprintf(fp, "\t%-.6e", partdiamfact * FracMom(0.0, 1.0/2.0, 0.0, moments[k-1]) / MAX(moments[k-1][0], 1.0e-30));
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");

  //Write Primary Particle Diameter
  fprintf(fp, "primpartdiam [nm]\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    fprintf(fp, "\t%-.6e", partdiamfact * FracMom(1.0, -1.0, 0.0, moments[k-1]) / MAX(moments[k-1][0], 1.0e-30));
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");

  //Write Scattering/Extinction Diameter
  fprintf(fp, "d63 [nm]\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    fprintf(fp, "\t%-.6e", partdiamfact * pow(FracMom(4.0,-3.0,0.0,moments[k-1]) / MAX(FracMom(1.0,0.0,0.0,moments[k-1]),1.0e-30), 1.0/3.0));
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");

#ifdef HMOM
  //Write Primary Particle Diameter of Large Particles
  fprintf(fp, "primpartdiam-l [nm]\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    fprintf(fp, "\t%-.6e", partdiamfact * FracMomLarge(1.0, -1.0, 0.0, moments[k-1]) / MAX(moments[k-1][0]-moments[k-1][fNSootMoments-1], 1.0e-30));
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");
#endif

  //Write Primary Particle Number
  fprintf(fp, "primpartnum\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    fprintf(fp, "\t%-.6e", FracMom(-2.0, 3.0, 0.0, moments[k-1]) / MAX(moments[k-1][0], 1.0e-30));
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");

#ifdef HMOM
  //Write Primary Particle Number of Large Particles
  fprintf(fp, "primpartnum-l\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    fprintf(fp, "\t%-.6e", FracMomLarge(-2.0, 3.0, 0.0, moments[k-1]) / MAX(moments[k-1][0]-moments[k-1][fNSootMoments-1], 1.0e-30));
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");
#endif

  //Write Hydrogen Site per Surface Ratio
  fprintf(fp, "Chi [1/m^2]\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    fprintf(fp, "\t%-.6e", chifact * FracMom(0.0, -1.0, 1.0, moments[k-1]) / MAX(moments[k-1][0], 1.0e-30));
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");

#ifdef HMOM
  //Write Hydrogen Site per Surface Ratio of Large Particles
  fprintf(fp, "Chi-l [1/m^2]\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    fprintf(fp, "\t%-.6e", chifact * FracMomLarge(0.0, -1.0, 1.0, moments[k-1]) / MAX(moments[k-1][0]-moments[k-1][fNSootMoments-1], 1.0e-30));
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");
#endif

  //Write Dimer Production Rate
  fprintf(fp, "wDimer\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    ComputeDimerParticules(species, temp[k-1], Y[k-1], density[k-1], 1.0, mixMolarMass[k-1], moments[k-1]);
    fprintf(fp, "\t%-.6e", dimer_rate1);
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");

  //Write Surface Growth Coefficient
  fprintf(fp, "SootSurfGrowthCoeff\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    fprintf(fp, "\t%-.6e", SurfaceGrowthCoeffFor(Y[k-1], temp[k-1], density[k-1], molarMass) - SurfaceGrowthCoeffBack(Y[k-1], temp[k-1], density[k-1], molarMass));
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");

  //Write Oxidation Coefficient
  fprintf(fp, "SootOxCoeff\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    fprintf(fp, "\t%-.6e", OxidationCoeff(Y[k-1], temp[k-1], density[k-1], molarMass));
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");

  //Write Soot Source Terms
  for (m = 0; m < fNSootMoments; ++m)
  {
    i = Geti(m);
    ii = double(Geti(m))/6.0;
    j = Getj(m);
    jj = double(Getj(m))/6.0;
    l = Getk(m);
    ll = double(Getk(m))/6.0;

    //Nucleation
    if (fNucleation)
    {
      fprintf(fp, "Nucleation%d-%d-%d\n", i, j, l);
      for (k = 0; k < gridPoints+2; ++k)
      {
	ComputeDimerParticules(species, temp[k-1], Y[k-1], density[k-1], 1.0, mixMolarMass[k-1], moments[k-1]);
	if (i==50)
	  fprintf(fp, "\t%-.6e", NucleationSource(0, 0, 0, temp[k-1]));
	else
	  fprintf(fp, "\t%-.6e", NucleationSource(ii, jj, ll, temp[k-1]));
	if ((k+1) % 5 == 0)
	  fprintf(fp, "\n");
      }
      if (k % 5)
	fprintf(fp, "\n");
    }

    //Coagulation
    if (fCoagulation)
    {
      fprintf(fp, "Coagulation%d-%d-%d\n", i, j, l);
      for (k = 0; k < gridPoints+2; ++k)
      {
	if (i==50)
	  fprintf(fp, "\t%-.6e", CoagulationSourceSmall(temp[k-1], moments[k-1]));
	else
	  fprintf(fp, "\t%-.6e", CoagulationSource(ii, jj, ll, temp[k-1], moments[k-1]));
	if ((k+1) % 5 == 0)
	  fprintf(fp, "\n");
      }
      if (k % 5)
	fprintf(fp, "\n");
    }

    //Condensation
    if (fCondensation)
    {
      fprintf(fp, "Condensation%d-%d-%d\n", i, j, l);
      for (k = 0; k < gridPoints+2; ++k)
      {
	ComputeDimerParticules(species, temp[k-1], Y[k-1], density[k-1], 1.0, mixMolarMass[k-1], moments[k-1]);

	if (i==50)
	  fprintf(fp, "\t%-.6e", CondensationSourceSmall(temp[k-1], moments[k-1]));
	else
	  fprintf(fp, "\t%-.6e", CondensationSource(ii, jj, ll, temp[k-1], moments[k-1]));
	if ((k+1) % 5 == 0)
	  fprintf(fp, "\n");
      }
      if (k % 5)
	fprintf(fp, "\n");
    }

    //Surface Growth
    if (fSurfaceGrowth)
    {
      fprintf(fp, "SurfGrowth%d-%d-%d\n", i, j, l);
      for (k = 0; k < gridPoints+2; ++k)
      {
	ComputeSootRateCoeffs(kSoot, temp[k-1], flame->GetReaction());
	ComputeCSootStar(kSoot, Y[k-1], density[k-1], molarMass, mixMolarMass[k-1]);
	if (i==50)
	  fprintf(fp, "\t%-.6e", SurfaceGrowthSourceSmall(moments[k-1], Y[k-1], temp[k-1], density[k-1], molarMass));
	else
	  fprintf(fp, "\t%-.6e", SurfaceGrowthSource(ii, jj, ll, moments[k-1], Y[k-1], temp[k-1], density[k-1], molarMass));
	if ((k+1) % 5 == 0)
	  fprintf(fp, "\n");
      }
      if (k % 5)
	fprintf(fp, "\n");
    }

    //Oxidation
    if (fSurfaceOxidation)
    {
      fprintf(fp, "SurfOx%d-%d-%d\n", i, j, l);
      for (k = 0; k < gridPoints+2; ++k)
      {
	ComputeSootRateCoeffs(kSoot, temp[k-1], flame->GetReaction());
	ComputeCSootStar(kSoot, Y[k-1], density[k-1], molarMass, mixMolarMass[k-1]);
	if (i==50)
	  fprintf(fp, "\t^-.63", OxidationSourceSmall(moments[k-1], Y[k-1], temp[k-1], density[k-1], molarMass));
	else
	  fprintf(fp, "\t%-.6e", OxidationSource(ii, jj, ll, moments[k-1], Y[k-1], temp[k-1], density[k-1], molarMass));
	if ((k+1) % 5 == 0)
	  fprintf(fp, "\n");
      }
      if (k % 5)
	fprintf( fp, "\n" );
    }
  }

  //Density Correction
  fprintf(fp, "RhoDot\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    fprintf(fp, "\t%-.6e", rhodot[k-1]);
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");
}
#endif
#ifdef VS
void T1DSoot::PrintFlameletFile(int gridPoints, T1DFlamePtr flame, FILE *fp)
{
  int i, j, k, l;
  double ii, jj;
  double ** moments = fMoments->mat;
  double fvFact = fMolarMassSoot / fSootDensity;
  double numdensfact = AVOGADRO * 1.0e-6; //kmole/m^3 -> parts/cm^3
  double partdiamfact = pow((6.0*fMolarMassSoot)/(PI*AVOGADRO*fSootDensity), 1.0/3.0) * 1.0e9; //nm
  double * molarMass = flame->GetSpecies()->GetMolarMass()->vec;
  double * density = flame->GetProperties()->GetDensity()->vec;
  double * mixMolarMass = flame->GetProperties()->GetMolarMass()->vec;
  double * mixViscosity = flame->GetProperties()->GetViscosity()->vec;
  double ** Y = flame->GetMassFracs()->mat;
  double * temp = flame->GetTemperature()->vec;
  double * kSoot = fSootRateCoeffs->vec;
  double * rhodot = fRhoDot->vec;
  double * enthdot = fEnthDot->vec;
  T1DSpeciesPtr species = flame->GetSpecies();

  //Write Moments
  for (l = 0; l < fNSootMoments; ++l) //
  {
    i = Geti(l);
    j = Getj(l);

    fprintf(fp, "conc-SOOT%d-%d\n", i, j);
    for (k = 0; k < gridPoints+2; ++k)
    {
      fprintf(fp, "\t%-.6e", moments[k-1][l]);
      if ((k+1) % 5 == 0)
	fprintf(fp, "\n");
    }
    if (k % 5)
      fprintf(fp, "\n");
  }

  //Write Soot Mass Fraction
  fprintf(fp, "MassFrac_Soot\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    fprintf(fp, "\t%-.6e", moments[k-1][1] * fMolarMassSoot / density[k-1]); 
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");

#ifdef HMOM
  //Alpha
  fprintf(fp, "alpha\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    fprintf(fp, "\t%-.6e", (moments[k-1][0]-moments[k-1][fNSootMoments-1]) / moments[k-1][0]); //
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");
#endif

  //Write Soot Volume Fraction
  fprintf(fp, "fv\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    fprintf(fp, "\t%-.6e", moments[k-1][1] * fvFact);
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");

  // Added by SD
  // Compute and write the integral of soot volume fraction
  // Double sootvolume = 0.0;
  // Double dsootvolume = 0.0;
  // Double * x=flame->GetSolver()->bt->GetGrid()->GetCurrentGrid()->GetX()->vec;
  // sootvolume = 0.5 * (moments[-1][1] * fvFact + moments[0][1] * fvFact) * (x[0] - flame->GetSolver()->bt->GetLeft());
  // for ( k = 0; k < gridPoints-1; ++k ) {
  //   dsootvolume = 0.5 * (moments[k+1-1][1] * fvFact + moments[k+1][1] * fvFact) * (x[k+1] - x[k]);
  //   sootvolume = sootvolume + dsootvolume;
  // }
  // dsootvolume = 0.5 * (moments[gridPoints-1][1] * fvFact + moments[gridPoints][1] * fvFact) * (flame->GetSolver()->bt->GetRight() - x[gridPoints-1]);
  // sootvolume = sootvolume + dsootvolume;
  // fprintf(fp, "SootVolume\n");
  // for ( k = 0; k < gridPoints+2; ++k){
  //   fprintf(fp, "\t%-.6e", sootvolume);
  //   if ((k+1) % 5 == 0)
  //     fprintf(fp, "\n");
  // }
  // if (k % 5)
  //   fprintf(fp, "\n");
 
  // Compute and write the integral of soot volume fraction
  Double sootvolume = 0.0;
  Double dsootvolume = 0.0;
  Double * x=flame->GetSolver()->bt->GetGrid()->GetCurrentGrid()->GetX()->vec;
  Double radiation = 0.0;
  Double dradiation = 0.0;
  sootvolume = 0.5 * (moments[-1][1] * fvFact + moments[0][1] * fvFact) * (x[0] - flame->GetSolver()->bt->GetLeft());
  radiation = sootvolume * pow(0.5 * (temp[-1] + temp[0]), 4.0);
  for ( k = 0; k < gridPoints-1; ++k ) {
    dsootvolume = 0.5 * (moments[k+1-1][1] * fvFact + moments[k+1][1] * fvFact) * (x[k+1] - x[k]);
    sootvolume = sootvolume + dsootvolume;
    dradiation = dsootvolume * pow(0.5 * (temp[k] + temp[k+1]), 4.0);
    radiation = radiation + dradiation;
  }
  dsootvolume = 0.5 * (moments[gridPoints-1][1] * fvFact + moments[gridPoints][1] * fvFact) * (flame->GetSolver()->bt->GetRight() - x[gridPoints-1]);
  sootvolume = sootvolume + dsootvolume;
  dradiation = dsootvolume * pow(0.5 * (temp[gridPoints-1] + temp[gridPoints]), 4.0);
  radiation = radiation + dradiation;

  fprintf(fp, "SootVolume\n");
  for ( k = 0; k < gridPoints+2; ++k){
    fprintf(fp, "\t%-.6e", sootvolume);
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");

  // Write Radiation
  fprintf(fp, "Radiation\n");
  for ( k = 0; k < gridPoints+2; ++k){
    fprintf(fp, "\t%-.6e", radiation);
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");
  // End SD

  //Write Number Density
  fprintf(fp, "numdens [cm^-3]\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    fprintf(fp, "\t%-.6e", moments[k-1][0] * numdensfact);
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");

  //Write Particle "Diameter" Based on Volume
  fprintf(fp, "partdiamV [nm]\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    fprintf(fp, "\t%-.6e", partdiamfact * FracMom(1.0/3.0, 0.0, moments[k-1]) / MAX(moments[k-1][0], 1.0e-30));
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");

  //Write Particle "Diameter" Based on Surface
  fprintf(fp, "partdiamS [nm]\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    fprintf(fp, "\t%-.6e", partdiamfact * FracMom(0.0, 1.0/2.0, moments[k-1]) / MAX(moments[k-1][0], 1.0e-30));
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");

  //Write Primary Particle Diameter
  fprintf(fp, "primpartdiam [nm]\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    fprintf(fp, "\t%-.6e", partdiamfact * FracMom(1.0, -1.0, moments[k-1]) / MAX(moments[k-1][0], 1.0e-30));
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");

  //Write Scattering/Extinction Diameter
  fprintf(fp, "d63 [nm]\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    fprintf(fp, "\t%-.6e", partdiamfact * pow(FracMom(4.0,-3.0,moments[k-1]) / MAX(FracMom(1.0,0.0,moments[k-1]),1.0e-30), 1.0/3.0));
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");

#ifdef HMOM
  //Write Primary Particle Diameter of Large Particles
  fprintf(fp, "primpartdiam-l [nm]\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    fprintf(fp, "\t%-.6e", partdiamfact * FracMomLarge(1.0, -1.0, moments[k-1]) / MAX(moments[k-1][0]-moments[k-1][fNSootMoments-1], 1.0e-30)); //
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");
#endif

  //Write Primary Particle Number
  fprintf(fp, "primpartnum\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    fprintf(fp, "\t%-.6e", FracMom(-2.0, 3.0, moments[k-1]) / MAX(moments[k-1][0], 1.0e-30));
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");

#ifdef HMOM
  //Write Primary Particle Number of Large Particles
  fprintf(fp, "primpartnum-l\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    fprintf(fp, "\t%-.6e", FracMomLarge(-2.0, 3.0, moments[k-1]) / MAX(moments[k-1][0]-moments[k-1][fNSootMoments-1], 1.0e-30)); //
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");
#endif

  double Df = 1.8;

  //Write Collision/Mobility Diameter
  fprintf(fp, "colldiam [nm]\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    fprintf(fp, "\t%-.6e", partdiamfact * FracMom(1.0-2.0/Df, 3.0/Df-1.0, moments[k-1]) / MAX(moments[k-1][0], 1.0e-30));
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");

 #ifdef HMOM
  //Write Collision/Mobility Diameter of Large Particles
  fprintf(fp, "colldiam-l [nm]\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    fprintf(fp, "\t%-.6e", partdiamfact * FracMomLarge(1.0-2.0/Df, 3.0/Df-1.0, moments[k-1]) / MAX(moments[k-1][0]-moments[k-1][fNSootMoments-1], 1.0e-30)); //
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");
#endif

  double sigma;

  //Write "Sigma" of Volume
  fprintf(fp, "sigmaV\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    sigma = sqrt( log( MAX((FracMom(2.0, 0.0, moments[k-1]) * FracMom(0.0, 0.0, moments[k-1])) /
		           (FracMom(1.0, 0.0, moments[k-1]) * FracMom(1.0, 0.0, moments[k-1])), 1.0) ) );
    fprintf(fp, "\t%-.6e", sigma);
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");

  //Write "Sigma" of Surface Area
  fprintf(fp, "sigmaS\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    sigma = sqrt( log( MAX((FracMom(0.0, 2.0, moments[k-1]) * FracMom(0.0, 0.0, moments[k-1])) /
		           (FracMom(0.0, 1.0, moments[k-1]) * FracMom(0.0, 1.0, moments[k-1])), 1.0) ) );
    fprintf(fp, "\t%-.6e", sigma);
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");

  //Write "Sigma" of Primary Particle Diameter
  fprintf(fp, "sigmadp\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    sigma = sqrt( log( MAX((FracMom(2.0, -2.0, moments[k-1]) * FracMom(0.0,  0.0, moments[k-1])) /
		           (FracMom(1.0, -1.0, moments[k-1]) * FracMom(1.0, -1.0, moments[k-1])), 1.0) ) );
    fprintf(fp, "\t%-.6e", sigma);
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");

  //Write "Sigma" of Primary Particle Number
  fprintf(fp, "sigmanp\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    sigma = sqrt( log( MAX((FracMom(-4.0, 6.0, moments[k-1]) * FracMom( 0.0, 0.0, moments[k-1])) /
		           (FracMom(-2.0, 3.0, moments[k-1]) * FracMom(-2.0, 3.0, moments[k-1])), 1.0) ) );
    fprintf(fp, "\t%-.6e", sigma);
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");

#ifdef HMOM
  //Write "Sigma" of Volume of Large Particles
  fprintf(fp, "sigmaV-l\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    sigma = sqrt( MAX(log( MAX((FracMomLarge(2.0, 0.0, moments[k-1]) * FracMomLarge(0.0, 0.0, moments[k-1])) /
			       MAX(FracMomLarge(1.0, 0.0, moments[k-1]) * FracMomLarge(1.0, 0.0, moments[k-1]), 1.0e-60), 1.0)), 1.0e-60) );
    fprintf(fp, "\t%-.6e", sigma);
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");

  //Write "Sigma" of Surface Area of Large Particles
  fprintf(fp, "sigmaS-l\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    sigma = sqrt( MAX(log( MAX((FracMomLarge(0.0, 2.0, moments[k-1]) * FracMomLarge(0.0, 0.0, moments[k-1])) /
		           MAX(FracMomLarge(0.0, 1.0, moments[k-1]) * FracMomLarge(0.0, 1.0, moments[k-1]), 1.0e-60), 1.0)), 1.0e-60));
    fprintf(fp, "\t%-.6e", sigma);
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");

  //Write "Sigma" of Primary Particle Diameter of Large Particles
  fprintf(fp, "sigmadp-l\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    sigma = sqrt( MAX(log( MAX((FracMomLarge(2.0, -2.0, moments[k-1]) * FracMomLarge(0.0,  0.0, moments[k-1])) /
			       MAX(FracMomLarge(1.0, -1.0, moments[k-1]) * FracMomLarge(1.0, -1.0, moments[k-1]), 1.0e-60), 1.0)), 1.0e-60));
    fprintf(fp, "\t%-.6e", sigma);
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");

  //Write "Sigma" of Primary Particle Number of Large Particles
  fprintf(fp, "sigmanp-l\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    sigma = sqrt( MAX(log( MAX((FracMomLarge(-4.0, 6.0, moments[k-1]) * FracMomLarge( 0.0, 0.0, moments[k-1])) /
			       MAX(FracMomLarge(-2.0, 3.0, moments[k-1]) * FracMomLarge(-2.0, 3.0, moments[k-1]), 1.0e-60), 1.0)), 1.0e-60));
    fprintf(fp, "\t%-.6e", sigma);
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");
#endif

  //Write Dimer Production Rate
  fprintf(fp, "Dimer_ProdRate [kmol/m3-2]\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    ComputeDimerParticules(species, temp[k-1], Y[k-1], density[k-1], 1.0, mixMolarMass[k-1], moments[k-1]);
    fprintf(fp, "\t%-.6e", dimer_rate1);
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");

  //Write Dimer Carbon Content
  fprintf(fp, "Dimer_nbrC\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    ComputeDimerParticules(species, temp[k-1], Y[k-1], density[k-1], 1.0, mixMolarMass[k-1], moments[k-1]);
    fprintf(fp, "\t%-.6e", 2.0*dimer_nbrC2);
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");

  //Write Dimer Hydrogen Content
  fprintf(fp, "Dimer_nbrH\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    ComputeDimerParticules(species, temp[k-1], Y[k-1], density[k-1], 1.0, mixMolarMass[k-1], moments[k-1]);
    fprintf(fp, "\t%-.6e", dimer_nbrH);
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");

  //Write Surface Growth Coefficient
  fprintf(fp, "SootSurfGrowthCoeff\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    ComputeSootRateCoeffs(kSoot, temp[k-1], flame->GetReaction());
    ComputeCSootStar(kSoot, Y[k-1], density[k-1], molarMass, mixMolarMass[k-1]);
    fprintf(fp, "\t%-.6e", SurfaceGrowthCoeffFor(Y[k-1], temp[k-1], density[k-1], molarMass) - SurfaceGrowthCoeffBack(Y[k-1], temp[k-1], density[k-1], molarMass));
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");

  //Write Oxidation Coefficient
  fprintf(fp, "SootOxCoeff\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    ComputeSootRateCoeffs(kSoot, temp[k-1], flame->GetReaction());
    ComputeCSootStar(kSoot, Y[k-1], density[k-1], molarMass, mixMolarMass[k-1]);
    fprintf(fp, "\t%-.6e", OxidationCoeff(Y[k-1], temp[k-1], density[k-1], molarMass));
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");

  //Write O2 Oxidation Coefficient (for fragmentation)
  fprintf(fp, "SootOxCoeff_O2\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    ComputeSootRateCoeffs(kSoot, temp[k-1], flame->GetReaction());
    ComputeCSootStar(kSoot, Y[k-1], density[k-1], molarMass, mixMolarMass[k-1]);
    fprintf(fp, "\t%-.6e", OxidationCoeff_O2(Y[k-1], temp[k-1], density[k-1], molarMass));
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");

  //Write Individual Oxidation Rates
  fprintf(fp, "SootOxRateO2 [1/s]\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    ComputeSootRateCoeffs(kSoot, temp[k-1], flame->GetReaction());
    ComputeCSootStar(kSoot, Y[k-1], density[k-1], molarMass, mixMolarMass[k-1]);
    fprintf(fp, "\t%-.6e", OxidationRate_O2(moments[k-1], Y[k-1], temp[k-1], density[k-1], molarMass));
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");

  fprintf(fp, "SootOxRateOH [1/s]\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    ComputeSootRateCoeffs(kSoot, temp[k-1], flame->GetReaction());
    ComputeCSootStar(kSoot, Y[k-1], density[k-1], molarMass, mixMolarMass[k-1]);
    fprintf(fp, "\t%-.6e", OxidationRate_OH(moments[k-1], Y[k-1], temp[k-1], density[k-1], molarMass));
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");


  //Write Soot Source Terms
  for (l = 0; l < fNSootMoments; ++l) //
  {
    i = Geti(l);
    ii = double(Geti(l))/6.0;
    j = Getj(l);
    jj = double(Getj(l))/6.0;

    //Nucleation
    if (fNucleation)
    {
      fprintf(fp, "Nucleation%d-%d\n", i, j);
      for (k = 0; k < gridPoints+2; ++k)
      {
	ComputeDimerParticules(species, temp[k-1], Y[k-1], density[k-1], 1.0, mixMolarMass[k-1], moments[k-1]);
	if (i==50)
	  fprintf(fp, "\t%-.6e", NucleationSource(0, 0, temp[k-1]));
	else
	  fprintf(fp, "\t%-.6e", NucleationSource(ii, jj, temp[k-1]));
	if ((k+1) % 5 == 0)
	  fprintf(fp, "\n");
      }
      if (k % 5)
	fprintf(fp, "\n");
    }

    //Coagulation
    if (fCoagulation)
    {
      fprintf(fp, "Coagulation%d-%d\n", i, j);
      for (k = 0; k < gridPoints+2; ++k)
      {
	if (i==50)
	  fprintf(fp, "\t%-.6e", CoagulationSourceSmall(temp[k-1], moments[k-1]));
	else
	  fprintf(fp, "\t%-.6e", CoagulationSource(ii, jj, temp[k-1], moments[k-1]));
	if ((k+1) % 5 == 0)
	  fprintf(fp, "\n");
      }
      if (k % 5)
	fprintf(fp, "\n");
    }

    //Condensation
    if (fCondensation)
    {
      fprintf(fp, "Condensation%d-%d\n", i, j);
      for (k = 0; k < gridPoints+2; ++k)
      {
	ComputeDimerParticules(species, temp[k-1], Y[k-1], density[k-1], 1.0, mixMolarMass[k-1], moments[k-1]);
	if (i==50)
	  fprintf(fp, "\t%-.6e", CondensationSourceSmall(temp[k-1], moments[k-1]));
	else
	  fprintf(fp, "\t%-.6e", CondensationSource(ii, jj, temp[k-1], moments[k-1]));
	if ((k+1) % 5 == 0)
	  fprintf(fp, "\n");
      }
      if (k % 5)
	fprintf(fp, "\n");
    }

    //Surface Growth
    if (fSurfaceGrowth)
    {
      fprintf(fp, "SurfGrowth%d-%d\n", i, j);
      for (k = 0; k < gridPoints+2; ++k)
      {
	ComputeSootRateCoeffs(kSoot, temp[k-1], flame->GetReaction());
	ComputeCSootStar(kSoot, Y[k-1], density[k-1], molarMass, mixMolarMass[k-1]);
	if (i==50)
	  fprintf(fp, "\t%-.6e", SurfaceGrowthSourceSmall(moments[k-1], Y[k-1], temp[k-1], density[k-1], molarMass));
	else
	  fprintf(fp, "\t%-.6e", SurfaceGrowthSource(ii, jj, moments[k-1], Y[k-1], temp[k-1], density[k-1], molarMass));
	if ((k+1) % 5 == 0)
	  fprintf(fp, "\n");
      }
      if (k % 5)
	fprintf(fp, "\n");
    }

    //Oxidation
    if (fSurfaceOxidation)
    {
      fprintf(fp, "SurfOx%d-%d\n", i, j);
      for (k = 0; k < gridPoints+2; ++k)
      {
	ComputeSootRateCoeffs(kSoot, temp[k-1], flame->GetReaction());
	ComputeCSootStar(kSoot, Y[k-1], density[k-1], molarMass, mixMolarMass[k-1]);
	if (i==50)
	  fprintf(fp, "\t%-.6e", OxidationSourceSmall(moments[k-1], Y[k-1], temp[k-1], density[k-1], molarMass));
	else
	  fprintf(fp, "\t%-.6e", OxidationSource(ii, jj, moments[k-1], Y[k-1], temp[k-1], density[k-1], molarMass));
	if ((k+1) % 5 == 0)
	  fprintf(fp, "\n");
      }
      if (k % 5)
	fprintf(fp, "\n");
    }

    //Fragmentation
    if (fFragmentation)
    {
      fprintf(fp, "Frag%d-%d\n", i, j);
      for (k = 0; k < gridPoints+2; ++k)
      {
	if (i==50)
	  fprintf(fp, "\t%-.6e", FragmentationSourceSmall(moments[k-1]));
	else
	  fprintf(fp, "\t%-.6e", FragmentationSource(ii, jj, moments[k-1]));
	if ((k+1) % 5 == 0)
	  fprintf(fp, "\n");
      }
      if (k % 5)
	fprintf(fp, "\n");
    }	  
  }

  //Density Correction
  fprintf(fp, "RhoDot [kg/m3-s]\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    fprintf(fp, "\t%-.6e", rhodot[k-1]);
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");

  //Enthalpy Correction
  fprintf(fp, "EnthDot [J/m3-s]\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    fprintf(fp, "\t%-.6e", enthdot[k-1]);
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");

  // Quantities Needed for Flamelets
    // Free Molecular Collisions: Square Root of Temperature
    fprintf(fp, "sqrtT\n");
    for (k = 0; k < gridPoints+2; ++k)
    {
      fprintf(fp, "\t%-.6e", sqrt(temp[k-1]));
      if ((k+1) % 5 == 0)
	fprintf(fp, "\n");
    }
    if (k % 5)
      fprintf(fp, "\n");

    // Continuum Collisions: Temperature Over Viscosity
    fprintf(fp, "T_MU\n");
    for (k = 0; k < gridPoints+2; ++k)
    {
      fprintf(fp, "\t%-.6e", temp[k-1] / mixViscosity[k-1]);
      if ((k+1) % 5 == 0)
	fprintf(fp, "\n");
    }
    if (k % 5)
      fprintf(fp, "\n");

    // Slip Factor: Viscosity Over Density Times Square Root of MolarMass Over Temperature
    fprintf(fp, "MUsqrtW_RHOsqrtT\n");
    for (k = 0; k < gridPoints+2; ++k)
    {
      fprintf(fp, "\t%-.6e", mixViscosity[k-1] / density[k-1] * sqrt((mixMolarMass[k-1]/1000.0) / temp[k-1]));
      if ((k+1) % 5 == 0)
	fprintf(fp, "\n");
    }
    if (k % 5)
      fprintf(fp, "\n");

    // Full Velocity
    fprintf(fp, "fullvel [m/s]\n");
    for (k = 0; k < gridPoints+2; ++k)
    {
      if (k!=0 && k!=gridPoints+1)
      {
	flame->SetFlameNode(k-1);
	fprintf(fp, "\t%-.6e", GetFullVelocity(flame, k-1));
      }
      else
	fprintf(fp, "\t%-.6e", 0.0);
      
      if ((k+1) % 5 == 0)
	fprintf(fp, "\n");
    }
    if (k % 5)
      fprintf(fp, "\n");
}
#endif
#ifdef V
void T1DSoot::PrintFlameletFile(int gridPoints, T1DFlamePtr flame, FILE *fp)
{
  int i, k, l;
  double ii;
  double ** moments = fMoments->mat;
  double fvFact = fMolarMassSoot / fSootDensity;
  double numdensfact = AVOGADRO * 1.0e-6; //kmole/m^3 -> parts/cm^3
  double partdiamfact = pow((6.0*fMolarMassSoot)/(PI*AVOGADRO*fSootDensity), 1.0/3.0) * 1.0e9; //nm
  double partsurffact = PI * pow((6.0*fMolarMassSoot)/(PI*AVOGADRO*fSootDensity), 2.0/3.0) * 1.0e12; //um^2
  double * molarMass = flame->GetSpecies()->GetMolarMass()->vec;
  double * density = flame->GetProperties()->GetDensity()->vec;
  double * mixMolarMass = flame->GetProperties()->GetMolarMass()->vec;
  double ** Y = flame->GetMassFracs()->mat;
  double * temp = flame->GetTemperature()->vec;
  double * kSoot = fSootRateCoeffs->vec;
  double * rhodot = fRhoDot->vec;
  T1DSpeciesPtr species = flame->GetSpecies();

  //Write Moments
  for (l = 0; l < fNSootMoments; ++l)
  {
    i = Geti(l);

    fprintf(fp, "conc-SOOT%d\n", i);
			
    for (k = 0; k < gridPoints+2; ++k)
    {
      fprintf(fp, "\t%-.6e", moments[k-1][l]);
      if ((k+1) % 5 == 0)
	fprintf(fp, "\n");
    }
    if (k % 5)
      fprintf(fp, "\n");
  }

  //Write Soot Mass Fraction
  fprintf(fp, "MassFrac_Soot\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    fprintf(fp, "\t%-.6e", moments[k-1][1] * fMolarMassSoot / density[k-1]);
    if ((k+1) % 5 == 0)
    {
      fprintf(fp, "\n");
    }
  }
  if (k % 5)
  {
    fprintf(fp, "\n");
  }

#ifdef HMOM
  //Alpha
  fprintf(fp, "alpha\n");
  for (k = 0; k < gridPoints+2; ++k)
  {
    fprintf(fp, "\t%-.6e", (moments[k-1][0]-moments[k-1][fNSootMoments-1]) / moments[k-1][0]);
    if ((k+1) % 5 == 0)
      fprintf(fp, "\n");
  }
  if (k % 5)
    fprintf(fp, "\n");
#endif

	//Write Soot Volume Fraction
	fprintf(fp, "fv\n");
	for (k = 0; k < gridPoints+2; ++k)
	{
		fprintf(fp, "\t%-.6e", moments[k-1][1] * fvFact);
		if ((k+1) % 5 == 0)
		{
			fprintf(fp, "\n");
		}
	}
	if (k % 5)
	{
		fprintf(fp, "\n");
	}

	//Write Number Density
	fprintf(fp, "numdens [cm^-3]\n");
	for (k = 0; k < gridPoints+2; ++k)
	{
		fprintf(fp, "\t%-.6e", moments[k-1][0] * numdensfact);
		if ((k+1) % 5 == 0)
		{
			fprintf(fp, "\n");
		}
	}
	if (k % 5)
	{
		fprintf(fp, "\n");
	}

	//Write Particle Diameter
	fprintf(fp, "partdiam [nm]\n");
	for (k = 0; k < gridPoints+2; ++k)
	{
		fprintf(fp, "\t%-.6e", partdiamfact * FracMom(1.0/3.0, moments[k-1]) / MAX(moments[k-1][0], 1.0e-30));
		if ((k+1) % 5 == 0)
		{
			fprintf(fp, "\n");
		}
	}
	if (k % 5)
	{
		fprintf(fp, "\n");
	}

	//Write Particle Diameter based on 
	fprintf(fp, "partdiam12 [nm]\n");
	for (k = 0; k < gridPoints+2; ++k)
	{
		fprintf(fp, "\t%-.6e", partdiamfact * pow(moments[k-1][2] / MAX(moments[k-1][1], 1.0e-30), 1.0/3.0));
		if ((k+1) % 5 == 0)
		{
			fprintf(fp, "\n");
		}
	}
	if (k % 5)
	{
		fprintf(fp, "\n");
	}

#ifdef HMOM
	//Write Particle Diameter of Large Particles
	fprintf(fp, "partdiam-l [nm]\n");
	for (k = 0; k < gridPoints+2; ++k)
	{
	  fprintf(fp, "\t%-.6e", partdiamfact * FracMomLarge(1.0/3.0, moments[k-1]) / MAX(FracMomLarge(0.0, moments[k-1]), 1.0e-30));
		if ((k+1) % 5 == 0)
		{
			fprintf(fp, "\n");
		}
	}
	if (k % 5)
	{
		fprintf(fp, "\n");
	}
#endif

	//Write Particle Surface
	fprintf(fp, "partsurf [um^2]\n");
	for (k = 0; k < gridPoints+2; ++k)
	{
		fprintf(fp, "\t%-.6e", partsurffact * FracMom(2.0/3.0, moments[k-1]) / MAX(moments[k-1][0], 1.0e-30));
		if ((k+1) % 5 == 0)
		{
			fprintf(fp, "\n");
		}
	}
	if (k % 5)
	{
		fprintf(fp, "\n");
	}

#ifdef HMOM
	//Write Particle Surface of Large Particles
	fprintf(fp, "partsurf-l [um^2]\n");
	for (k = 0; k < gridPoints+2; ++k)
	{
	  fprintf(fp, "\t%-.6e", partsurffact * FracMomLarge(2.0/3.0, moments[k-1]) / MAX(FracMomLarge(0.0, moments[k-1]), 1.0e-30));
		if ((k+1) % 5 == 0)
		{
			fprintf(fp, "\n");
		}
	}
	if (k % 5)
	{
		fprintf(fp, "\n");
	}
#endif

#ifndef LINDSTEDT
	//Write Surface Growth Coefficient
	fprintf(fp, "SootSurfGrowthCoeff\n");
	for (k = 0; k < gridPoints+2; ++k)
	{
	  fprintf(fp, "\t%-.6e", SurfaceGrowthCoeffFor(Y[k-1], temp[k-1], density[k-1], molarMass) - SurfaceGrowthCoeffBack(Y[k-1], temp[k-1], density[k-1], molarMass));
		if ((k+1) % 5 == 0)
		{
			fprintf(fp, "\n");
		}
	}
	if (k % 5)
	{
		fprintf(fp, "\n");
	}

	//Write Oxidation Coefficient
	fprintf(fp, "SootOxCoeff\n");
	for (k = 0; k < gridPoints+2; ++k)
	{
	  fprintf(fp, "\t%-.6e", OxidationCoeff(Y[k-1], temp[k-1], density[k-1], molarMass));
		if ((k+1) % 5 == 0)
		{
			fprintf(fp, "\n");
		}
	}
	if (k % 5)
	{
		fprintf(fp, "\n");
	}

	//Write Soot Source Terms
	for (l = 0; l < fNSootMoments; ++l)
	{
		i = Geti(l);
		ii = double(Geti(l))/6.0;

		//Nucleation
		if (fNucleation)
		{
			fprintf(fp, "Nucleation%d\n", i);
			for (k = 0; k < gridPoints+2; ++k)
			{
				ComputeDimerParticules(species, temp[k-1], Y[k-1], density[k-1], 1.0, mixMolarMass[k-1], moments[k-1]);
				if (i == 50)
				  fprintf(fp, "\t%-.6e", NucleationSource(0, temp[k-1]));
				else
				  fprintf(fp, "\t%-.6e", NucleationSource(ii, temp[k-1]));
				if ((k+1) % 5 == 0)
				{
					fprintf(fp, "\n");
				}
			}
			if (k % 5)
			{
				fprintf(fp, "\n");
			}
		}

		//Coagulation
		if (fCoagulation)
		{
			fprintf(fp, "Coagulation%d\n", i);
			for (k = 0; k < gridPoints+2; ++k)
			{
			  if (i == 50)
			    fprintf(fp, "\t%-.6e", CoagulationSourceSmall(temp[k-1], moments[k-1]));
			  else
			    fprintf(fp, "\t%-.6e", CoagulationSource(ii, temp[k-1], moments[k-1]));
				if ((k+1) % 5 == 0)
				{
					fprintf(fp, "\n");
				}
			}
			if (k % 5)
			{
				fprintf(fp, "\n");
			}
		}

		//Condensation
		if (fCondensation)
		{
			fprintf(fp, "Condensation%d\n", i);
			for (k = 0; k < gridPoints+2; ++k)
			{
				ComputeDimerParticules(species, temp[k-1], Y[k-1], density[k-1], 1.0, mixMolarMass[k-1], moments[k-1]);
				if (i == 50)
				  fprintf(fp, "\t%-.6e", CondensationSourceSmall(temp[k-1], moments[k-1]));
				else
				  fprintf(fp, "\t%-.6e", CondensationSource(ii, temp[k-1], moments[k-1]));
				if ((k+1) % 5 == 0)
				{
					fprintf(fp, "\n");
				}
			}
			if (k % 5)
			{
				fprintf(fp, "\n");
			}
		}

		//Surface Growth
		if (fSurfaceGrowth)
		{
			fprintf(fp, "SurfGrowth%d\n", i);
			for (k = 0; k < gridPoints+2; ++k)
			{
				ComputeSootRateCoeffs(kSoot, temp[k-1], flame->GetReaction());
				ComputeCSootStar(kSoot, Y[k-1], density[k-1], molarMass, mixMolarMass[k-1]);
				if (i == 50)
				  fprintf(fp, "\t%-.6e", SurfaceGrowthSourceSmall(moments[k-1], Y[k-1], temp[k-1], density[k-1], molarMass));
				else
				  fprintf(fp, "\t%-.6e", SurfaceGrowthSource(ii, moments[k-1], Y[k-1], temp[k-1], density[k-1], molarMass));
				if ((k+1) % 5 == 0)
				{
					fprintf(fp, "\n");
				}
			}
			if (k % 5)
			{
				fprintf(fp, "\n");
			}
		}

		//Oxidation
		if (fSurfaceOxidation)
		{
			fprintf(fp, "SurfOx%d\n", i);
			for (k = 0; k < gridPoints+2; ++k)
			{
				ComputeSootRateCoeffs(kSoot, temp[k-1], flame->GetReaction());
				ComputeCSootStar(kSoot, Y[k-1], density[k-1], molarMass, mixMolarMass[k-1]);
				if (i == 50)
				  fprintf(fp, "\t%-.6e", OxidationSourceSmall(moments[k-1], Y[k-1], temp[k-1], density[k-1], molarMass));
				else
				  fprintf(fp, "\t%-.6e", OxidationSource(ii, moments[k-1], Y[k-1], temp[k-1], density[k-1], molarMass));
				if ((k+1) % 5 == 0)
				{
					fprintf(fp, "\n");
				}
			}
			if (k % 5)
			{
				fprintf( fp, "\n" );
			}
		}
	}

	//Density Correction
	fprintf(fp, "RhoDot\n");
	for (k = 0; k < gridPoints+2; ++k)
	{
		fprintf(fp, "\t%-.6e", rhodot[k-1]);
		if ((k+1) % 5 == 0)
		{
			fprintf(fp, "\n");
		}
	}
	if (k % 5)
	{
		fprintf(fp, "\n");
	}

	//PAH Production Rate
	double prodRate[3];

	//PAH Production Rate
	fprintf (fp, "Dimer_ProdRate [mol/m^3s]\n");
	for (k = 0; k < gridPoints+2; ++k)
	{
		ComputeDimerProdRate(Y[k-1], temp[k-1], density[k-1], molarMass, prodRate);
		fprintf (fp, "\t%-.6e", 1000.0 * prodRate[0]);
		if ((k+1) % 5 == 0)
		{
			fprintf (fp, "\n");
		}
	}
	if ((k+1) % 5)
	{
		fprintf (fp, "\n");
	}

	//PAH nbrC2
	fprintf (fp, "Dimer_nbrC2\n");
	for (k = 0; k < gridPoints+2; ++k)
	{
		ComputeDimerProdRate(Y[k-1], temp[k-1], density[k-1], molarMass, prodRate);
		fprintf (fp, "\t%-.6e", prodRate[1] / prodRate[0]);
		if ((k+1) % 5 == 0)
		{
			fprintf (fp, "\n");
		}
	}
	if ((k+1) % 5)
	{
		fprintf (fp, "\n");
	}

	//PAH nbrH
	fprintf (fp, "Dimer_nbrH\n");
	for (k = 0; k < gridPoints+2; ++k)
	{
		ComputeDimerProdRate(Y[k-1], temp[k-1], density[k-1], molarMass, prodRate);
		fprintf (fp, "\t%-.6e", prodRate[2] / prodRate[0]);
		if ((k+1) % 5 == 0)
		{
			fprintf (fp, "\n");
		}
	}
	if ((k+1) % 5)
	{
		fprintf (fp, "\n");
	}
#else
	if (fNucleation)
	{
	  fprintf (fp, "Nucleation\n");
	  for (k = 0; k < gridPoints+2; ++k)
	  {
	    fprintf (fp, "\t%-.6e", LindstedtNucleation(temp[k-1], Y[k-1], density[k-1], molarMass));
	    if ((k+1) % 5 == 0)
	    {
	      fprintf (fp, "\n");
	    }
	  }
	  if ((k+1) % 5)
	  {
	    fprintf (fp, "\n");
	  }
	}
	 
	if (fCoagulation)
	{
	  fprintf (fp, "Coagulation\n");
	  for (k = 0; k < gridPoints+2; ++k)
	  {
	    fprintf (fp, "\t%-.6e", LindstedtCoagulation(moments[k-1], temp[k-1]));
	    if ((k+1) % 5 == 0)
	    {
	      fprintf (fp, "\n");
	    }
	  }
	  if ((k+1) % 5)
	  {
	    fprintf (fp, "\n");
	  }
	}

	if (fSurfaceGrowth)
	{
	  fprintf (fp, "SurfaceGrowth\n");
	  for (k = 0; k < gridPoints+2; ++k)
	  {
	    fprintf (fp, "\t%-.6e", LindstedtGrowth(moments[k-1], temp[k-1], Y[k-1], density[k-1], molarMass));
	    if ((k+1) % 5 == 0)
	    {
	      fprintf (fp, "\n");
	    }
	  }
	  if ((k+1) % 5)
	  {
	    fprintf (fp, "\n");
	  }
	}

	if (fSurfaceOxidation)
	{
	  fprintf (fp, "Oxidation\n");
	  for (k = 0; k < gridPoints+2; ++k)
	  {
	    fprintf (fp, "\t%-.6e", LindstedtOxidation(moments[k-1], temp[k-1], Y[k-1], density[k-1], molarMass));
	    if ((k+1) % 5 == 0)
	    {
	      fprintf (fp, "\n");
	    }
	  }
	  if ((k+1) % 5)
	  {
	    fprintf (fp, "\n");
	  }
	}
#endif
}
#endif

void T1DSoot::PrintRHSSoot( TNewtonPtr bt, T1DFlamePtr flame )
{
  TAdaptiveGridPtr	grid = bt->GetGrid();
  TGridPtr		currentGrid = grid->GetCurrentGrid();
  NodeInfoPtr		nodeInfo = bt->GetNodeInfo();
  int			i, k;
  int			N = currentGrid->GetNGridPoints();
  FILE			*fp = NULL;
  double			dummy;
  int			counter;
	
  counter = ( int ) ( modf( bt->GetNIter()/10.0, &dummy ) * 10.0 );
  sprintf( flame->GetOutFileBuff(), "%srhsSoot%d.dout", flame->GetOutputPath(), counter );
  if ( !( fp = fopen( flame->GetOutFileBuff(), "w") ) ) { 
    cerr << "#warning: unable to open file " << flame->GetOutFileBuff() << NEWL;
    exit(2);
  }
  
  fprintf( fp, "*\n%-12s", "eta" );
  for ( i = 0; i < fNSootMoments; ++i ) {
    fprintf( fp, "\tConv_%-d\tDiff_%-d\tThPh_%-d\tNucl_%-d\tCoag_%-d\tCoaN_%-d\tCond_%-d\tSuGr_%-d\tSuOx_%-d", i, i, i, i, i, i, i, i, i );
  }
  fprintf( fp, "\n" );
  
  for ( k = 0; k < N; ++k ){
    bt->SetNodeInfo( flame, k );
    PrintRHSSoot( flame, nodeInfo, fp );
  }
  fclose( fp );
}

void T1DSoot::PrintRHSSoot( T1DFlamePtr flame, NodeInfoPtr nodeInfo, FILE *fp )
{
  int 			i;
  TFlameNodePtr	flameNode = flame->GetFlameNode();
  Double			temp = flameNode->temp[kCurr];
//   Double			*pahMoments = flameNode->pahMoments;
  Double			*moments = flameNode->moments;
  Double			*Y = flameNode->Y[kCurr];
  Double			*kSoot = fSootRateCoeffs->vec;
  Double			*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
  Double			density = flameNode->mixDensity[kCurr];
  
//   ComputeFractionalMoments( moments );
//   ComputePhi( temp );
  ComputeSootRateCoeffs( kSoot, temp, flame->GetReaction() );
  ComputeCSootStar( kSoot, Y, density, molarMass, flameNode->mixMolarMass[kCurr] );
  
  fprintf( fp, "%-.6e", *nodeInfo->x );
  
  for ( i = 0; i < fNSootMoments; ++i ) {

//     fprintf( fp, "\t%-.6e", SootConvection( i, flame, nodeInfo, TRUE ) );
    
//     if ( fSizeDepDiff ) {
//       fprintf( fp, "\t%-.6e", -SootDiffusion( i, kPhysical, flame, nodeInfo ) );
//     }
//     else {
//       fprintf( fp, "\t%-.6e", -SootDiffusionNew( i, kPhysical, flame, nodeInfo ) );
//     }
//     fprintf( fp, "\t%-.6e", -SootThermoPhoresis( i, kPhysical, flame, nodeInfo ) );
    
//     fprintf( fp, "\t%-.6e", -NucleationNew( i, temp, pahMoments ) );
//     fprintf( fp, "\t%-.6e", -SourceCoagulationNew( i, temp, moments ) );
//     fprintf( fp, "\t%-.6e", -SourceSurfDepCoag( i, moments, Y, temp, density, molarMass ) );
//     fprintf( fp, "\t%-.6e", -SourceCondensationNew( i, temp, pahMoments, moments, Y, density, molarMass ) );
//     fprintf( fp, "\t%-.6e", -SourceSurfGrowthNew( i, moments, Y, density, molarMass ) );
//     fprintf( fp, "\t%-.6e", -SourceSootOxidationNew( i, moments, Y, density, molarMass ) );
  }
  fprintf( fp, "\n" );
}

void T1DSoot::ComputeDiffusivity( T1DFlamePtr flame )
{
  // omega_ij = 1.0

  TFlameNodePtr	flameNode = flame->GetFlameNode();
  Double			temp = flameNode->temp[kCurr];
  static Double	pi = 4.0 * atan( 1.0 );
  static Double	fact = 1.5 * sqrt( 0.5 / pi );
  static Double	dMin = pow( 6.0 / pi * fMolarMassSoot / ( fSootDensity * AVOGADRO ), 1.0 / 3.0 ); // [m]
  static Double	dMin2 = dMin * dMin; // [m^2]
  
  flameNode->diffSoot[kCurr] = fact / flameNode->mixDensity[kCurr] 
    * sqrt( flameNode->mixMolarMass[kCurr] * RGAS * temp ) / AVOGADRO / dMin2;
	
}
//#endif // ZEROD
