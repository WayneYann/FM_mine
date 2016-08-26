#include "FlameMaster.h"
#include "T0DIsoChor.h"

#undef V
#define VS
#undef VSH

#undef LINDSTEDT

//Initialize the part relevant to PAH formation
void TSoot::InitPAH(TInputDataPtr input)
{
  //List of species to participate in nucleation
  AromInd = new int[20];
  AromNbrC2 = new double[20];
  AromNbrH = new double[20];
  AromStick = new double[20];

  int i = 0;
  
  if ((AromInd[i] = input->FindSpecies ("A2-C10H8")) != -1)
  {
    AromNbrC2[i] = 5.0;
    AromNbrH[i] = 8.0;
    AromStick[i] = 2.0e-3;
    i++;
  }
  if ((AromInd[i] = input->FindSpecies ("A2R5-C12H8")) != -1)
  {
    AromNbrC2[i] = 6.0;
    AromNbrH[i] = 8.0;
    AromStick[i] = 4.0e-3;
    i++;
  }
  if ((AromInd[i] = input->FindSpecies ("P2-C12H10")) != -1)
  {
    AromNbrC2[i] = 6.0;
    AromNbrH[i] = 10.0;
    AromStick[i] = 8.5e-3;
    i++;
  }
  if ((AromInd[i] = input->FindSpecies ("A3-C14H10")) != -1)
  {
    AromNbrC2[i] = 7.0;
    AromNbrH[i] = 10.0;
    AromStick[i] = 1.5e-2;
    i++;
  }
  if ((AromInd[i] = input->FindSpecies ("A3R5-C16H10")) != -1)
  {
    AromNbrC2[i] = 8.0;
    AromNbrH[i] = 10.0;
    AromStick[i] = 2.5e-2;
    i++;
    }
  if ((AromInd[i] = input->FindSpecies ("A4-C16H10")) != -1)
  {
    AromNbrC2[i] = 8.0;
    AromNbrH[i] = 10.0;
    AromStick[i] = 2.5e-2;
    i++;
  }
  if ((AromInd[i] = input->FindSpecies ("FLTN-C16H10")) != -1)
  {
    AromNbrC2[i] = 8.0;
    AromNbrH[i] = 10.0;
    AromStick[i] = 2.5e-2;
    i++;
  }
  if ((AromInd[i] = input->FindSpecies ("A4R5-C18H10")) != -1)
  {
    AromNbrC2[i] = 9.0;
    AromNbrH[i] = 10.0;
    AromStick[i] = 3.9e-2;
    i++;
    }
  

  // SD for DLR2012 adding the following PAHs
  /*
  if ((AromInd[i] = input->FindSpecies ("BGHIF-H10C18")) != -1)
  {
    AromNbrC2[i] = 9.0;
    AromNbrH[i] = 10.0;
    AromStick[i] = 3.9e-2;
    i++;
  }
  if ((AromInd[i] = input->FindSpecies ("C18H12")) != -1)
  {
    AromNbrC2[i] = 9.0;
    AromNbrH[i] = 12.0;
    AromStick[i] = 3.9e-2;
    i++;
  }
  if ((AromInd[i] = input->FindSpecies ("BAPYR-H12C20")) != -1)
  {
    AromNbrC2[i] = 10.0;
    AromNbrH[i] = 12.0;
    AromStick[i] = 6.1e-2;
    i++;
  }
  */
  
  // SD for KAUST use the following PAHs
  /*
  if ((AromInd[i] = input->FindSpecies ("A2-H8C10")) != -1)
  {
    AromNbrC2[i] = 5.0;
    AromNbrH[i] = 8.0;
    AromStick[i] = 2.0e-3;
    i++;
  }
  if ((AromInd[i] = input->FindSpecies ("A2R5-H8C12")) != -1)
  {
    AromNbrC2[i] = 6.0;
    AromNbrH[i] = 8.0;
    AromStick[i] = 4.0e-3;
    i++;
  }
  if ((AromInd[i] = input->FindSpecies ("P2-H10C12")) != -1)
  {
    AromNbrC2[i] = 6.0;
    AromNbrH[i] = 10.0;
    AromStick[i] = 8.5e-3;
    i++;
  }
  if ((AromInd[i] = input->FindSpecies ("A3-H10C14")) != -1)
  {
    AromNbrC2[i] = 7.0;
    AromNbrH[i] = 10.0;
    AromStick[i] = 1.5e-2;
    i++;
  }
  if ((AromInd[i] = input->FindSpecies ("A3R5-H10C16")) != -1)
  {
    AromNbrC2[i] = 8.0;
    AromNbrH[i] = 10.0;
    AromStick[i] = 2.5e-2;
    i++;
  }
  if ((AromInd[i] = input->FindSpecies ("A4-H10C16")) != -1)
  {
    AromNbrC2[i] = 8.0;
    AromNbrH[i] = 10.0;
    AromStick[i] = 2.5e-2;
    i++;
  }
  if ((AromInd[i] = input->FindSpecies ("FLTN-H10C16")) != -1)
  {
    AromNbrC2[i] = 8.0;
    AromNbrH[i] = 10.0;
    AromStick[i] = 2.5e-2;
    i++;
  }
  if ((AromInd[i] = input->FindSpecies ("A4R5-H10C18")) != -1)
  {
    AromNbrC2[i] = 9.0;
    AromNbrH[i] = 10.0;
    AromStick[i] = 3.9e-2;
    i++;
  }
  // Different from TheSoot
  
  if ((AromInd[i] = input->FindSpecies ("CHRYSEN-C18H12")) != -1)
  {
    AromNbrC2[i] = 9.0;
    AromNbrH[i] = 12.0;
    AromStick[i] = 3.9e-2;
    i++;
  }
  if ((AromInd[i] = input->FindSpecies ("BAPYR-C20H12")) != -1)
  {
    AromNbrC2[i] = 10.0;
    AromNbrH[i] =12.0;
    AromStick[i] = 6.1e-2;
    i++;
  }
  if ((AromInd[i] = input->FindSpecies ("BEPYREN-C20H12")) != -1)
  {
    AromNbrC2[i] = 10.0;
    AromNbrH[i] = 12.0;
    AromStick[i] = 6.1e-2;
    i++;
  }  
  if ((AromInd[i] = input->FindSpecies ("PERYLEN-C20H12")) != -1)
  {
    AromNbrC2[i] = 10.0;
    AromNbrH[i] = 12.0;
    AromStick[i] = 6.1e-2;
    i++;
    }
  if ((AromInd[i] = input->FindSpecies ("BGHIPER-C22H12")) != -1)
  {
    AromNbrC2[i] = 11.0;
    AromNbrH[i] = 12.0;
    AromStick[i] = 8.7e-2;
    i++;
  }
  if ((AromInd[i] = input->FindSpecies ("ANTHAN-C22H12")) != -1)
  {
    AromNbrC2[i] = 11.0;
    AromNbrH[i] = 12.0;
    AromStick[i] = 8.7e-2;
    i++;
  }
  if ((AromInd[i] = input->FindSpecies ("CORONEN-C24H12")) != -1)
  {
    AromNbrC2[i] = 12.0;
    AromNbrH[i] = 12.0;
    AromStick[i] = 1.2e-1;
    i++;
    }*/
  //SD gradually change sticking coefficients
  
  for(int tempn = 0; tempn < i; tempn++)
    {AromStick[tempn] *= 1.0;
    }
  
   
  nArom = i;

#ifndef LINDSTEDT
  nucl_nbrC2 = 2.0*AromNbrC2[0];
#else
  nucl_nbrC2 = 25.0;
#endif
  nucl_nbrH = 2.0*AromNbrH[0];

  for (int z = 1; z < nArom; z++)
  {
    nucl_nbrC2 = MIN(nucl_nbrC2, 2.0*AromNbrC2[z]);
    nucl_nbrH = MIN(nucl_nbrH, 2.0*AromNbrH[z]);
  }

  nucl_surf = pow(nucl_nbrC2, 2.0/3.0);

  //for (int q = 0; q < nArom; q++)
  //  AromStick[q] *= 0.05;
}


void TSoot::ComputeDimerProdRate(double * Y, double temp, double density, double * molarMass, double * prodRate)
{
#ifndef LINDSTEDT
  //Reset
  prodRate[0] = 0.0;
  prodRate[1] = 0.0;
  prodRate[2] = 0.0;

  //Frenklach path
  for (int i = 0; i < nArom; i++)
  {
    double aromconc = density * MAX(Y[AromInd[i]],1.0e-60) / molarMass[AromInd[i]]; //Mueller
    double wDimer = 0.5 * GetBetaDimer(temp,AromNbrC2[i]) * aromconc * aromconc;

    if (wDimer < 0.0)
      wDimer = 0.0;
    
    //Sticky coeffs
    wDimer *= AromStick[i];
    
    prodRate[0] += wDimer; //Number density
    prodRate[1] += 2.0 * AromNbrC2[i] * wDimer; //nbrC2 density
    prodRate[2] += 2.0 * AromNbrH[i] * wDimer; //nbrH density
  }

  if (prodRate[0] == 0.0)
    prodRate[0] = 1.0e-60;
#endif
}


//---------- Compute PAH  ---------------

void T0DSoot::ComputeDimerParticules(T0DFlamePtr object, double * moments)
{
//   T0DIsoChorPtr	flame = (T0DIsoChorPtr) object;
//   T0DPropertiesPtr props = flame->GetProperties();
//   T0DSpeciesPtr	species = flame->GetSpecies();

//   int firstSpec = flame->GetOffsetFirstSpecies();
//   int fTemp = flame->GetOffsetTemperature(); 
//   int fSootVariables = flame->GetOffsetSoot();
//   double * y = flame->GetSolution()->vec;
//   double * Y = &y[firstSpec];
//   double temp = y[fTemp];

//   double * molarMass = species->GetMolarMass()->vec;
//   double density = props->GetDensity();
//   double viscosity = props->GetMixViscosity();
//   double mixMolarMass = props->GetMixMolarMass();
  
//   TSoot::ComputeDimerParticules(species, temp, Y, density, viscosity, mixMolarMass, moments);
}

void T1DSoot::ComputeDimerParticules(T1DFlamePtr flame, double * moments)
{
  TFlameNodePtr flameNode = flame->GetFlameNode();
  T1DSpeciesPtr species = flame->GetSpecies();

  double * Y = flameNode->Y[kCurr];
  double density = flameNode->mixDensity[kCurr];
  double viscosity = flameNode->mixViscosity[kCurr];
  double mixMolarMass = flameNode->mixMolarMass[kCurr];
  double temp = flameNode->temp[kCurr];

  TSoot::ComputeDimerParticules(species, temp, Y, density, viscosity, mixMolarMass, moments);
}

void T1DSoot::ComputeDimerParticules(T1DFlamePtr flame, double * moments, int gridpoint)
{
  double ** Y = flame->GetMassFracs()->mat;
  double * density = flame->GetProperties()->GetDensity()->vec;
  double * mixMolarMass = flame->GetProperties()->GetMolarMass()->vec;
  double * temp = flame->GetTemperature()->vec;

  T1DSpeciesPtr species = flame->GetSpecies();

  TSoot::ComputeDimerParticules(species, temp[gridpoint], Y[gridpoint], density[gridpoint], 1.0, mixMolarMass[gridpoint], moments);
}

void T1DSoot::ComputeDimerParticules(TSpeciesPtr species, double temp, double * Y, double density, double viscosity, double mixMolarMass, double * moments)
{
  TSoot::ComputeDimerParticules(species, temp, Y, density, viscosity, mixMolarMass, moments);
}

void TSoot::ComputeDimerParticules(TSpeciesPtr species, double temp, double * Y, double density, double viscosity, double mixMolarMass, double * moments)
{
#ifndef LINDSTEDT
  double * molarMass = species->GetMolarMass ()->vec;
  double prodRate[3];

  //Compute production rate of dimers of PAHs
  ComputeDimerProdRate(Y, temp, density, molarMass, prodRate);

  //Set number of carbon atoms
  dimer_rate = prodRate[0];
  if (prodRate[0] != 0.0)
  {
    //Constant nucleation size
    //Conserve dimer mass
    dimer_nbrC2 = nucl_nbrC2;
    dimer_nbrH = nucl_nbrH;
    prodRate[0] = prodRate[1] / dimer_nbrC2;
  }
  else
  {
    dimer_nbrC2 = nucl_nbrC2;
    dimer_nbrH = nucl_nbrH;
  }
  dimer_rate1 = prodRate[0];
  
  //Sink term due to "nucleation"
  double betaN = 0.0;
  if (fNucleation)
    betaN = GetBetaNucl(temp, dimer_nbrC2);

  //Sink term due to "condensation"
  double betaC = 0.0;
  if (fCondensation)
    betaC = GetBetaCond(temp, dimer_nbrC2, moments);

  //Solve quadratic equation for fictive dimer species
  double Delta = betaC * betaC + 4.0 * betaN * prodRate[0];
  if (Delta >= 0.0)
    dimer_conc = (sqrt(Delta) - betaC) / (2.0 * betaN);
  else
  {
    cout << prodRate[1] << "\t" << prodRate[2] << NEWL;
    cout << dimer_nbrC2 << "\t" << dimer_nbrH << NEWL;
    cerr << "betaN: " << betaN << NEWL;
    cerr << "betaC: " << betaC << NEWL;
    cerr << "prodRate: " << prodRate[0] << NEWL;
    cerr << "error: negative Delta: " << Delta << NEWL;
    cerr << moments[0] << endl;
    cerr << moments[1] << endl;
    cerr << moments[2] << endl;
    cerr << moments[3] << endl;
    cerr << moments[4] << endl;
    cerr << moments[5] << endl;
    exit(2);
  }
#else
  dimer_nbrC2 = 25.0;
#endif
}

//Update the production rates of the additional species like (H, OH, ...)
//to account from consumption in the PAH and Soot reactions
void TSoot::UpdateProductionRates(TSpeciesPtr species, TReactionPtr reaction, double * prodRate, double density, double * Y, double temp, double * moments)
{
  if (fSootUpdateProdRates)
  {
    double * molarMass = species->GetMolarMass()->vec;
    double mixMolarMass = 28.0;
#ifndef LINDSTEDT
    if (fNucleation)
    {
      //Frenklach Path
      for (int i = 0; i < nArom; i++)
      {
	double aromconc = density * MAX(Y[AromInd[i]],1.0e-60) / molarMass[AromInd[i]]; //Mueller
	double wDimer = 0.5 * GetBetaDimer(temp, AromNbrC2[i]) * aromconc * aromconc;
	
	//Sticking Coefficients
	wDimer *= AromStick[i];
	
	//Rate of Dimerization
	prodRate[AromInd[i]] -= 2.0 * molarMass[AromInd[i]] * wDimer;
      }
    }
    
//     //Surface Reaction and Oxidation
//     if (fSurfaceGrowth || fSurfaceOxidation)
//     {
//       double * kSoot = fSootRateCoeffs->vec;
//       double * wSoot = fSootReactionRate->vec;

//       ComputeSootRateCoeffs(kSoot, temp, reaction);
//       ComputeCSootStar(kSoot, Y, density, molarMass, mixMolarMass);
//       ComputeSootReactionRates(wSoot, Y, moments, density, molarMass);

//       if (fSurfaceGrowth)
//       {
// 	prodRate[f_OH] += molarMass[f_OH] * (-wSoot[ks1f] + wSoot[ks1b] - wSoot[ks6]);
// 	prodRate[f_H2O] += molarMass[f_H2O] * (+wSoot[ks1f] - wSoot[ks1b]);      
// 	prodRate[f_H] += molarMass[f_H] * (-wSoot[ks2f] + wSoot[ks2b] + wSoot[ks3f] - wSoot[ks3b]);
// 	prodRate[f_H2] += molarMass[f_H2] * (+ wSoot[ks2f] - wSoot[ks2b]);
// 	prodRate[f_C2H2] += molarMass[f_C2H2] * (-wSoot[ks4]);
//       }
//       if (fSurfaceOxidation)
//       {
// 	prodRate[f_O2] += molarMass[f_O2] * (-wSoot[ks5]);
// 	prodRate[f_CO] += molarMass[f_CO] * (+2.0*wSoot[ks5] + wSoot[ks6]);
//       }
//     }
#else
    if (fNucleation)
    {
      prodRate[f_C2H2] -= molarMass[f_C2H2] * (1.0e4*exp(-21100.0/temp)) * (density*Y[f_C2H2]/molarMass[f_C2H2]);
      prodRate[f_H2]   += molarMass[f_H2]   * (1.0e4*exp(-21100.0/temp)) * (density*Y[f_C2H2]/molarMass[f_C2H2]);
    }
      
    if (fSurfaceGrowth)
    {
      double f_S = sqrt(PI*pow(6.0*fMolarMassSoot/(PI*fSootDensity),2.0/3.0)) * pow(moments[1], 1.0/3.0) * pow(moments[0]*AVOGADRO, 1.0/6.0);

      prodRate[f_C2H2] -= molarMass[f_C2H2] * (6.0e3*exp(-12100.0/temp)) * f_S * (density*Y[f_C2H2]/molarMass[f_C2H2]);
      prodRate[f_H2]   += molarMass[f_H2]   * (6.0e3*exp(-12100.0/temp)) * f_S * (density*Y[f_C2H2]/molarMass[f_C2H2]);
    }

    if (fSurfaceOxidation)
    {
      double surf = PI*pow(6.0*moments[1]*fMolarMassSoot/(PI*fSootDensity*moments[0]*AVOGADRO),2.0/3.0)*moments[0]*AVOGADRO;

      prodRate[f_O2] -= 0.5 * molarMass[f_O2] * (1.0e4*sqrt(temp)*exp(-19680.0/temp)) * surf * (density*Y[f_O2]/molarMass[f_O2]);
      prodRate[f_CO] +=       molarMass[f_CO] * (1.0e4*sqrt(temp)*exp(-19680.0/temp)) * surf * (density*Y[f_O2]/molarMass[f_O2]);
    }
#endif
  }
  else
    return;
}


//Density Correction
void TSoot::CalcRhoDot(TSpeciesPtr species, TReactionPtr reaction, double * rhodot, double density, double * Y, double temp, double * moments)
{
  *rhodot = 0.0;

  if (fSootUpdateProdRates)
  {
    double * molarMass = species->GetMolarMass()->vec;

#ifndef LINDSTEDT
    if (fNucleation)
    {
      //Frenklach path
      for (int i = 0; i < nArom; i++)
      {
	double aromconc = density * MAX(Y[AromInd[i]],1.0e-60) / molarMass[AromInd[i]]; //Mueller
	double wDimer = 0.5 * GetBetaDimer(temp,AromNbrC2[i]) * aromconc * aromconc;
	
	//Sticky coeffs
	wDimer *= AromStick[i];
	
	//Rate of Dimerization
	*rhodot -= 2.0 * molarMass[AromInd[i]] * wDimer;
      }
    }
    
    //Surface Reaction and Oxidation
//     if (fSurfaceGrowth || fSurfaceOxidation)
//     {
//       double * kSoot = fSootRateCoeffs->vec;
//       double * wSoot = fSootReactionRate->vec;

//       ComputeSootRateCoeffs(kSoot, temp, reaction);
//       ComputeCSootStar(kSoot, Y, density, molarMass, mixMolarMass);
//       ComputeSootReactionRates(wSoot, Y, moments, density, molarMass);
      
//       if (fSurfaceGrowth)
//       {
// 	*rhodot += molarMass[f_OH] * (-wSoot[ks1f] + wSoot[ks1b] - wSoot[ks6]);
// 	*rhodot += molarMass[f_H2O] * (+wSoot[ks1f] - wSoot[ks1b]);      
// 	*rhodot += molarMass[f_H] * (-wSoot[ks2f] + wSoot[ks2b] + wSoot[ks3f] - wSoot[ks3b]);
// 	*rhodot += molarMass[f_H2] * (+wSoot[ks2f] - wSoot[ks2b]);
// 	*rhodot += molarMass[f_C2H2] * (-wSoot[ks4]);
//       }
      
//       if (fSurfaceOxidation)
//       {
// 	*rhodot += molarMass[f_O2] * (-wSoot[ks5]);
// 	*rhodot += molarMass[f_CO] * (+2.0*wSoot[ks5] + wSoot[ks6]);
//       }
//     }
#else
    if (fNucleation)
    {
      *rhodot -= molarMass[f_C2H2] * (1.0e4*exp(-21100.0/temp)) * (density*Y[f_C2H2]/molarMass[f_C2H2]);
      *rhodot += molarMass[f_H2]   * (1.0e4*exp(-21100.0/temp)) * (density*Y[f_C2H2]/molarMass[f_C2H2]);
    }
      
    if (fSurfaceGrowth)
    {
      double f_S = sqrt(PI*pow(6.0*fMolarMassSoot/(PI*fSootDensity),2.0/3.0)) * pow(moments[1], 1.0/3.0) * pow(moments[0]*AVOGADRO, 1.0/6.0);

      *rhodot -= molarMass[f_C2H2] * (6.0e3*exp(-12100.0/temp)) * f_S * (density*Y[f_C2H2]/molarMass[f_C2H2]);
      *rhodot += molarMass[f_H2]   * (6.0e3*exp(-12100.0/temp)) * f_S * (density*Y[f_C2H2]/molarMass[f_C2H2]);
    }

    if (fSurfaceOxidation)
    {
      double surf = PI*pow(6.0*moments[1]*fMolarMassSoot/(PI*fSootDensity*moments[0]*AVOGADRO),2.0/3.0)*moments[0]*AVOGADRO;

      *rhodot -= 0.5 * molarMass[f_O2] * (1.0e4*sqrt(temp)*exp(-19680.0/temp)) * surf * (density*Y[f_O2]/molarMass[f_O2]);
      *rhodot +=       molarMass[f_CO] * (1.0e4*sqrt(temp)*exp(-19680.0/temp)) * surf * (density*Y[f_O2]/molarMass[f_O2]);
    }
#endif
  }
  else
    return;
}

void TSoot::CalcEnthDot(TSpeciesPtr species, TReactionPtr reaction, double * enthdot, double density, double * Y, double temp, double * moments, double * enthalpy)
{
  *enthdot = 0.0;

  if (fSootUpdateProdRates)
  {
    double * molarMass = species->GetMolarMass()->vec;
    
#ifndef LINDSTEDT
    if (fNucleation)
    {
      //Frenklach path
      for (int i = 0; i < nArom; i++)
      {
	double aromconc = density * MAX(Y[AromInd[i]],1.0e-60) / molarMass[AromInd[i]]; //Mueller
	double wDimer = 0.5 * GetBetaDimer(temp,AromNbrC2[i]) * aromconc * aromconc;
	
	//Sticky coeffs
	wDimer *= AromStick[i];
	
	//Rate of Dimerization
	*enthdot -= 2.0 * molarMass[AromInd[i]] * wDimer * enthalpy[AromInd[i]];
      }
    }

    //Surface Reaction and Oxidation
//     if (fSurfaceGrowth || fSurfaceOxidation)
//     {
//       double * kSoot = fSootRateCoeffs->vec;
//       double * wSoot = fSootReactionRate->vec;

//       ComputeSootRateCoeffs(kSoot, temp, reaction);
//       ComputeCSootStar(kSoot, Y, density, molarMass, mixMolarMass);
//       ComputeSootReactionRates(wSoot, Y, moments, density, molarMass);
      
//       if (fSurfaceGrowth)
//       {
// 	*enthdot += molarMass[f_OH] * (-wSoot[ks1f] + wSoot[ks1b] - wSoot[ks6]) * enthalpy[f_OH];
// 	*enthdot += molarMass[f_H2O] * (+wSoot[ks1f] - wSoot[ks1b]) * enthalpy[f_H2O];      
// 	*enthdot += molarMass[f_H] * (-wSoot[ks2f] + wSoot[ks2b] + wSoot[ks3f] - wSoot[ks3b]) * enthalpy[f_H];
// 	*enthdot += molarMass[f_H2] * (+wSoot[ks2f] - wSoot[ks2b]) * enthalpy[f_H2];
// 	*enthdot += molarMass[f_C2H2] * (-wSoot[ks4]) * enthalpy[f_C2H2];
//       }
      
//       if (fSurfaceOxidation)
//       {
// 	*enthdot += molarMass[f_O2] * (-wSoot[ks5]) * enthalpy[f_O2];
// 	*enthdot += molarMass[f_CO] * (+2.0*wSoot[ks5] + wSoot[ks6]) * enthalpy[f_CO];
//       }
//     }
#else
    if (fNucleation)
    {
      *enthdot -= enthalpy[f_C2H2] * molarMass[f_C2H2] * (1.0e4*exp(-21100.0/temp)) * (density*Y[f_C2H2]/molarMass[f_C2H2]);
      *enthdot += enthalpy[f_H2]   * molarMass[f_H2]   * (1.0e4*exp(-21100.0/temp)) * (density*Y[f_C2H2]/molarMass[f_C2H2]);
    }

    if (fSurfaceGrowth)
    {
      double f_S = sqrt(PI*pow(6.0*fMolarMassSoot/(PI*fSootDensity),2.0/3.0)) * pow(moments[1], 1.0/3.0) * pow(moments[0]*AVOGADRO, 1.0/6.0);

      *enthdot -= enthalpy[f_C2H2] * molarMass[f_C2H2] * (6.0e3*exp(-12100.0/temp)) * f_S * (density*Y[f_C2H2]/molarMass[f_C2H2]);
      *enthdot += enthalpy[f_H2]   * molarMass[f_H2]   * (6.0e3*exp(-12100.0/temp)) * f_S * (density*Y[f_C2H2]/molarMass[f_C2H2]);
    }

    if (fSurfaceOxidation)
    {
      double surf = PI*pow(6.0*moments[1]*fMolarMassSoot/(PI*fSootDensity*moments[0]*AVOGADRO),2.0/3.0)*moments[0]*AVOGADRO;

      *enthdot -= enthalpy[f_O2] * 0.5 * molarMass[f_O2] * (1.0e4*sqrt(temp)*exp(-19680.0/temp)) * surf * (density*Y[f_O2]/molarMass[f_O2]);
      *enthdot += enthalpy[f_CO] *       molarMass[f_CO] * (1.0e4*sqrt(temp)*exp(-19680.0/temp)) * surf * (density*Y[f_O2]/molarMass[f_O2]);
    }
#endif
  }
  else
    return;
}

void TSoot::CalcDimerProductionRate(TSpeciesPtr species, TReactionPtr reaction, double * dimerProdRate, double * dimerNbrC, double * dimerNbrH, double density, double * Y, double temp, double * moments)
{
  *dimerProdRate = 0.0;
  *dimerNbrC = 2.0*dimer_nbrC2;
  *dimerNbrH = dimer_nbrH;

  double * molarMass = species->GetMolarMass()->vec;

  for (int i = 0; i < nArom; i++)
  {
    double aromconc = density * MAX(Y[AromInd[i]],1.0e-60) / molarMass[AromInd[i]]; //Mueller
    double wDimer = 0.5 * GetBetaDimer(temp,AromNbrC2[i]) * aromconc * aromconc;

    if (wDimer < 0.0)
      wDimer = 0.0;

    wDimer *= AromStick[i];

    *dimerProdRate += 2.0 * AromNbrC2[i] * wDimer;
  }

  *dimerProdRate /= nucl_nbrC2;
}


// Compute the term: sqrt(1/i+1/j) * (i^1/3+j^1/3)^2
double TSoot::GetBetaDimer(double temp, double i)
{
  double C = GetC(temp) / 2.2;
  //No van der Waals enhancement factor for PAH

  return C * (4.0 * sqrt(2.0) * pow(i,1.0/6.0));
}

//Compute the term: sqrt(1/i+1/j) * (i^1/3+j^1/3)^2
double TSoot::GetBetaNucl(double temp, double i)
{
  double C = GetC(temp);

  return C * (4.0 * sqrt(2.0) * pow(i,1.0/6.0));
}

//Compute the term: sqrt(1/i+1/j) * (i^1/3+j^1/3)^2
double TSoot::GetBetaCond(double temp, double i, double * moments)
{
#ifdef V
  double C = GetC(temp) / 2.2;

  double S0 =		FracMom( 2.0/3.0, moments) * pow(i, - 1.0/2.0) +
		2.0 * 	FracMom( 1.0/3.0, moments) * pow(i, - 1.0/6.0) + 
			FracMom(     0.0, moments) * pow(i,   1.0/6.0) + 
		0.5 *	FracMom(-1.0/3.0, moments) * pow(i,   1.0/2.0) +
			FracMom(-2.0/3.0, moments) * pow(i,   5.0/6.0) +
		0.5 *	FracMom(-    1.0, moments) * pow(i,   7.0/6.0);
  return C * S0;
#endif
#ifdef VS
  double C_fm = GetC(temp) / 2.2;
  double C_cont = 8.0 * RGAS/AVOGADRO * temp / (3.0*mu) * AVOGADRO;
  double lambda = 3.0*mu/rho*sqrt(PI*Wmix/(8.0*RGAS*temp)) / pow(6.0*fMolarMassSoot/(PI*fSootDensity*AVOGADRO),1.0/3.0);

  double fm;
  double cont;

  double Df = 1.8;
  double av = 1.0 - (2.0 / Df);
  double as = (3.0 / Df) - 1.0;

  double S0v_fm =	FracMom(2.0*av    , 2.0*as, moments) * pow(i, -3.0/6.0) +
		  2.0 * FracMom(    av    ,     as, moments) * pow(i, -1.0/6.0) +
			FracMom(       0.0,    0.0, moments) * pow(i,  1.0/6.0) +
                  0.5 * FracMom(2.0*av-1.0, 2.0*as, moments) * pow(i,  3.0/6.0) +
			FracMom(    av-1.0,     as, moments) * pow(i,  5.0/6.0) +
                  0.5 * FracMom(      -1.0,    0.0, moments) * pow(i,  7.0/6.0);
  return C_fm * S0v_fm;

//   double S0v_cont = 2.000 *          FracMom( 0.0   ,  0.0   , moments) * pow(i,  0.0/3.0) +
//                                      FracMom(     av,      as, moments) * pow(i, -1.0/3.0) +
//                                      FracMom(-    av, -    as, moments) * pow(i,  1.0/3.0) +
//                     1.257 * lambda * FracMom(-    av, -    as, moments) * pow(i,  0.0/3.0) +
//                     1.257 * lambda * FracMom( 0.0   ,  0.0   , moments) * pow(i, -1.0/3.0) +
//                     1.257 * lambda * FracMom(-2.0*av, -2.0*as, moments) * pow(i,  1.0/3.0) +
//                     1.257 * lambda * FracMom(     av,      as, moments) * pow(i, -2.0/3.0);

//   cont = C_cont * S0v_cont;

//   return (fm * cont) / (fm + cont);
#endif
#ifdef VSH
  double C_fm = GetC(temp) / 2.2;

  double Df = 1.8;
  double av = 1.0 - (2.0 / Df);
  double as = (3.0 / Df) - 1.0;

  double S0v_fm =	FracMom(2.0*av    , 2.0*as, 0.0, moments) * pow(i, -3.0/6.0) +
                  2.0 * FracMom(    av    ,     as, 0.0, moments) * pow(i, -1.0/6.0) +
                        FracMom(       0.0,    0.0, 0.0, moments) * pow(i,  1.0/6.0) +
                  0.5 * FracMom(2.0*av-1.0, 2.0*as, 0.0, moments) * pow(i,  3.0/6.0) +
                        FracMom(    av-1.0,     as, 0.0, moments) * pow(i,  5.0/6.0) +
                  0.5 * FracMom(      -1.0,    0.0, 0.0, moments) * pow(i,  7.0/6.0);

  return C_fm * S0v_fm;
#endif
}	

double TSoot::GetC(double temp)
{
  //Units of C are [m^3 / (kmole-s)]

  double V1 = fMolarMassSoot / (AVOGADRO * fSootDensity);
  double A = pow(6.0/PI, 1.0/3.0);

  double C = 2.2 * A * A * sqrt((PI * RGAS * temp) / (2.0 * AVOGADRO * fSootDensity)) * pow(V1, 1.0/6.0);

  return C * AVOGADRO;
}
