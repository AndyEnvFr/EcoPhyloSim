/*
 * main.cpp
 *
 *  Created on: 24.06.2014
 *      Author: Paul
 *              Betim Musa <musab@informatik.uni-freiburg.de>
 */

#include "./RandomGen.h"
#include "Grid.h"
#include "PhylSimModel.h"
// #include "Utils/RandomGen.h"
// #include "Parameters.h"
#include "iostream"
#include "stdexcept"
#include "string"

int main() {
  /*    unsigned int runs = 100; // number of Generations
      bool neutral = true; // neutral model or not
      bool dd = false; // Density dependent or independent model
      bool env = false; // environmentally dependent or independent
  model bool mort = true; // mortality fitness bool repro = false; //
  reproductive fitness int dispersal = 1; // 1 = global dispersal, 2 =
  nearest neighbor dispersal, 3= kernel dispersal int xDim = 30; //
  Number of grid cells in x-direction int yDim = 30; // number of grid
  cells in y-direction
  //	srand(1500);
      double specrate = 1.0;
      int densityCut = 2;
      int dispersalCut = 2;
      int mortalityStrength = 10;
      std::string saveLoc = "C:/Users/Stefan/Documents/";
      RandomGen ran;
      ran.seedrand(1000);

      double envStrength = 1.0;
      double compStrength = 1.0;
      int fission = 1;
      double redQueen = 0.0;
      double redQueenStrength = 0.0;
      int protracted = 0;
      std::vector<double> airmat(xDim * yDim);
      std::vector<double> soilmat(xDim * yDim);*/

  // Just to test
  // Parameters* pa = new Parameters();
  // std::cout << "Value of numberOfRuns: " <<
  // pa->getParameterValue<int>(std::string("numberOfRuns")) <<
  // std::endl;

  int x = 30;
  int y = 30;
  int dispersal = 1;
  int runs = 100;
  double specRate = 1;
  bool negativeDens = false; // 0
  bool positiveDens = false; // 0
  bool env = false;          // 0
  bool neutral = true;       // 1
  bool mort = true;          // 1
  int mortStrength = 10;
  bool repro = false;
  int dispersalCutoff = 1;
  int densityCutoff = 1;
  int seed = 975;
  double envStrength = 1;
  double pDDStrength = 1;
  double nDDStrength = 1;
  int fission = 0;
  double redQueen = 0;
  double redQueenStrength = 0;
  int protracted = 0;
  std::vector<double> airmat(1);
  std::vector<double> soilmat(1);
  double envNicheWidth = 0.03659906; // initial value, when nicheWidth was fixed
  double nDDNicheWidth = 0.1;
  double pDDNicheWidth = 0.1;
  bool prunePhylogeny = 1;

  RandomGen ran;
  ran.seedrand(seed); // seed is int while seedrand expects unsigned int

  // avoid inconsistent input
  if (neutral && negativeDens)
    throw std::runtime_error("A neutral model can't be density dependent!");
  if (neutral && positiveDens)
    throw std::runtime_error("A neutral model can't be (positive) density dependent!");
  if (neutral && env)
    throw std::runtime_error("A neutral model can't be depending on the environment!");

  std::string saveLoc = "C:/Users/Tankr_000/Documents/";

  // Running the model

  PhylSimModel Model(x, y, dispersal, runs, specRate, negativeDens, positiveDens, env, neutral, mort,
                     mortStrength, repro, dispersalCutoff, densityCutoff, saveLoc, nDDStrength, pDDStrength,
                     envStrength, fission, redQueen, redQueenStrength, protracted, airmat, soilmat, nDDNicheWidth,
                     pDDNicheWidth, envNicheWidth);
  Model.update(runs);

  if (prunePhylogeny) {
    Model.m_Global->m_Phylogeny.prunePhylogeny(runs);
    Model.m_Global->m_Phylogeny.writePhylogenyR(1, Model.m_Global->m_Phylogeny.m_PrunedPhylo);
  }

  Model.m_Global->m_Phylogeny.writePhylogenyR(1, Model.m_Global->m_Phylogeny.m_FullPhylogeny);

  //	Model.get_data();
  //	Model.getclimate();
  //	Model.gettraits();
  /* if (dispersal == 1) {

       Model.m_Global->m_Phylogeny.prunePhylogeny(runs);
       Model.m_Global->m_Phylogeny.writePhylogeny(1,
   Model.m_Global->m_Phylogeny.m_PrunedPhylo, 'P');
       Model.m_Global->m_Phylogeny.writePhylogeny(1,
   Model.m_Global->m_Phylogeny.m_FullPhylogeny, 'F');
       Model.m_Global->m_Phylogeny.writeSpeciesData();
   } else {

       Model.m_Local->m_Phylogeny.prunePhylogeny(runs);
       Model.m_Local->m_Phylogeny.writePhylogeny(1,
   Model.m_Local->m_Phylogeny.m_PrunedPhylo, 'P');
       Model.m_Local->m_Phylogeny.writePhylogeny(1,
   Model.m_Local->m_Phylogeny.m_FullPhylogeny, 'F');
       Model.m_Local->m_Phylogeny.writeSpeciesData();
   }*/

  return 0;
}
