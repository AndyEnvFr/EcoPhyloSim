/*
 * Author: Betim Musa
 * This is the main interface between the phylogeny simulation (written
 * in C++) and R.
 *
 * renamed variables for positive density dependen implementation:
 * compStrength -> nDDStrength
 * NicheWidth -> envNicheWidth
 */

#include <Rcpp.h>

#include "./RandomGen.h"
#include "Grid.h"
#include "Individual.h"
#include "PhylSimModel.h"
#include "Species.h"
#include "debug.h"
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

// basic file operations for DEBUG, remove when done with debugging
#include <fstream>
#include <iostream>

using namespace Rcpp;

//' Core phylosim model
//' @export
// [[Rcpp::export]]
List callModel(int x, int y, int dispersal, IntegerVector runs, double specRate, bool negativeDens, bool positiveDens,
               bool env, bool neutral, bool mort, int mortStrength, bool repro, int dispersalCutoff, int nDensCutoff,
               int pDensCutoff, int seed, double envStrength, double nDDStrength, double pDDStrength, int fission,
               double redQueen, double redQueenStrength, int protracted, NumericVector airmatR, NumericVector soilmatR,
               bool prunePhylogeny, double envNicheWidth, double nDDNicheWidth, double pDDNicheWidth) {
#ifdef DEBUG
  std::ofstream debugFile;
  debugFile.open("debug.txt");
  std::cout << "Writing parameters to debug.txt..." << ";" << std::endl;
  debugFile << "int x = " << x << ";" << std::endl;
  debugFile << "int y = " << y << ";" << std::endl;
  debugFile << "int dispersal = " << dispersal << ";" << std::endl;
  debugFile << "IntegerVector runs = " << runs << ";" << std::endl;
  debugFile << "double specRate = " << specRate << ";" << std::endl;
  debugFile << "bool negativeDens = " << negativeDens << ";" << std::endl;
  debugFile << "bool positiveDens = " << positiveDens << ";" << std::endl;
  debugFile << "bool env = " << env << ";" << std::endl;
  debugFile << "bool neutral = " << neutral << ";" << std::endl;
  debugFile << "bool mort = " << mort << ";" << std::endl;
  debugFile << "int mortStrength = " << mortStrength << ";" << std::endl;
  debugFile << "bool repro = " << repro << ";" << std::endl;
  debugFile << "int dispersalCutoff = " << dispersalCutoff << ";" << std::endl;
  debugFile << "int negative densityCutoff = " << nDensCutoff << ";" << std::endl;
  debugFile << "int positive densityCutoff = " << pDensCutoff << ";" << std::endl;
  debugFile << "int seed = " << seed << ";" << std::endl;
  debugFile << "double envStrength = " << envStrength << ";" << std::endl;
  debugFile << "double nDDStrength = " << nDDStrength << ";" << std::endl;
  debugFile << "double pDDStrength = " << pDDStrength << ";" << std::endl;
  debugFile << "int fission = " << fission << ";" << std::endl;
  debugFile << "double redQueen = " << redQueen << ";" << std::endl;
  debugFile << "double redQueenStrength = " << redQueenStrength << ";" << std::endl;
  debugFile << "int protracted = " << protracted << ";" << std::endl;
  debugFile << "NumericVector airmatR = " << airmatR << ";" << std::endl;
  debugFile << "airmatR = " << airmatR.begin() << ", " << airmatR.end() << ";" << std::endl;
  debugFile << "NumericVector airmatR = " << soilmatR << ";" << std::endl;
  debugFile << "soilmatR = " << soilmatR.begin() << ", " << soilmatR.end() << ";" << std::endl;
  debugFile << "double envNicheWidth = " << envNicheWidth << ";" << std::endl;
  debugFile << "double pDDNicheWidth = " << pDDNicheWidth << ";" << std::endl;
  debugFile << "double nDDNicheWidth = " << nDDNicheWidth << ";" << std::endl;
  debugFile.close();
#endif

  RandomGen ran;
  ran.seedrand(seed); // seed is int while seedrand expects unsigned int

  std::string tempSaveLoc = "./";

  int nRuns = runs.length(); // number of consecutive simulations that
                             // are executed sequentially
  int prevNGen = 0;

  std::vector<double> airmat(airmatR.begin(), airmatR.end());
  std::vector<double> soilmat(soilmatR.begin(), soilmatR.end());

  Rcpp::List outList = Rcpp::List::create();

  PhylSimModel phylSimModel(x, y, dispersal, runs[nRuns - 1], specRate, negativeDens, positiveDens, env, neutral, mort,
                            mortStrength, repro, dispersalCutoff, nDensCutoff, pDensCutoff, tempSaveLoc, nDDStrength,
                            pDDStrength, envStrength, fission, redQueen, redQueenStrength, protracted, airmat, soilmat,
                            nDDNicheWidth, pDDNicheWidth, envNicheWidth);

  for (int step = 0; step < nRuns; step++) {

    Rcpp::IntegerVector specOut(x * y);
    Rcpp::NumericVector traitOut(x * y);
    Rcpp::NumericVector envOut(x * y);
    Rcpp::NumericVector neutralOut(x * y);
    Rcpp::NumericVector compOut(x * y);
    Rcpp::CharacterVector phyloOut = "";

    int curNGen = runs[step] - prevNGen;

#ifdef DEBUG
    std::cout << "Run model for" << curNGen << " generations\n";
#endif
    phylSimModel.update(curNGen);

    prevNGen = runs[step];

    int indCounter = 0;

    if (dispersal == 1) {
      for (int i = 0; i < x; i++) {
        for (int j = 0; j < y; j++) {
          specOut[indCounter] = phylSimModel.m_Global->m_Individuals[i][j].m_Species->m_ID;
          traitOut[indCounter] = phylSimModel.m_Global->m_Individuals[i][j].m_Mean;
          envOut[indCounter] = phylSimModel.m_Global->m_Environment[i * y + j].first;
          neutralOut[indCounter] = phylSimModel.m_Global->m_Individuals[i][j].m_NeutralMarker;
          compOut[indCounter] = phylSimModel.m_Global->m_Individuals[i][j].m_CompetitionMarker;
          indCounter++;
        }
      }
#ifndef PHYL_OFF
      phylSimModel.m_Global->m_Phylogeny.prunePhylogeny(prevNGen);
      std::string phyloPass("\0");
#endif

#ifdef DEBUG
      std::cout << "before writePhylogenyR" << std::endl;
      std::cout << prunePhylogeny << std::endl;
      // consistency checks -> length, etc
#endif

#ifndef PHYL_OFF
      if (prunePhylogeny) {
        phyloPass =
            phylSimModel.m_Global->m_Phylogeny.writePhylogenyR(1, phylSimModel.m_Global->m_Phylogeny.m_PrunedPhylo);
      } else {
        phyloPass =
            phylSimModel.m_Global->m_Phylogeny.writePhylogenyR(1, phylSimModel.m_Global->m_Phylogeny.m_FullPhylogeny);
      }
#endif

#ifdef DEBUG
      std::cout << "after writePhylogenyR" << std::endl;
#endif

#ifndef PHYL_OFF
      char *cstr = new char[phyloPass.length() + 1];
      std::strcpy(cstr, phyloPass.c_str());
      phyloOut[0] = cstr;
      delete[] cstr;
#endif

    } else {
      for (int i = 0; i < x; i++) {
        for (int j = 0; j < y; j++) {
          specOut[indCounter] = phylSimModel.m_Local->m_Individuals[i][j].m_Species->m_ID;
          traitOut[indCounter] = phylSimModel.m_Local->m_Individuals[i][j].m_Mean;
          envOut[indCounter] = phylSimModel.m_Local->m_Environment[i * y + j].first;
          neutralOut[indCounter] = phylSimModel.m_Local->m_Individuals[i][j].m_NeutralMarker;
          compOut[indCounter] = phylSimModel.m_Local->m_Individuals[i][j].m_CompetitionMarker;
          indCounter++;
        }
      }

#ifndef PHYL_OFF
      phylSimModel.m_Local->m_Phylogeny.prunePhylogeny(prevNGen);
      std::string phyloPass("\0");
#endif

#ifndef PHYL_OFF
      if (prunePhylogeny) {
        phyloPass =
            phylSimModel.m_Local->m_Phylogeny.writePhylogenyR(1, phylSimModel.m_Local->m_Phylogeny.m_PrunedPhylo);
      } else {
        phyloPass =
            phylSimModel.m_Local->m_Phylogeny.writePhylogenyR(1, phylSimModel.m_Local->m_Phylogeny.m_FullPhylogeny);
      }
#endif

#ifndef PHYL_OFF
      char *cstr = new char[phyloPass.length() + 1];
      std::strcpy(cstr, phyloPass.c_str());
      phyloOut[0] = cstr;
      delete[] cstr;
#endif
    }

    Rcpp::List listResults =
        Rcpp::List::create(Rcpp::Named("Species") = specOut, Rcpp::Named("EnvTrait") = traitOut,
                           Rcpp::Named("NeutralTrait") = neutralOut, Rcpp::Named("CompetitionTrait") = compOut,
                           Rcpp::Named("Environment") = envOut, Rcpp::Named("Phylogeny") = phyloOut);

    outList.push_back(listResults);
  }

  return outList;
}
