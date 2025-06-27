/*
 * model.c++
 *
 *  Created on: 20.06.2014
 *      Author: Paul
 *              Betim Musa <musab@informatik.uni-freiburg.de>
 */

#include "./RandomGen.h"
#include <fstream>
#include <iostream>
#include <string>
// #include "Utils/RandomGen.h"
#include "Grid.h"
#include "Individual.h"
#include "PhylSimModel.h"
#include "Species.h"

PhylSimModel::PhylSimModel(int X, int Y, int dispersal, int simulationEnd, double specRate, bool negativeDens,
                           bool positiveDens, bool env, bool neutral, bool mort, int mortalityStrength, bool repro,
                           int dispersalCutoff, int nDensCutoff, int pDensCutoff, std::string saveLocation,
                           double nDDStrength, double pDDStrength, double envStrength, int fission, double redQueen,
                           double redQueenStrength, int protracted, std::vector<double> airmat,
                           std::vector<double> soilmat, double nDDNicheWidth, double pDDNicheWidth,
                           double envNicheWidth) {

#ifdef DEBUG
  std::cout << "Running simulation with \n";
  std::cout << "Dispersal " << dispersal << "; cutoff " << dispersalCutoff << "\n";
  std::cout << "Competition " << negativeDens << "; strength " << nDDStrength << " cutoff " << densityCutoff << "\n";
  std::cout << "Mutualism " << positiveDens << "; strength " << pDDStrength << " cutoff " << densityCutoff << "\n";
  std::cout << "Environment " << env << "; strength " << envStrength << "\n";
  std::cout << "\n---- debug message for development purposes, remove "
               "debug switch in debug.h for turning this off \n\n";
#endif

  if (dispersal == 1) {
    m_Global = new GlobalEnvironment(
        X, Y, dispersal, neutral, negativeDens, positiveDens, env, mort, repro, simulationEnd, specRate,
        dispersalCutoff, nDensCutoff, pDensCutoff, mortalityStrength, nDDStrength, pDDStrength, envStrength, fission,
        redQueen, redQueenStrength, protracted, airmat, soilmat, nDDNicheWidth, pDDNicheWidth, envNicheWidth);
    m_Local = NULL;
  } else if (dispersal == 2 || dispersal == 3) {
    m_Global = NULL;
    m_Local = new LocalEnvironment(
        X, Y, dispersal, neutral, negativeDens, positiveDens, env, mort, repro, simulationEnd, specRate,
        dispersalCutoff, nDensCutoff, pDensCutoff, mortalityStrength, nDDStrength, pDDStrength, envStrength, fission,
        redQueen, redQueenStrength, protracted, airmat, soilmat, nDDNicheWidth, pDDNicheWidth, envNicheWidth);
  }

  timeStep = 0;
  m_Dispersal = dispersal;
  m_X_coordinate = X;
  m_Y_coordinate = Y;
}

PhylSimModel::~PhylSimModel() {
  delete m_Global;
  delete m_Local;
}

void PhylSimModel::update(unsigned int runs) {
  for (unsigned int generation = 1; generation < runs + 1; generation++) {

    // std::cout << "generation :" <<  generation << '/' << runs <<
    // '\n';

#ifdef DEBUG
    if (generation % 1000 == 0) {
      std::cout << "generation :" << generation << '/' << runs << '\n';
    }
#endif

    if (m_Dispersal == 1) {
      m_Global->increaseAge(timeStep + generation);
      m_Global->reproduce(timeStep + generation);
      m_Global->speciation(timeStep + generation);
      //  std::cout << "Disp generation :" <<  generation << '/' << runs
      //  << '\n';
    } else if (m_Dispersal == 3) {
      m_Local->increaseAge(timeStep + generation);
      m_Local->reproduce(timeStep + generation);
      m_Local->speciation(timeStep + generation);
      //   std::cout << "Disp generation :" <<  generation << '/' <<
      //   runs << '\n';
    }
  }
  timeStep += runs;
}

// void model::get_data(){
//
//		std::ofstream test_matrix("..\\test_out.txt");
//		for(int ba=0;ba<Landscape.get_dimensions().first;ba++){
//			test_matrix << '\n';
//
//		for(int bu=0;bu<Landscape.get_dimensions().second;bu++){
//			test_matrix <<
// Landscape.individuals[ba][bu].Species->ID <<',';
//		}}
//
// }

void PhylSimModel::getclimate() {
  int x = m_Local->get_dimensions().first;
  int y = m_Local->get_dimensions().second;
  std::ofstream temperature_matrix("..\\temperature_out.txt");
  for (int ba = 0; ba < x; ba++) {
    temperature_matrix << '\n';

    for (int bu = 0; bu < y; bu++) {
      temperature_matrix << m_Local->m_Environment[ba * y + bu].first << ',';
    }
  }
}

// void model::gettraits(){
//	std::ofstream trait_matrix("..\\trait_out.txt");
//		for(int ba=0;ba<Landscape.get_dimensions().first;ba++){
//			trait_matrix << '\n';
//
//		for(int bu=0;bu<Landscape.get_dimensions().second;bu++){
//			trait_matrix <<
// Landscape.individuals[ba][bu].mean <<',';
//		}}
//
// }
