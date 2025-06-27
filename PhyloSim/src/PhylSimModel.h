/*
 * model.h
 *
 *  Created on: 20.06.2014
 *      Author: Paul
 *
 * renamed variables after implementation of positive density dependence [Andy]
 * [m_]compStrength -> [m_]nDDStrength
 * [m_]nicheWidth -> [m_]envNicheWidth
 * dens -> negativeDens
 */

#ifndef PHYLSIMMODEL_H_
#define PHYLSIMMODEL_H_

#include "debug.h"

#include <string>

#include "Grid.h"

class PhylSimModel {
public:
  int m_Dispersal;
  int m_X_coordinate;
  int m_Y_coordinate;
  int timeStep;
  // TODO FH I don't see a compelling reason to implement this as
  // pointers, move to normal class fields?
  GlobalEnvironment *m_Global;
  LocalEnvironment *m_Local;

  PhylSimModel(int x, int y, int dispersal, int simulationEnd, double specRate, bool negativeDens, bool positiveDens,
               bool env, bool neutral, bool mort, int mortStrength, bool repro, int dispersalCutoff, int densityCutoff,
               std::string saveLocation, double nDDStrength, double pDDStrength, double envStrength, int fission,
               double redQueen, double redQueenStrength, int protracted, std::vector<double> airmat,
               std::vector<double> soilmat, double nDDNicheWidth, double pDDNicheWidth, double envNicheWidth);

  ~PhylSimModel();

  void get_data();

  void getclimate();

  void update(unsigned int runs);

  void gettraits();

private:
};

#endif /* PHYLSIMMODEL_H_ */
