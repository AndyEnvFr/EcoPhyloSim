/*
 * grid.h
 *
 *  Created on: 20.06.2014
 *      Author: Paul
 *              Betim Musa <musab@informatik.uni-freiburg.de>
 *              Andrea Ingrosso
 *
 * renamed variables after implementation of positive density dependence [Andy]
 * [m_]compStrength -> [m_]nDDStrength
 * [m_]nicheWidth -> [m_]envNicheWidth
 * m_DD -> m_nDD
 * dd -> ndd
 */

#ifndef GRID_H_
#define GRID_H_

#include <utility> // needed for std::pair
#include <vector>

#include "Individual.h"
#include "Phylogeny.h"
#include "RandomGen.h"

/**
 * The landscape is the main component of our world. It is separated
 * into global environment and local environment. When the species are
 * globally independently (their reproduction doesn't depend on the
 * other species) distributed then we choose the global environment. In
 * the other case, when the distribution of the species is dependent to
 * the local environmental conditions, then we choose a local
 * environment.
 */
class Landscape {

public:
int m_Cutoff;
int m_nDensCutoff;
int m_pDensCutoff;
int m_Dispersal_type;
// The size of the landscape.
int m_Xdimensions, m_Ydimensions;
// Counter for all the species that are born.
unsigned long long m_Global_Species_Counter;
unsigned int m_mortalityStrength;

// Indicates whether the species dispersion depends on local
// conditions or not.
bool m_Neutral;
bool m_nDD;
bool m_pDD;
bool m_Env;
// switches to determine wheather fitness should affect mortality or
// reproduction
bool m_mortality;
bool m_reproduction;

// The size of the landscape is defined as m_Xdimensions *
// m_Ydimensions
unsigned int m_LandscapeSize;
int m_KernelSize;

double m_AirTemperature;
double m_SoilMoistureRange;
double m_GradientStep;
double m_Speciation_Rate;

double cellsWithin_N_DensCutoff;
double cellsWithin_P_DensCutoff;
double m_envNicheWidth;
double m_nDDNicheWidth;
double m_pDDNicheWidth;
double m_envStrength;
double m_nDDStrength;
double m_pDDStrength;
int m_fission;
double m_redQueen;
double m_redQueenStrength;
int m_protracted;
std::vector<double> airmat;
std::vector<double> soilmat;

// Change the temperature in the environment by the given magnitude.
void tempChange(int sign, double magnitude);

// Change the moisture  in the environment by the given magnitude.
void moistChange(int sign, double magnitude);

// accessbile for subclasses. Are overwritten there
virtual double calculateRelatedness(int, int, int) { return 0.0; }
virtual void densityUpdate(int, int) {}

Landscape();

Landscape(int xsize, int ysize, int type, bool neutral, bool ndd, bool pdd, bool env, bool mort, bool repro,
  unsigned int simulationEnd, double specRate, int dispersalCutoff, int nDensCutoff, int pDensCutoff,
  unsigned int mortalityStrength, double envStrength, double nDDStrength, double pDDStrength, int fission,
  double redQueen, double redQueenStrength, int protracted, std::vector<double> airmat,
  std::vector<double> soilmat, double envNicheWidth, double nDDNicheWidth, double pDDNicheWidth);
  
  virtual ~Landscape();
  
  // TODO(Betim): Should it really be public?
  Phylogeny m_Phylogeny;

  // A 2-D array of individuals.
  Individual **m_Individuals;
  std::vector<std::pair<double, double>> m_Environment;
  RandomGen m_RandomGenerator;

  unsigned int m_SimulationEnd;

  // Start the reproduction of the species.
  virtual void reproduce(unsigned int generation);

  // Increase the age of each individual by 1.
  void increaseAge(unsigned int generation);

  void speciation(unsigned int generation);

  std::pair<int, int> get_dimensions();
};

class GlobalEnvironment : public Landscape {
public:
  GlobalEnvironment();

  GlobalEnvironment(int xsize, int ysize, int type, bool neutral, bool ndd, bool pdd, bool env, bool mort, bool repro,
                    unsigned int runs, double specRate, int dispersalCutoff, int nDensCutoff, int pDensCutoff,
                    unsigned int mortalityStrength, double envStrength, double nDDStrength, double pDDStrength,
                    int fission, double redQueen, double redQueenStrength, int protracted, std::vector<double> airmat,
                    std::vector<double> soilmat, double envNicheWidth, double nDDNicheWidth, double pDDNicheWidth);

  virtual ~GlobalEnvironment();

  void reproduce(unsigned int generation);
};

class LocalEnvironment : public Landscape {
public:
  LocalEnvironment();

  // for local environment: get densities from neighbors
  double calculateRelatedness(int focus_x, int focus_y, int cutoff) override;
  void densityUpdate(int x, int y) override;


  LocalEnvironment(int xsize, int ysize, int type, bool neutral, bool ndd, bool pdd, bool env, bool mort, bool repro,
                   unsigned int runs, double specRate, int dispersalCutoff, int nDensCutoff, int pDensCutoff,
                   unsigned int mortalityStrength, double envStrength, double nDDStrength, double pDDStrength,
                   int fission, double redQueen, double redQueenStrength, int protracted, std::vector<double> airmat,
                   std::vector<double> soilmat, double envNicheWidth, double nDDNicheWidth, double pDDNicheWidth);

  virtual ~LocalEnvironment();

  void reproduce(unsigned int generation);

  // used to get LocalDensities for positive and negative density dependence
};

#endif /* GRID_H_ */
