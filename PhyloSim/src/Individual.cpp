/*
 * individual.cpp
 *
 *  Created on: 20.06.2014
 *      Author: Paul
 *              Betim Musa
 *              Florian Hartig
 *
 *
 * renamed variables after implementation of positive density dependence [Andy]
 * [m_]compStrength -> [m_]nDDStrength
 * [m_]nicheWidth -> [m_]envNicheWidth
 * dd -> ndd
 */

#include "Individual.h"

#include <cmath>
#include <float.h>
#include <iostream>
#include <stdexcept>

Individual::Individual() {
  this->m_Species = NULL;
  this->m_X_coordinate = -999;
  this->m_Y_coordinate = -999;

  this->m_nLocalDensity = 0.0; // density experienced around itself, will
  // be updated automatically
  this->m_pLocalDensity = 0.0;

  this->m_Age = 0;
  this->m_incip_Age = -99999999;

  // THESE SEEM OBSOLTE ???
  // this ->m_FitnessWeight = 0.5;
  // this ->m_DensityStrength = 0.4;
  // this ->m_Weight = 1.0;
  // END OBSOLETE

  this->m_envNicheWidth = 0.03659906; // environmental niche width, see getFitness
  this->m_Mean = 0.5;                 // environmental trait
  this->m_CompetitionMarker = 0.5;    // competition trait
  this->m_NeutralMarker = 0.5;        // neutral trait

  this->m_envStrength = 1;
  this->m_nDDStrength = 1;
  this->m_pDDStrength = 1;

  this->m_dispersalDistance = 0.0; // parameter for dispersal kernel
}

// COPY CONSTRUCTOR
// better change the overload below to deep copy
Individual::Individual(const Individual &ind) {

  std::cout << "CHECK IF THIS WORKS" << std::endl;

  this->m_Species = ind.m_Species;
  this->m_X_coordinate = ind.m_X_coordinate;
  this->m_Y_coordinate = ind.m_Y_coordinate;
  this->m_nLocalDensity = ind.m_nLocalDensity;
  this->m_pLocalDensity = ind.m_pLocalDensity;
  this->m_Age = 0;
  this->m_incip_Age = -99999999;
  //	this -> m_FitnessWeight = ind.m_FitnessWeight;
  //	this -> m_DensityStrength = ind.m_DensityStrength;
  //	this -> m_Weight = ind.m_Weight;
  this->m_envNicheWidth = ind.m_envNicheWidth;
  this->m_Mean = ind.m_Mean;

  this->m_CompetitionMarker = ind.m_CompetitionMarker;
  this->m_NeutralMarker = ind.m_NeutralMarker;
  this->m_dispersalDistance = ind.m_dispersalDistance;

  this->m_envStrength = ind.m_envStrength;
  this->m_nDDStrength = ind.m_nDDStrength;
  this->m_pDDStrength = ind.m_pDDStrength;
}

Individual::~Individual() {}

// FH - This is a weird construction ... this operator is used
// for creating a new individual, and applying evolution?
// why not create a child function or something like that
// that would be far more logical
void Individual::operator=(const Individual &ind) {
  this->m_Species = ind.m_Species;
  this->m_X_coordinate = -999;
  this->m_Y_coordinate = -999;
  this->m_nLocalDensity = ind.m_nLocalDensity;
  this->m_pLocalDensity = ind.m_pLocalDensity;
  this->m_Age = 0;
  this->m_incip_Age = -99999999;

  //	this -> m_FitnessWeight = ind.m_FitnessWeight;
  //	this -> m_DensityStrength = ind.m_DensityStrength;
  //	this -> m_Weight = ind.m_Weight;

  this->m_envNicheWidth = ind.m_envNicheWidth;

  this->m_Mean = ind.m_Mean;
  this->m_CompetitionMarker = ind.m_CompetitionMarker;
  this->m_NeutralMarker = ind.m_NeutralMarker;
  this->m_dispersalDistance = ind.m_dispersalDistance;

  this->m_dispersalDistance = ind.m_dispersalDistance;

  this->m_envStrength = ind.m_envStrength;
  this->m_nDDStrength = ind.m_nDDStrength;
  this->m_pDDStrength = ind.m_pDDStrength;
}

// Dispersal Kernel
// One would think it's easier to calculate this directly in the
// dispersal Function, but for some reason it seems faster that way

double Individual::kernel(double distance) { return exp(-distance / m_dispersalDistance); }

double Individual::dispersal(int dispersal_type, double distance) {
  if (dispersal_type == 3) // kernel
  {
    // return exp(-distance / cutoff / 2.0) ;
    return kernel(distance);      // for some weird reason, this option is
                                  // considerably faster!!!
  } else if (dispersal_type == 2) // nearest neighbor
  {
    if (distance > 1.0)
      return 0.0;
    else
      return 1.0;
  } else
    throw std::invalid_argument("Problem in the parameters of dispersal function");
}

double Individual::getSeedsTo(int rel_x, int rel_y, int dispersal_type, double temp, bool env, bool ndd, bool pdd,
                              int generation, double redQueenStrength, double redQueen) {
  double dispersal_weight = 0.0;
  dispersal_weight = dispersal(dispersal_type, euclidian_distance(rel_x, rel_y)); // Kernel or NN

  if (env || ndd || pdd) {
    double fitness_weight = getFitness(temp, env, ndd, pdd, generation, redQueenStrength, redQueen);
    return (dispersal_weight * fitness_weight);
  } else {
    return (dispersal_weight);
  }
}

/**
 * gets the fitness of the current individual. Fitness has basline 1 and can range from [0,3]
 * @param temp environmental parameter
 * @param env environment acting
 * @param ndd density acting (favors dissimilar traits)
 * @param pdd positive density dependence (favors similar traits)
 * @param envNicheWidth variance for the environmental fitness kernel
 * @return Fitness
 */
double Individual::getFitness(double temp, bool env, bool ndd, bool pdd, int generation, double redQueenStrength,
                              double redQueen) {
  double out = 1.0; // baseline fitness allows subtraction via ndd

  // changed by andy
  if (env)
    out += m_envStrength * exp(-0.5 * pow((temp - m_Mean) / m_envNicheWidth, 2.0)); // environmental niche

  // localdensity is transformed via gaussian kernel in grid.cpp LocalEnvironment::densityUpdate
  if (ndd)
    out -= m_nDDStrength * m_nLocalDensity;

#ifdef DEBUG_ANDY
  std::cout << "m_nLocalDensity: " << m_nLocalDensity << "  m_nDDStrength: " << m_nDDStrength << "\n";
#endif

  // pdd and ndd are defined likewise. Just pdd is added and ndd is subtracted from fitness
  if (pdd)
    out += m_pDDStrength * m_pLocalDensity;

#ifdef DEBUG_ANDY
  std::cout << "m_pLocalDensity: " << m_pLocalDensity << "  m_pDDStrength: " << m_pDDStrength << "\n";
#endif

  // Implementation of the redQueen Mechanism                      !!! IGNORED FOR CDD AND PDD IMPLEMENTATION !!! needs
  // to be adjusted
  if ((redQueenStrength != 0) || (redQueen != 0)) {

    // The new fitness value is calculated as a function of the specie's
    // age
    out = out +
          (out * redQueenStrength * std::pow(2.71828, (-redQueen * (generation - 1 - m_Species->m_Date_of_Emergence))));
  }
  return out;
}

double Individual::euclidian_distance(int x, int y) { return sqrt((x * x) + (y * y)); }

void Individual::evolve() {

  // if (m_X_coordinate == 0 && m_Y_coordinate == 0) printInfo();

  double width = 0.01;

  double upperBound = 1.0;
  double lowerBound = 0.0;

  double weightSpecies = 0.2;

  // Environment

  m_Mean = (1.0 - weightSpecies) * m_Mean + weightSpecies * m_Species->m_Mean +
           m_RandomGenerator.randomDouble(-width, width);
  if (m_Mean > upperBound)
    m_Mean = upperBound - (m_Mean - upperBound);
  else if (m_Mean < lowerBound)
    m_Mean = lowerBound + std::abs(m_Mean);

  // Competition

  m_CompetitionMarker = (1.0 - weightSpecies) * m_CompetitionMarker + weightSpecies * m_Species->m_CompetitionMean +
                        m_RandomGenerator.randomDouble(-width, width);
  if (m_CompetitionMarker > upperBound)
    m_CompetitionMarker = upperBound - (m_CompetitionMarker - upperBound);
  else if (m_CompetitionMarker < lowerBound)
    m_CompetitionMarker = lowerBound + std::abs(m_CompetitionMarker);

  // Neutral

  m_NeutralMarker = (1.0 - weightSpecies) * m_NeutralMarker + weightSpecies * m_Species->m_NeutralMean +
                    m_RandomGenerator.randomDouble(-width, width);
  if (m_NeutralMarker > upperBound)
    m_NeutralMarker = upperBound - (m_NeutralMarker - upperBound);
  else if (m_NeutralMarker < lowerBound)
    m_NeutralMarker = lowerBound + std::abs(m_NeutralMarker);

  reportBirth(); // ATTENTION: reportBirth is called here!
}

void Individual::evolveDuringSpeciation() {

  m_Age = 0;

  // EVOLUTION DURING SPECIATION

  //    if (false) {
  //
  //        double width = 0.01;
  //
  //        double upperBound = 1.0;
  //        double lowerBound = 0.0;
  //
  //        // Environment
  //
  //        m_Mean += m_RandomGenerator.randomDouble(-width, width);
  //        if (m_Mean > upperBound) m_Mean = upperBound - (m_Mean -
  //        upperBound); else if (m_Mean < lowerBound) m_Mean =
  //        lowerBound + std::abs(m_Mean);
  //
  //        // Competition
  //
  //        m_CompetitionMarker +=
  //        m_RandomGenerator.randomDouble(-width, width); if
  //        (m_CompetitionMarker > upperBound) m_CompetitionMarker =
  //        upperBound - (m_CompetitionMarker - upperBound); else if
  //        (m_CompetitionMarker < lowerBound) m_CompetitionMarker =
  //        lowerBound + std::abs(m_CompetitionMarker);
  //
  //        //Neutral
  //
  //        m_NeutralMarker += m_RandomGenerator.randomDouble(-width,
  //        width); if (m_NeutralMarker > upperBound) m_NeutralMarker =
  //        upperBound - (m_NeutralMarker - upperBound); else if
  //        (m_NeutralMarker < lowerBound) m_NeutralMarker = lowerBound
  //        + std::abs(m_NeutralMarker);
  //
  //    }

  // END EVOLUTION

  m_Species->m_Mean = m_Mean;
  m_Species->m_CompetitionMean = m_CompetitionMarker;
  m_Species->m_NeutralMean = m_NeutralMarker;

  m_Species->m_FirstMean = m_Mean;
  m_Species->m_FirstComp = m_CompetitionMarker;
  m_Species->m_FirstNeutral = m_NeutralMarker;

  reportBirth(); // ATTENTION: reportBirth is called here!
}

// TODO move this in the species class
void Individual::reportDeath(int generation) {
  m_Species->removeIndividual(m_Mean, m_CompetitionMarker, m_NeutralMarker, generation);
}

void Individual::reportBirth() { m_Species->addIndividual(m_Mean, m_CompetitionMarker, m_NeutralMarker); }

void Individual::printInfo() {
  std::cout << "Location: " << m_X_coordinate << m_Y_coordinate << " EnvTrait:" << m_Mean << " Spec " << m_Species->m_ID
            << " mean" << m_Species->m_Mean << " ... \n";
}
