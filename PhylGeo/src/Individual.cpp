/*
 * individual.cpp
 *
 *  Created on: 20.06.2014
 *      Author: Paul
 *              Betim Musa
 *              Florian Hartig
 */

#include "Individual.h"

#include <float.h>
#include <cmath>
#include <stdexcept>
#include <iostream>



Individual::Individual()
{
  this ->m_Species = new Species();
	this ->m_X_coordinate = -999;
	this ->m_Y_coordinate = -999;
  
	this ->m_LocalDensity = 0.0; // density experienced around itself, will be updated automatically 
	this ->m_Age = 0;
  
  // THESE SEEM OBSOLTE ???
	//this ->m_FitnessWeight = 0.5;
	//this ->m_DensityStrength = 0.4;
	//this ->m_Weight = 1.0;
  // END OBSOLETE
  
	this -> m_Variance = 0.03659906;
	this -> m_Mean = 0.5;
	this -> m_CompetitionMarker = 0.5;
	this -> m_NeutralMarker = 0.5;
  
  this -> m_envStrength = 1;
  this -> m_compStrength = 1;
  
	this -> m_dispersalDistance = 0.0; // parameter for dispersal kernel 

}

// COPY CONSTRUCTOR
// better change the overload below to deep copy
Individual::Individual(const Individual &ind)
{
  
  std::cout << "CHEKC IF THIS WORKS";

	this ->m_Species = ind.m_Species;
	this -> m_X_coordinate = ind.m_X_coordinate;
	this -> m_Y_coordinate = ind.m_Y_coordinate;
	this -> m_LocalDensity = ind.m_LocalDensity;
	this -> m_Age = 0;
//	this -> m_FitnessWeight = ind.m_FitnessWeight;
//	this -> m_DensityStrength = ind.m_DensityStrength;
//	this -> m_Weight = ind.m_Weight;
	this -> m_Variance = ind.m_Variance;
	this -> m_Mean = ind.m_Mean;
  
	this -> m_CompetitionMarker = ind.m_CompetitionMarker; 
	this -> m_NeutralMarker = ind.m_NeutralMarker; 
	this -> m_dispersalDistance = ind.m_dispersalDistance;
  
  this -> m_envStrength = ind.m_envStrength;
  this -> m_compStrength = ind.m_compStrength;
}

Individual::~Individual(){

}

// FH - This is a weird construction ... this operator is used 
// for creating a new individual, and applying evolution?
// why not create a child function or something like that
// that would be far more logical
void Individual::operator=(const Individual &ind)
{
	this ->m_Species = ind.m_Species;
	this -> m_X_coordinate = -999;
	this -> m_Y_coordinate = -999;
	this -> m_LocalDensity = ind.m_LocalDensity;
	this -> m_Age = 0;
//	this -> m_FitnessWeight = ind.m_FitnessWeight;
//	this -> m_DensityStrength = ind.m_DensityStrength;
//	this -> m_Weight = ind.m_Weight;

	this -> m_Variance = ind.m_Variance;
  
  this -> m_Mean = ind.m_Mean;
  this -> m_CompetitionMarker = ind.m_CompetitionMarker; 
	this -> m_NeutralMarker = ind.m_NeutralMarker; 
	this -> m_dispersalDistance = ind.m_dispersalDistance;

	this -> m_dispersalDistance = ind.m_dispersalDistance;
  
  this -> m_envStrength = ind.m_envStrength;
  this -> m_compStrength = ind.m_compStrength;
}


  // Dispersal Kernel
  // One would think it's easier to calculate this directly in the dispersal
  // Function, but for some reason it seems faster that way

	double Individual::kernel(double distance){
		return exp(-distance/m_dispersalDistance);
	}



	double Individual::dispersal(int dispersal_type, double distance)
	{
		if (dispersal_type == 3) // kernel
		{
      //return exp(-distance / cutoff / 2.0) ; 
			return kernel(distance); // for some weird reason, this option is considerably faster!!!
		}
    else if (dispersal_type == 2) //nearest neighbor
		{
			if(distance > 1.0) return 0.0;
			else return 1.0;
		}
		else throw std::invalid_argument( "Problem in the parameters of dispersal function" );
	}


	double Individual::getSeedsTo(int rel_x, int rel_y, int dispersal_type, double temp, bool env, bool dd)
	{
    double dispersal_weight = 0.0;
		dispersal_weight = dispersal(dispersal_type, euclidian_distance(rel_x, rel_y)); // Kernel or NN
    
		if(env || dd) {
      double fitness_weight = getFitness(temp, env, dd);
      return(dispersal_weight * fitness_weight);
		}else{
      return(dispersal_weight);
		}
	}
  
   /**
   * gets the fitness of the current individual 
   * @param temp environmental parameter
   * @param env environment acting
   * @param dd density acting 
   * @return Fitness
   */
  double Individual::getFitness(double temp, bool env, bool dd)
	{
    double out = (DBL_MIN*100.0); 
		if(env) out += m_envStrength * 1.2 * exp(-0.5 * pow((temp - m_Mean) / m_Variance, 2.0)) + 1-m_envStrength; // environmental niche
    if(dd) out += m_compStrength * m_LocalDensity + 1- m_compStrength + (DBL_MIN*100.0);
    return out;
	}


	double Individual::euclidian_distance(int x, int y)
	{
		return sqrt((x*x) + (y*y));
	}

	void Individual::evolve()
	{
    
			double width = 0.05;
      
			double upperBound = 1.0;
			double lowerBound = 0.0;
      
      double weightSpecies = 0.0;
      
      // Environment

      std::cout << m_X_coordinate << m_Y_coordinate << m_Mean << "mean" << m_Species->m_Mean << " ... ";
    
      
      m_Mean = (1.0 - weightSpecies) * m_Mean +  weightSpecies * m_Species->m_Mean + m_RandomGenerator.randomDouble(-width, width);
      std::cout << m_Mean;
      if(m_Mean > upperBound) m_Mean = upperBound - (m_Mean - upperBound);
		  else if(m_Mean < lowerBound) m_Mean = lowerBound + std::abs(m_Mean);
      std::cout << m_Mean << "\n";
      
      // Competition
      
      m_CompetitionMarker = (1.0 - weightSpecies) * m_CompetitionMarker +  weightSpecies * m_Species->m_CompetitionMean + m_RandomGenerator.randomDouble(-width, width);
		  if(m_CompetitionMarker > upperBound) m_CompetitionMarker = upperBound - (m_CompetitionMarker - upperBound);
		  else if(m_CompetitionMarker < lowerBound) m_CompetitionMarker = lowerBound + std::abs(m_CompetitionMarker);
      
      //Neutral
        
      m_NeutralMarker = (1.0 - weightSpecies) * m_NeutralMarker +  weightSpecies * m_Species->m_NeutralMean + m_RandomGenerator.randomDouble(-width, width);
		  if(m_NeutralMarker > upperBound) m_NeutralMarker = upperBound - (m_NeutralMarker - upperBound);
		  else if(m_NeutralMarker < lowerBound) m_NeutralMarker = lowerBound + std::abs(m_NeutralMarker);
      
      reportBirth();      
	}
  
  
  void Individual::evolveDuringSpeciation()
	{
    
      m_Age = 0;
    
      // EVOLUTION DURING SPECIATION
    
			double width = 0.05;
      
			double upperBound = 1.0;
			double lowerBound = 0.0;
      
      // Environment
      
      //m_Mean +=  m_RandomGenerator.randomDouble(-width, width);
      if(m_Mean > upperBound) m_Mean = upperBound - (m_Mean - upperBound);
		  else if(m_Mean < lowerBound) m_Mean = lowerBound + std::abs(m_Mean);
      
      // Competition
      
      //m_CompetitionMarker += m_RandomGenerator.randomDouble(-width, width);
		  if(m_CompetitionMarker > upperBound) m_CompetitionMarker = upperBound - (m_CompetitionMarker - upperBound);
		  else if(m_CompetitionMarker < lowerBound) m_CompetitionMarker = lowerBound + std::abs(m_CompetitionMarker);
      
      //Neutral
              
      //m_NeutralMarker += m_RandomGenerator.randomDouble(-width, width);
		  if(m_NeutralMarker > upperBound) m_NeutralMarker = upperBound - (m_NeutralMarker - upperBound);
		  else if(m_NeutralMarker < lowerBound) m_NeutralMarker = lowerBound + std::abs(m_NeutralMarker);
      
      m_Species->m_Mean = m_Mean;
      m_Species->m_CompetitionMean = m_CompetitionMarker;
      m_Species->m_NeutralMean = m_NeutralMarker;
      
      m_Species->m_FirstMean = m_Mean;
      m_Species->m_FirstComp = m_CompetitionMarker;
      m_Species->m_FirstNeutral = m_NeutralMarker;
      
      reportBirth();

	}
  
  // TODO move this in the species class
   // TODO move this in the species class
  void Individual::reportDeath(int generation){
      m_Species->removeIndividual(m_Mean, m_CompetitionMarker, m_NeutralMarker, generation);
  }
  
  void Individual::reportBirth(){
      m_Species->addIndividual(m_Mean, m_CompetitionMarker, m_NeutralMarker);
  } 


