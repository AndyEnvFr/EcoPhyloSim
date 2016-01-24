/*
 * gridMod.cpp
 *
 *  Created on: 12.01.2015
 *      Author: Paul
 *              Betim Musa <musab@informatik.uni-freiburg.de>
 */


#include <cfloat>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>
#include <numeric>
#include <math.h>
#include<algorithm>

#include "Grid.h"
#include "Individual.h"
#include "Species.h"
#include "debug.h"
//#include "./Utils/RandomGen.h"
#include "./RandomGen.h"


//Landscape::Landscape()
//{
//  m_Cutoff= 4;
//  m_Neutral = false;
//	m_DD = true;
//	m_Env = true;
//	m_DensityCutoff = 2;
//	m_Runs= 100;
//	m_Dispersal_type = 3;
//	m_Global_Species_Counter = 1;
//	m_Speciation_Rate = 2.0;
//	m_Xdimensions = 256;
//	m_Ydimensions = 256;
//	m_LandscapeSize = m_Xdimensions * m_Ydimensions;
//	m_RandomGenerator.seedrand(1500);
//
//
//
//
//	this->m_Individuals = new Individual * [m_Xdimensions];
//	for(int cols = 0; cols < m_Xdimensions;cols++)
//	{
//	   m_Individuals[cols] = new Individual [m_Ydimensions];
//	}
//
//	initialize(m_Xdimensions,m_Ydimensions, m_Runs);
//
//	m_AirTemperature =  0.0; // Celsius
//	m_GradientStep = 0.0078125;
//	m_SoilMoistureRange = 101.0; // percent
//
//	for(int i = 0; i < m_Xdimensions; i++)
//	{
//		for(int j = 0; j < m_Ydimensions; j++)
//		{
//			doubleCell envi;
//			envi.first = m_AirTemperature;
//			envi.second = 1.0;
//			m_Environment.push_back(envi);
//		}
//		if(i < m_Xdimensions / 2) m_AirTemperature += m_GradientStep;
//		else m_AirTemperature -= m_GradientStep;
//	}
//}


Landscape::Landscape(int xsize, int ysize, int type, bool neutral, bool dd, bool env, 
bool mort, bool repro, unsigned int runs, double specRate, int dispersalCutoff, 
int DensityCutoff, unsigned int mortalityStrength, 
double envStrength, double compStrength, int fission, double redQueen, double redQueenStrength, int protracted,
std::vector<double> airmat, std::vector<double> soilmat)
{

   m_Cutoff= dispersalCutoff;
   m_Neutral = neutral;
   m_DD = dd;
   m_Env = env;
   m_mortality = mort;
   m_reproduction = repro;
   m_SimulationEnd = runs;
   m_DensityCutoff = DensityCutoff;
   m_Dispersal_type = type;
   m_Global_Species_Counter = 1;
   m_Speciation_Rate = specRate;
   m_LandscapeSize = xsize * ysize;
   m_Xdimensions = xsize;
   m_Ydimensions = ysize;
   m_mortalityStrength = mortalityStrength;
   m_envStrength = envStrength;
   m_compStrength = compStrength;
   m_fission = fission;
   m_redQueen = redQueen;
   m_redQueenStrength = redQueenStrength;
   m_protracted = protracted;

   //	func.seedrand(1500);

   // Construct grid of individuals
   this ->m_Individuals = new Individual * [m_Xdimensions];
   for(int cols = 0; cols < m_Xdimensions;cols++)
   {
      m_Individuals[cols] = new Individual[m_Ydimensions];
   }

   // Initialization
   
   // Create new species 
   Species * spec = new Species(1,1, 0, std::make_pair<int,int>(0,0), m_SimulationEnd);

   // Add individuas of this new species across the grid
   for(int cols = 0; cols < m_Xdimensions; cols++)
   {
      for(int rows=0; rows < m_Ydimensions; rows++)
      {
         this ->m_Individuals[cols][rows].m_Species = spec;
         this ->m_Individuals[cols][rows].m_X_coordinate = cols;
         this->m_Individuals[cols][rows].m_Y_coordinate =rows;
         this->m_Individuals[cols][rows].reportBirth();
         this->m_Individuals[cols][rows].m_dispersalDistance = m_Cutoff / 2.0;       
         this->m_Individuals[cols][rows].m_envStrength = m_envStrength;
         this->m_Individuals[cols][rows].m_compStrength = m_compStrength;      
         
         //this->individuals[cols][rows].Species->date_of_extinction = runs;
      }
   }
   m_Phylogeny.updatePhylogeny(m_Individuals[xsize-1][ysize-1].m_Species);


  // Set up Environment

   if((airmat.size() == 1) && (soilmat.size() == 1)){
   m_AirTemperature =  0.0; // Celsius
   m_GradientStep = (1.0/(double)m_Xdimensions)*2.0;
   m_SoilMoistureRange = 101.0; // percent

   for(int i = 0; i < m_Xdimensions; i++)
   {
      for(int j = 0; j < m_Ydimensions; j++)
      {
         std::pair<double, double> envi;
         envi.first = m_AirTemperature;
         envi.second = 1.0;
         m_Environment.push_back(envi);
      }
      if(i < m_Xdimensions / 2.0) m_AirTemperature += m_GradientStep;
      else m_AirTemperature -= m_GradientStep;
   }
   }

   if(airmat.size() > 1){
	   for(int i = 0; i < m_Xdimensions; i++)
	      {
	         for(int j = 0; j < m_Ydimensions; j++)
	         {
	            std::pair<double, double> envi;
	            envi.first = airmat[i*m_Xdimensions + j];
	            if(soilmat.size() > 1){
	            envi.second = soilmat[i*m_Xdimensions + j];
	            } else envi.second = 1.0;
	            m_Environment.push_back(envi);
	         }
	      }
   }
   
   // Grid Geometry calculations 
   
   cellsWithinDensityCutoff = 0.0 ;

   for(int i= - m_DensityCutoff; i <= m_DensityCutoff; i++)
   {
      int yLim = floor(sqrt(m_DensityCutoff* m_DensityCutoff - i*i )); // avoid diagonal bias
      for(int j = - yLim ; j <=  yLim ; j++)
      {
        cellsWithinDensityCutoff += 1.0 ;
      }
   } 
   
}

Landscape::~Landscape()
{
   for(int i = 0; i < m_Xdimensions; i++)
   {
      delete [] m_Individuals[i];
   }
   delete[] m_Individuals;
}

void Landscape::reproduce(unsigned int generation) {
   // (Betim): This method is needed in order to use polymorphism method calls.
   // (e.g. Model->updateEnvironment, the argument Landscape needs the method "reproduce"
   // to be defined, although it is void.)
}

std::pair<int, int> Landscape::get_dimensions()
{
   std::pair<int,int>(dimensions);
   dimensions = std::make_pair(m_Xdimensions, m_Ydimensions);
   return dimensions;
}


void Landscape::increaseAge(unsigned int generation)
{
   for (int rows = 0; rows < this->m_Xdimensions; rows++)
   {
      for (int cols =0; cols < this->m_Ydimensions; cols++)
      {
         this->m_Individuals[rows][cols].m_Age += 1;
         this->m_Individuals[rows][cols].m_incip_Age += 1;

      }
   }
}

void Landscape::tempChange(int sign, double magnitude)
{
   for(int i = 0; i < m_Xdimensions; i++)
   {
      for(int j = 0; j < m_Ydimensions; j++)
      {
         m_Environment[i * m_Ydimensions + j].first = m_Environment[i * m_Ydimensions + j].first * sign * magnitude;
      }
   }
}

void Landscape::moistChange(int sign, double magnitude)
{
   for(int i = 0; i < m_Xdimensions; i++)
   {
      for(int j = 0; j < m_Ydimensions; j++)
      {
         m_Environment[i * m_Ydimensions + j].second = m_Environment[i * m_Ydimensions + j].second * sign * magnitude;
      }
   }
}


GlobalEnvironment::GlobalEnvironment(int xsize, int ysize, int type, bool neutral, bool dd, bool env, bool mort, bool repro, unsigned int runs, double specRate, int dispersalCutoff, int densityCutoff, unsigned int mortalityStrength,double envStrength, double compStrength, int fission, double redQueen, double redQueenStrength, int protracted, std::vector<double> airmat, std::vector<double> soilmat) :
	                              Landscape(xsize,  ysize,  type,  neutral,  dd,  env, mort, repro,  runs, specRate, dispersalCutoff, densityCutoff, mortalityStrength, envStrength, compStrength, fission,  redQueen, redQueenStrength, protracted, airmat, soilmat)
{

}

GlobalEnvironment::~GlobalEnvironment() {
   // TODO: Implement destructor!
}

void GlobalEnvironment::reproduce(unsigned int generation)
{
  
  /////////////////////////////////////////////////////
  // NEUTRAL CASE
  
   if(m_Neutral && (m_redQueen==0) && (m_redQueenStrength==0))
   {
      #ifdef DEBUG
      std::cout<<"In global neutral \n";
      #endif
     
      int x_coordinate, y_coordinate, x_parent, y_parent;

      //			srand(time(0));
      for(unsigned int event = 0; event < m_LandscapeSize; event++)
      {

         x_coordinate =  m_RandomGenerator.randomInt(0,m_Xdimensions-1); 
         y_coordinate = m_RandomGenerator.randomInt(0,m_Ydimensions-1); 

         x_parent = m_RandomGenerator.randomInt(0,m_Xdimensions-1); 
         y_parent = m_RandomGenerator.randomInt(0,m_Ydimensions-1); 
         
         //				std::cout << x_coordinate << "::" << x_parent << "::" << y_coordinate << "::" << y_parent << '\n';
         
         m_Individuals[x_coordinate][y_coordinate].reportDeath(generation);

         m_Individuals[x_coordinate][y_coordinate] = m_Individuals[x_parent][y_parent]; // overloaded operator, deep copy ! // report birth and evolve is automatic
         m_Individuals[x_coordinate][y_coordinate].m_X_coordinate = x_coordinate;
         m_Individuals[x_coordinate][y_coordinate].m_Y_coordinate = y_coordinate;
         m_Individuals[x_coordinate][y_coordinate].evolve();
      }
   }
   
  /////////////////////////////////////////////////////
  // NON-NEUTRAL CASE
   
   else // Density dependence and / or Environmental dependence
   {

      #ifdef DEBUG
      std::cout<<"In global non-neutral \n";
      #endif
     
      int x_parent = 0;
      int y_parent = 0;
      
      unsigned int numberOfRuns = 0;
      double seedSum = 0.0;
      int new_parent = 0;
      unsigned int array_length;
      
//      double averageCompetitionTrait = 0;
//      double spread = 0;
      
      std::vector<std::pair<int, int> > parents (m_LandscapeSize);

      double weights[m_LandscapeSize];
      double cumWeights[m_LandscapeSize] ;


      int count = 0;
      for(int x = 0; x < m_Xdimensions; x++)
      {
         for(int y = 0 ; y <  m_Ydimensions ; y++)
         {
            parents[count].first = x;
            parents[count].second = y;
            count ++;
         }
      }
 
      // IF FITNESS ACTS ON REPRODUCTION, calculate average fitness for lottery competition
      if(m_reproduction){
        
         #ifdef DEBUG
         std::cout<<"In global non-neutral, reproduction fitness \n";
         #endif
        
        array_length = 0;
        for(int kernel_x = 0; kernel_x < m_Xdimensions; kernel_x++)
        {
           for(int kernel_y = 0 ; kernel_y <  m_Ydimensions ; kernel_y++)
           {
              weights[array_length] = m_Individuals[kernel_x][kernel_y].getFitness(m_Environment[kernel_x * m_Ydimensions + kernel_y].first, m_Env, m_DD, generation, m_redQueenStrength, m_redQueen);
              seedSum += weights[array_length];
              array_length ++;
           }
        }
        
        // TODO this seems inefficient, can just calculate cumweights directly
                   
        cumWeights[0] = weights[0] ;
        for(unsigned int i =1; i< m_LandscapeSize; i++)
        {
           cumWeights[i] = weights[i] + cumWeights[i-1];

        }
        
      }
      
      //if(m_mortality) numberOfRuns = m_LandscapeSize*2;
      //else numberOfRuns = m_LandscapeSize;
      // TODO could make numberofruns a parameter 
      
      numberOfRuns = m_LandscapeSize;

      int event = 0;
      int numberDeath = 0;
      
      // TODO - seems to me the next for loop is not yet corrected if fitness goes on reproduction, 
      
      while(numberDeath < numberOfRuns)
      {
         event++;

         int x_coordinate = m_RandomGenerator.randomInt(0,m_Xdimensions-1);
         int y_coordinate = m_RandomGenerator.randomInt(0,m_Ydimensions-1);

         if(event % m_mortalityStrength != 0 && m_mortality){
            double weight = m_Individuals[x_coordinate][y_coordinate].getFitness(m_Environment[x_coordinate * m_Ydimensions + y_coordinate].first, m_Env, m_DD,  generation, m_redQueenStrength, m_redQueen);
            // important!! the frequency in relation to the base mortality controls the intensity of the mechanisms
        	  double chanceOfDeath = m_RandomGenerator.randomDouble(0.0,1.0);
            //std::cout<<weight<<"\n"; // DEBUG
            if(weight > chanceOfDeath){
              continue;
            }
          }
          numberDeath++;

         m_Individuals[x_coordinate][y_coordinate].reportDeath(generation);

         if(m_reproduction){
            new_parent = m_RandomGenerator.multinomialDraw(cumWeights, m_LandscapeSize-1, seedSum); 
         }else{
           new_parent = m_RandomGenerator.randomInt(0, m_LandscapeSize-1);
         } 
         
         x_parent = parents[new_parent].first;
         y_parent = parents[new_parent].second;

         m_Individuals[x_coordinate][y_coordinate] = m_Individuals[x_parent][y_parent];  // overloaded operator, deep copy ! // report birth and evolve is automatic
         m_Individuals[x_coordinate][y_coordinate].m_X_coordinate = x_coordinate;
         m_Individuals[x_coordinate][y_coordinate].m_Y_coordinate = y_coordinate;
         
         m_Individuals[x_coordinate][y_coordinate].evolve();
         
         parents[x_coordinate*m_Ydimensions + y_coordinate].first = x_coordinate;
         parents[x_coordinate*m_Ydimensions + y_coordinate].second = y_coordinate;
         
  

         // start reproduction updates
         // this seems to be suboptimally programmed at the moment
         // below, the individual class functions should be used for calculating values
         // need to chek if other speedups are possible
         // However, I left it for now because we dont use it for the mortality 
         
         
      // UPDATE RELATEDNESS
      // texperimental global relatedness calculations
      
      
            
//         m_Individuals[x_coordinate][y_coordinate].m_LocalDensity = spread + std::abs(averageCompetitionTrait - m_Individuals[x_coordinate][y_coordinate].m_CompetitionMarker);


//      if(!m_reproduction && m_DD && event % 1000 == 0 )
//      {
//        averageCompetitionTrait = 0;
//
//        for(int x = 0; x < m_Xdimensions; x++)
//        {
//           for(int y = 0 ; y <  m_Ydimensions ; y++)
//           {
//              averageCompetitionTrait += m_Individuals[x][y].m_CompetitionMarker;
//           }
//        }
//        averageCompetitionTrait = averageCompetitionTrait / (double) m_LandscapeSize ;
//
//        spread = 0;
//        for(int x = 0; x < m_Xdimensions; x++)
//        {
//           for(int y = 0 ; y <  m_Ydimensions ; y++)
//           {
//              spread += std::abs(m_Individuals[x][y].m_CompetitionMarker -averageCompetitionTrait);
//           }
//        }
//        spread = spread / averageCompetitionTrait;
//
//        for(int x = 0; x < m_Xdimensions; x++)
//        {
//           for(int y = 0 ; y <  m_Ydimensions ; y++)
//           {
//              m_Individuals[x][y].m_LocalDensity = spread + std::abs(averageCompetitionTrait - m_Individuals[x][y].m_CompetitionMarker);
//           }
//        }
//      }
      
      
      if(!m_reproduction && m_DD) densityUpdate(x_coordinate,y_coordinate);

      
         
        // only if env and dens affect reproduction
        if(m_reproduction){
          
          std::cout<<"REPRODUCTION AT THE MOMENT NOT IMPLEMENTED";
          
//        #ifdef DEBUG
//        std::cout<<"In global non-neutral, reproduction fitness \n";
//        #endif
//
//           if(!(m_DD))
//           {
//              double newWeight = 0.0;
//              
//              double envFitnessParent = 1.2 * exp(-0.5 * pow((m_Environment[x_coordinate * m_Ydimensions + y_coordinate].first - m_Individuals[x_coordinate][y_coordinate].m_Mean) / m_Individuals[x_coordinate][y_coordinate].m_Variance, 2.0));
//              //						double envFitnessPropagule = (1.0 / (individuals[kernel_x][kernel_y].variance * sqrt(2.0 * 3.147))) * exp(-0.5 * pow((environment[x_coordinate][y_coordinate].first - individuals[kernel_x][kernel_y].mean) / individuals[kernel_x][kernel_y].variance, 2.0)); // environmental influence !
//              newWeight = envFitnessParent  + (DBL_MIN*100.0) ; //weights plus base value
//  
//              //						double envFitnessParent = (1.0 / (individuals[x_coordinate][y_coordinate].variance * sqrt(2.0 * 3.147))) * exp(-0.5 * pow((environment[x_coordinate][y_coordinate].first - individuals[x_coordinate][y_coordinate].mean) / individuals[x_coordinate][y_coordinate].variance, 2.0)); // environmental influence !
//              //						double envFitnessPropagule = (1.0 / (individuals[x_coordinate][y_coordinate].variance * sqrt(2.0 * 3.147))) * exp(-0.5 * pow((environment[x_coordinate][y_coordinate].first - individuals[x_coordinate][y_coordinate].mean) / individuals[x_coordinate][y_coordinate].variance, 2.0)); // environmental influence !
//              //						newWeight = envFitnessParent * envFitnessPropagule + (DBL_MIN*100.0);
//  
//              unsigned int vecPos = x_coordinate*m_Ydimensions + y_coordinate;
//              double oldWeight = weights[vecPos];
//              weights[vecPos] = newWeight;
//  
//              double diff = newWeight - oldWeight; // recalculate accumulated weights
//              seedSum += diff;
//              
//              // TODO make debug around this 
//              //if (oldWeight < 0.)std::cout << oldWeight << '\n';
//              //if (newWeight < 0.)std::cout << newWeight << '\n';
//  
//              if(vecPos == 0)
//              {
//                 cumWeights[0] = weights[0];
//              }
//              else
//              {
//                 cumWeights[vecPos] = cumWeights[vecPos -1] + newWeight;
//              }
//  
//              for (unsigned int i = vecPos+1; i < m_LandscapeSize; i++)
//              {
//                 cumWeights[i] += diff;
//                 //						std::cout << cumWeights[i] << '\n';
//              }
//              cumWeights[0] = weights[0] ;
//              for(unsigned int i =1; i< m_LandscapeSize; i++)
//              {
//                 cumWeights[i] = weights[i] + cumWeights[i-1];
//              }
//           }
//  
//           else if(m_DD)
//           {
//              int  densityKernel_x, densityKernel_y, focus_x, focus_y;
//              double relatedness = 0.0;
//              double cells = double((m_DensityCutoff * 2 +1) * (m_DensityCutoff * 2 +1) - 1) ;
//  
//              for(int relativeX= - m_DensityCutoff; relativeX <= m_DensityCutoff; relativeX++)
//              {
//                 for(int relativeY = - m_DensityCutoff ; relativeY <=  m_DensityCutoff ; relativeY++)
//                 {
//                    // get focus cell
//                    focus_x = ((x_coordinate + relativeX + m_Xdimensions) % m_Xdimensions);
//                    focus_y = ((y_coordinate + relativeY + m_Ydimensions) % m_Ydimensions);
//  
//                    for(int relativeX2= - m_DensityCutoff; relativeX2 <= m_DensityCutoff; relativeX2++)
//                    {
//                       for(int relativeY2 = - m_DensityCutoff ; relativeY2 <=  m_DensityCutoff ; relativeY2++)
//                       {
//                          densityKernel_x = ((focus_x + relativeX2 + m_Xdimensions) % m_Xdimensions);
//                          densityKernel_y = ((focus_y + relativeY2 + m_Ydimensions) % m_Ydimensions);
//  
//                          if (!(densityKernel_x == focus_x && densityKernel_y == focus_y))
//                          {
//                             relatedness += std::abs(m_Individuals[focus_x][focus_y].m_CompetitionMarker - m_Individuals[densityKernel_x][densityKernel_y].m_CompetitionMarker);
//                          }
//                       }
//                    }
//                    m_Individuals[focus_x][focus_y].m_LocalDensity = relatedness / cells;
//                    //std::cout<<m_Individuals[focus_x][focus_y].m_LocalDensity<<"\n"; // DEBUG
//  
//                    // TODO  THIS IS PROBALY USELESS IF MORTALITY FITNESS 
//                    if(m_Env)
//                    {
//                       double envFitnessParent = 1.2 * exp((-0.5 * (m_Environment[focus_x * m_Ydimensions + focus_y].first - m_Individuals[focus_x][focus_y].m_Mean) / m_Individuals[focus_x][focus_y].m_Variance * (m_Environment[focus_x * m_Ydimensions + focus_y].first - m_Individuals[focus_x][focus_y].m_Mean) / m_Individuals[focus_x][focus_y].m_Variance)); // environmental influence !
//                       //								double envFitnessPropagule = (1.0 / (individuals[focus_x][focus_y].variance * sqrt(2.0 * 3.147))) * exp((-0.5 * (environment[x_coordinate][y_coordinate].first - individuals[focus_x][focus_y].mean) / individuals[focus_x][focus_y].variance * (environment[x_coordinate][y_coordinate].first - individuals[focus_x][focus_y].mean) / individuals[focus_x][focus_y].variance)); // environmental influence !
//  
//                       seedSum -= weights[focus_x*m_Ydimensions + focus_y ];
//                       weights[focus_x*m_Ydimensions + focus_y ] = envFitnessParent  * m_Individuals[focus_x][focus_y].m_LocalDensity + (DBL_MIN*100.0);
//                       seedSum += weights[focus_x*m_Ydimensions + focus_y ];
//                    }
//                    else
//                    {
//                       seedSum -= weights[focus_x*m_Ydimensions + focus_y ];
//                       weights[focus_x*m_Ydimensions + focus_y ] = m_Individuals[focus_x][focus_y].m_LocalDensity + (DBL_MIN*100.0);
//                       seedSum +=  weights[focus_x*m_Ydimensions + focus_y ];
//                    }
//                    // end useless                  
//                 }
//              }
//  
//              unsigned int start =  ((x_coordinate -m_DensityCutoff + m_Xdimensions) % m_Xdimensions) * m_Ydimensions + ((y_coordinate - m_DensityCutoff + m_Ydimensions) % m_Ydimensions) ;
//              unsigned int end = m_LandscapeSize;
//  
//              if(start == 0) cumWeights[0] = weights[0] ;
//              else cumWeights[start] = cumWeights[start-1] + weights[start] ;
//  
//              for(unsigned int kk = start+1; kk < end; kk++)
//              {
//                 cumWeights[kk] = weights[kk] + cumWeights[kk-1];
//                 //						std::cout << seedSum << " : " << cumWeights[kk] << " : " << kk << " : " << array_length-1 <<  '\n';
//              }
//              //					std::cout << seedSum << " : " << cumWeights[array_length-1] << " : " << array_length-1  << " : " << start  << " : " << end <<  '\n';
//           }
//         
        }
        // end reproduction calculations

         //				std::cout << seedSum << " : " << cumWeights[array_length-1] << " : " << array_length-1 <<  '\n';

      }
   }
}


LocalEnvironment::LocalEnvironment(int xsize, int ysize, int type, bool neutral, bool dd, bool env,bool mort, bool repro, unsigned int runs, double specRate, int dispersalCutoff, int densityCutoff, unsigned int mortalityStrength,double envStrength, double compStrength, int fission, double redQueen, double redQueenStrength, int protracted, std::vector<double> airmat, std::vector<double> soilmat) :
	                              Landscape(xsize,  ysize,  type,  neutral,  dd,  env, mort, repro,  runs, specRate, dispersalCutoff, densityCutoff, mortalityStrength, envStrength, compStrength, fission, redQueen, redQueenStrength, protracted, airmat, soilmat)
{

}

LocalEnvironment::~LocalEnvironment() {
   // TODO: Implement destructor!
}

void LocalEnvironment::reproduce(unsigned int generation)
{
  
   #ifdef DEBUG
   std::cout<<"In local \n";
   #endif
  
   int x_coordinate;
   int y_coordinate;
   int x_parent = 0, kernel_x = 0;
   int y_parent = 0, kernel_y = 0;
   unsigned int array_length;
   double seedSum = 0.0;
   int new_parent;
   unsigned int numberOfRuns = 0;
   m_KernelSize = ((2*m_Cutoff+1)*(2*m_Cutoff+1)-1);

   //if(m_mortality) numberOfRuns =  m_LandscapeSize*2;
   //else numberOfRuns =  m_LandscapeSize;
   
   numberOfRuns =  m_LandscapeSize;
  
  // TODO 

   int event = 0;
   int numberDeath = 0;

   while(numberDeath < numberOfRuns)
   {
      event++;
      
      // Take a new coordinate
      
      x_coordinate = m_RandomGenerator.randomInt(0,m_Xdimensions-1);
      y_coordinate = m_RandomGenerator.randomInt(0,m_Ydimensions-1);

      // If fitness acts on mortality, break here with a chance ~ fitness
      
      if(m_mortality && event % m_mortalityStrength != 0){ // important!! the frequency in relation to the base mortality controls the intensity of the mechanisms
    	  double weight = m_Individuals[x_coordinate][y_coordinate].getFitness(m_Environment[x_coordinate * m_Ydimensions + y_coordinate].first, m_Env, m_DD, generation, m_redQueenStrength, m_redQueen);
        double chanceOfDeath = m_RandomGenerator.randomDouble(0.0,1.0);
		  if(weight > chanceOfDeath) continue;
	    }
      
      // If we continue, count up number of deaths and perform replacement
      
      numberDeath++;
      
      m_Individuals[x_coordinate][y_coordinate].reportDeath(generation);

      ////////////////////////////////////////////
      // DISPERSAL 
      
      std::vector<std::pair<int, int> > parents (m_KernelSize);

      double weights[m_KernelSize];
      array_length = 0;
      seedSum = 0.0;

      for(int relativeX = - m_Cutoff; relativeX <= m_Cutoff; relativeX++)
      {
         int yLims = floor(sqrt(m_Cutoff* m_Cutoff - relativeX*relativeX )); // avoid diagonal bias
         for(int relativeY = - yLims ; relativeY <=  yLims ; relativeY++)
         {
            kernel_x = (x_coordinate + relativeX + m_Xdimensions) % m_Xdimensions;
            kernel_y = (y_coordinate + relativeY + m_Ydimensions) % m_Ydimensions;
            if (!(kernel_x == x_coordinate && kernel_y == y_coordinate)) 
            {
               parents[array_length].first = kernel_x;
               parents[array_length].second =  kernel_y;
               if(m_reproduction){
                 weights[array_length] = m_Individuals[kernel_x][kernel_y].getSeedsTo(relativeX,relativeY, m_Dispersal_type, m_Environment[kernel_x * m_Ydimensions + kernel_y].first, m_Env, m_DD, generation, m_redQueenStrength, m_redQueen);
               }else if(m_mortality){
                 weights[array_length] = m_Individuals[kernel_x][kernel_y].dispersal(m_Dispersal_type, m_Individuals[kernel_x][kernel_y].euclidian_distance(relativeX, relativeY));
               }else throw 11;
               seedSum += weights[array_length];
               array_length +=1;
            }
         }
      }

      // calculate cummulative weights for the multinomial function
      // TODO could be moved in the upper loop
      // TODO check multinomial implementation correct, should cumweights start with 0?

      double cumWeights[array_length] ;
      cumWeights[0] = weights[0];
      for(unsigned int i =1; i< array_length; i++)
      {
         cumWeights[i] = weights[i] + cumWeights[i-1];
      }

      new_parent = m_RandomGenerator.multinomialDraw(cumWeights, array_length-1, seedSum);
      x_parent = parents[new_parent].first;
      y_parent = parents[new_parent].second;

      m_Individuals[x_coordinate][y_coordinate] = m_Individuals[x_parent][y_parent];  // overloaded operator, deep copy ! // report birth and evolve is automatic
      m_Individuals[x_coordinate][y_coordinate].m_X_coordinate = x_coordinate;
      m_Individuals[x_coordinate][y_coordinate].m_Y_coordinate = y_coordinate;
      
      m_Individuals[x_coordinate][y_coordinate].evolve();
 
      // UPDATE RELATEDNESS

      if(m_DD) densityUpdate(x_coordinate,y_coordinate);
   }
}


// update the density around an individual with 
void Landscape::densityUpdate(int x, int y){

   int focus_x, focus_y, neighborX, neighborY;
   double relatedness;

   for(int X1= - m_DensityCutoff; X1 <= m_DensityCutoff; X1++)
   {
      int yLims = floor(sqrt(m_DensityCutoff* m_DensityCutoff - X1*X1 )); 
      for(int Y1 = - yLims ; Y1 <=  yLims ; Y1++)
      {
         // get focus cell
         focus_x = ((x + X1 + m_Xdimensions) % m_Xdimensions);
         focus_y = ((y + Y1 + m_Ydimensions) % m_Ydimensions);

         relatedness = 0.0;

         for(int X2= - m_DensityCutoff; X2 <= m_DensityCutoff; X2++)
         {
            int yLims2 = floor(sqrt(m_DensityCutoff* m_DensityCutoff - X2*X2 )); // avoid diagonal bias
            for(int Y2 = - yLims2 ; Y2 <=  yLims2 ; Y2++)
            {
               neighborX = ((focus_x + X2 + m_Xdimensions) % m_Xdimensions);
               neighborY = ((focus_y + Y2 + m_Ydimensions) % m_Ydimensions);

               if (!(neighborX == focus_x && neighborY == focus_y))
               {
                 double a = m_Individuals[focus_x][focus_y].m_CompetitionMarker;
                 double b = m_Individuals[neighborX][neighborY].m_CompetitionMarker;
                 double diff = std::abs(a -b);
                 relatedness += diff;
                 //if (diff < 0.5) relatedness += diff / 0.5;
                 //else relatedness += 1;
               }
            }
         }
         m_Individuals[focus_x][focus_y].m_LocalDensity = relatedness / cellsWithinDensityCutoff;
         //std::cout<<m_Individuals[focus_x][focus_y].m_LocalDensity<<"\n"; //DEBUG
      }
   } 
}



void Landscape::speciation (unsigned int generation)
{

	// In this function is the main part of the implementation of the speciation mechanisms.

	  int specRate = m_RandomGenerator.randomPoisson(m_Speciation_Rate);


	  for (int i = 0; i < specRate; i++)
	  {

		// std::cout << specRate << std::endl;
	    int x = m_RandomGenerator.randomInt(0,m_Xdimensions-1); // rand() % xdimensions;
	    int y = m_RandomGenerator.randomInt(0,m_Ydimensions-1); // rand() % ydimensions;



  // In every case the counter for the incipient ages set to zero.
 // This value is updated every time step. Which individuals are chosen is dependent on the
 // other speciation mechanisms (fission or no fission).

	    if(m_fission ==0){
	       m_Individuals[x][y].m_incip_Age = 0;
	    }

	    if(m_fission != 0){
	    	std::vector<int> xvec;
	        std::vector<int> yvec;
	      for(int k = 0; k < m_Xdimensions; k++){
	        for(int j = 0; j < m_Ydimensions; j++){
	          if(m_Individuals[k][j].m_Species == m_Individuals[x][y].m_Species)
	        	xvec.push_back(k);
	            yvec.push_back(j);
	        }}


	        	if(m_fission == 1){
	        		for(unsigned int z =0; z< xvec.size(); z+=2){
	        			m_Individuals[xvec[z]][yvec[z]].m_incip_Age = 0;
	        		}
	        	}

	        	if(m_fission ==2){
	    	        int sum = std::accumulate(xvec.begin(), xvec.end(), 0);
	    	        int mean = sum / xvec.size();
	    	        if(sum == 0) mean=-1;

	    	        for(unsigned int z = 0; z < xvec.size(); z++){
	    	          if(xvec[z] < mean){
	    	           m_Individuals[xvec[z]][yvec[z]].m_incip_Age = 0;
	        	}
	          }
	          }

	}}



	    // In the same timestep new species arise.
	   // Which individuals become part of the new species is dependent on the incipient age of the
	  // individual and the other speciation mechanism (fission or no fission).

	if(m_fission == 0){
		for(int k = 0; k < m_Xdimensions; k++){
	           for(int j = 0; j < m_Ydimensions; j++){
	        	if(m_Individuals[k][j].m_incip_Age == m_protracted){
	        		m_Global_Species_Counter+=1;

	        		  m_Individuals[k][j].m_Species->m_Children.push_back(m_Global_Species_Counter);
	        		  m_Individuals[k][j].reportDeath(generation);
                      m_Individuals[k][j].m_Species = new Species(m_Global_Species_Counter , m_Individuals[k][j].m_Species->get_species_ID(),
                    		  generation, std::make_pair(k, j), m_SimulationEnd);
                      m_Individuals[k][j].evolveDuringSpeciation();
	        	      m_Phylogeny.updatePhylogeny(m_Individuals[k][j].m_Species);

	        	      if(m_DD) densityUpdate(k,j);
	        	}
	          }
	         }
	       }

	if(m_fission != 0){
		  int oldspec;
		for(int k = 0; k < m_Xdimensions; k++){
	       for(int j = 0; j < m_Ydimensions; j++){
			  if(m_Individuals[k][j].m_incip_Age == m_protracted){
				      m_Global_Species_Counter+=1;
					  m_Individuals[k][j].m_Species->m_Children.push_back(m_Global_Species_Counter);
					  m_Individuals[k][j].reportDeath(generation);

					  oldspec =  m_Individuals[k][j].m_Species->m_ID;
					  m_Individuals[k][j].m_Species = new Species(m_Global_Species_Counter , m_Individuals[k][j].m_Species->get_species_ID(),
					                    		  generation, std::make_pair(k, j), m_SimulationEnd);

					  m_Individuals[k][j].evolveDuringSpeciation();
				      m_Phylogeny.updatePhylogeny(m_Individuals[k][j].m_Species);
					  if(m_DD) densityUpdate(k,j);
					  m_Individuals[k][j].m_incip_Age = -99999999999999;

					  for(int z = 0; z< m_Xdimensions; z++){
						  for(int q=0; q< m_Ydimensions; q++){
							  if((m_Individuals[z][q].m_incip_Age == m_protracted) && (m_Individuals[z][q].m_Species->m_ID ==  oldspec)){
								  m_Individuals[z][q].reportDeath(generation);
								  m_Individuals[z][q].m_Species = m_Individuals[k][j].m_Species;
								  m_Individuals[z][q].evolveDuringSpeciation();
								  if(m_DD) densityUpdate(z,q);
								  m_Individuals[z][q].m_incip_Age = -99999999999999;
							  }
						  }
					  }
				  }
		 	}}}


}
        



