/*
 * species.cpp
 *
 *  Created on: 26.06.2014
 *      Author: Paul
 */

#include "Species.h"

Species::Species()
{
	this->m_ID = 0;
	this->m_Ancestor = 0;
	this->m_Count = 0;
	this->m_Date_of_Emergence = 0;
	this->m_Origin = std::make_pair(0,0);
	this->m_Variance = 0.04659906;
	this->m_Mean = 1.0;
	this->m_FirstMean = 0.5;
	this->m_FirstComp = 0.5;
	this->m_FirstNeutral = 0.5;
	this->m_MeanSum = 0.0;
	this->m_CompetitionSum = 0.0;
	this->m_NeutralSum = 0.0;
	this->m_Date_of_Extinction = 100;
	this -> m_CompetitionMean = 0.5;
	this -> m_NeutralMean = 0.5;
}

Species::Species(const Species & spec)
{
	this->m_ID = spec.m_ID;
	this->m_Ancestor = spec.m_Ancestor;
	this->m_Count = spec.m_Count;
	this->m_Date_of_Emergence = spec.m_Date_of_Emergence;
	this->m_Origin = spec.m_Origin;
	this->m_Variance = spec.m_Variance;
	this->m_Mean =spec.m_Mean;
	this->m_FirstMean = spec.m_FirstMean;
	this->m_FirstComp = spec.m_FirstComp;
	this->m_FirstNeutral =spec.m_FirstNeutral;
	this->m_MeanSum = spec.m_MeanSum;
	this->m_CompetitionSum = spec.m_CompetitionSum;
	this->m_NeutralSum = spec.m_NeutralSum;
	this->m_Date_of_Extinction = spec.m_Date_of_Extinction;
	this -> m_CompetitionMean = spec.m_CompetitionMean;
	this -> m_NeutralMean = spec.m_NeutralMean;
	this ->m_Children = spec.m_Children;
}

void Species::operator=(const Species & spec)
{
	this->m_ID = spec.m_ID;
	this->m_Ancestor = spec.m_Ancestor;
	this->m_Count = spec.m_Count;
	this->m_Date_of_Emergence = spec.m_Date_of_Emergence;
	this->m_Origin = spec.m_Origin;
	this->m_Variance = spec.m_Variance;
	this->m_Mean =spec.m_Mean;
	this->m_FirstMean = spec.m_FirstMean;
	this->m_FirstComp = spec.m_FirstComp;
	this->m_FirstNeutral =spec.m_FirstNeutral;
	this->m_MeanSum = spec.m_MeanSum;
	this->m_CompetitionSum = spec.m_CompetitionSum;
	this->m_NeutralSum = spec.m_NeutralSum;
	this->m_Date_of_Extinction = spec.m_Date_of_Extinction;
	this -> m_CompetitionMean = spec.m_CompetitionMean;
	this -> m_NeutralMean = spec.m_NeutralMean;
	this ->m_Children = spec.m_Children;
}

Species::Species(unsigned long long id, unsigned long long father, unsigned long long date, std::pair<int, int> birthplace, unsigned int simulationEnd)
{
	this->m_ID = id;
	this->m_Ancestor = father;
	this->m_Count = 0;
	this->m_Date_of_Emergence = date;
	this->m_Origin = birthplace;
	this->m_Variance = 0.04659906;
	this->m_Mean = 1.0;
	this->m_FirstMean = 1.0;
	this->m_FirstComp = 0.5;
	this->m_FirstNeutral = 0.5;
	this->m_MeanSum = 0.0;
	this->m_CompetitionSum = 0.0;
	this->m_NeutralSum = 0.0;
	this->m_Date_of_Extinction = simulationEnd;
	this -> m_CompetitionMean = 0.5;
	this -> m_NeutralMean = 0.5;
}
Species::~Species(){

}

//	void species::increase_count()
//	{
//		count++ ;
//	}
//	void  species::decrease_count()
//	{
//		count-- ;
//	}

	void Species::addIndividual(double env, double comp, double neutral){
    m_Count += 1;
		m_MeanSum += env;
		m_CompetitionSum += comp;
		m_NeutralSum += neutral;
    updateMean();
	}

	void Species::removeIndividual(double env, double comp, double neutral, int generation){ 
    
    m_Count -=1;
  
    if(m_Count < 0){
       m_Date_of_Extinction = generation;
    }
    else{
       m_MeanSum -= env;
       m_CompetitionSum -= comp;
    	 m_NeutralSum -= neutral;
       updateMean();
    }
	}

	void Species::updateMean(){
		m_Mean = m_MeanSum / (double(m_Count));
		m_CompetitionMean = m_CompetitionSum / (double(m_Count));
		m_NeutralMean = m_NeutralSum / (double(m_Count));
	}

	unsigned long long Species::get_species_count(){
		return m_Count;
	}

	unsigned long long Species::get_species_ID(){
		return m_ID;
	}



