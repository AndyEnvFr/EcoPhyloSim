/*
 * model.c++
 *
 *  Created on: 20.06.2014
 *      Author: Paul
 *              Betim Musa <musab@informatik.uni-freiburg.de>
 */

#include <iostream>
#include <fstream>
#include <string>
#include "./RandomGen.h"
//#include "Utils/RandomGen.h"
#include "Individual.h"
#include "Grid.h"
#include "PhylSimModel.h"
#include "Species.h"


PhylSimModel::PhylSimModel(int X, int Y, int dispersal, int simulationEnd, double specRate, bool dens,
                           bool env, bool neutral, bool mort, int mortalityStrength, bool repro, int dispersalCutoff,
                           int densityCutoff, std::string saveLocation, double envStrength, double compStrength,
                           int fission,
                           double redQueen, double redQueenStrength, int protracted, std::vector<double> airmat,
                           std::vector<double> soilmat, double nicheWidth) {

#ifdef DEBUG
    std::cout<<"Running simulation with \n";
    std::cout<<"Dispersal "<< dispersal << "; cutoff " << dispersalCutoff <<  "\n";
    std::cout<<"Competition "<< dens << "; strength " << compStrength <<  " cutoff " << densityCutoff << "\n";
    std::cout<<"Environment "<< env << "; strength " << envStrength<<  "\n";
    std::cout<<"\n---- debug message for development purposes, remove debug switch in debug.h for turning this off \n\n";
#endif

    if (dispersal == 1) {
        m_Global = new GlobalEnvironment(X, Y, dispersal, neutral, dens, env, mort, repro, simulationEnd, specRate,
                                         dispersalCutoff, densityCutoff, mortalityStrength, envStrength, compStrength,
                                         fission, redQueen, redQueenStrength, protracted, airmat, soilmat, nicheWidth);
        m_Local = NULL;
    } else if (dispersal == 2 || dispersal == 3) {
        m_Global = NULL;
        m_Local = new LocalEnvironment(X, Y, dispersal, neutral, dens, env, mort, repro, simulationEnd, specRate,
                                       dispersalCutoff, densityCutoff, mortalityStrength, envStrength, compStrength,
                                       fission, redQueen, redQueenStrength, protracted, airmat, soilmat, nicheWidth);
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

        // std::cout << "generation :" <<  generation << '/' << runs << '\n';

#ifdef DEBUG
        if(generation % 1000 == 0){
          std::cout << "generation :" <<  generation << '/' << runs << '\n';
        }
#endif

        if (m_Dispersal == 1) {
            m_Global->increaseAge(timeStep + generation);
            m_Global->reproduce(timeStep + generation);
            m_Global->speciation(timeStep + generation);
            //  std::cout << "Disp generation :" <<  generation << '/' << runs << '\n';
        } else if (m_Dispersal == 3) {
            m_Local->increaseAge(timeStep + generation);
            m_Local->reproduce(timeStep + generation);
            m_Local->speciation(timeStep + generation);
            //   std::cout << "Disp generation :" <<  generation << '/' << runs << '\n';
        }
    }
    timeStep += runs;
}

//void model::get_data(){
//
//		std::ofstream test_matrix("..\\test_out.txt");
//		for(int ba=0;ba<Landscape.get_dimensions().first;ba++){
//			test_matrix << '\n';
//
//		for(int bu=0;bu<Landscape.get_dimensions().second;bu++){
//			test_matrix << Landscape.individuals[ba][bu].Species->ID <<',';
//		}}
//
//}

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

//void model::gettraits(){
//	std::ofstream trait_matrix("..\\trait_out.txt");
//		for(int ba=0;ba<Landscape.get_dimensions().first;ba++){
//			trait_matrix << '\n';
//
//		for(int bu=0;bu<Landscape.get_dimensions().second;bu++){
//			trait_matrix << Landscape.individuals[ba][bu].mean <<',';
//		}}
//
//}
