/*
 * paramters.cpp
 *
 *  Created on: 05.03.2015
 *      Author: Paul
 *              Betim Musa <musab@informatik.uni-freiburg.de>
 */
#include "Parameters.h"
#include "StringConversion.h"
#include "pugixml.hpp"
#include <ctime>
#include <iostream>
#include <map>
#include <string>
#include <vector>
//
using namespace std;
using namespace pugi;
//
// Define static variables
string Parameters::s_XMLFileName;
string Parameters::s_PathToXmlFile;
//
bool Parameters::s_WriteCurrentParameterToXMlFile;
bool Parameters::s_CreateParameters;
bool Parameters::s_ReadParametersFromXMLFile;
bool Parameters::s_LogInfos;
//
//
map<string, int> Parameters::m_IntParameters;
map<string, double> Parameters::m_DoubleParameters;
map<string, bool> Parameters::m_BoolParameters;

//
Parameters::Parameters() {
  // Init static member
  Parameters::s_XMLFileName = string("PhylogenyParameters.xml");
  Parameters::s_PathToXmlFile = string("XML/Data/");
  //  Parameters::m_WriteCurrentParameterToXMlFile = false;
  //  Parameters::m_ReadParametersFromXMLFile = !Parameters::m_WriteCurrentParameterToXMlFile;
  Parameters::s_WriteCurrentParameterToXMlFile = true;
  Parameters::s_ReadParametersFromXMLFile = true;
  Parameters::s_CreateParameters = true;
  Parameters::s_LogInfos = true;
  //
  if (s_CreateParameters) {
    createParameters();
  }
  //  if (Parameters::m_WriteCurrentParameterToXMlFile && !Parameters::m_ReadParametersFromXMLFile) {
  //     writeParametersToXMLFile();
  //  } else if (!Parameters::m_WriteCurrentParameterToXMlFile && Parameters::m_ReadParametersFromXMLFile) {
  //     readParameterFromXMLFile();
  //  }
  if (Parameters::s_ReadParametersFromXMLFile) {
    readParameterFromXMLFile();
  }
  if (Parameters::s_WriteCurrentParameterToXMlFile) {
    writeParametersToXMLFile();
  }
}

//
Parameters::~Parameters() {
  m_IntParameters.clear();
  m_BoolParameters.clear();
  m_DoubleParameters.clear();
}

//
void Parameters::readParameterFromXMLFile(string xmlFileName) {
  xml_document parameterDocument;
  if (s_LogInfos) {
    cout << "Filename of the XML document to load parameters from: " << xmlFileName.c_str() << endl;
  }
  // Load xml document
  xml_parse_result pars_res = parameterDocument.load_file(xmlFileName.c_str(), format_default, encoding_utf8);
  if (s_LogInfos) {
    // Print loading result
    cout << "Loading result: " << (pars_res == 1 ? "SUCCESS" : "FAIL") << endl;
  }
  xml_node params = parameterDocument.child("Phylogeny_Parameter_Document");
  // Extract parameter from the XML document.
  for (auto itr = params.children().begin(); itr != params.children().end(); ++itr) {
    string paramName = string((*itr).attribute("Name").value());
    string paramType = string((*itr).attribute("Type").value());
    string paramValue = string((*itr).attribute("Value").value());
    //
    if (s_LogInfos) {
      cout << "ChildName: " << (*itr).name() << endl;
      cout << "Parametername: " << paramName << endl;
      cout << "Parametertype: " << paramType << endl;
      cout << "Parameter value: " << paramValue << endl;
    }
    // Put parameter in according map
    // TODO(Betim): Catch conversion exceptions!
    if (paramType == string("int")) {
      int value = str_to_int(paramValue);
      m_IntParameters.insert(pair<string, int>(paramName, value));
    } else if (paramType == string("double")) {
      double value = str_to_double(paramValue);
      m_DoubleParameters.insert(pair<string, double>(paramName, value));
    } else if (paramType == string("bool")) {
      bool value = false;
      // Bool has only 2 possibilities, so we can do conversion by hand.
      if (paramValue == string("true"))
        value = true;
      m_BoolParameters.insert(pair<string, bool>(paramName, value));
    }
  }
}

//
string Parameters::getCurrentTimeStamp() {
  // Get current time and date
  time_t theTime = time(NULL);
  struct tm *localTime = localtime(&theTime);
  int min = localTime->tm_min;
  int hour = localTime->tm_hour;
  int day = localTime->tm_mday;
  int month = localTime->tm_mon + 1;
  int year = localTime->tm_year + 1900;
  string creationTime =
      to_string(hour) + ":" + to_string(min) + "__" + to_string(day) + "." + to_string(month) + "." + to_string(year);
  return creationTime;
}

//
void Parameters::writeParametersToXMLFile() {
  xml_document parameterDocument;
  xml_node creationTimeNode = parameterDocument.append_child("time_stamp");
  // Get current time stamp
  string creationTime = Parameters::getCurrentTimeStamp();
  creationTimeNode.append_attribute("created_at") = creationTime.c_str();
  //
  xml_node rootNode = parameterDocument.append_child("Phylogeny_Parameter_Document");
  //
  // Traverse backwards the parameter maps
  for (auto itr = m_IntParameters.rbegin(); itr != m_IntParameters.rend(); ++itr) {
    xml_node c = rootNode.append_child("Parameter");
    c.append_attribute("Name") = (*itr).first.c_str();
    c.append_attribute("Type") = "int";
    c.append_attribute("Value") = (*itr).second;
  }
  //
  for (auto itr = m_DoubleParameters.rbegin(); itr != m_DoubleParameters.rend(); ++itr) {
    xml_node c = rootNode.append_child("Parameter");
    c.append_attribute("Name") = (*itr).first.c_str();
    c.append_attribute("Type") = "double";
    c.append_attribute("Value") = (*itr).second;
  }
  //
  for (auto itr = m_BoolParameters.rbegin(); itr != m_BoolParameters.rend(); ++itr) {
    xml_node c = rootNode.append_child("Parameter");
    c.append_attribute("Name") = (*itr).first.c_str();
    c.append_attribute("Type") = "bool";
    c.append_attribute("Value") = (*itr).second;
  }
  //
  if (s_LogInfos) {
    cout << "Printing XML parameter document to stdout" << endl;
    parameterDocument.save(cout);
  }
  //
  string fileName = Parameters::s_PathToXmlFile + Parameters::s_XMLFileName;
  const char *file = fileName.c_str();
  //  const char* file =Parameters::m_XMLFileName.c_str();
  if (s_LogInfos) {
    cout << "Creating XML document!" << endl;
    cout << "Filename of the XML document: " << file << endl;
  }
  bool saveSuccessful = parameterDocument.save_file(file, "\t", format_default, encoding_utf8);
  if (s_LogInfos) {
    cout << "Saving XML file..." << (saveSuccessful ? "SUCCESS" : "FAIL") << endl;
  }
}

//
void Parameters::readParameterFromXMLFile() {
  string xmlFileName = Parameters::s_PathToXmlFile + Parameters::s_XMLFileName;
  Parameters::readParameterFromXMLFile(xmlFileName);
}

//
//
//
void Parameters::createParameters() {
  vector<pair<string, int>> intParameters;
  vector<pair<string, double>> doubleParameters;
  vector<pair<string, bool>> boolParameters;
  //
  // Parameters for Landscape
  intParameters.push_back(pair<string, int>(string("LandscapeSizeX"), 100));
  intParameters.push_back(pair<string, int>(string("LandscapeSizeY"), 200));
  intParameters.push_back(pair<string, int>(string("numberOfRuns"), 1000));
  //
  //
  // Parameters for Dispersion
  intParameters.push_back(pair<string, int>(string("DispersalMode"), 2));
  intParameters.push_back(pair<string, int>(string("DispersalCutOff"), 2));
  intParameters.push_back(pair<string, int>(string("nDensCutoff"), 2));
  intParameters.push_back(pair<string, int>(string("pDensCutoff"), 2));
  //
  doubleParameters.push_back(pair<string, double>(string("SpeciationRate"), 0.7));
  doubleParameters.push_back(pair<string, double>(string("DensityStrength"), 0.9));
  //
  boolParameters.push_back(pair<string, bool>(string("NeutralEnvironment"), false));
  boolParameters.push_back(pair<string, bool>(string("DensityDependence"), false));
  boolParameters.push_back(pair<string, bool>(string("EnvironmentalDependence"), true));
  //
  //
  // Parameters for Evolution
  doubleParameters.push_back(pair<string, double>(string("SpeciesInfluence"), 0.7));
  doubleParameters.push_back(pair<string, double>(string("RandomInfluence"), 0.7));
  //
  //
  // Put parameters into container
  m_IntParameters.insert(intParameters.begin(), intParameters.end());
  m_DoubleParameters.insert(doubleParameters.begin(), doubleParameters.end());
  m_BoolParameters.insert(boolParameters.begin(), boolParameters.end());
  //
  //
}
//
