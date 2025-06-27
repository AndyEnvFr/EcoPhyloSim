/*
 * Parameters.h
 *
 *  Created on: 29.05.2015
 *      Author: Betim Musa <musab@informatik.uni-freiburg.de>
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <iostream>
#include <map>
#include <string>

using namespace std;

class Parameters {

public:
  // Landscape
  int X;
  int Y;
  int runs;

  // Dispersion
  int dispersalMode;
  int dispersalCutoff;
  int nDensCutoff;
  int pDensCutoff;
  double speciationRate;
  double densityStrength;
  bool neutral;
  bool densityDependence;
  bool environmentalDependence;

  // Evolution
  double speciesInfluence;
  double randomInfluence;

  // Constructor / Destructor
  Parameters();

  virtual ~Parameters();

  /**
   * Writes the parameters from maps m_IntParameters, m_DoubleParameters, m_BoolParameters into an XML document.
   * The document has the following structure.
   * <?xml version="1.0"?>
   *  <time_stamp created_at="16:46__5.6.2015" />
   *  <Phylogeny_Parameter_Document>
   *     <Parameter Name="nameOfParameter" Type="typeOfParameter" Value="valueOfParameter" />
   *     ....
   *  </Phylogeny_Parameter_Document>
   */
  static void writeParametersToXMLFile();
  //
  /**
   * Reads/Loads parameters from the XMl document specified by m_XMLFileName. The parameters are then put into the
   * according parameter map (e.g. m_IntParameters, or m_DoubleParameters or m_BoolParameters).
   */
  static void readParameterFromXMLFile();

  //
  static void readParameterFromXMLFile(string xmlFileName);
  //
  /**
   * Template method to get a parameter specified by "parameterName".
   * See the specialized methods below for more information.
   */
  template <typename T> static T getParameterValue(string parameterName);
  //
private:
  /**
   * Specifies the filename of the XML document.
   */
  static string s_XMLFileName;
  /**
   * Specifies the path to the XML document file.
   */
  static string s_PathToXmlFile;
  //
  static bool s_WriteCurrentParameterToXMlFile;
  static bool s_CreateParameters;
  static bool s_ReadParametersFromXMLFile;
  static bool s_LogInfos;
  //
  static map<string, int> m_IntParameters;
  static map<string, double> m_DoubleParameters;
  static map<string, bool> m_BoolParameters;

  //
  void createParameters();

  static string getCurrentTimeStamp();
};
//
/**
 * Returns the integer parameter from the map m_IntParameters, if it contains the parameter specified by
 * "parameterName". Otherwise it returns the default value 0.
 */
template <> inline int Parameters::getParameterValue<int>(string parameterName) {
  auto itr = m_IntParameters.find(parameterName);
  if (itr != m_IntParameters.end()) {
    return itr->second;
  } else {
    cout << "Key is not in integer map. " << endl;
    cout << "Default value 0 is returned. " << endl;
  }
  return 0;
}

/**
 * Returns the double parameter from the map m_DoubleParameters, if it contains the parameter specified by
 * "parameterName". Otherwise it returns the default value 0.0.
 */
template <> inline double Parameters::getParameterValue<double>(string parameterName) {
  auto itr = m_DoubleParameters.find(parameterName);
  if (itr != m_DoubleParameters.end()) {
    return itr->second;
  } else {
    cout << "Key is not in double map." << endl;
    cout << "Default value 0.0 is returned. " << endl;
  }
  return 0.0;
}

/**
 * Returns the boolean parameter from the map m_BoolParameters, if it contains the parameter specified by
 * "parameterName". Otherwise it returns the default value 'false'.
 */
template <> inline bool Parameters::getParameterValue<bool>(string parameterName) {
  auto itr = m_BoolParameters.find(parameterName);
  if (itr != m_BoolParameters.end()) {
    return itr->second;
  } else {
    cout << "Key is not in boolean map." << endl;
    cout << "Default value 'false' is returned. " << endl;
  }
  return false;
}

#endif /* PARAMETERS_H_ */
