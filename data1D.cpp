/*
 * This is the implementation file of class "Data1D".
 * Used for:bin energy spectrum, bin angular distribution,
 *          store intensity dependence.
 * Author:Sui Luo. 
 * Version:First written in Feb 5, 2013. 
 *         Last modified in Oct 1, 2013.
 *         Version 0.3
*/

#include "data1D.h"

//constructor
Data1D::Data1D(double d1, double d2) {
  key = d1;
  weight = d2;
}

//modify the data
void Data1D::modify(double d1, double d2) {
  key = d1;
  weight = d2;
}

//return key
double Data1D::getKey() const {
  return key;
}

//return weight
double Data1D::getWeight() const {
  return weight;
}

//print data1D
string Data1D::toString() const {
  stringstream s;
  s << key << " " << weight << endl;
  return s.str();
}

