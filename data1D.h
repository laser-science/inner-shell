/*
 * This is the interface file of class "Data1D".
 * Used for:bin energy spectrum, bin angular distribution,
 *          store intensity dependence.
 * Author:Sui Luo. 
 * Version:First written in Feb 5, 2013. 
 *         Last modified in Oct 1, 2013.
 *         Version 0.3
*/

#ifndef DATA1D_H
#define DATA1D_H

#include <iostream>
#include <sstream>
#include <string>
using namespace std;

class Data1D {

 public:
  Data1D(double, double);
  void modify(double, double);
  double getKey() const;
  double getWeight() const;
  string toString() const;
  friend class Collection1D;//not necessary
  friend struct by_key;

 protected:
  double key, weight;
};

//define the operator to be used in sort
struct by_key{
  bool operator() (Data1D const &D1, Data1D const &D2){
    return D1.key < D2.key;
  }
};

#endif /*DATA1D_H*/
