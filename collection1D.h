/*
 * This is the interface file of class "Collection1D".
 * Used for:bin energy spectrum, bin angular distribution,
 *          store intensity dependence.
 * Author:Sui Luo.
 * Version:First written in Feb 05, 2013.
 *         Last modified in May 04, 2015.
 *         Version 0.5
*/

/**********************************************************************/
/*header*/
/**********************************************************************/
#ifndef COLLECTION1D_H
#define COLLECTION1D_H

#include "data1D.h"
#include "logger.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace std;

/**********************************************************************/
/*class definition*/
/**********************************************************************/
class Collection1D {

  //public functions
 public:
  //constructor
  Collection1D();

  //desctructor
  ~Collection1D();

  //get min, max, dimension, limit and # of bin
  void getParameter(double,
		    double,
		    int,
		    double,
		    double);

  //get min, max, dimension, limit and # of bin - overloading
  void getParameter(int,
	       double,
	       double);

  //get logger input
  void getLog(Logger);

  //get output file name
  void getOutputName(string);

  //add new Data1D
  void addData(Data1D);

  //output unsorted data container
  void output();

  //sort data, bin data, convolute data and output - customerize
  void run();

  bool double_equals(double a, double b, double epsilon);

  //private functions and memebers
 private:
  //declare data container for input
  vector<Data1D>* collection;

  //declare data container for output
  vector<Data1D>* collection_out;

  //member of min/max/convolution
  double min, max, cov1, cov2;

  //member of dimension
  int dim;

  //member of output file name
  string fileInfo;

  //member of if external min/max flag
  bool externalMinMax;

  //member of logger
  Logger log;

  //sorc data container
  void sortC();

  //bin sorted data container
  void binC();

  //output binned+sorted data container
  void outputbinC();

  //convolute data using +- fixed percent method
  void convolution1C();

  //convolute data using +- fixed value method
  void convolution2C();

  //customize
  void customize();

};

#endif /*COLLECTION1D_H*/
