/*
 * This is the interface file of class "Logger".
 * Used for:pass ofstream file from main program.
 * Author:Sui Luo. 
 * Version:First written in Feb 23, 2013. 
 *         Last modified in Oct 01, 2013.
 *         Version 0.3
*/

#ifndef LOGGER_H
#define LOGGER_H

#include <iostream>
#include <fstream>
#include <string>
using namespace std;

class Logger {

 public:
  Logger(ofstream&);
  Logger();
  void operator= (Logger);//operator overloading
  void writeLog(string);
  void closeLog();
  friend class collection1D;
  friend class collection2D;
 private:
  ofstream* filePointer;

};

#endif /*LOGGER_H*/
