/*
 * This is the implementation file of class "Logger".
 * Used for:pass ofstream file from main program.
 * Author:Sui Luo. 
 * Version:First written in Feb 23, 2013. 
 *         Last modified in Oct 01, 2013.
 *         Version 0.3
*/

#include "logger.h"

//constructor
Logger::Logger(ofstream& ofs1)
  :filePointer(&ofs1){}

//constructor
Logger::Logger() {
  filePointer = NULL;
}

//operator overloading
void Logger::operator= (Logger log2) {
  filePointer = log2.filePointer;
}

//write to log file
void Logger::writeLog(string s1) {
  *filePointer << s1;
}

//close log file
void Logger::closeLog() {
  filePointer->close();
}
