/*
 * IonProbLi.cpp
 * Calculates Final Ion Populations for an Intensity and Wavelength
 * Makes a Table of Final Ionization Population for a given Wavelength and Intensity Range
 * Also files for Population over time
 * For Lithium
 * Based on carbon_pop.f and 1e_v1.f
 * Zach Germain, Department of Physics and Astronomy, University of Delaware
 * Last Updated: 8/19/19
 */

/* In order to change the intensity range:
 * ************Change IntensitySIMinMax; input minimum and magnitudes
 * ************Change nInt; the number of intensities you wish to compute
 * In order to change the wavelength:
 * ************Change wavelength in atomic units
 * In order to change Element:
 * ************You need to change all the atomic parameters
 * ************Change File Names with Symbol
 * ************Add/Eliminate some of the prints to files that are based on no. elements
 * ************Add/Eliminate some of the population rate equations
 */

/* How the Calculation is Done
 * For some user-specified range of intensities, ArraysCalc calculates the
 * *****intensities in SI and au as well as electric field in atomic units
 * Loops over all intensities
 * For each intensity, it computes populations by using RKSuite
 * The differential equation for each ion looks something like this:
 * *****dpop_i/dt=-rate_i*pop_i+rate_(i-1)*pop_(i-1)
 * The rate is the tunneling ionization rate and is calculated by using
 * *****the ADK model
 * Outputs populations for each intensity to separate files
 * Outputs final ion populations for all intensities to a main file
 */

/****************************************************************************/
/****************************************************************************/
#include<iostream>
#include<fstream>
#include<string>
#include<sstream>

#include<cmath>
#include<ctime>

#include "rksuite.h"

using namespace std;
/****************************************************************************/
/****************************************************************************/


/*********************************Constants**********************************/
/****************************************************************************/
//Speed of Light in au
 double c = 137.03545;
//Pi
 double pi = 4*atan(1);
//Atomic Unit of Intensity
 double intAu = 6.43640931e15;
 /****************************************************************************/
 /**********************************END***************************************/


 /****************************Atomic Parameters*******************************/
 /****************************************************************************/
//Number of Ions, used for indexing
const int nIon = 4; //Z+1: Last element is 1 as placeholder
//Ionization Potential
 double Ip[nIon] =
		{
		0.1988535,
		2.77972,
		4.5001122,
		1
		};
//Nstar: effective principle quantum number
 double nstar[nIon] =
		{
		1.588535,
		0.848232,
		0.999988,
		1
		};
//flm: function of l,m
 double flm[nIon] =
		{
		1,
		1,
		1,
		1
		};
//C2nl: coefficient for ADK
 double C2nl[nIon] =
		{
		2.40312,
		4.2081,
		4.0000,
		1
		};
//Ion Numbers
 double ionNum[nIon] =
		{
		0,
		1,
		2,
		3,
		};
 //Orbital quantum numbers
 int lNum[nIon] =
 		{
		0,
		0,
		0,
 		1
 		};
 /****************************************************************************/
 /**********************************END***************************************/


 /*******************************Laser Parameters*****************************/
 /****************************************************************************/
//Wavelength of Light: Changeable by user inputted as aus
//si~au
//200nm~3779, 400nm~7558, 800nm~15118,
//1.6um~30235, 3.2um~60471, 6.4um~120942
 double wavelength = 7558,
//Frequency of light: calculated from wavelength
			 freq = (2*pi*c)/wavelength,
//Period of the Pulse
			 period = wavelength/c,
//Intensity of Light: Minimum and Maximum Magnitudes
//Stored as W/cm^2
			intensitySIMinMax[2] = {log10(1e12),log10(1e18)},
//intensitySIMinMax converted to au
//Not in magnitudes
			intensityAuMinMax[2] = {pow(10.0,intensitySIMinMax[0])/intAu, pow(10.0,intensitySIMinMax[1])/intAu};
//Desired Number of Intensities
const int nInt = 1000;
//Array of Intensity in SI
double intSIArr[nInt],
//Array of Intensity in atomic units
	   intAuArr[nInt],
//Amplitude of Electric Field for an Intensity
       eAmp[nInt],
//Electric Field at some Point t
	   eField[nInt];
/****************************************************************************/
/*************************************END************************************/

/*******************Exterior Loop Variables for Problem**********************/
/****************************************************************************/
//These current values are for working inside the loop for some current
//iterative or time
double currIntSI,
	   currIntAu,
	   currEAmp,
	   eFieldt,
//Time Variables For Integration
	   numSteps = 1000000*(wavelength/7558), //Number of steps of integration, changeable
	   tStart = 0,
	   tFinal = period, //Last time its integrated too
	   tDelta = period/numSteps, //How much endpoint changes after integration
	   tIntEnd = tStart + tDelta, //Changing endpoint for integration
	   tCurr;
int iInt; //Iterative for intensities
/****************************************************************************/
/*************************************END************************************/


/**************************Function Declarations*****************************/
/****************************************************************************/
double GetEField(double tNow, int iEAmp);
void ArraysCalc();
double ADKRate(int iIon);
string getStringFromNumber(int n);
void popDiffEqs(double tHave, double popHave[], double dpopHave[]);
/****************************************************************************/
/************************************END*************************************/

//Begins Main Program
int main()
{
	//Begins time recorded
	float elapsed1 = (float)clock();

	//Initialization of Main Arrays for Problem
	ArraysCalc();

	/****************************Variable Declaration****************************/
	/****************************************************************************/
	//RKSuite Variable Intialization
	int meth = 2,
	//CFlag initialization. 1 works, >1 Error
		cFlag;
	//Tolerance
	double tol = pow(10,-4),
	//Threshold array
		thres[nIon],
	//This value gives starting point of integration to program
		hStart = 0;
	//Does program give messages?
	bool mesage = false,
	//Does assess error?
		errAss = false;
	for (int i = 0; i < nIon; i++)
		{
			thres[i]=10e-8;
		}

	//Population Declaration
	double rateStart[nIon];
	/****************************************************************************/
	/*************************************END************************************/


	/***************************File Creation and Headers************************/
	/****************************************************************************/
	//Creates And Outputs to Log Files
	ofstream outLog;
	outLog.open("IonFinalPopLogLi.txt");
	outLog << "Wavelength(au): " << wavelength << endl
		   << "Frequency(au): " << freq << endl
		   << "Period(au): " << period << endl
		   << "Minimum Intensity(si): " << pow(10.0, intensitySIMinMax[0]) << endl
		   << "Maximum Intensity(si): " << pow(10.0, intensitySIMinMax[1]) << endl
		   << "Minimum Intensity(au): " << intensityAuMinMax[0] << endl
		   << "Maximum Intensity(au): " << intensityAuMinMax[1] << endl;

	//Creates and Writes Headers for Main Data File
	ofstream intProb;
	intProb.open("IntensityvsFinalPopLi.dat");
	intProb << "Intensity(SI)" << " "
			<< "Intensity(AU)" << " "
			<< "Final_Ion_Population_for_Li" << " "
			<< "Final_Ion_Population_for_Li1+" << " "
			<< "Final_Ion_Population_for_Li2+" << " "
			<< "Final_Ion_Population_for_Li3+"
			<< endl;
	//Flag file
	ofstream flagFile("Flag_File_Li.txt");
	/****************************************************************************/
	/*************************************END************************************/


	//This for loop iterates over all intensities
	for (iInt = 0; iInt < nInt; iInt++)
	{
		/*****************************Variable Assignation****************************/
		/****************************************************************************/
		//Assignation for Current Values for Intensity
		currIntSI = intSIArr[iInt];
		currIntAu = intAuArr[iInt];
		currEAmp = eAmp[iInt];
		eFieldt = 0;

		//Reinitialization of Variables
		tStart = 0,
	    tDelta = period/numSteps,
	    tIntEnd = tStart + tDelta;
		rateStart[0] = 1;
		for (int i = 1; i < nIon; i++)
		{
			rateStart[i] = 0;
		}
		double popNow[nIon],
			   dpopNow[nIon];
		for (int i = 1; i < nIon; i++)
		{
			dpopNow[i] = 0;
		}
		popNow[0] = 1;
		for (int i = 1; i < nIon; i++)
		{
			popNow[i] = 0;
		}
		double popAnte[nIon];
		popAnte[0] = 1;
		for (int i = 1; i < nIon; i++)
			{
				popAnte[i] = 0;
			}
		/****************************************************************************/
		/************************************END*************************************/


		/***************************File for This Intensity**************************/
		/****************************************************************************/
		/*//Intensity Output File Creation and Headers
		string outIntName = "int_" + getStringFromNumber(iInt) + "Li.dat";
		const char* fileNameInt = outIntName.c_str();
		ofstream outInt(fileNameInt);
		outInt << "Intensity(au)" << " "
			   << "Time(au)" << " "
			   << "Electric_Field(au)" << " "
			   << "Population_Li" << " "
			   << "Population_Li+" << " "
			   << "Population_Li2+" << " "
			   << "Population_Li3+"
			   << endl;*/
		/****************************************************************************/
		/*************************************END************************************/


		//RKSuite Setup
		RKSUITE rkSuite;
		rkSuite.setup(nIon,tStart,rateStart,tIntEnd,tol,thres, meth,
				"C",errAss,hStart,mesage);

		//This loop iterates over all time
		while (tIntEnd <= tFinal)
		{
			eFieldt = GetEField(tIntEnd-tDelta, iInt);
			rkSuite.ct(popDiffEqs, tCurr, popNow, dpopNow, cFlag);
			//Loop sets ion pops to zero if they are negative
			for (int i = 0; i < nIon; i++)
			{
				if (popNow[i] < 0)
				{
					popNow[i] = 0;
				}
			}
			//Outputs to intensity file
/*			outInt << currIntAu << " "
				   << tIntEnd << " "
				   << eFieldt << " "
				   << popNow[0] << " "
				   << popNow[1] << " "
				   << popNow[2] << " "
				   << popNow[3]
				   << endl;*/
			//Flag File Tests
			for (int f = 0; f < nIon; f++)
			{
				if (popAnte[f] != 0 && popAnte[f]>1e-30 && popAnte[0]!=1)
				{
					double fracDiff = abs(popNow[f]/popAnte[f]);
					if (fracDiff >= 2)
					{
						flagFile << "For intensity " << intSIArr[iInt]
								 << " there has been a flag at time "
								 << tIntEnd << "." << endl
								 << "The difference between "
								 << "the current and previous population of ion "
								 << f << " is greater than 2." << endl
								 << "It is " << fracDiff << endl;
						flagFile << endl << endl << endl;
					}
					else if (fracDiff >= 1.5)
					{
						flagFile << "For intensity " << intSIArr[iInt]
								 << " there has been a flag at time "
								 << tIntEnd << "." << endl
								 << "The difference between "
								 << "the current and previous population of ion "
								 << f << " is greater than 1.5." << endl
								 << "It is " << fracDiff << endl;
						flagFile << endl << endl << endl;
					}
				}
				popAnte[f] = popNow[f];
			}
			//Adds to and resets the endpoint for the integration
			tIntEnd+=tDelta;
			rkSuite.reset(tIntEnd);


		} //END time iteration

		//Outputs to whole file
		intProb << currIntSI << " "
				<< currIntAu << " "
				<< popNow[0] << " "
				<< popNow[1] << " "
				<< popNow[2] << " "
				<< popNow[3]
				<< endl;

		/*outInt.close();*/
		tCurr = 0;
	} //END iteration over Intensities

	//Writes to log the elapsed time
	float elapsed2 = (double)clock();
	float tElapse = (elapsed2-elapsed1)/CLOCKS_PER_SEC;
	outLog << "Time Taken: " << tElapse << endl;

	//Closes file
	outLog.close();
	intProb.close();

}//END main

/*********************************Functions**********************************/
/****************************************************************************/

//Calculates the Electric Field for a previously described eAmp
//Takes in the current time and index for eAmp for it
//Returns the eFIeld at the time as for that eAmp
double GetEField(double tNow, int iEAmp)
{
	//Sine model for theory
	return eAmp[iEAmp]*sin(freq*tNow);
}
//End GetEField


//Calculates the Intensity and Amplitude of Electric Field Arrays in SI and Au units
void ArraysCalc()
{
	//Performs these Intensity Calculations using Magnitudes
	double siDiff = (intensitySIMinMax[1]-intensitySIMinMax[0])/(nInt-1);
	intSIArr[0] = pow(10, intensitySIMinMax[0]);
	intAuArr[0] = intSIArr[0]/intAu;
	eAmp[0] = sqrt((8.0*pi*intAuArr[0])/c);

	for (int i = 1; i < nInt; i++)
	{
		intSIArr[i] = pow(10,intensitySIMinMax[0]+i*siDiff);
		intAuArr[i] = intSIArr[i]/intAu;
		eAmp[i] = sqrt((8.0*pi*intAuArr[i])/c);
	}
}
//End ArraysCalc


//Calculates ADK for indexed Ion and Intensity
//iIon is the index of the Ion number desired
//iintAuArr is the current intensity desired in au
double ADKRate(int iIon)
{
	if (eFieldt == 0.0)
	{
		return 0.0;
	}
	else
	{
		double epsilon = pow(2*Ip[iIon], 1.5);
		double factor = epsilon/fabs(eFieldt);
		double nmpower1 = 2*nstar[iIon] - 1; //l=0
		double nmpower2 = 2*nstar[iIon] - 2; //l=1
		double adk1 = C2nl[iIon]*Ip[iIon]*flm[iIon]*sqrt(3.0/(pi*factor))*
				pow(2*factor,nmpower1)*exp(-(2.0*factor)/3.0);
		double adk2 = C2nl[iIon]*Ip[iIon]*flm[iIon]*sqrt(3.0/(pi*factor))*
						pow(2*factor,nmpower2)*exp(-(2.0*factor)/3.0);
		if (lNum[nIon] == 0)
		{
			return adk1;
		}
		else if (lNum[nIon] == 1)
		{
			return (adk1 + 2*adk2)/3;
		}
	}
}
//End ADKRate


//This function calculates the differential rate equations
//The equations follow this general formula:
//popi'=-ADKi*popi+ADK(i-1)*pop(i-1)
//The first ion can't be added to
//The last ion can't be subtracted from
void popDiffEqs(double tHave, double popHave[], double dpopHave[])
{
	dpopHave[0] = -ADKRate(0)*popHave[0];
	dpopHave[1] = -ADKRate(1)*popHave[1] + ADKRate(0)*popHave[0];
	dpopHave[2] = -ADKRate(2)*popHave[2] + ADKRate(1)*popHave[1];
	dpopHave[3] = 						   ADKRate(2)*popHave[2];
}
//Ends popDiffEqs


/*this function return string of input number*/
//Taken From cutoff_v3_Ne.cpp by Sui Luo
string getStringFromNumber(int n){

  stringstream ss;
  ss << n;
  return ss.str();

  /*end of function*/
}
/****************************************************************************/
/****************************************************************************/
