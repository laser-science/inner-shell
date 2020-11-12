/****************************************************************************/
/*
  Aim to see adk-weighted return flux vs. return energy spectrum.
  On 10% ionization species-intensity combination.
  Use 3D-trajectory-ensemble for rescatter calculation.
  Sui Luo. 
  Department of Physics and Astronomy, University of Delaware
  Commented and Edited Last by Zach Germain and David Milliken
  Last Updated: 10:00 08/16/2019
*/
/****************************************************************************/

/****************************************************************************/
/*Headers*/
/****************************************************************************/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <algorithm>
#include <vector>

/*rksuite*/
#include "rksuite.h"

/*distribution generation*/
#include "nr3.h"
#include "ran.h"
#include "deviates.h"

using namespace std;

/****************************************************************************/
/*constants*/
/****************************************************************************/
/*speed of light*/
const double c = 137.03545;

/*PI*/
const double pi = 4.0*atan(1.0);

/*atomic unit of intensity*/
const double int_au = FILL_IN;

/*full width half maximum of pulse*/
const double fwhm = FILL_IN;

/*wavelength of laser, written to be hard coded in nm and then converted to au */
const double wavelength_nm = FILL_IN;
const double wavelength = wavelength_nm*18.897;

/*frequency of pulse*/
const double freq = (2.0*pi*c)/wavelength;

/*period of pulse*/
const double period = wavelength/c;

/*sigma*/
const double sigma = fwhm/(2.0*sqrt(2.0*fabs(log(2.0))));

/*number of equations of motion*/
const int neq = 6;

/*number of time steps*/
//NEVER USED?
const int n_time = 30000;

/*number of rescatter angle*/
const int n_rescatter_angle = 180;

/*declare atomic parameters*/
/*Every vector is a list of parameters
 *The list elements correspond to the ions of the same index as in the ionNum list
 *Change n_charge to the number of elements in the lists, generally the number of ions you're running
 */
/*declare atomic parameters*/
const int n_charge = 1;
const double ip[n_charge] = {
		FILL_IN
};

const double nstar[n_charge] = {
		FILL_IN
};

const double c2nl[n_charge] = {
		FILL_IN
};

const double flm[n_charge] = {
		FILL_IN
};

const double intensityVect[n_charge] = {
		FILL_IN
};

const int ionNum[n_charge] = {
  FILL_IN
};

/*These are random seeds
 * When adding more ions, add more values to seedVect
 * The value doesn't matter, as it is a random seed and any input should give a random output when called
 */
const int seedVect[n_charge] =
{
	0
};

/*ADK threshold*/
//NEVER USED?
const double adk_thresh = 1.0e-300;
const double prob_e_thresh = 1.0e-6;

/****************************************************************************/
/*global*/
/****************************************************************************/
/*input intensity*/
double int_si;
/*intensity in atomic unit*/
double intensity;
/*EM peak amplitude from laser in atomic units*/
double eAmpMax;
/*initial phase of pulse*/
double iniphase;
/*job number*/
int n_job;

/****************************************************************************/
/*subroutine and function declare*/
/****************************************************************************/
double GetADK(const double,
              const double,
              const double,
              const double,
              const double[],
              const double);

void Derivs(double,
            double[],
            double[]);

void GetEmField(const double,
                const double,
                double*,
                double&);

double GetExitPoint(const double, 
		    const double, 
		    const double, 
		    const double);

void GetInitial(const double,
		const double,
		const double,
		const int,
		const int,
		vector<double>*,
		vector<double>*,
		vector<double>*,
		vector<double>*);

int start(int flag);

string getStringFromNumber(const int);

/****************************************************************************/
/*main program*/
/****************************************************************************/
/*This program runs through all ions and runs start for each of them
 * It only runs at one wavelength at a time - changing the wavelength can be done at line 53
 * Inputing an intensity as 0 will tell the program to skip that index
 */
int main(){
	for(int i = 0; i < n_charge; i++){
		if(intensityVect[i]!=0){
			start(i);
			cout << "Ion number " << ionNum[i] << " finished" << endl;
		}
		else
		{
			cout << "Ion number " << ionNum[i] << " had intensity 0" << endl;
		}
	}
	return 0;
}
/****************************************************************************/
/*start program*/
/****************************************************************************/
/*This program produces data only if rescattering happens
 * It runs all the trajectories for one ion, determined by flag */
int start(int flag){

	/*record running time*/
	clock_t elapsed1,elapsed2;
	elapsed1 = clock();

	//150 is an upper limit
	//int stepTime = 100; //v3 uses this
	int stepTime = 100;

	/*declare population evolution variables*/
	double
	rate_adk;
	/*double These variables are used when the GetEmField is used
	el_cpn[neq], around line 360
	elfield;*/
	int target_ion; //The ion number

	/*declare trajectory integration variables*/
	double
	y_start[neq], //Holds initial momentums and positions for a trajectory
	x_pre,
	x_now,
	traj_exit, //x-position of trajectory exit
	rescatter_time,
	rescatter_kin;
	int nSample = 1000;

	/*initialize rksuite parameters*/
	RKSUITE rksuite;
	double y[neq],yp[neq],ymax[neq];
	double twant, tnow;
	double tol = 1e-6; // =1e-6 for method 2
	double thres[neq] = {1e-10,1e-10,1e-10,1e-10,1e-10,1e-10}; //1e-10 for method 2
	double hstart = 0.0; //Gives control of choosing initial integration point
	int method = 2; //Method 2 is most efficient, medium-error
	int uflag; //For error flagging
	bool mesage = false;
	bool errass = false;

	/*job number*/
	/*The job number is to choose an ion based on the flag currently being run
	*/
	n_job = flag;

	/*target ion*/
	target_ion = ionNum[n_job];

	/*generate initial phase of pulse*/
	iniphase = (double(0)/stepTime)*pi;

	/*set intensity and EM parameters*/
	int_si = intensityVect[n_job]; //Selects intensity specific to ion
	intensity = int_si/int_au; //Conversion to au
	eAmpMax = sqrt((8.0*pi*intensity)/c); //Max amplitude for EM from Laser

	/*output declare*/
	/*This names the files based on their index, ion, and wavelength
	 *When running a different element, change "U" to whatever element being run to avoid confusing different files*/
	string file_str_Log = "logg_FILL_IN_" + getStringFromNumber(n_job) + "+" + getStringFromNumber(ionNum[n_job]) + "_" + getStringFromNumber(wavelength_nm) + "nm.txt";
	string file_str_Data = "data_FILL_IN_" + getStringFromNumber(n_job)  + "+" + getStringFromNumber(ionNum[n_job]) + "_" + getStringFromNumber(wavelength_nm) + "nm.dat";
	const char *file_char_Log = file_str_Log.c_str();
	const char *file_char_Data = file_str_Data.c_str();
	ofstream outLog(file_char_Log, ios::out);
	ofstream outData(file_char_Data, ios::out);

	/*set pulse time domain parameters*/
	double t_delta = period/4.0/stepTime;
	double t_integ_delta = period/stepTime;

	double omega = (2.0*pi*137)/wavelength;
	double a0 = eAmpMax/(c*omega);
	double gamma = (sqrt(2.0*ip[n_job]*pow(c,2.0))*pow(a0,3.0))/(16.0*omega);

	/*write log of project parameters*/
	outLog << ">>> wavelength = " << wavelength << endl;
	outLog << ">>> period(a.u.) = " << period << endl;
	outLog << ">>> freq(a.u.) = " << freq << endl;
	outLog << ">>> intensity(s.i.) = " << int_si << endl;
	outLog << ">>> intensity(a.u.) = " << intensity << endl;
	outLog << ">>> eAmpMax(a.u.) = " << eAmpMax << endl;
	outLog << ">>> t_delta(a.u.) = " << t_delta << endl;
	outLog << ">>> t_integ_delta(a.u.) = " << t_integ_delta << endl;
	outLog << ">>> range = 0 - 1/4Period" << endl;
	outLog << ">>> Ion Number: " << target_ion << endl;
	outLog << ">>> Gamma_r: " << gamma << endl;

	/*write title of rescatter data*/
	outData << "birth_phase" << " "
	<< "adk_rate_dt" << " "
	<< "return_phase" << " "
	<< "ini_X" << " "
	<< "ini_Y" << " "
	<< "ini_Z" << " "
	<< "res_X" << " "
	<< "res_Y" << " "
	<< "res_Z" << " "
	<< "return_kin" << endl;

	/*declare time domain variable*/
	double t_start = period/4.0-t_delta; //Different for every phase
	double t_final, t_integ_final;

	/*phase average*/
	for(int pp = 0; pp < stepTime; pp++)
	{

		/*current start*/
		t_start += t_delta;
		t_final  = t_start + 1.01*period;
		//Larger for safety in RKSuite
		t_integ_final = t_final + t_integ_delta;

		/*initialize trajectory*/
		for (int i = 0; i < neq; i++)
		{
			y_start[i] = 0.0;
		}

		/*get ADK rate*/
		rate_adk = GetADK(ip[n_job],
						  nstar[n_job],
						  c2nl[n_job],
						  flm[n_job],
						  y_start,
						  t_start);

		/*get field*/
		/*
		GetEmField(t_start,
		y_start[5],
		el_cpn,
		elfield);
		*/

		/*initialize vector to store initial position and momentum*/
		vector<double>* iniY = new vector<double>();
		vector<double>* iniZ = new vector<double>();
		vector<double>* iniPy = new vector<double>();
		vector<double>* iniPz = new vector<double>();

		/*calculate exit point*/
		traj_exit = GetExitPoint(y_start[5],
		t_start,
		ip[n_job],
		target_ion);

		/*get initial conditions for entire set of trajectories*/
		//Pancake-like model of electron distribution. 2D CLOUD?
		GetInitial(t_start,
		y_start[5],
		ip[n_job],
		nSample,
		seedVect[n_job],
		iniY,
		iniZ,
		iniPy,
		iniPz);

		/*loop through all trajectories*/
		for (int i = 0; i < nSample; i++)
		{

			/*get initial conditions for momentum and position*/
			y_start[0] = 0.0;
			y_start[1] = (*iniPy)[i];
			y_start[2] = (*iniPz)[i];
			y_start[3] = traj_exit;
			y_start[4] = (*iniY)[i];
			y_start[5] = (*iniZ)[i];

			/*pre-set rksuite.cpp ut*/
			rksuite.setup(neq,
			t_start,
			y_start,
			t_integ_final,
			tol,
			thres,
			method,
			"U",
			errass,
			hstart,
			mesage);

			/*initialize twant to the current phase start*/
			twant = t_start;

			/*initialize x_pre*/
			//x position at beginning of integration, assigned to trajectory exit
			x_pre = y_start[3];

			/*trajectory integration*/
			/*This while loop integrates until the break is reached or
			* the electron isn't rescattered by t_final
			*/
			while (twant <= t_final)
			{
				/*time step increment*/
				twant += t_integ_delta/6;

				/*rksuite integrate*/
				rksuite.ut(Derivs,twant,tnow,y,yp,ymax,uflag);

				/*if error*/
				if (uflag > 3)
				{
					cout << "UT uflag > 3 happens: twant = " << twant << endl;
				}

				/*update x_now*/
				x_now = y[3];

				/*check if rescatter happens*/
				if ((x_pre*x_now) < 0.0)
				{
					rescatter_time = tnow-t_start;
					/*This will print out the data if the electron rescatters
					* Then, it breaks the trajectory loop
					*/
					if (rescatter_time <= period)
					{

						rescatter_kin = sqrt(pow(y[0],2)+
						pow(y[1],2)+
						pow(y[2],2)+
						pow(c,2))*c-pow(c,2);

						outData << t_start/period*2.0*pi << " "
						<< rate_adk*t_delta << " "
						<< rescatter_time/period*2.0*pi << " "
						<< y_start[3] << " "
						<< y_start[4] << " "
						<< y_start[5] << " "
						<< y[3] << " "
						<< y[4] << " "
						<< y[5] << " "
						<< rescatter_kin << endl;

						break;
					}
				}

				/*update x_pre*/
				x_pre = x_now;

				/*end of trajectory integration*/
			}

		/*end of loop through all trajecotries*/
		}


	/*end of phase averaging*/
	}

	/*time*/
	elapsed2 = clock();
	float elapsed_diff = ((float)elapsed2-(float)elapsed1);
	outLog << "The code elapsed time is = "
	<< elapsed_diff/CLOCKS_PER_SEC << "seconds";

	/*close output*/
	outLog.close();
	outData.close();

	/*end of main program*/
	return 0;


}

/****************************************************************************/
/*function adk
  this function returns the adk rate*/
/****************************************************************************/
double GetADK(const double ip,
              const double nstar,
              const double c2nl,
              const double flm,
              const double y[neq],
              const double t)
{
	/*calculate ADK rate and ionization probability.
	The following derivation follows page 25~27
	in David Neal Fittinghoff's Ph.D. Thesis Dec 1993.*/

	/*variable declare*/
	double rate_adk, epsilon, factor, nmpower1;
	double e_cpn[neq], efield;

	/*get EM components and amplitude*/
	GetEmField(t, y[5], e_cpn, efield);

	/*check if EM == 0*/
	if (efield == 0.0)
	{
		rate_adk = 0.0;
	}
	else
	{
		/*calculate factors*/
		epsilon = pow((2.0*ip),1.5);
		factor = epsilon/fabs(efield);
		nmpower1 = 2.0*nstar-1.0; /*assume m=0 for all charge states*/

		/*calculate rate*/
		rate_adk = c2nl*sqrt(3.0/(pi*factor))*ip*flm*
		  (pow((2.0*factor),nmpower1))*exp(-(2.0*factor)/3.0);
	}
	return rate_adk;
/*end of function*/
}

/****************************************************************************/
/*function derivs
  this subroutine defines the equation of motion*/
/****************************************************************************/
void Derivs(double tgot,
            double ygot[neq],
            double ypgot[neq])
{
	/*declare EM field variable*/
	double e_cpn[neq];
	double eff;

	/*get EM field components*/
	GetEmField(tgot, ygot[5], e_cpn, eff);

	/*calculate gamma factor*/
	double gammasqrt = sqrt(pow(ygot[0],2)+
						  pow(ygot[1],2)+
						  pow(ygot[2],2)+
						  pow(c,2));

	/*define equation of motion - relativistic*/

	ypgot[3] = ygot[0]*c/gammasqrt;
	ypgot[4] = ygot[1]*c/gammasqrt;
	ypgot[5] = ygot[2]*c/gammasqrt;
	ypgot[0] = -e_cpn[0]-(ygot[1]*e_cpn[5]-ygot[2]*e_cpn[4])*c/gammasqrt;
	ypgot[1] = -e_cpn[1]-(ygot[2]*e_cpn[3]-ygot[0]*e_cpn[5])*c/gammasqrt;
	ypgot[2] = -e_cpn[2]-(ygot[0]*e_cpn[4]-ygot[1]*e_cpn[3])*c/gammasqrt;

	/*define equation of motion - non-relativistic*/
	/*
	ypgot[3] = ygot[0];
	ypgot[4] = ygot[1];
	ypgot[5] = ygot[2];
	ypgot[0] = -e_cpn[0]-(ygot[1]*e_cpn[5]-ygot[2]*e_cpn[4]);
	ypgot[1] = -e_cpn[1]-(ygot[2]*e_cpn[3]-ygot[0]*e_cpn[5]);
	ypgot[2] = -e_cpn[2]-(ygot[0]*e_cpn[4]-ygot[1]*e_cpn[3]);
	*/
/*end of subroutine*/
}

/****************************************************************************/
/*function emfield
  this subroutine calculates the EM field*/
/****************************************************************************/
void GetEmField(const double t,
                const double zz,
                double* e_cpn,
                double& eff)
{
	/*in this project we use homogeneous field, so E=E(t)*/
	/*calculate EM field amplitude*/
	double wn = 2*pi/wavelength;
	//eff = eAmpMax*sin(freq*t+iniphase-wn*zz)*exp(-pow(((t-zz/c)/(2*sigma)),2));
	eff = eAmpMax*sin(freq*t-wn*zz);
	/*calculate EM components*/
	e_cpn[0] = eff;
	e_cpn[1] = 0.0;
	e_cpn[2] = 0.0;
	e_cpn[3] = 0.0;
	/*e_cpn[4] = 0.0;*/
	e_cpn[4] = e_cpn[0]/c;
	e_cpn[5] = 0.0;
/*end of subroutine*/
}

/****************************************************************************/
/*function exitpoint
  this function returns the exit point*/
/****************************************************************************/
double GetExitPoint(const double z, 
		    const double t, 
		    const double ip, 
		    const double charge)
//CURRENTLY RETURNS 0?
{
	/*calculate EM field*/
	double e_cpn[neq], eff;
	GetEmField(t, z, e_cpn, eff);

	/*the electron ionized at opposite direction to the field*/
	double deterTerm = pow(ip,2)-4.0*abs(e_cpn[0])*charge;

	/*declare exitX*/
	double exitX;

	/*if EM field == 0*/
	if (e_cpn[0] == 0.0)
	{
		exitX = 0.0;
	}
	/*if EM field != 0*/
	else
	{
		/*if ATI regime*/
		if (deterTerm <= 0.0)
		{
			/*determine field direction*/
			if (e_cpn[0] > 0.0)
			{
				exitX = -ip/(2.0*abs(e_cpn[0]));
			}
			else
			{
				exitX =  ip/(2.0*abs(e_cpn[0]));
			}
		}
		/*if tunnel regime*/
		else
		{
			/*determine field direction*/
			if (e_cpn[0] > 0.0)
			{
				exitX = -(ip+sqrt(deterTerm))/(2.0*abs(e_cpn[0]));
			}
			else
			{
				exitX =  (ip+sqrt(deterTerm))/(2.0*abs(e_cpn[0]));
			}
		}
	}
	exitX = 0.0;
	return exitX;

/*end of function*/
}

/****************************************************************************/
/*function initial condition
  this function returns the initial position and momentum*/
/****************************************************************************/
void GetInitial(const double t,
		const double zz,
		const double ip,
		const int nSample,
		const int seed,
		vector<double>* iniY,
		vector<double>* iniZ,
		vector<double>* iniPy,
		vector<double>* iniPz){

  /*calculate EM field*/
  double e_cpn[neq], eff;
  GetEmField(t, zz, e_cpn, eff);

  /*calculate the spatial uncertainty width*/
  double yz_width = (pow(2*ip,0.25))/sqrt(2.0*abs(e_cpn[0]));

  /*declare normal distribution generator*/
  Normaldev ng1(0.0,yz_width,seed);
  Normaldev ng2(0.0,0.5/yz_width,seed);

  /*declare random generator*/
  Ran myran(seed);

  /*declare intermediate variables*/
  double angRange = 360.0;
  double del,delP,angle;

  for (int i = 0; i < nSample; i++) {
    del = abs(ng1.dev());
    angle = angRange*myran.doub();
    iniY->push_back(del*cos(angle/180.0*pi));
    iniZ->push_back(del*sin(angle/180.0*pi));

    delP = abs(ng2.dev());
    angle = angRange*myran.doub();
    iniPy->push_back(delP*cos(angle/180.0*pi));
    iniPz->push_back(delP*sin(angle/180.0*pi));
  }

  /*end of function*/
}

/****************************************************************************/
/*this function return string of input number*/
/****************************************************************************/
string getStringFromNumber(const int n){

  stringstream ss;
  ss << n;
  return ss.str();

  /*end of function*/
}
