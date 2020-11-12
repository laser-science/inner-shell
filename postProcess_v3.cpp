/*
 * This program reads batch electron rescatter files
 * and generate the mean and std of wavepacket.
 * bin and normalize.
 * Author:Sui Luo. 
 * Last modified: 8/16/2019.
 * Commented and Edited last by David Milliken
 * We used data from cutoff_V5.cpp to generate the data processed here
 */


/*********************** Header ***********************/
#include "data1D.h"
#include "collection1D.h"
#include "logger.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <cmath>
#include <cstdlib>

using namespace std;

/***************** Declare Subroutines ****************/
void procMean(string, string, string, string);
void getMean(vector<double>, 
		vector<double>,
		vector<double>,
		vector<double>,
		vector<double>,
		vector<double>,
		double&, double&,
		double&, double&,
		double&, double&);
double getFlux(double, double);
vector<vector<double> > binEnergy(vector<vector<double> >, double); 
void binAndNormalize(string, vector<vector<double> >, string, string);
//void binAndNormalize2(string, vector<vector<double> >); /* binAndNormalize2 wasn't used by us and was commented out upon finding the file */

/******************** Main Program ********************/
/*
 * When run, hard code the wavelength and what range the for loop should go through (each of the flags of sets being run)
 */

int main()
{
	/*generate input file name&path*/
	string path_str;
	path_str="FILL_IN"; /* Put here the path to the folder your input data is in */

	/*generate input file*/
	string file_str;
	for (int i = 0; i <= 0; i++) //the range which i runs through should be the range of indexes from your input files
								 //this is the first integer which appears in the file name after the element
	{
		stringstream ind,ion,wavelength;
		ind << i; //ind will be the index value currently being called
		ion << FILL_IN;
		wavelength << FILL_IN; //input the wavelength you're running
		file_str = path_str + "data_FILL_IN_" + ind.str() + "+" + ion.str() + "_" +  wavelength.str() + "nm.dat";
		procMean(file_str,ind.str(),ion.str(),wavelength.str());
	}

	/*end of the program*/
	return 0;
}

/******************** Subroutines *********************/

/******************************************************/

void procMean(string file_str, string flag, string ion, string wavelength) {

	/*generate output file and title columns */
	string out_str = "out_FILL_IN_" + flag + "+" + ion + "_" + wavelength + "nm.dat";
	const char *output_file_char = out_str.c_str();
	ofstream output(output_file_char, ios::out);
	output << "release_phase" << " "
			<< "return_phase" << " "
			<< "adk_rate_dt" << " "
			<< "ini_deflection" << " "
			<< "ini_width" << " "
			<< "Wavepacket_deflection" << " "
			<< "Wavepacket_width" << " "
			<< "Gamma_r" << " "
			<< "Return_kinetic" << " "
			<< "Return_flux_ADKPopu" << " "
			<< "Return_kinetic/dE" << endl;

	/*generate target file name with full path&name*/
	const char *input_file_char = file_str.c_str();

	/*open file*/
	ifstream input;
	input.open(input_file_char);

	/*check existance of file*/
	if(!input)
	{
		cout << "Error: file could not be opened >> " + file_str << endl;
		exit(1);
	}

	/* create container to hold our data */
	vector<vector<double> > datas;

	/*temperoral variable*/
	string t0, t1, t2, t3, t4, t5, t6, t7, t8, t9;//title in each file
	double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9;//data in each file

	/*read in titles*/
	input >> t0 >> t1 >> t2 >> t3 >> t4 >> t5 >> t6 >> t7 >> t8 >> t9;

	/*variables to record*/
	double tstart = 0.0;
	double rate = 0.0; //here 'rate' input is actually adk_rate*dt
	vector<double> vt, vy, vz, vyi, vzi, vkin;
	bool firstRow = true;
	//int loops = 0; /* Loops was used when input data was corrupt to find where the function stopped */

	/*read in all data in one file*/
	while(input >> d0 >> d1 >> d2 >> d3 >> d4 >> d5 >> d6 >> d7 >> d8 >> d9)
	{
		//loops++;
		//cout << loops << endl;

		if (!firstRow && abs(d0-tstart) > 0.000000000001) /* This if statement ensures we do everything after this at least once before what's inside the if statment */
		{
			/* It is a complex version of a do, while loop */
			double mean, std, meani, stdi, kin, rePhase;
			getMean(vt,vy,vz,vyi,vzi,vkin,mean,std,meani,stdi,kin,rePhase);

			/*get new record*/
			vector<double> tmpdata;
			tmpdata.push_back(tstart); // release phase
			tmpdata.push_back(rePhase); // return phase
			tmpdata.push_back(rate); // adk rate dt
			tmpdata.push_back(meani); // initial deflection
			tmpdata.push_back(stdi); // initial width
			tmpdata.push_back(mean); // wavepacket deflection
			tmpdata.push_back(std); // wavepacket width
			tmpdata.push_back(mean*mean/(2*std*std)); //Gamma_r
			tmpdata.push_back(kin); // return kinetic energy
			tmpdata.push_back(getFlux(mean,std)*rate); // return flux ADK population
			datas.push_back(tmpdata); // places the vector of data into datas

			/*reset*/
			vt.clear();
			vyi.clear();
			vzi.clear();
			vy.clear();
			vz.clear();
			vkin.clear();
		} /* end of if */

		tstart = d0;
		if (firstRow) /* check to see if we should begin to go into the previous if statement for future iterations */
		{
			if (tstart > 0.0)
			{
				firstRow = false;
			}
		} /* end of double if */

		/* add data to our recorded variables */
		rate = d1; // adk rate dt
		vt.push_back(d2); // return phase
		vyi.push_back(d4); // initial y
		vzi.push_back(d5); // initial z
		vy.push_back(d7); // return y
		vz.push_back(d8); // return z
		vkin.push_back(d9); // return kinetic energy

	} /* end of while loop */


	/* close input */
	input.close();
	cout << "End of read in file >> " + file_str << endl;

	/*get /dE flux and add it to each vector in datas*/
	datas[0].push_back(datas[0][9]/abs(datas[1][8]-datas[0][8]));
	/* Note: now [i][8] is our energy and [i][9] is flux, unlike in our input data */
	for (int i = 1; i < datas.size(); i++)
	{
		datas[i].push_back(datas[i][9]/abs(datas[i][8]-datas[i-1][8]));
	}


	/*output*/
	for (int i = 0; i < datas.size(); i++)
	{
		for (int j = 0; j < datas[i].size(); j++)
		{
			output << datas[i][j] << " ";
		}
		output << endl;
	}


	output.close();

	/*bin energy/dE*/
  
	/* create and open second output */
	string out_str2 = "bin_FILL_IN_" + flag + "+" + ion + "_" + wavelength + "nm.dat";
	const char *output_file_char2 = out_str2.c_str();
	ofstream output2(output_file_char2, ios::out);
	output2 << "energy/dE_bin" << " "
			<< "Return_flux_ADKPopu" << endl;

	/* bin data and write it in second output */
	double binSize = 1.0;
	vector<vector<double>> bins = binEnergy(datas,binSize);
	for (int i = 0; i < bins.size(); i++)
	{
		output2 << bins[i][0] << " "
				<< bins[i][1] /*<< " " << bins[i][2]*/ << endl;
		/* the third value can be added to check to see how many points are put in each bin, but that was initially just used to bug check */
	}
	output2.close();

     
	/*bin and normalize*/
	binAndNormalize(flag, datas, ion, wavelength);
	//binAndNormalize2(flag,datas);

	/*end of the subroutine*/
	return;

}

/******************************************************/
// finds the means and stand deviations
void getMean(vector<double> vt,
		vector<double> vy,
		vector<double> vz,
		vector<double> vyi,
		vector<double> vzi,
		vector<double> vkin,
		double &mean, double &std,
		double &meani, double &stdi,
		double &kin, double &rePhase) {

	/*mean time*/
	double mt = 0.0;
	for (int i = 0; i < vt.size(); i++)
	{
		mt += vt[i];
	}
	mt /= vt.size();
	rePhase = mt;

	/*mean return y*/
	double my = 0.0;
	for (int i = 0; i < vy.size(); i++)
	{
		my += vy[i];
	}
	my /= vy.size();

	/*mean return z*/
	double mz = 0.0;
	for (int i = 0; i < vz.size(); i++)
	{
		mz += vz[i];
	}
	mz /= vz.size();

	mean = sqrt(my*my + mz*mz);

	/*std return y*/
	double sqry = 0.0;
	for (int i = 0; i < vy.size(); i++)
	{
		sqry += (vy[i]-my)*(vy[i]-my);
	}

	/*std return z*/
	double sqrz = 0.0;
	for (int i = 0; i < vz.size(); i++)
	{
		sqrz += (vz[i]-mz)*(vz[i]-mz);
	}

	std = sqrt((sqry+sqrz)/vy.size());

	/*mean initial y*/
	double myi = 0.0;
	for (int i = 0; i < vyi.size(); i++)
	{
		myi += vyi[i];
	}
	myi /= vyi.size();

	/*mean initial z*/
	double mzi = 0.0;
	for (int i = 0; i < vzi.size(); i++)
	{
		mzi += vzi[i];
	}
	mzi /= vzi.size();

	meani = sqrt(myi*myi + mzi*mzi);

	/*std initial y*/
	double sqryi = 0.0;
	for (int i = 0; i < vyi.size(); i++)
	{
		sqryi += (vyi[i]-myi)*(vyi[i]-myi);
	}

	/*std initial z*/
	double sqrzi = 0.0;
	for (int i = 0; i < vzi.size(); i++)
	{
		sqrzi += (vzi[i]-mzi)*(vzi[i]-mzi);
	}

	stdi = sqrt((sqryi+sqrzi)/vyi.size());

	/*mean return kin*/
	double mkin = 0.0;
	for (int i = 0; i < vkin.size(); i++)
	{
		mkin += vkin[i];
	}
	mkin /= vkin.size();
	kin = mkin;
}

/******************************************************/
double getFlux(double mean, double std) {
	double pi = 4.0*atan(1.0);
	double flux;
	/*2D Gaussian wave function*/
	flux = 1.0/(2.0*pi*std*std)*exp(-mean*mean/(2.0*std*std));

	return flux;
}

/******************************************************/
vector<vector<double> > binEnergy(vector<vector<double> > datas, 
		double binSize) {
	/* Find max energy */
	double maxEner = 0.0;
	for (int i = 0; i < datas.size(); i++)
	{
		if (datas[i][8] > maxEner)
		{
			maxEner = datas[i][8];
		}
	}

	int dim = maxEner/binSize + 1;

	/* Create bins */
	vector<vector<double> > bins;
	for (int i = 0; i < dim; i++)
	{
		vector<double> temp = {(double)i+binSize/2.0,0.0};
		bins.push_back(temp);
	}

	/*map bin*/
	for (int i = 0; i < datas.size(); i++)
	{
		int binNum = datas[i][8]/binSize;
		bins[binNum][1] += datas[i][10];
	}

	return bins;
}

/******************************************************/
void binAndNormalize(string flag,
		vector<vector<double> > datas, string ion, string wavelength) {


	/*container*/
	Collection1D bound1;
	Collection1D bound2;

	/*variables*/
	bool toBound2 = false;
	double preEner = 0.0;
	double minFlux = 1.0e-500;

	/*get the sum of adk_rate_td*/
	double adk_rate_dt_sum = 0.0;
	for (int i = 0; i < datas.size(); i++)
	{
		adk_rate_dt_sum += datas[i][2];
	}

	/*divide and conquer*/
	for (int i = 0; i < datas.size(); i++)
	{
		if (datas[i][10] > minFlux)
		{
			if (datas[i][8] < preEner && !toBound2) //checks to see if the data is a short or long scattering
			{
				toBound2 = true;
			}

			/* adds to the appropriate data set the values */
			if (toBound2)
			{
				Data1D newdata0(datas[i][8],datas[i][10]/adk_rate_dt_sum);
				bound2.addData(newdata0);
			}
			else
			{
				Data1D newdata0(datas[i][8],datas[i][10]/adk_rate_dt_sum);
				bound1.addData(newdata0);
			}
			preEner = datas[i][8];
		}
	} /* end for loop */

	/*open log file*/
	string log_file_str="binNorm_FILL_IN_log_" + flag + "+" + ion + "_" + wavelength + "nm_Bound.txt";
	const char *log_file_char=log_file_str.c_str();
	ofstream logfile(log_file_char);
	Logger logger1(logfile);

	/*output string*/
	string out_str_1 = "binNorm_FILL_IN_" + flag + "+" + ion + "_" + wavelength + "nm_bound1.dat";
	string out_str_2 = "binNorm_FILL_IN_" + flag + "+" + ion + "_" + wavelength + "nm_bound2.dat";

	/* bin, normalize, and then print to the output files our data */
	bound1.getOutputName(out_str_1);
	bound1.getLog(logger1);
	bound1.run();

	bound2.getOutputName(out_str_2);
	bound2.getLog(logger1);
	bound2.run();
}

/******************************************************/
//void binAndNormalize2(string flag,
//    vector<vector<double> > datas) {
//  /*container*/
//  Collection1D bound;
//
//  /*variables*/
//  double minFlux = 1.0e-30;
//
//  /*divide and conquer*/
//  for (int i = 0; i < datas.size(); i++) {
//    if (datas[i][10] > minFlux) {
//      Data1D newdata0(datas[i][8],datas[i][10]);
//      bound.addData(newdata0);
//    }
//  }
//
//  /*open log file*/
//  string log_file_str="binNorm_Ne_log_Accum.txt";
//  const char *log_file_char=log_file_str.c_str();
//  ofstream logfile(log_file_char);
//  Logger logger1(logfile);
//
//  /*output string*/
//  string out_str_1 = "binNorm_Ne_" + flag + "_accum.dat";
//
//  bound.getOutputName(out_str_1);
//  bound.getLog(logger1);
//  bound.run();
//}
