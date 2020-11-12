/*
 * This is the implementation file of class "Collection1D".
 * Used for:bin energy spectrum, bin angular distribution,
 *          convoluate distribution
 *          store intensity dependence.
 * Author:Sui Luo.
 * Version:First written in Feb 05, 2013.
 *         Last modified in June 24, 2019.
 *         Version 0.8
 * Extra comments have been added by David Milliken on 6/24/2019 to add readability and ease of understanding
 */

#include "collection1D.h"

//constructor
Collection1D::Collection1D()
{
	collection = new vector<Data1D>;
	collection_out = new vector<Data1D>;
}

//destructor
Collection1D::~Collection1D()
{
	collection->clear();
	delete collection;
	collection_out->clear();
	delete collection_out;
}

//get min, max, dimension, limit and # of bin
void Collection1D::getParameter(double d1,
				double d2,
				int i1,
				double d3,
				double d4)
{
	//set parameters
	externalMinMax = true;
	min = d1;
	max = d2;
	dim = i1;
	cov1 = d3;
	cov2 = d4;
}

//get min, max, dimension, limit and # of bin - overloading
void Collection1D::getParameter(int i1,
				double d3,
				double d4)
{
	//set parameters
	externalMinMax = false;
	dim = i1;
	cov1 = d3;
	cov2 = d4;
}

//get the log file ofstream reference
void Collection1D::getLog(Logger l1)
{
	log = l1;
}

//get the output file name
void Collection1D::getOutputName(string s1)
{
	fileInfo = s1;
}

//add data to collection
void Collection1D::addData(Data1D D1)
{
	collection->push_back(D1);
}

//sort collection according the key
void Collection1D::sortC()
{
	/* sorts the data based on the key (energy) */
	sort(collection->begin(), collection->end(), by_key());
	log.writeLog("sort is done.\n");

	//find min/max if external min/max not given
	if (!externalMinMax)
	{
		min = collection->front().getKey();
		max = collection->back().getKey();
	}
}

//discretize the data in collection
void Collection1D::binC()
{
	//set bin size
	double binSize = (max - min)/((double)dim);
	double current = min;
	double binKey, binWeight;

	//initialize iterator
	vector<Data1D>::iterator iv = collection->begin();

	//loop starts
	while ((current + binSize) <= max)
	{
		current += binSize;
		binKey = current - binSize/2.0;
		binWeight = 0.0;

		if (current > max)
		{
			log.writeLog("boundary overflow.\n");
			exit(1);
		}

		while ((iv->getKey() <= current) && (iv < collection->end()))
		{
			binWeight += iv->getWeight();
			iv++;

		}

		//record bin data
		Data1D tempData(binKey, binWeight);
		collection_out->push_back(tempData);

	} /* end of outer while loop */

	log.writeLog("bin is done.\n");
}

//convolute the container based on *percentage*
void Collection1D::convolution1C()
{
	//set bin size
	double binSize = (max - min)/((double)dim);

	//loop starts
	for (int i = 1; i <= dim-1; i++)
	{
		double current = min + i*binSize;
		double binWeight = 0.0;

		for (vector<Data1D>::iterator iv = collection->begin();
				iv < collection->end();
				iv++)
		{
			if (iv->getKey() >= (1.0-cov1)*current && iv->getKey() <= (1.0+cov1)*current)
			{
				binWeight += iv->getWeight();
			}
		}

		//add new data
		Data1D tempData(current,binWeight);
		collection_out->push_back(tempData);
	} /* end of outer for loop */

	log.writeLog("convolution_1 is done.\n");
}

//convolute the container based on *fix-value*
void Collection1D::convolution2C()
{

	//set bin size
	double binSize = (max - min)/((double)dim);

	//loop starts
	for (int i = 1; i <= dim-1; i++)
	{
		double current = min + i*binSize;
		double binWeight = 0.0;

		for (vector<Data1D>::iterator iv = collection->begin();
				iv < collection->end();
				iv++)
		{
			if (iv->getKey() >= (current-cov2) && iv->getKey() <= (current+cov2))
			{
				binWeight += iv->getWeight();
			}
		}

		//add new data
		Data1D tempData(current,binWeight);
		collection_out->push_back(tempData);
	} /* end of outer for loop */

	log.writeLog("convolution_2 is done.\n");
}

//output the unbinned collection
void Collection1D::output()
{
	const char *filename_char = fileInfo.c_str();
	ofstream output(filename_char, ios::out);

	output << "Energy" << " " << "Flux" << endl;

	for (vector<Data1D>::iterator iv = collection->begin();
			iv < collection->end();
			iv++)
	{
		output << iv->toString();
	}

	output.close();
	log.writeLog("output unbinned is done.\n");
}

//output the binned collection
void Collection1D::outputbinC()
{
	const char *filename_char = fileInfo.c_str();
	ofstream output(filename_char, ios::out);

	output << "Energy" << " " << "Flux" << endl;

	for (vector<Data1D>::iterator iv = collection_out->begin();
			iv < collection_out->end();
			iv++)
	{
		output << iv->toString();
	}

	output.close();
	log.writeLog("output binned is done.\n");
}

//run the file
void Collection1D::run()
{
	if (collection->size() == 0)
	{
		cout << "collection size = 0. exit" << endl;
		return;
	}
	sortC();
	//binC();
	//convolution1C();//for energy spectrum
	//convolution2C();//for azimuthal distribution
	customize();
	//output();
	outputbinC();
}

//customize
void Collection1D::customize()
{
	/*hard-code base*/
	double base = 1.08; //1.08

	/*temporal vectors*/
	vector<int> tmpKey;
	vector<double> tmpWeight;

	/*convert to int*/
	for (vector<Data1D>::iterator iv = collection->begin();
			iv < collection->end();
			iv++)
	{
		/* Here iv is of the type iterator, meaning it iterates through each value in Key when the ++ function is used */
		tmpKey.push_back((int)(log10(iv->getKey())/log10(base)+0.5)); //takes the log base 1.4 of the key (Energy) and adds it to tmpKey
		tmpWeight.push_back(iv->getWeight()); // adds the weight (Flux) to tmpWeight




	} /* end of for loop */

	int dim = tmpKey[tmpKey.size()-1]-tmpKey[0]+1;

	/*temporal vectors*/
	vector<int> conKey(dim, 0);
	for (int i = 0; i < conKey.size(); i++)
	{
		conKey[i] = tmpKey[0] + i;
	}
	vector<double> conWeight(conKey.size(), 0.0);

	/*concatenate*/
	int itr = 0;
	double minWeight = 100000.0;
	for (int i = 0; i < conKey.size(); i++)
	{

		while (itr < tmpKey.size() && conKey[i]==tmpKey[itr])
		{
			conWeight[i] += tmpWeight[itr];
			itr++;
			/* This gathers all of the data that fits within the same key */
		}
		/*update min weight*/
		if (conWeight[i] > 0.0 && conWeight[i] < minWeight)
		{
			minWeight = conWeight[i];
		}

	}


	/*interpolate*/
	itr = 0;
	int preItr;
	while (itr < conKey.size())
	{
		if (conWeight[itr] < minWeight) /* checks if our flux value is less than the minimum (should happen when it's 0) */
		{
			preItr = itr - 1;
			while (conWeight[itr] < minWeight) /* checks to the next energy bins until one isn't 0 */
			{
				itr++;
			}

			double del = (conWeight[itr]-conWeight[preItr])/(itr-preItr); /* find the delta from before and after the line of 0s */
			preItr++;

			while (preItr < itr)
			{
				conWeight[preItr] += conWeight[preItr-1] + del; /* fills the empty bins with flux at steps of size delta */
				preItr++;
			}
		}
		else
		{
			itr++;
		}
	} /* end of outer while loop */



	/*record in output container*/
	for (int i = 0; i < conKey.size(); i++)
	{
		Data1D tempData(pow(base,conKey[i]),
				conWeight[i]/(pow(base,conKey[i])-pow(base,conKey[i]-1))); /* returns the bin values to their initial values from the 1.4 scaled values */
		collection_out->push_back(tempData);
		stringstream ss, rr, re;
		ss << conKey[i];
		rr << conWeight[i];
		re << pow(base,conKey[i]);
		log.writeLog(ss.str());
		log.writeLog("\t");
		log.writeLog(re.str());
		log.writeLog("\t");
		log.writeLog(rr.str());
		log.writeLog("\n");
	}
}

bool Collection1D::double_equals(double a, double b, double epsilon)
{
	if(abs(a-b) > epsilon)
	{
		return false;
	}
	else
	{
		return true;
	}
}

