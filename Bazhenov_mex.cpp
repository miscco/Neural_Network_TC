/*
*	Copyright (c) 2016 Michael Schellenberger Costa mschellenbergercosta@gmail.com
*
*	Permission is hereby granted, free of charge, to any person obtaining a copy
*	of this software and associated documentation files (the "Software"), to deal
*	in the Software without restriction, including without limitation the rights
*	to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*	copies of the Software, and to permit persons to whom the Software is
*	furnished to do so, subject to the following conditions:
*
*	The above copyright notice and this permission notice shall be included in
*	all copies or substantial portions of the Software.
*
*	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
*	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
*	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
*	THE SOFTWARE.
*/

/****************************************************************************************************/
/* 		Implementation of the Bazhenov2002 model as MATLAB routine (mex compiler)					*/
/* 		mex command is given by:																	*/
/* 		mex CXXFLAGS="\$CXXFLAGS -std=c++11 -O3 -fopenmp"											*/
/*			Bazhenov.cpp Pyramidal_Neuron.cpp Inhibitory_Neuron.cpp									*/
/*		The Simulation requires the following boost libraries:	Random								*/
/****************************************************************************************************/
#include "mex.h"
#include "matrix.h"

#include "Data_Storage.h"
#include "Initialize_Neurons.h"
#include "Iterate_ODE.h"

/****************************************************************************************************/
/*										Fixed simulation settings									*/
/****************************************************************************************************/
extern const int T		= 10;								/* Simulation length in s				*/
extern const int res 	= 5E4;								/* Number of iteration steps per s		*/
extern const int red 	= 1E1;								/* Fraction of stime steps saved		*/
extern const double dt 	= 1E3/res;							/* Duration of a timestep in ms			*/
extern const int N_e	= 128;								/* Number of pyramidal  cells			*/
extern const int N_i	= 32;								/* Number of inhibitory cells			*/
extern const int N_t	= 128;								/* Number of thalamocortical cells		*/
extern const int N_r	= 32;								/* Number of thalamic reticular cells	*/
extern const int N_Cores= 7;								/* Number of CPU cores					*/
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/
mxArray* SetMexArray(int N, int M);

/****************************************************************************************************/
/*										Simulation routine	 										*/
/*										lhs defines outputs											*/
/*										rhs defines inputs											*/
/****************************************************************************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	/* Set the seed */
	srand(time(NULL));

	/* Initialize the populations */
	vector<Pyramidal_Neuron> PY(0);
	vector<Inhibitory_Neuron> IN(0);
	vector<Thalamocortical_Neuron> TC(0);
	vector<Reticular_Neuron> RE(0);
	setupNetwork(PY, IN, TC, RE);

	/* Data container in MATLAB format */
	vector<mxArray*> Data;
	Data.push_back(SetMexArray(N_e, T*res/red));	// Ve
	Data.push_back(SetMexArray(N_i, T*res/red));	// Vi
	Data.push_back(SetMexArray(N_e, T*res/red));	// Ca

	/* Pointer to the data blocks */
	vector<double*> pData(Data.size(), NULL);
	for(unsigned i=0; i<Data.size(); ++i)
		pData[i] = mxGetPr(Data[i]);

	/* Simulation */
	int count = 0;
	for (int t=0; t< T*res; ++t) {
		Iterate_ODE(PY, IN);
		if(t%red==0){
			get_data(count, PY, IN, pData);
			count++;
		}

	}

	/* Return the data containers */
	for(unsigned i=0; i<Data.size(); ++i)
		plhs[i] = Data[i];

	return;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*									Create MATLAB data container									*/
/****************************************************************************************************/
mxArray* SetMexArray(int N, int M) {
	mxArray* Array	= mxCreateDoubleMatrix(0, 0, mxREAL);
	mxSetM(Array, N);
	mxSetN(Array, M);
	mxSetData(Array, mxMalloc(sizeof(double)*M*N));
	return Array;
}
/****************************************************************************************************/
/*										 		end													*/
/****************************************************************************************************/
