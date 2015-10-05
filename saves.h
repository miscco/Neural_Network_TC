/*
*	Copyright (c) 2015 Michael Schellenberger Costa mschellenbergercosta@gmail.com
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
/*									Functions for data storage										*/
/****************************************************************************************************/
#pragma once
#include <vector>
#include "Pyramidal_Neuron.h"
#include "Inhibitory_Neuron.h"
using std::vector;
/****************************************************************************************************/
/*											Save data												*/
/****************************************************************************************************/
inline void get_data(int counter, vector<Pyramidal_Neuron>& PY, vector<Inhibitory_Neuron>& IN,
					 vector<double*> pData) {
	/* Parameters for the parallelization */
	extern const int N_Cores;
	extern const int N_e, N_i;
	/* NOTE As C++ and Matlab have a different storage order (Row-major vs Column-major), the index
	 * has to be adapted! For an NxM matrix A, element A(i,j) is accessed by A(j+i*M) rather than
	 * the usual A(i+j*N)
	 */
	#pragma omp parallel for num_threads(N_Cores) schedule(static)
	for(int i=0; i<N_e; ++i)
		pData[0][i+N_e*counter] = PY[i].Vs[0];
	#pragma omp parallel for num_threads(N_Cores) schedule(static)
	for(int i=0; i<N_i; ++i)
		pData[1][i+N_i*counter] = IN[i].V [0];
	#pragma omp parallel for num_threads(N_Cores) schedule(static)
	for(int i=0; i<N_e; ++i)
		pData[2][i+N_e*counter] = PY[i].Ca[0];
}
/****************************************************************************************************/
/*										 		end													*/
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
