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
#ifndef ODE_H
#define ODE_H
#include <cmath>
#include <vector>
#include "Inhibitory_Neuron.h"
#include "Pyramidal_Neuron.h"
#include "Reticular_Neuron.h"
#include "Thalamocortical_Neuron.h"

using std::vector;

void Iterate_ODE(vector<Pyramidal_Neuron>& PY,
				 vector<Inhibitory_Neuron>& IN,
				 vector<Thalamocortical_Neuron>& TC,
				 vector<Reticular_Neuron>& RE) {
	/* Parameters for the parallelization */
	extern const int N_Cores;

	/* First get all the RK terms */
	for (int i=0; i < 4; i++) {
		#pragma omp parallel for num_threads(N_Cores) schedule(dynamic)
		for(unsigned j=0; j < PY.size(); j++)
			PY[j].set_RK(i);

		#pragma omp parallel for num_threads(N_Cores) schedule(dynamic)
		for(unsigned j=0; j < IN.size(); j++)
			IN[j].set_RK(i);

		#pragma omp parallel for num_threads(N_Cores) schedule(dynamic)
		for(unsigned j=0; j < TC.size(); j++)
			TC[j].set_RK(i);

		#pragma omp parallel for num_threads(N_Cores) schedule(dynamic)
		for(unsigned j=0; j < RE.size(); j++)
			RE[j].set_RK(i);
	}

	/* Add the RK terms up*/
	#pragma omp parallel for num_threads(N_Cores) schedule(dynamic)
	for(unsigned j=0; j < PY.size(); j++)
		PY[j].add_RK();

	#pragma omp parallel for num_threads(N_Cores) schedule(dynamic)
	for(unsigned j=0; j < IN.size(); j++)
		IN[j].add_RK();

	#pragma omp parallel for num_threads(N_Cores) schedule(dynamic)
	for(unsigned j=0; j < TC.size(); j++)
		TC[j].add_RK();

	#pragma omp parallel for num_threads(N_Cores) schedule(dynamic)
	for(unsigned j=0; j < RE.size(); j++)
		RE[j].add_RK();
}

#endif // ODE_H
