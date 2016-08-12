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
#ifndef INITIALIZE_Neurons_H
#define INITIALIZE_Neurons_H
# define M_PI           3.14159265358979323846  /* pi */
#include <cmath>
#include <exception>
#include <vector>

#include "Random_Stream.h"

#include "Inhibitory_Neuron.h"
#include "Pyramidal_Neuron.h"
#include "Reticular_Neuron.h"
#include "Thalamocortical_Neuron.h"

using std::vector;

enum neuronType {
	PYRAMIDAL = 0,
	INHIBITORY,
	THALAMOCORTICAL,
	RETICULAR
};

static vector<vector<double>> getParameterDistribution(neuronType Type) {
	srand(time(NULL));

	extern const vector<int> NumCells;
	vector<vector<double>> param;
	vector<double> mP, dP;

	/* Get the distributions of the initial parameters */
	switch(Type) {
	case PYRAMIDAL:
		param = vector<vector<double>>(NumCells[PYRAMIDAL], vector<double>(3));
		/* E_L, g_L , g_sd */
		mP	= {-60.95, 66.7E-3, 1.75E-3};
		dP	= {0.3, 6.7E-3, 0.1E-3};
		break;
	case INHIBITORY:
		param= vector<vector<double>>(NumCells[INHIBITORY], vector<double>(2));
		/* E_L, g_L */
		mP	= {-63.8, 102.5E-3};
		dP	= {0.15, 2.5E-3};
		break;
	case THALAMOCORTICAL:
		param = vector<vector<double>>(NumCells[THALAMOCORTICAL], vector<double>(2));
		/* E_L, g_L */
		mP	= {-60.95, 66.7E-3};
		dP	= {0.3, 6.7E-3};
		break;
	case RETICULAR:
		param= vector<vector<double>>(NumCells[RETICULAR], vector<double>(2));
		/* E_L, g_L */
		mP	= {-63.8, 102.5E-3};
		dP	= {0.15, 2.5E-3};
		break;
	default:
		throw std::runtime_error("Unknown neuron type!");
	}

	/* Initialize the random number generators */
	vector<random_stream_normal> MTRand;
	for(unsigned i=0; i<mP.size(); i++) {
		MTRand.push_back(random_stream_normal(mP[i], dP[i]));
	}

	/* Get the randomly distributed parameters */
	for(unsigned i=0; i<param.size(); i++)
		for(unsigned j=0; j<mP.size(); ++j)
			param[i][j] = MTRand[j]();

	return param;
}

template<class NEURON>
static vector<NEURON> initializeNeurons(neuronType type) {
	extern const vector<int> NumCells;

	/* Generate the randomly distributed parameter values */
	vector<vector<double>> param = getParameterDistribution(type);

	/* Initialize the neurons */
	vector<NEURON> neurons;
	for (int i=0; i<NumCells[type]; i++)
		neurons.push_back(NEURON(param[i]));

	return neurons;
}


static vector<vector<int>> getConnectivity(neuronType post, neuronType pre) {
	extern const vector<int> NumCells;
	double length = 5*NumCells[PYRAMIDAL];
	/* Sigma for the normal distribution */
	double sigma;
	switch (pre) {
	case PYRAMIDAL:
		sigma = 250/length*NumCells[post];
		break;
	case INHIBITORY:
		sigma = 125/length*NumCells[post];
		break;
	case THALAMOCORTICAL:
		sigma = 125/length*NumCells[post];
		break;
	case RETICULAR:
		sigma = 125/length*NumCells[post];
		break;
	default:
		throw std::runtime_error("Unknown connection type!");
	}

	/* Generate the random number generator for the number of connections */
	random_stream_normal MTRand_N  = random_stream_normal(20, 5);

	/* Generate the random number generator for the target neurons */
	random_stream_normal MTRand_T = random_stream_normal(0, sigma);

	vector<vector<int>> connectivity(NumCells[post], vector<int>(0));
	for (int i=0; i < NumCells[pre]; i++) {
		unsigned N_con = (abs((int)MTRand_N()));
		for(unsigned j=0; j < N_con; j++) {
			int Target;
			/* Self connections are not allowed */
			do {
				Target = ((int)MTRand_T()+i)%NumCells[post];
				if(Target < 0)
					Target += NumCells[post];
			} while (Target == i && pre == post);
			connectivity[Target].push_back(i);
		}
	}
	return connectivity;
}

void connectNeurons(vector<Pyramidal_Neuron>& PY,
					vector<Inhibitory_Neuron>& IN,
					vector<Thalamocortical_Neuron>& TC,
					vector<Reticular_Neuron>& RE) {
	/* Generate random connectivity matrices. For every Neuron[i] they store
	 * the index of all neurons it RECEIVES input from
	 */
	vector<vector<int>> conPP = getConnectivity(PYRAMIDAL, PYRAMIDAL);
	vector<vector<int>> conPI = getConnectivity(PYRAMIDAL, INHIBITORY);
	vector<vector<int>> conPT = getConnectivity(PYRAMIDAL, THALAMOCORTICAL);
	for (unsigned i=0; i < PY.size(); i++) {
		for (unsigned j=0; j < conPP[i].size(); j++)
			PY[i].PY_Con.push_back(&PY[conPP[i][j]]);
		for (unsigned j=0; j < conPI[i].size(); j++)
			PY[i].IN_Con.push_back(&IN[conPI[i][j]]);
		for (unsigned j=0; j < conPT[i].size(); j++)
			PY[i].TC_Con.push_back(&TC[conPT[i][j]]);
	}

	vector<vector<int>> conIP = getConnectivity(INHIBITORY, PYRAMIDAL);
	vector<vector<int>> conII = getConnectivity(INHIBITORY, INHIBITORY);
	vector<vector<int>> conIT = getConnectivity(INHIBITORY, THALAMOCORTICAL);
	for (unsigned i=0; i < IN.size(); i++) {
		for (unsigned j=0; j < conIP[i].size(); j++)
			IN[i].PY_Con.push_back(&PY[conIP[i][j]]);
		for (unsigned j=0; j < conII[i].size(); j++)
			IN[i].IN_Con.push_back(&IN[conII[i][j]]);
		for (unsigned j=0; j < conIT[i].size(); j++)
			IN[i].TC_Con.push_back(&TC[conIT[i][j]]);
	}

	vector<vector<int>> conTR = getConnectivity(THALAMOCORTICAL, RETICULAR);
	vector<vector<int>> conTP = getConnectivity(THALAMOCORTICAL, PYRAMIDAL);
	for (unsigned i=0; i < TC.size(); i++) {
		for (unsigned j=0; j < conTR[i].size(); j++)
			TC[i].RE_Con.push_back(&RE[conTR[i][j]]);
		for (unsigned j=0; j < conTP[i].size(); j++)
			TC[i].PY_Con.push_back(&PY[conTP[i][j]]);
	}

	vector<vector<int>> conRT = getConnectivity(RETICULAR, THALAMOCORTICAL);
	vector<vector<int>> conRR = getConnectivity(RETICULAR, RETICULAR);
	vector<vector<int>> conRP = getConnectivity(RETICULAR, PYRAMIDAL);
	for (unsigned i=0; i < RE.size(); i++) {
		for (unsigned j=0; j < conRT[i].size(); j++)
			RE[i].TC_Con.push_back(&TC[conRT[i][j]]);
		for (unsigned j=0; j < conRR[i].size(); j++)
			RE[i].RE_Con.push_back(&RE[conRR[i][j]]);
		for (unsigned j=0; j < conRP[i].size(); j++)
			RE[i].PY_Con.push_back(&PY[conRP[i][j]]);
	}
}


void setupNetwork(vector<Pyramidal_Neuron>& PY,
				  vector<Inhibitory_Neuron>& IN,
				  vector<Thalamocortical_Neuron>& TC,
				  vector<Reticular_Neuron>& RE) {
	/* Seed the random number generator */
	srand(time(NULL));

	/* Initialize the individual neurons */
	PY = initializeNeurons<Pyramidal_Neuron>(PYRAMIDAL);
	IN = initializeNeurons<Inhibitory_Neuron>(INHIBITORY);
	TC = initializeNeurons<Thalamocortical_Neuron>(THALAMOCORTICAL);
	RE = initializeNeurons<Reticular_Neuron>(RETICULAR);

	connectNeurons(PY, IN, TC, RE);
}

#endif // INITIALIZE_Neurons_H
