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
#ifndef INITIALIZE_Neurons_H
#define INITIALIZE_Neurons_H
# define M_PI           3.14159265358979323846  /* pi */
#include <vector>
#include <cmath>
#include "Pyramidal_Neuron.h"
#include "Inhibitory_Neuron.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/variate_generator.hpp>
using std::vector;

/****************************************************************************************************/
/*										Typedefs for RNG											*/
/****************************************************************************************************/
typedef boost::random::mt11213b							ENG;	  /* Mersenne Twister		*/
typedef boost::random::normal_distribution<double>		DIST;	  /* Discrete distribution	*/
typedef boost::random::variate_generator<ENG,DIST>		GEN;      /* Variate generator		*/
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/
vector<vector<double>> Get_Param_Values(int Type) {
	extern const int N_e, N_i;
	vector<vector<double>> temp;
	vector<double> mP, dP;
	/* Get the right number of the initial parameter */
	switch(Type) {
	case 1:
		temp = vector<vector<double>>(N_e, vector<double>(3));
		/* E_L, g_L , g_sd */
		mP	= {-60.95, 66.7E-3, 1.75E-3};
		dP	= {0.3, 6.7E-3, 0.1E-3};
		break;
	case 2:
		temp= vector<vector<double>>(N_i, vector<double>(2));
		/* E_L, g_L */
		mP	= {-63.8, 102.5E-3};
		dP	= {0.15, 2.5E-3};
		break;
	}

	/* Set the random number generators */
	vector<GEN> MTRand;
	for(unsigned i=0; i<mP.size(); ++i) {
		MTRand.push_back(GEN(ENG(rand()), DIST(mP[i], dP[i])));
	}

	/* Get the randomly distributed parameters */
	for(unsigned i=0; i<temp.size(); ++i)
		for(unsigned j=0; j<mP.size(); ++j)
			temp[i][j] = MTRand[j]();

	return temp;
}


vector<vector<int>> Get_Connectivity(int Type) {
	extern const int N_e, N_i;
	double length = 5*N_e;
	/* Sigma for the normal distribution */
	double sigma;
	int N_con, Target;
	int N1 = 0 , N2 = 0;
	switch (Type) {
	case 1: /* Pyramidal to Pyramidal */
		N1 = N_e;
		N2 = N_e;
		sigma = 250/length*N2;
		break;
	case 2: /* Pyramidal to Inhibitory */
		N1 = N_e;
		N2 = N_i;
		sigma = 250/length*N2;
		break;
	case 3: /* Inhibitory to Pyramidal */
		N1 = N_i;
		N2 = N_e;
		sigma = 125/length*N2;
		break;
	case 4: /* Inhibitory to Inhibitory */
		N1 = N_i;
		N2 = N_i;
		sigma = 125/length*N2;
		break;
	}
	/* Generate the randomg number generator for the number of connections */
	GEN MTRand_N  = GEN(ENG(rand()), DIST(20, 5));
	/* Generate the random number generator for the target neurons */
	GEN MTRand_T = GEN(ENG(rand()), DIST(0, sigma));

	vector<vector<int>> Connectivity(N2, vector<int>(0));
	for (int i=0; i< N1; ++i) {
		N_con = (abs((int)MTRand_N()));
		for(int j=0; j<N_con; ++j) {
			/* Self connections are not allowed */
			do {
				Target = ((int)MTRand_T()+i)%N2;
				if(Target < 0)
					Target += N2;
			} while (Target == i && ((Type ==1) | (Type==4)));
			Connectivity[Target].push_back(i);
		}
	}
	return Connectivity;
}

void Initialize_Neurons(vector<Pyramidal_Neuron>& PY, vector<Inhibitory_Neuron>& IN) {
	extern const int N_e, N_i;
	/* seed rand() */
	srand(time(NULL));

	/* Generate the randomly distributed parameter values */
	vector<vector<double>> Param_PY = Get_Param_Values(1);
	vector<vector<double>> Param_IN = Get_Param_Values(2);


	/* Generate random connectivity matrices */
	/* For every Neuron [i] they store the index of all neurons it RECEIVES input from */
	vector<vector<int>> Con_ee = Get_Connectivity(1);
	vector<vector<int>> Con_ei = Get_Connectivity(2);
	vector<vector<int>> Con_ie = Get_Connectivity(3);
	vector<vector<int>> Con_ii = Get_Connectivity(4);

	/* Initialize the neurons */
	for (int i=0; i<N_e; ++i)
		PY.push_back(Pyramidal_Neuron (Param_PY[i], Con_ee[i].size(),Con_ie[i].size()));
	for (int i=0; i<N_i; ++i)
		IN.push_back(Inhibitory_Neuron(Param_IN[i], Con_ei[i].size(),Con_ii[i].size()));



	/* Connect the neurons */
	for (int i=0; i<N_e; ++i) {
		/* Get the pointers to the Source neurons */
		vector<Pyramidal_Neuron*> Temp_pPY;
		for(unsigned j=0; j<Con_ee[i].size(); ++j)
			Temp_pPY.push_back(&PY[Con_ee[i][j]]);

		vector<Inhibitory_Neuron*> Temp_pIN;
		for(unsigned j=0; j<Con_ie[i].size(); ++j)
			Temp_pIN.push_back(&IN[Con_ie[i][j]]);

		/* Connect the neurons */
		PY[i].set_Connections(Temp_pPY, Temp_pIN);
	}

	for (int i=0; i<N_i; ++i) {
		/* Get the pointer to the Source neurons */
		vector<Pyramidal_Neuron*> Temp_pPY(0);
		for(unsigned j=0; j<Con_ei[i].size(); ++j)
			Temp_pPY.push_back(&PY[Con_ei[i][j]]);

		vector<Inhibitory_Neuron*> Temp_pIN(0);
		for(unsigned j=0; j<Con_ii[i].size(); ++j)
			Temp_pIN.push_back(&IN[Con_ii[i][j]]);

		/* Connect the neurons */
		IN[i].set_Connections(Temp_pPY, Temp_pIN);
	}
}

#endif // INITIALIZE_Neurons_H
