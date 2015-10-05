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
#ifndef INHIBITORY_NEURON_H
#define INHIBITORY_NEURON_H
#pragma once
#include <cmath>
#include <vector>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include "Pyramidal_Neuron.h"
#include "Thalamocortical_Neuron.h"
using std::vector;
class Pyramidal_Neuron;
class Thalamocortical_Neuron;

/****************************************************************************************************/
/*										Typedefs for RNG											*/
/****************************************************************************************************/
typedef boost::random::mt11213b                    	ENG;    /* Mersenne Twister		*/
typedef boost::random::normal_distribution<double>	DIST;   /* Normal Distribution	*/
typedef boost::random::variate_generator<ENG,DIST> 	GEN;    /* Variate generator	*/
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*									Macro for vector initialization									*/
/****************************************************************************************************/
#ifndef _INIT
#define _INIT(x)	{x, 0.0, 0.0, 0.0, 0.0}
#endif
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*						Implementation of the pyramidal neuron after Bazhenov2002 					*/
/****************************************************************************************************/

class Inhibitory_Neuron {
public:
	Inhibitory_Neuron(void) {}

    Inhibitory_Neuron(vector<double> Param, vector<double>Connectivity)
    : PY_Con(Connectivity[0], NULL),
      IN_Con(Connectivity[1], NULL),
      TC_Con(Connectivity[2], NULL),
      E_L(Param[0]), g_L(Param[1])
	{}

	/* Connect the Neuron with the target neurons */
	void    set_Connections(vector<Pyramidal_Neuron*>&, vector<Inhibitory_Neuron*>&);

	/* Set strength of input current */
	void	set_Input(double I) {Input = I;}

	/* ODE functions */
	void 	set_RK		(int);
	void 	add_RK	 	(void);

private:
    /* Initializes the random number generators */
    void	set_RNG(void);

    /* Noise function */
    double 	noise_xRK 	(int, int) const;
    double 	noise_aRK 	(int) const;

    /* Current functions */
    double 	I_L     (int) const;
    double 	I_Na    (int) const;
    double 	I_K     (int) const;

    /* Synaptic currents */
    double 	I_AMPA  (int) const;
    double 	I_NMDA  (int) const;
    double 	I_GABA  (int) const;

    /* Gating functions */
    double 	alpha_h_Na(int) const;
    double 	alpha_n_K (int) const;
    double 	beta_h_Na (int) const;
    double 	beta_n_K  (int) const;

	/* Random number generators */
	vector<GEN>		MTRands;

	/* Container for noise */
	vector<double>	Rand_vars;

	/* Connections (Neurons that target THIS Neuron*/
	vector<Pyramidal_Neuron*>  PY_Con;
	vector<Inhibitory_Neuron*> IN_Con;
    vector<Thalamocortical_Neuron*> TC_Con;

	/* Paramter constants */
	/* Membrane conductivity */
	const int		C_m		= 1;

	/* Averaged membrane area */
	const double	A_i		= 20E-5;

	/* Reversal potentials */
	const double	E_L		= -63.8;
	const int		E_K		= -90;
	const int		E_Na	= 55;

	const int		E_AMPA	= 0;
	const int		E_NMDA	= 0;
	const int		E_GABA  = -70;

	/* Channel conductivities */
	const double	g_L		= 102.5E-3;
	const double	g_Na	= 35;
	const double	g_K		= 9;

	const double	g_AMPA	= 2.25E-6;
	const double	g_NMDA	= 0.5E-6;
	const double	g_GABA	= 0.165E-6;

	/* Synapse time constants */
	const int		tau_GABA= 10;

	/* Noise parameters */
	const double	dphi	= 5E-3;

	/* Input current */
	double			Input	= 0.0;

    /* Parameters for the RK iteration */
    const vector<double> A = {0.5, 0.5, 1.0, 1.0};
    const vector<double> B = {0.75, 0.75, 0.0, 0.0};

    /* Variables of the neuron */
    vector<double> 	V		= _INIT(E_L),		/* Dendritic membrane voltage			*/
                    h_Na	= _INIT(0.0),		/* inactivation of Na channel			*/
                    n_K		= _INIT(0.0),   	/* activation 	of K  channel			*/
                    s_GABA	= _INIT(0.0);   	/* Fraction of open AMPA channels		*/

    friend class Pyramidal_Neuron;
    friend class Thalamocortical_Neuron;
};
#endif // INHIBITORY_NEURON_H
